"""Sinc-DVR vibrational solver (port of DVR.f, Soares Neto 1995).

Pipeline: CSV (r, V) -> fit extended Rydberg order 6 -> DVR Hamiltonian
(Colbert-Miller) -> eigenvalues -> spectrum, level differences and
spectroscopic constants (we, wexe, weye, alfae, gamae) for J=0 and J=1, then
Dunham's expansion for arbitrary J and the rovibrational partition function
with U, Cv, S and G over a temperature grid.

Everything but the reduced mass is derived from the CSV: grid = r-range,
level count = all bound states below the fitted De. Each run also builds an
article-style LaTeX PDF (tables + figures) unless --no-pdf.

Units: bohr / hartree internally; output in cm-1 (legacy conversion factor).
"""
import argparse
import os
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

CM = 219474.631      # hartree -> cm-1 (legacy literal, keep for regression)
AMU_TO_ME = 1822.88839   # amu -> electron masses (legacy literal)
NPOINTS = 500        # DVR grid size; dense enough, matches legacy reference run
JMAX = 5             # highest J tabulated from Dunham's expansion (--jmax)
C2 = 1.4387768775    # hc/k, cm*K (second radiation constant): eps[cm-1]/T -> beta
RGAS = 8.314462618   # J/(mol K)
TEMPS = (100.0, 200.0, 298.15, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
         900.0, 1000.0)   # default temperature grid for the thermo table

# ------------------------------------------------------------------- units
# The DVR runs in bohr/hartree; a CSV may arrive in anything. Factors -> a.u.
R_TO_BOHR = {"bohr": 1.0, "angstrom": 1.8897261254578281}
E_TO_HARTREE = {"hartree": 1.0, "cm-1": 1.0 / CM, "eV": 1.0 / 27.211386245988,
                "kcal/mol": 1.0 / 627.509474, "kJ/mol": 1.0 / 2625.4996394799}
# Substrings that name a unit in a column header. "au"/"a.u." is bohr for r and
# hartree for E -- both are factor 1, so the collision is harmless.
_R_ALIASES = {"angstrom": ("angstrom", "angst", "ang", "å"),
              "bohr": ("bohr", "a.u.", "au", "a0", "atomic")}
_E_ALIASES = {"cm-1": ("cm-1", "cm^-1", "cm**-1", "cm1", "wavenumber"),
              "kcal/mol": ("kcal",), "kJ/mol": ("kj",), "eV": ("ev",),
              "hartree": ("hartree", "e_h", "eh", "a.u.", "au", "atomic")}


# ---------------------------------------------------------------- potential
def rydberg6(r, De, Re, xe, c1, c2, c3, c4, c5, c6):
    """Extended Rydberg: -De*(1 + sum ck*xx^k)*exp(-c1*xx) + De + xe."""
    xx = r - Re
    poly = 1.0 + c1*xx + c2*xx**2 + c3*xx**3 + c4*xx**4 + c5*xx**5 + c6*xx**6
    return -De * poly * np.exp(-c1 * xx) + De + xe


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _match(field, aliases):
    f = field.lower()
    for unit, keys in aliases.items():
        if any(k in f for k in keys):
            return unit
    return None


def _guess_e_unit(span):
    """From the energy range of a bound diatomic. Decisive: a well is never
    ~1 hartree deep, and never as shallow as ~100 cm-1 with 4+ levels.
    kcal/mol and kJ/mol overlap both and are only honoured from a header."""
    if span < 1.0:
        return "hartree"
    return "cm-1" if span > 100.0 else "eV"


def _guess_r_unit(Re):
    """From the position of the well. Only the extremes decide: 1.2-5.0 reads
    as a real bond either way (3.4 bohr = 1.8 A, or 3.4 A = 6.4 bohr)."""
    if Re < 1.2:            # 1.2 bohr = 0.63 A, shorter than any real bond
        return "angstrom"
    if Re > 5.0:            # 5 A, longer than all but the heaviest vdW dimers
        return "bohr"
    return None


def resolve_units(fields, r, v, r_unit=None, e_unit=None):
    """Units of a CSV: caller's override > header name > magnitude of the data.

    `fields` are the header column names ([] if the file has no header).
    Returns (r_unit, e_unit, how), how mapping "r"/"E" -> how it was decided.
    """
    how = {"r": "given", "E": "given"}
    if r_unit is None and fields:
        r_unit, how["r"] = _match(fields[0], _R_ALIASES), "header"
    if r_unit is None:
        r_unit, how["r"] = _guess_r_unit(r[np.argmin(v)]), "guess"
    if r_unit is None:
        r_unit, how["r"] = "bohr", "ambiguous"
    if e_unit is None and fields:
        e_unit, how["E"] = _match(fields[1], _E_ALIASES), "header"
    if e_unit is None:
        e_unit, how["E"] = _guess_e_unit(np.ptp(v)), "guess"
    return r_unit, e_unit, how


def load_csv(path, r_unit=None, e_unit=None):
    """Two columns r, V -> (r in bohr, V in hartree, units).

    Comma or whitespace separated. Units are auto-detected (see resolve_units)
    and converted; `units` is the (r_unit, e_unit, how) triple that was applied.
    """
    with open(path) as f:
        first = f.readline()
    delim = "," if "," in first else None
    fields = [f.strip() for f in first.strip().lstrip("#").split(delim or None)]
    header = bool(fields) and not _is_number(fields[0])
    data = np.loadtxt(path, delimiter=delim, skiprows=1 if header else 0)
    r, v = data[:, 0], data[:, 1]
    ru, eu, how = resolve_units(fields if header else [], r, v, r_unit, e_unit)
    return r * R_TO_BOHR[ru], v * E_TO_HARTREE[eu], (ru, eu, how)


def fit_rydberg6(r, v):
    """Fit the 9 parameters; returns dict with De, Re, xe, c1..c6.

    Variable projection: for fixed (Re, c1) the model is linear in
    (xe, De, De*c2..De*c6), solved by lstsq; outer 2D least squares over
    (Re, c1). Plain curve_fit on all 9 stalls in a local minimum here.
    """
    def inner(Re, c1):
        xx = r - Re
        e = np.exp(-c1 * xx)
        M = np.column_stack([np.ones_like(xx), 1.0 - (1.0 + c1 * xx) * e]
                            + [-(xx**k) * e for k in range(2, 7)])
        coef, *_ = np.linalg.lstsq(M, v, rcond=None)
        return coef, M @ coef - v

    imin = np.argmin(v)
    xe0, re0 = v[imin], r[imin]
    de0 = max(v[-1], v[0]) - xe0            # asymptote minus minimum
    # harmonic guess: V''(Re) = De*(c1^2 - 2 c2) with c2..c6 = 0
    mask = np.abs(r - re0) < 5 * (r[1] - r[0])
    k = 2.0 * np.polyfit(r[mask] - re0, v[mask] - xe0, 2)[0]
    c10 = np.sqrt(abs(k) / de0) if de0 > 0 else 1.0
    with np.errstate(over="ignore"):
        sol = least_squares(lambda q: inner(*q)[1], [re0, c10],
                            xtol=1e-15, ftol=1e-15, gtol=1e-15)
    Re, c1 = sol.x
    coef, res = inner(Re, c1)
    xe, De = coef[0], coef[1]
    p = {"De": De, "Re": Re, "xe": xe, "c1": c1}
    for k, b in zip(range(2, 7), coef[2:]):
        p[f"c{k}"] = b / De
    p["rms"] = float(np.sqrt(np.mean(res**2)))
    return p


# ---------------------------------------------------------------- DVR core
def build_hamiltonian(A, B, N, mu, Vfunc):
    """Colbert-Miller particle-in-a-box sinc-DVR (BUILD, DVR.f:106-121)."""
    i, j = np.indices((N, N)) + 1
    fct1 = np.pi**2 / (4.0 * mu * (B - A) ** 2)
    with np.errstate(divide="ignore"):
        T = (-1.0) ** (i - j) * fct1 * (
            1.0 / np.sin((i - j) * np.pi / (2 * (N + 1))) ** 2
            - 1.0 / np.sin((i + j) * np.pi / (2 * (N + 1))) ** 2)
    d = np.arange(1, N + 1)
    x = A + (B - A) * d / (N + 1)
    T[np.diag_indices(N)] = fct1 * ((2.0 * (N + 1) ** 2 + 1.0) / 3.0
                                    - 1.0 / np.sin(d * np.pi / (N + 1)) ** 2)
    return T + np.diag(Vfunc(x)), x


def solve(A, B, N, mu, p, J, elow, ehigh):
    """Eigenvalues (hartree, relative to potential minimum) for rotation J.

    Returns every eigenvalue inside (elow, ehigh). With ehigh = De this is the
    set of bound vibrational levels (level count is automatic, not a knob).
    """
    def Vtot(r):  # fitted potential re-zeroed at its minimum + centrifugal
        return (rydberg6(r, p["De"], p["Re"], 0.0, p["c1"], p["c2"], p["c3"],
                         p["c4"], p["c5"], p["c6"])
                + J * (J + 1) / (2.0 * mu * r**2))
    H, _ = build_hamiltonian(A, B, N, mu, Vtot)
    return eigh(H, eigvals_only=True, subset_by_value=(elow, ehigh))


# ------------------------------------------------------ spectroscopy/output
def spectro_constants(E):
    """we, wexe, weye from first 4 levels (EIGCLL, DVR.f:356-364)."""
    if len(E) < 4:
        raise ValueError(f"need >=4 bound levels for we/wexe/weye, got {len(E)}")
    d1, d2, d3 = E[1] - E[0], E[2] - E[0], E[3] - E[0]
    we = (141 * d1 - 93 * d2 + 23 * d3) / 24 * CM
    wexe = (13 * d1 - 11 * d2 + 3 * d3) / 4 * CM
    weye = (3 * d1 - 3 * d2 + d3) / 6 * CM
    return we, wexe, weye


def alfae_gamae(E, we, wexe, weye):
    """alfae, gamae from a J=1 ladder plus the we/wexe/weye of the J=0 ladder.

    This is the legacy two-run methodology: run J=0, read we/wexe/weye off its
    fort.4, paste them into DVR.f:342-346, switch V to J=1 and run again
    (DVR.f:347-353 then combines the pasted J=0 constants with the J=1
    spacings). Here both runs happen in one process, so the paste is just the
    argument. Only defined for J != 0: fed the J=0 ladder itself the formulas
    degenerate into a self-consistency residual (~0), not a constant.

    With E(v,J) = G(v) + [Be - alfae(v+1/2) + gamae(v+1/2)^2] J(J+1) and J=1,
    eliminating G(v) with the J=0 constants leaves two equations for the two
    unknowns -- these closed forms.
    """
    Ecm = E * CM
    d1, d2 = Ecm[1] - Ecm[0], Ecm[2] - Ecm[0]
    alfae = (-12 * d1 + 4 * d2 + 4 * we - 23 * weye) / 8
    gamae = (-2 * d1 + d2 + 2 * wexe - 9 * weye) / 4
    return alfae, gamae


def rot_constants(E0, E1):
    """alfae, gamae from Bv = (E_v(J=1)-E_v(J=0))/2, v = 0,1,2."""
    Bv = (E1[:3] - E0[:3]) / 2.0 * CM
    gamae = (Bv[0] - 2 * Bv[1] + Bv[2]) / 2.0
    alfae = (Bv[0] - Bv[1]) + 2.0 * gamae
    Be = Bv[0] + alfae / 2.0 - gamae / 4.0
    return alfae, gamae, Be, Bv


def Bv_dunham(v, c):
    """B_v = Be - alfae(v+1/2) + gamae(v+1/2)^2, the J(J+1) bracket of eq. 5."""
    x = np.asarray(v) + 0.5
    return c["Be"] - c["alfae"] * x + c["gamae"] * x**2


def dunham(v, J, c):
    """eps(v,J) in cm-1 from Dunham's expansion (Baggio 2017, eq. 5).

    Truncated at exactly the constants the two-run methodology yields, so no J
    beyond 1 costs a diagonalization: DVR runs at J=0 and J=1, this extrapolates
    the rest. The printed eq. 5 drops the (v+1/2)^2 off gamae in its first line;
    its own expanded form and eqs. 9-10 keep it, so follow those.
    """
    x = np.asarray(v) + 0.5
    G = c["we"][0] * x - c["wexe"][0] * x**2 + c["weye"][0] * x**3
    return G + Bv_dunham(v, c) * J * (J + 1)


def vmax_cutoffs(c, De_cm, vcap=999):
    """(vmax, vmax_ok): the two vibrational summation limits of Baggio 2017.

    vmax (eq. 6) is the largest v with eps(v,0) < De. The cubic G(v) turns over
    somewhere and then descends back under De forever, so scan upward and stop
    at the first v that is either unbound or past the turning point -- "largest
    integer for which the inequality holds" only means anything on the rising
    branch.

    vmax_ok (eq. 10) caps that further at the last v with B_v > 0: beyond it the
    Euler-Maclaurin rotational integral diverges, so a partition function must
    stop earlier than dissociation allows (24 -> 17 for NiPc in the paper).
    """
    v = np.arange(vcap + 1)
    e = dunham(v, 0, c)
    ok = e < De_cm
    ok[1:] &= np.diff(e) > 0
    vmax = int(np.argmin(ok)) - 1 if not ok.all() else vcap
    pos = Bv_dunham(v, c) > 0
    vpos = int(np.argmin(pos)) - 1 if not pos.all() else vcap
    return vmax, min(vmax, vpos)


def thermo(T, c, vmax, jcap=20000):
    """Rovibrational partition function and thermodynamics of the diatomic.

        Q(T) = sum_{v=0}^{vmax} sum_J (2J+1) exp(-[eps(v,J)-eps(0,0)]*C2/T)

    The double sum is explicit over the Dunham levels -- no Euler-Maclaurin, so
    nothing here assumes T is high enough for the rotational integral to stand
    in for the sum. `vmax` must still be the eq.-10 cutoff of vmax_cutoffs():
    past it B_v < 0, eps falls with J, and the J sum diverges rather than
    converging (the same divergence Euler-Maclaurin runs into, just visible one
    step earlier).

    Zero of energy is the lowest level eps(0,0), so Q -> 1, U -> 0 and S -> 0 as
    T -> 0. Everything returned is internal (rovibrational) and per mole: no
    translation, no electronic degeneracy, no pV term, so H == U and G == A.

    Returns arrays over T: Q, U (J/mol), Cv (J/mol/K), S (J/mol/K), G (J/mol).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if T.min() <= 0:
        raise ValueError("temperatures must be > 0 K")
    v = np.arange(vmax + 1)
    Bmin = Bv_dunham(v, c).min()
    if Bmin <= 0:
        raise ValueError("B_v <= 0 within v <= vmax; pass the eq.-10 cutoff "
                         "(vmax_cutoffs()[1]), not the dissociation one")
    # Truncate J where the slowest-decaying term (smallest B_v, highest T) is
    # dead: B_v*J(J+1)*C2/T = 40 puts exp(-40) ~ 4e-18 below the v=0,J=0 term.
    Jmax = min(int(np.sqrt(40.0 * T.max() / (C2 * Bmin))) + 10, jcap)
    J = np.arange(Jmax + 1)
    eps = dunham(v[:, None], J[None, :], c) - dunham(0, 0, c)     # (v, J)
    w = (2 * J + 1.0)[None, :] * np.exp(-(C2 / T[:, None, None]) * eps)
    Q = w.sum(axis=(1, 2))
    # <eps> and <eps^2> give U and Cv analytically: U/R = C2*<eps>,
    # Cv/R = (C2/T)^2 * var(eps). No finite differences in T.
    e1 = (w * eps).sum(axis=(1, 2)) / Q
    e2 = (w * eps**2).sum(axis=(1, 2)) / Q
    return {"T": T, "Q": Q, "Jmax": Jmax,
            "U": RGAS * C2 * e1,
            "Cv": RGAS * (C2 / T) ** 2 * (e2 - e1**2),
            "S": RGAS * (np.log(Q) + C2 * e1 / T),
            "G": -RGAS * T * np.log(Q)}


def all_constants(results):
    """Every spectroscopic constant (cm-1) as a dict, for text + PDF report."""
    we0, wexe0, weye0 = spectro_constants(results[0])
    we1, wexe1, weye1 = spectro_constants(results[1])
    # methodology: J=0 constants + J=1 spacings. Bv only supplies Be, which the
    # methodology does not produce.
    alfae, gamae = alfae_gamae(results[1], we0, wexe0, weye0)
    _, _, Be, Bv = rot_constants(results[0], results[1])
    return {"we": (we0, we1), "wexe": (wexe0, wexe1), "weye": (weye0, weye1),
            "Be": Be, "alfae": alfae, "gamae": gamae, "Bv": Bv}


def report(out, p, results, units=None, jmax=JMAX, temps=TEMPS):
    w = out.write
    if units:
        ru, eu, how = units
        w("UNIDADES DO CSV (convertidas para bohr/hartree)\n" + "-" * 46 + "\n")
        w(f" r = {ru}  ({how['r']})\n V = {eu}  ({how['E']})\n\n")
    w("AJUSTE RYDBERG 6 (unidades atomicas)\n" + "-" * 46 + "\n")
    w(f" Re (dist. equilibrio, bohr) = {p['Re']:.12f}\n")
    w(f" xe (energia equilibrio, hartree) = {p['xe']:.12e}\n")
    w(f" De (hartree) = {p['De']:.12e}   ({p['De']*CM:.6f} cm-1)\n")
    for k in ("c1", "c2", "c3", "c4", "c5", "c6"):
        w(f" {k} = {p[k]:.12e}\n")
    w(f" rms do ajuste = {p['rms']:.3e} hartree\n\n")

    ref = spectro_constants(results[0])   # J=0 constants feed alfae/gamae of J!=0
    for J, E in results.items():
        Ecm = E * CM
        w(f"J = {J}\n")
        w("ENERGIA VIBRACIONAL (cm-1)\n" + "-" * 46 + "\n")
        for i, e in enumerate(Ecm, 1):
            w(f"{i:6d}   {e:.13f}\n")
        w("\nESPECTRO VIBRACIONAL (cm-1)\n" + "-" * 46 + "\n")
        for i in range(1, len(Ecm)):
            w(f"{i:6d}   {Ecm[i] - Ecm[0]:.13f}\n")
        w("\nDIFERENCA ENTRE OS NIVEIS VIBRACIONAIS (cm-1)\n" + "-" * 46 + "\n")
        for i in range(len(Ecm) - 1):
            w(f"{i + 1:6d}   {Ecm[i + 1] - Ecm[i]:.13f}\n")
        we, wexe, weye = spectro_constants(E)
        w("\nCONSTANTES ESPECTROSCOPICAS\n" + "-" * 46 + "\n")
        w(f" WE(cm-1)   = {we:.13f}\n")
        w(f" WEXE(cm-1) = {wexe:.13f}\n")
        w(f" WEYE(cm-1) = {weye:.13e}\n")
        alfae, gamae = alfae_gamae(E, *ref)   # sempre com as constantes do J=0
        w(f" ALFAE(cm-1)= {alfae:.13e}\n")
        w(f" GAMAE(cm-1)= {gamae:.13e}\n")
        if J:
            w(" (das energias deste run com WE/WEXE/WEYE do J=0: a metodologia)\n")
        else:
            w(" (residuo de consistencia, nao constantes: em J=0 o fator J(J+1)\n"
              "  zera e o +4WE-23WEYE da formula cancela a parte vibracional,\n"
              "  entao sobra ~0. Longe de zero = constantes e niveis de runs\n"
              "  diferentes. Os alfae/gamae fisicos sao os do J=1.)\n")
        w("\n")

    alfae, gamae, Be, Bv = rot_constants(results[0], results[1])
    w("CONSTANTES ROTACIONAIS (de Bv = (E_v(J=1)-E_v(J=0))/2)\n" + "-" * 46 + "\n")
    for v, b in enumerate(Bv):
        w(f" B{v}(cm-1) = {b:.13e}\n")
    w(f" BE(cm-1)    = {Be:.13e}\n")
    w(f" ALFAE(cm-1) = {alfae:.13e}\n")
    w(f" GAMAE(cm-1) = {gamae:.13e}\n")

    c = all_constants(results)
    De_cm = p["De"] * CM
    vmax, vmax_ok = vmax_cutoffs(c, De_cm)
    n = len(results[0])
    dev = dunham(np.arange(n), 0, c) - results[0] * CM
    i = int(np.argmax(np.abs(dev)))
    w("\nEXPANSAO DE DUNHAM (Baggio 2017, eq. 5)\n" + "-" * 46 + "\n")
    w(" eps(v,J) = we(v+1/2) - wexe(v+1/2)^2 + weye(v+1/2)^3\n")
    w("            + [Be - alfae(v+1/2) + gamae(v+1/2)^2] J(J+1)\n")
    w(" com as constantes do J=0 e o alfae/gamae da metodologia dos dois runs\n")
    w(f" De (cm-1) = {De_cm:.6f}\n")
    w(f" vmax  (eq. 6, eps(v,0) < De)     = {vmax}\n")
    w(f" vmax  (eq. 10, Bv > 0)           = {vmax_ok}\n")
    w(f" desvio Dunham-DVR em J=0 sobre os {n} niveis do DVR: max {dev[i]:+.4f}"
      f" cm-1 em v={i}, rms {np.sqrt(np.mean(dev**2)):.4f}\n")
    jj = range(jmax + 1)
    w("\n eps(v,J) (cm-1)\n")
    w("     v" + "".join(f"{'J=' + str(j):>17}" for j in jj) + "\n")
    for vq in range(vmax + 1):
        w(f"{vq:6d}" + "".join(f"{dunham(vq, j, c):17.6f}" for j in jj) + "\n")

    if temps is not None and len(temps):
        t = thermo(temps, c, vmax_ok)
        w("\nFUNCAO DE PARTICAO E TERMODINAMICA ROVIBRACIONAL\n" + "-" * 46
          + "\n")
        w(" Q(T) = soma_v soma_J (2J+1) exp(-[eps(v,J)-eps(0,0)]*hc/kT),\n")
        w(f" soma dupla explicita ate v = {vmax_ok} (eq. 10, Bv > 0) e "
          f"J = {t['Jmax']}.\n")
        w(" Zero em eps(0,0): Q->1, U->0, S->0 quando T->0. Contribuicao\n")
        w(" rovibracional apenas -- sem translacao, sem eletronica, sem pV\n")
        w(" (entao H = U e G = A).\n\n")
        w(f"{'T (K)':>9}{'Q':>15}{'U (kJ/mol)':>14}{'Cv (J/mol/K)':>15}"
          f"{'S (J/mol/K)':>14}{'G (kJ/mol)':>14}\n")
        for i, tt in enumerate(t["T"]):
            w(f"{tt:9.2f}{t['Q'][i]:15.6e}{t['U'][i] / 1000:14.6f}"
              f"{t['Cv'][i]:15.6f}{t['S'][i]:14.6f}"
              f"{t['G'][i] / 1000:14.6f}\n")


# --------------------------------------------------------------------- main
def _drop_box_states(E):
    """Trim spurious near-threshold levels. On a finite grid the outer wall
    quantizes the continuum edge, adding fake levels just below De whose
    spacing turns back up. For a single-well potential the true spacing
    decreases monotonically toward dissociation, so keep levels only up to the
    first point where the spacing stops shrinking."""
    if len(E) < 3:
        return E
    dE = np.diff(E)
    rising = np.nonzero(dE[1:] > dE[:-1])[0]     # first spacing increase
    return E if rising.size == 0 else E[:rising[0] + 2]


def run_analysis(csv_path, mu, r_unit=None, e_unit=None):
    """Fit + diagonalize a (r, V) CSV. Units, grid (= r-range) and level count
    (= bound states below De) are all automatic.
    Returns (r, v, p, results, A, B, units) with r/v already in bohr/hartree."""
    r, v, units = load_csv(csv_path, r_unit, e_unit)
    A, B = r[0], r[-1]
    p = fit_rydberg6(r, v)
    results = {J: _drop_box_states(solve(A, B, NPOINTS, mu, p, J, 0.0, p["De"]))
               for J in (0, 1)}
    return r, v, p, results, A, B, units


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", help="CSV com colunas r (bohr), V (hartree)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--mass", type=float,
                   help="massa reduzida em massas de eletron")
    g.add_argument("--mass-amu", type=float, dest="mass_amu",
                   help="massa reduzida em amu (convertida x1822.88839)")
    ap.add_argument("--r-unit", choices=sorted(R_TO_BOHR), dest="r_unit",
                    help="unidade da coluna r (padrao: detecta do CSV)")
    ap.add_argument("--e-unit", choices=sorted(E_TO_HARTREE), dest="e_unit",
                    help="unidade da coluna V (padrao: detecta do CSV)")
    ap.add_argument("--jmax", type=int, default=JMAX,
                    help=f"maior J tabulado por Dunham (padrao: {JMAX}). O DVR "
                         "so roda em J=0 e 1; o resto e extrapolado pela eq. 5")
    ap.add_argument("--temps", type=str, default=None,
                    help="temperaturas (K) da tabela termodinamica, separadas "
                         f"por virgula (padrao: {TEMPS[0]:g}..{TEMPS[-1]:g}). "
                         "Vazio (--temps '') pula a tabela")
    ap.add_argument("--no-pdf", action="store_true",
                    help="pula geracao do relatorio LaTeX/PDF")
    args = ap.parse_args(argv)

    temps = TEMPS if args.temps is None else tuple(
        float(x) for x in args.temps.replace(",", " ").split())

    if args.mass_amu is not None:
        mu = args.mass_amu * AMU_TO_ME
    elif args.mass is not None:
        mu = args.mass
    else:
        ap.error("informe a massa reduzida: --mass (massas de eletron) ou "
                 "--mass-amu (amu)")

    r, v, p, results, A, B, units = run_analysis(args.csv, mu,
                                                 r_unit=args.r_unit,
                                                 e_unit=args.e_unit)

    base = os.path.splitext(os.path.basename(args.csv))[0]
    txt = f"{base}_out.txt"
    with open(txt, "w") as f:
        report(f, p, results, units, args.jmax, temps)
    report(sys.stdout, p, results, units, args.jmax, temps)
    ru, eu, how = units
    print(f"\n[dvr] unidades do CSV: r = {ru} ({how['r']}), "
          f"V = {eu} ({how['E']})")
    if how["r"] == "ambiguous":
        print(f"[dvr] AVISO: nada no CSV diz a unidade de r, e Re = "
              f"{r[np.argmin(v)]:.2f} bohr e plausivel tanto em bohr quanto em "
              "angstrom. Assumi bohr. Se estiver errado, rode com "
              "--r-unit angstrom (ou nomeie a coluna, ex. 'radii_bohr').",
              file=sys.stderr)
    if how["E"] == "guess":
        print(f"[dvr] AVISO: nada no CSV nomeia a unidade de V; deduzi {eu} da "
              "profundidade do poco. O criterio erra em poco raso (um poco de "
              "<1 eV dado em eV le como hartree, 27x fundo demais). Se estiver "
              "errado, rode com --e-unit (ou nomeie a coluna, ex. 'E (eV)').",
              file=sys.stderr)
    print(f"[dvr] texto: {txt}   niveis: J=0 {len(results[0])}, "
          f"J=1 {len(results[1])}")

    if not args.no_pdf:
        import dvrreport
        dvrreport.build_report(base, r, v, p, results, A, B, mu, units,
                               args.jmax, temps)

    return p, results


if __name__ == "__main__":
    main()
