"""Sinc-DVR vibrational solver (port of DVR.f, Soares Neto 1995).

Pipeline: CSV (r, V) -> fit extended Rydberg order 6 -> DVR Hamiltonian
(Colbert-Miller) -> eigenvalues -> spectrum, level differences and
spectroscopic constants (we, wexe, weye, alfae, gamae) for J=0 and J=1.

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


# ---------------------------------------------------------------- potential
def rydberg6(r, De, Re, xe, c1, c2, c3, c4, c5, c6):
    """Extended Rydberg: -De*(1 + sum ck*xx^k)*exp(-c1*xx) + De + xe."""
    xx = r - Re
    poly = 1.0 + c1*xx + c2*xx**2 + c3*xx**3 + c4*xx**4 + c5*xx**5 + c6*xx**6
    return -De * poly * np.exp(-c1 * xx) + De + xe


def load_csv(path):
    """Two columns r (bohr), V (hartree). Comma or whitespace separated."""
    with open(path) as f:
        first = f.readline()
    delim = "," if "," in first else None
    skip = 0 if any(c.isdigit() for c in first.split(delim or None)[0]) else 1
    data = np.loadtxt(path, delimiter=delim, skiprows=skip)
    return data[:, 0], data[:, 1]


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


def rot_constants(E0, E1):
    """alfae, gamae from Bv = (E_v(J=1)-E_v(J=0))/2, v = 0,1,2."""
    Bv = (E1[:3] - E0[:3]) / 2.0 * CM
    gamae = (Bv[0] - 2 * Bv[1] + Bv[2]) / 2.0
    alfae = (Bv[0] - Bv[1]) + 2.0 * gamae
    Be = Bv[0] + alfae / 2.0 - gamae / 4.0
    return alfae, gamae, Be, Bv


def all_constants(results):
    """Every spectroscopic constant (cm-1) as a dict, for text + PDF report."""
    we0, wexe0, weye0 = spectro_constants(results[0])
    we1, wexe1, weye1 = spectro_constants(results[1])
    alfae, gamae, Be, Bv = rot_constants(results[0], results[1])
    return {"we": (we0, we1), "wexe": (wexe0, wexe1), "weye": (weye0, weye1),
            "Be": Be, "alfae": alfae, "gamae": gamae, "Bv": Bv}


def report(out, p, results, mu):
    w = out.write
    w("AJUSTE RYDBERG 6 (unidades atomicas)\n" + "-" * 46 + "\n")
    w(f" Re (dist. equilibrio, bohr) = {p['Re']:.12f}\n")
    w(f" xe (energia equilibrio, hartree) = {p['xe']:.12e}\n")
    w(f" De (hartree) = {p['De']:.12e}   ({p['De']*CM:.6f} cm-1)\n")
    for k in ("c1", "c2", "c3", "c4", "c5", "c6"):
        w(f" {k} = {p[k]:.12e}\n")
    w(f" rms do ajuste = {p['rms']:.3e} hartree\n\n")

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
        w(f" WEYE(cm-1) = {weye:.13e}\n\n")

    alfae, gamae, Be, Bv = rot_constants(results[0], results[1])
    w("CONSTANTES ROTACIONAIS (de Bv = (E_v(J=1)-E_v(J=0))/2)\n" + "-" * 46 + "\n")
    for v, b in enumerate(Bv):
        w(f" B{v}(cm-1) = {b:.13e}\n")
    w(f" BE(cm-1)    = {Be:.13e}\n")
    w(f" ALFAE(cm-1) = {alfae:.13e}\n")
    w(f" GAMAE(cm-1) = {gamae:.13e}\n")

    # Eq. 7 of Silva et al., J Mol Model 24:235 (2018) — same finite-diff
    # scheme the reference article uses (J=1 spacings + J=0 we/wexe/weye)
    we, wexe, weye = spectro_constants(results[0])
    E1 = results[1] * CM
    d11, d21 = E1[1] - E1[0], E1[2] - E1[0]
    alfae7 = (-12 * d11 + 4 * d21 + 4 * we - 23 * weye) / 8
    gamae7 = (-2 * d11 + d21 + 2 * wexe - 9 * weye) / 4
    w("\nCONSTANTES ROTACIONAIS (Eq. 7 do artigo, J Mol Model 24:235)\n" + "-" * 46 + "\n")
    w(f" ALFAE(cm-1) = {alfae7:.13e}\n")
    w(f" GAMAE(cm-1) = {gamae7:.13e}\n")


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


def run_analysis(csv_path, mu, npoints=NPOINTS):
    """Fit + diagonalize a (r, V) CSV. Grid = r-range, levels = bound states
    below De (both auto). Returns (r, v, p, results, A, B)."""
    r, v = load_csv(csv_path)
    A, B = r[0], r[-1]
    p = fit_rydberg6(r, v)
    results = {J: _drop_box_states(solve(A, B, npoints, mu, p, J, 0.0, p["De"]))
               for J in (0, 1)}
    return r, v, p, results, A, B


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", help="CSV com colunas r (bohr), V (hartree)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--mass", type=float,
                   help="massa reduzida em massas de eletron")
    g.add_argument("--mass-amu", type=float, dest="mass_amu",
                   help="massa reduzida em amu (convertida x1822.88839)")
    ap.add_argument("--no-pdf", action="store_true",
                    help="pula geracao do relatorio LaTeX/PDF")
    args = ap.parse_args(argv)

    if args.mass_amu is not None:
        mu = args.mass_amu * AMU_TO_ME
    elif args.mass is not None:
        mu = args.mass
    else:
        ap.error("informe a massa reduzida: --mass (massas de eletron) ou "
                 "--mass-amu (amu)")

    r, v, p, results, A, B = run_analysis(args.csv, mu)

    base = os.path.splitext(os.path.basename(args.csv))[0]
    txt = f"{base}_out.txt"
    with open(txt, "w") as f:
        report(f, p, results, mu)
    report(sys.stdout, p, results, mu)
    print(f"\n[dvr] texto: {txt}   niveis: J=0 {len(results[0])}, "
          f"J=1 {len(results[1])}")

    if not args.no_pdf:
        import dvrreport
        dvrreport.build_report(base, r, v, p, results, A, B, mu)

    return p, results


if __name__ == "__main__":
    main()
