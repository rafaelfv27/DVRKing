"""Regression gate vs legacy fort.4 / fort.97 (reference run of DVR.exe).

potential.csv holds the 500 (r, V) pairs logged by the legacy run in fort.98
(the fort.* files were deleted); they are exactly the grid points of the
reference run.

Part A (exact, 1e-6 cm-1): DVR core with V taken directly from those pairs.
Legacy bug: fort.4 levels were printed with the literal 219474.631 WITHOUT D0
(DVR.f:281), i.e. a REAL*4 factor = 219474.625; the constants use
219474.631D0. REF levels are rescaled by CM/float32(CM) before comparing.
Part B (fit-limited, 0.1 cm-1): full CSV -> Rydberg-6 fit -> levels.
The legacy data deviates from a perfect Rydberg-6 by rms 1.2e-7 hartree, so
the fitted potential cannot reproduce fort.4 tighter than ~0.05 cm-1.
"""
import os
import re

import numpy as np
from scipy.linalg import eigh

import dvr
import dvrreport

F4SCALE = dvr.CM / float(np.float32(dvr.CM))   # undo single-precision print
REF = np.array([9.7211221261708740, 28.870890234651771, 47.635691885033019,
                66.011545713639862, 83.994308177884406]) * F4SCALE
REF_WE = 19.530919626995409
A, B, N, MU = 13.04, 30.20, 500, 218947.07152


MU_OMEGA = 12522.30592
# test/J={0,1}.txt: legacy DVR.exe runs of Li_Omega.csv at MU_OMEGA.
REF_WE_OMEGA = {0: 397.34944605871459, 1: 397.34003799225286}
# ALFAE/GAMAE of J=1.txt only; J=0.txt prints garbage there (see below)
REF_ALFAE, REF_GAMAE = 4.70736802133647769e-3, -6.12243828309599891e-5


def _ref_levels(path):
    txt = open(path, encoding="latin-1").read().split("ESPECTRO")[0]
    return np.array([float(x) for x in re.findall(r"^\s+\d+\s+([\d.]+)\s*$",
                                                  txt, re.M)])


def test_li_omega_vs_reference():
    """Fit-limited like Part B, but looser: the Li_Omega points scatter ~3 cm-1
    around the best Rydberg-6 (fit rms 8.5e-6 hartree), so the error grows with
    v (0.04 cm-1 at v=1, 4.9 at v=20). Anything past ~5 means a real change.

    ALFAE/GAMAE are compared against J=1.txt only. DVR.f:347-353 builds them
    from the we/wexe/weye hardcoded at DVR.f:342-346, which the legacy workflow
    overwrote with the previous (J=0) run's values before switching V to J=1 and
    running again -- that two-run methodology is what dvr.report reproduces. The
    J=0.txt run never had them refreshed, so its ALFAE/GAMAE belong to another
    molecule; dvr.report no longer prints them for J=0.
    """
    r, v, _ = dvr.load_csv("Li_Omega.csv")
    p = dvr.fit_rydberg6(r, v)
    E = {}
    for J in (0, 1):
        ref = _ref_levels(f"test/J={J}.txt")
        E[J] = dvr.solve(r[0], r[-1], N, MU_OMEGA, p, J, 0.0, p["De"])
        err = np.abs(E[J][:len(ref)] * dvr.CM - ref)
        we = dvr.spectro_constants(E[J])[0]
        print(f"Li_Omega vs test/J={J}.txt: |dE| v=1 {err[0]:.3f}, max "
              f"{err.max():.3f}, |dwe| {abs(we - REF_WE_OMEGA[J]):.3f} cm-1")
        assert err[0] < 0.1 and err.max() < 5.0
        assert abs(we - REF_WE_OMEGA[J]) < 0.2

    # methodology: J=1 spacings + J=0 constants, the pairing the legacy runs used
    alfae, gamae = dvr.alfae_gamae(E[1], *dvr.spectro_constants(E[0]))
    print(f"  alfae {alfae:.6e} (ref {REF_ALFAE:.6e}), "
          f"gamae {gamae:.6e} (ref {REF_GAMAE:.6e})")
    assert abs(alfae - REF_ALFAE) < 5e-4 and abs(gamae - REF_GAMAE) < 5e-5


def test_units():
    """Detection + conversion. The two committed CSVs are both bohr/hartree and
    must survive untouched: potential.csv says so in its header, Li_Omega.csv
    says nothing and falls back to the bohr default."""
    r, v, (ru, eu, how) = dvr.load_csv("potential.csv")
    assert (ru, eu) == ("bohr", "hartree") and how == {"r": "header",
                                                       "E": "header"}
    assert r[0] == 13.07425149700598688 and v[0] == 5.543277100975916964e-03

    r, v, (ru, eu, how) = dvr.load_csv("Li_Omega.csv")
    assert (ru, eu) == ("bohr", "hartree")     # Re=3.4 is ambiguous -> default
    assert how == {"r": "ambiguous", "E": "guess"}
    assert r[0] == 1.511781 and v[0] == 0.043473447219

    # an angstrom/cm-1 file converts; same physics, different numbers on disk
    tmp = "test/_units.csv"
    ang, cm = dvr.R_TO_BOHR["angstrom"], dvr.CM
    np.savetxt(tmp, np.column_stack([r / ang, v * cm]), delimiter=",",
               header="r (Angstrom),E (cm-1)", comments="")
    r2, v2, (ru2, eu2, how2) = dvr.load_csv(tmp)
    assert (ru2, eu2) == ("angstrom", "cm-1") and how2["r"] == "header"
    assert np.allclose(r2, r) and np.allclose(v2, v)
    os.remove(tmp)

    # decisive magnitudes with no header at all
    assert dvr._guess_e_unit(0.05) == "hartree"
    assert dvr._guess_e_unit(11000.0) == "cm-1"
    assert dvr._guess_r_unit(0.74) == "angstrom"   # H2 in angstrom
    assert dvr._guess_r_unit(15.0) == "bohr"       # vdW grid in bohr
    assert dvr._guess_r_unit(3.4) is None          # honestly undecidable

    # an explicit override beats both header and data
    _, _, (ru3, _, how3) = dvr.load_csv("potential.csv", r_unit="angstrom")
    assert ru3 == "angstrom" and how3["r"] == "given"
    print("units: header, guess, override and conversion all OK")


def test_dunham():
    """Dunham (Baggio 2017 eq. 5) must reproduce what it was built from.

    we/wexe/weye come from finite differences of the first four J=0 levels, so
    the cubic G(v) has to return those four spacings exactly; Bv(v) is the
    parabola through B_0..B_2, so eps(v,1)-eps(v,0) has to return 2*B_v. Both
    hold to roundoff -- anything worse means the expansion or the constants
    drifted apart. Extrapolation beyond v=3 is checked only for sanity: the
    cubic is fit to 4 levels and cannot track a real anharmonic ladder to
    dissociation.
    """
    r, v, _ = dvr.load_csv("Li_Omega.csv")
    p = dvr.fit_rydberg6(r, v)
    results = {J: dvr.solve(r[0], r[-1], N, MU_OMEGA, p, J, 0.0, p["De"])
               for J in (0, 1)}
    c = dvr.all_constants(results)
    E0 = results[0] * dvr.CM

    g = dvr.dunham(np.arange(4), 0, c) - dvr.dunham(0, 0, c)
    assert np.abs(g - (E0[:4] - E0[0])).max() < 1e-8

    Bv = (results[1][:3] - results[0][:3]) / 2.0 * dvr.CM
    eps = dvr.dunham(np.arange(3), 1, c) - dvr.dunham(np.arange(3), 0, c)
    print(f"dunham: G(v) exact to {np.abs(g - (E0[:4] - E0[0])).max():.2e}, "
          f"2Bv to {np.abs(eps - 2 * Bv).max():.2e} cm-1")
    # the two routes to alfae/gamae -- the J=0-constants methodology and the
    # parabola through B_0..B_2 -- turn out to be the same identity, agreeing to
    # roundoff, so Bv_dunham reproduces Bv and not merely approximates it
    assert np.abs(eps - 2 * Bv).max() < 1e-9

    vmax, vmax_ok = dvr.vmax_cutoffs(c, p["De"] * dvr.CM)
    print(f"        vmax {vmax}, vmax(Bv>0) {vmax_ok}, "
          f"DVR bound levels {len(E0)}")
    assert 0 < vmax < 999 and vmax_ok <= vmax
    assert dvr.dunham(vmax, 0, c) < p["De"] * dvr.CM
    assert dvr.Bv_dunham(np.arange(vmax_ok + 1), c).min() > 0


def test_report_without_figures():
    """make_tex must omit every figure whose path is missing.

    An empty \\includegraphics{} is a fatal LaTeX error. When plotting fails
    build_report passes figs={}, so emitting the blocks anyway would take the
    tables down with them and end the run without any PDF -- which is what the
    old empty-string fallback did, despite promising the opposite.
    """
    r, v, _ = dvr.load_csv("Li_Omega.csv")
    p = dvr.fit_rydberg6(r, v)
    results = {J: dvr.solve(r[0], r[-1], N, MU_OMEGA, p, J, 0.0, p["De"])
               for J in (0, 1)}
    tex = dvrreport.make_tex("test/_figless", p, results, {}, MU_OMEGA)
    body = open(tex, encoding="utf-8").read()
    os.remove(tex)
    assert "includegraphics" not in body
    # the tables are the point: they must survive a plotting failure
    assert r"\begin{tabular}" in body and r"\begin{longtable}" in body
    print("report: figs={} -> no includegraphics, tables intact")


def test_esc():
    """base comes from the CSV filename, so it can carry any LaTeX special."""
    assert dvrreport._esc("run_50%_scan") == r"run\_50\%\_scan"
    assert dvrreport._esc("a$b#c") == r"a\$b\#c"
    assert dvrreport._esc("x{y}") == r"x\{y\}"
    assert dvrreport._esc("back\\slash") == r"back\textbackslash{}slash"
    print("esc: LaTeX specials escaped")


def test_drop_box_states():
    # monotone-decreasing spacing (physical) is kept whole
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0])
    assert len(dvr._drop_box_states(E)) == 5
    # spacing that turns back up (box artifact) is trimmed at the upturn
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0, 42.0])  # last gap 8>7 rises
    assert len(dvr._drop_box_states(E)) == 5


if __name__ == "__main__":
    test_units()
    test_esc()
    test_drop_box_states()
    test_report_without_figures()
    test_li_omega_vs_reference()
    test_dunham()
    # Part A: kinetic matrix + eigh, exact potential values
    rv = np.loadtxt("potential.csv", delimiter=",")
    H, x = dvr.build_hamiltonian(A, B, N, MU, lambda _: rv[:, 1])
    assert np.allclose(x, rv[:, 0])
    E = eigh(H, eigvals_only=True, subset_by_index=(0, 9)) * dvr.CM
    err_core = np.abs(E[:5] - REF).max()
    we_core = dvr.spectro_constants(E / dvr.CM)[0]
    print(f"core: max|E-fort.4| = {err_core:.3e}  |we-fort.97| = {abs(we_core-REF_WE):.3e} cm-1")
    assert err_core < 1e-6 and abs(we_core - REF_WE) < 1e-6

    # Part B: full pipeline through the Rydberg-6 fit
    p, results = dvr.main(["potential.csv", "--mass", str(MU), "--no-pdf"])
    err_fit = np.abs(results[0][:5] * dvr.CM - REF).max()
    print(f"\nfit pipeline: max|E-fort.4| = {err_fit:.3e} cm-1 (fit rms {p['rms']:.2e} hartree)")
    assert err_fit < 0.1
    print("REGRESSION OK")
