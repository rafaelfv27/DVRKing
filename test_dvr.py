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


def test_drop_box_states():
    # monotone-decreasing spacing (physical) is kept whole
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0])
    assert len(dvr._drop_box_states(E)) == 5
    # spacing that turns back up (box artifact) is trimmed at the upturn
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0, 42.0])  # last gap 8>7 rises
    assert len(dvr._drop_box_states(E)) == 5


if __name__ == "__main__":
    test_units()
    test_drop_box_states()
    test_li_omega_vs_reference()
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
