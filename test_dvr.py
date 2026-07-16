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
    overwrote with the previous (J=0) run's values before the J=1 run -- that is
    Eq. 7 of J Mol Model 24:235 and it is what dvr.report reproduces. The J=0.txt
    run never had them refreshed, so its ALFAE/GAMAE belong to another molecule.
    """
    r, v = dvr.load_csv("Li_Omega.csv")
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

    # Eq. 7: J=1 spacings + J=0 constants, the pairing the legacy runs used
    alfae, gamae = dvr.eq7_constants(E[1], *dvr.spectro_constants(E[0]))
    print(f"  alfae {alfae:.6e} (ref {REF_ALFAE:.6e}), "
          f"gamae {gamae:.6e} (ref {REF_GAMAE:.6e})")
    assert abs(alfae - REF_ALFAE) < 5e-4 and abs(gamae - REF_GAMAE) < 5e-5


def test_drop_box_states():
    # monotone-decreasing spacing (physical) is kept whole
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0])
    assert len(dvr._drop_box_states(E)) == 5
    # spacing that turns back up (box artifact) is trimmed at the upturn
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0, 42.0])  # last gap 8>7 rises
    assert len(dvr._drop_box_states(E)) == 5


if __name__ == "__main__":
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
