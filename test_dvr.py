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
import numpy as np
from scipy.linalg import eigh

import dvr

F4SCALE = dvr.CM / float(np.float32(dvr.CM))   # undo single-precision print
REF = np.array([9.7211221261708740, 28.870890234651771, 47.635691885033019,
                66.011545713639862, 83.994308177884406]) * F4SCALE
REF_WE = 19.530919626995409
A, B, N, MU = 13.04, 30.20, 500, 218947.07152


def test_drop_box_states():
    # monotone-decreasing spacing (physical) is kept whole
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0])
    assert len(dvr._drop_box_states(E)) == 5
    # spacing that turns back up (box artifact) is trimmed at the upturn
    E = np.array([0.0, 10.0, 19.0, 27.0, 34.0, 42.0])  # last gap 8>7 rises
    assert len(dvr._drop_box_states(E)) == 5


if __name__ == "__main__":
    test_drop_box_states()
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
