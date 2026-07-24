"""Checks for dvr.thermo(): the double sum and the analytic U/Cv/S.

Run: python3 test_thermo.py
"""
import numpy as np

import dvr


def _consts(we, wexe=0.0, weye=0.0, Be=1.0, alfae=0.0, gamae=0.0):
    """The constants dict dunham()/Bv_dunham() read, with knobs set to zero."""
    return {"we": (we, we), "wexe": (wexe, wexe), "weye": (weye, weye),
            "Be": Be, "alfae": alfae, "gamae": gamae}


def test_factorizes_when_harmonic():
    """With every anharmonic/coupling term zero, eps(v,J)-eps(0,0) is
    we*v + Be*J(J+1), so Q must be the product of the two 1D sums. Catches a
    wrong degeneracy, a wrong energy zero or a botched broadcast."""
    we, Be, vmax = 400.0, 1.5, 30
    c = _consts(we, Be=Be)
    T = np.array([100.0, 298.15, 1000.0])
    t = dvr.thermo(T, c, vmax)

    v, J = np.arange(vmax + 1), np.arange(t["Jmax"] + 1)
    b = dvr.C2 / T[:, None]
    qv = np.exp(-b * we * v).sum(axis=1)
    qr = ((2 * J + 1) * np.exp(-b * Be * J * (J + 1))).sum(axis=1)
    assert np.allclose(t["Q"], qv * qr, rtol=1e-12), t["Q"] / (qv * qr)


def test_cv_is_dU_dT():
    """Cv comes from the variance of eps, U from its mean -- two different
    moments of the same weights. If either is wrong they stop agreeing."""
    c = _consts(350.0, wexe=2.5, weye=-4e-3, Be=0.8, alfae=5e-3, gamae=-6e-5)
    T = np.array([150.0, 298.15, 700.0])
    h = 1e-3 * T
    up = dvr.thermo(T + h, c, 25)["U"]
    dn = dvr.thermo(T - h, c, 25)["U"]
    assert np.allclose(dvr.thermo(T, c, 25)["Cv"], (up - dn) / (2 * h),
                       rtol=1e-6)


def test_limits():
    """Zero at eps(0,0): only the ground level survives as T -> 0."""
    c = _consts(400.0, Be=1.5)
    # 0.1 K, not 1 K: the first rotational quantum is only 3*Be = 4.5 cm-1, and
    # kT at 1 K is 0.7 cm-1 -- J=1 is still worth 1% of Q there.
    t = dvr.thermo([0.1, 100.0, 500.0, 1000.0], c, 30)
    assert abs(t["Q"][0] - 1.0) < 1e-9          # Q -> 1
    assert abs(t["S"][0]) < 1e-6 and abs(t["U"][0]) < 1e-6
    assert np.all(np.diff(t["Q"]) > 0)          # Q rises with T
    assert np.all(t["Cv"] > 0) and np.all(t["S"] >= 0)
    # G = U - TS holds by construction; check it survived the arithmetic
    assert np.allclose(t["G"], t["U"] - t["T"] * t["S"], atol=1e-9)


def test_rejects_negative_Bv():
    """Past the eq.-10 cutoff B_v < 0 and the J sum diverges. Must refuse."""
    c = _consts(400.0, Be=1.0, alfae=0.5)       # B_v < 0 from v = 2 on
    try:
        dvr.thermo([300.0], c, 10)
    except ValueError as e:
        assert "B_v" in str(e)
    else:
        raise AssertionError("accepted a vmax with B_v <= 0")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print("THERMO OK")
