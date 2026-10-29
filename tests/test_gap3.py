"""Gap 3 at d = 0.78: the ingredients of the computer-assisted exclusion of a common minimal
surface (research/calc/gap3_neck.py, research/calc/gap3_sweep.py).

Fast checks only: the algebra behind the thin-neck and approach lemmas (SymPy), one piece of
the neck certificate, and one z0 box per exclusion event. The full certificates are
`python -m research.calc.gap3_neck` (15-25 min on 4 cores) and `python -m research.calc.gap3_sweep`
(hours); the slow test below reruns the first."""

import math

import pytest
import sympy as sp

from research.calc import gap3_neck as gn
from research.calc import gap3_sweep as gs


def test_rescaled_neck_system():
    """rho = lam R, z = z_w + lam Z, u = lam U, v = lam V turn the sigma-system into the same
    system for (R, Z, U, V) with the same P = 4 (u psi_z - v psi_rho)/psi; at lam = 0
    (P = 0) the catenoid R = cosh(sigma), Z = s sigma, U = sinh(sigma), V = s solves it."""
    sig, lam, s = sp.symbols("sigma lambda s")
    R, Z, U, V = (sp.Function(n)(sig) for n in "RZUV")
    P = sp.Function("P")(sig)
    rho, z, u, v = lam * R, sp.Symbol("z_w") + lam * Z, lam * U, lam * V
    # the sigma-system for (rho, z, u, v), divided by lam, is the system for (R, Z, U, V)
    sys_full = [sp.diff(rho, sig) - u, sp.diff(z, sig) - v,
                sp.diff(u, sig) - (rho - v * P), sp.diff(v, sig) - u * P]
    sys_scaled = [R.diff(sig) - U, Z.diff(sig) - V, U.diff(sig) - (R - V * P), V.diff(sig) - U * P]
    for a, b in zip(sys_full, sys_scaled):
        assert sp.simplify(a - lam * b) == 0
    cat = {R: sp.cosh(sig), Z: s * sig, U: sp.sinh(sig), V: s}
    for e in sys_scaled:
        e0 = e.subs(P, 0).subs(cat).doit()
        assert sp.simplify(e0.subs(s, 1)) == 0 and sp.simplify(e0.subs(s, -1)) == 0


def test_catenoid_dilation_jacobi_field():
    """At lam = 0 the lam-variation is the dilation of the catenoid. Its normal component
    f = phi / R = -s (1 - sigma tanh(sigma)) solves the catenoid's rotational Jacobi equation
    f'' + 2 sech^2(sigma) f = 0 (conformal coordinates), vanishes at sigma = +-1.1997, and has
    the sign pattern the neck lemma certifies: phi(0) = -s, s phi(+-1.5) > 0."""
    sig, s = sp.symbols("sigma s")
    phi = -s * (sp.cosh(sig) - sig * sp.sinh(sig))  # -V R + U Z on the catenoid
    f = sp.simplify(phi / sp.cosh(sig))
    assert sp.simplify(sp.diff(f, sig, 2) + 2 * f / sp.cosh(sig) ** 2) == 0
    root = sp.nsolve(sig * sp.tanh(sig) - 1, sig, 1.2)
    assert abs(root - 1.19967864) < 1e-7 and root < gn.SIGMA_1
    for sg in (-gn.SIGMA_1, gn.SIGMA_1):
        assert float(phi.subs({s: 1, sig: sg})) > 0.8
    assert phi.subs(sig, 0) == -s


def test_approach_lemma_inequalities():
    """The inequalities behind event N (research/calc/gap3_sweep.py docstring):
    rho^2 - (rho_w + K x/2)^2 = x (1 - K rho_w - K^2 x / 4) with x = rho^2 - rho_w^2;
    arccosh(t) <= log(2 t) for t >= 1; w log(2 r / w) <= 2 r / e for 0 < w <= r."""
    rho, rw, K = sp.symbols("rho rho_w K", positive=True)
    x = rho**2 - rw**2
    assert sp.expand(rho**2 - (rw + K * x / 2) ** 2 - x * (1 - K * rw - K**2 * x / 4)) == 0
    for t in (1.0, 1.001, 2.0, 1e3, 1e8):
        assert math.acosh(t) <= math.log(2 * t) + 1e-15
    w = sp.symbols("w", positive=True)
    g = w * sp.log(2 / w)  # r = 1 by scaling
    crit = sp.solve(sp.diff(g, w), w)
    assert crit == [2 / sp.E] and sp.simplify(g.subs(w, 2 / sp.E) - 2 / sp.E) == 0
    assert sp.diff(g, w, 2).subs(w, 2 / sp.E) < 0


def test_k_bound_dominates_samples():
    """K_N bounds 4 (|psi_z| + |psi_rho|) / psi on [0, 0.02] x [ZW_LO, ZW_HI] (sampled check
    of the interval bound, and its size)."""
    K = float(gs.K_N)
    worst = 0.0
    for i in range(21):
        rho = 0.02 * i / 20
        for k in range(101):
            z = gn.ZW_LO + (gn.ZW_HI - gn.ZW_LO) * k / 100
            psi, pr, pz = 1.0, 0.0, 0.0
            for c in (0.39, -0.39):
                r = math.hypot(rho, z - c)
                psi += 0.25 / r
                pr -= 0.25 * rho / r**3
                pz -= 0.25 * (z - c) / r**3
            worst = max(worst, 4 * (abs(pz) + abs(pr)) / psi)
    assert worst <= K < 1.1 * worst


def test_neck_lemma_on_one_piece():
    """The widest-neck corner of the certified region, near the bottom of the individual
    horizon (z_w ~ 0.2126), both directions."""
    for s in (-1, 1):
        piece = (gn.LAM0 - 2.0**-11, gn.LAM0, 0.2109375, 0.2109375 + 2.0**-9, s)
        assert gn.check_piece(piece) > 0.5


def test_neck_lemma_fails_before_the_catenoid_zero():
    """Falsification: at sigma = 1.0 < 1.1997 the Jacobi field has not changed sign yet, so the
    same piece must NOT be certified there."""
    piece = (gn.LAM0 - 2.0**-11, gn.LAM0, 0.2109375, 0.2109375 + 2.0**-9, -1)
    assert gn.check_piece(piece, sigma1=1.0) is None


@pytest.mark.parametrize("box, event", [
    ((0.570312500, 0.570556640625), "E3"),        # below z_top: dives into the throat
    ((0.58966064453125, 0.58978271484375), "N"),  # contains z_top: thin neck
    ((0.640625, 0.640869140625), "E1"),           # near the minimum of F
])
def test_sweep_box(box, event):
    lo, hi, ev = gs.run_box(box)[:3]
    assert ev == event


@pytest.mark.slow
def test_neck_certificate():
    n, worst = gn.certify_region()
    assert n > 0 and worst > 0


@pytest.mark.parametrize("z", [0.6218058259214554, 0.6377814681180715])
def test_sweep_keeps_true_common_surfaces(z):
    """Falsification: at d = 0.76 < d_c common minimal surfaces exist (axis heights from
    research/calc/bl_fold.py). The z0 box of width 2^-14 containing each must not be excluded,
    while a box between them (F = -0.0077 there) is."""
    from flint import arb
    d, w = arb("0.76"), 2.0**-14
    lo = math.floor(z / w) * w
    assert gs.run_box((lo, lo + w), d)[2] not in gs.EXCLUDED
    if z < 0.63:
        assert gs.run_box((0.62896728515625, 0.62896728515625 + w), d)[2] == "E1"
