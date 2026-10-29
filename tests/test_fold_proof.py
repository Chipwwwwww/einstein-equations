"""Building blocks of the (incomplete) validated integrator in research/calc/fold_proof.py.

These tests are the ones that caught an earlier, invalid "certificate": python-flint
silently truncates power series at ctx.cap (default 10), which dropped Taylor terms and the
remainder. Every enclosure below must contain the exact solution, for several orders.
"""

import pytest
import sympy as sp
from flint import arb, ctx

from research.calc import fold_proof as fp


def jet(value):
    return fp.J([arb(value)] + [arb(0)] * 5)


def test_series_cap_is_large_enough():
    assert ctx.cap >= 64


@pytest.mark.parametrize("L,h", [(10, 0.25), (18, 0.25), (26, 0.5)])
@pytest.mark.parametrize(
    "name",
    ["exp", "riccati", "nonautonomous"],
)
def test_validated_steps_enclose_exact_solutions(name, L, h):
    rhs, exact = {
        "exp": (lambda y: [y[0] * 0 + 1, y[1]], arb(1).exp()),                 # y' = y
        "riccati": (lambda y: [y[0] * 0 + 1, -(y[1] * y[1])], arb(1) / 2),     # y' = -y^2
        "nonautonomous": (lambda y: [y[0] * 0 + 1, y[0] * y[1]], (arb(1) / 2).exp()),  # y' = t y
    }[name]
    y, _ = fp.integrate(rhs, [jet(0), jet(1)], 0.0, 1.0, h, L)
    assert y[1].c[0].contains(exact)
    assert float(y[1].c[0].rad()) < 1e-4


def test_jet_arithmetic_matches_sympy_derivatives():
    """f(z, d, t) = sin(z d) / sqrt(1 + z^2 + t) at (z, d, t) = (0.3, 0.7, 0): the jet
    coefficients are f, f_z, f_d, f_t, f_zz / 2, f_zd."""
    Z, D, T = sp.symbols("z d t")
    f = sp.sin(Z * D) / sp.sqrt(1 + Z**2 + T)
    at = {Z: sp.Rational(3, 10), D: sp.Rational(7, 10), T: 0}
    expected = [f, f.diff(Z), f.diff(D), f.diff(T), f.diff(Z, 2) / 2, f.diff(Z, D)]
    zj = fp.J([arb("0.3"), arb(1), arb(0), arb(0), arb(0), arb(0)])
    dj = fp.J([arb("0.7"), arb(0), arb(1), arb(0), arb(0), arb(0)])
    tj = fp.J([arb(0), arb(0), arb(0), arb(1), arb(0), arb(0)])
    s, _ = fp.sin_cos(zj * dj)
    val = s * fp.inv(fp.sqrt(1 + zj * zj + tj))
    for c, e in zip(val.c, expected):
        # sp.N(., 50) is the exact value rounded to 50 digits: enlarge it by that rounding
        assert c.overlaps(arb(str(sp.N(e.subs(at), 50)), "1e-45"))
        assert float(c.rad()) < 1e-40


def test_axis_enclosure_contains_the_regular_start():
    """The axis fixed-point boxes are verified (T(X) inside int X) and A(0) = G(0, z0)/2 lies in
    the tightened box for A."""
    zj = fp.J([arb("0.62"), arb(1), arb(0), arb(0), arb(0), arb(0)])
    dj = fp.J([arb("0.7662"), arb(0), arb(1), arb(0), arb(0), arb(0)])
    A, Q, Z = fp.axis_enclosure(zj, dj, fp.S1)
    zero = fp.J([arb(0)] * 6)
    A0 = fp.G(zero + 0, zj, zero + 0, zero + 1, dj) * 0.5
    assert all(a.overlaps(b) for a, b in zip(A.c, A0.c))


@pytest.mark.slow
def test_fold_certificate():
    """The computer-assisted part of (F) (notes/09, H4): Krawczyk verifies a unique fold point
    of the shooting function in the box FOLD_CENTRE +- BOX_RADIUS, with d^2F/dz0^2 > 0,
    dF/dd > 0, dF/dt < 0 there, a positive Jacobi field and an embedded curve. ~10 min."""
    res = fp.certify_fold()
    assert res["verified"]
    F = res["at_fold"]
    assert F.c[4] > 0 and F.c[2] > 0 and F.c[3] < 0  # F_zz/2, F_d, F_t
    for chk in (res["checks_X"], res["checks_K"]):
        assert chk.jacobi_min > 0 and chk.sin_max < 0 and chk.rho_min > 0 and chk.axis_A_max < 0
        assert chk.z_switch_min > 0
    assert res["slope"].overlaps(arb("2.385183", "1e-6"))  # (*) and the direct solve
