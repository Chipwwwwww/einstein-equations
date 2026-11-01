"""Gap 3 lemmas (research/calc/bl_global.py): barrier certificates for Brill-Lindquist data
at d = 0.78, and the desingularised form of the meridian ODE near the axis."""

import sympy as sp

from research.calc import bl_global as bg


def test_throat_barrier_is_certified():
    """Lemma G1: spheres of radius s <= 0.155 about a puncture have H < 0 (normal away from
    the puncture), by interval arithmetic over (s, theta) boxes including s = 0."""
    assert bg.certify_throat(d="0.78", s_star="0.155") < 0


def test_far_barrier_is_certified():
    """Lemma G2: spheres |x| = R >= 1.05 have H > 0 (boxes up to R = 2, analytic beyond)."""
    assert bg.certify_far(d="0.78", R_star="1.05") > 0


def test_throat_quantity_limit_at_the_puncture():
    """s psi^2 H -> 2 - 4 = -2 as s -> 0, for every direction."""
    for th in (0.0, 1.0, 3.0):
        assert abs(bg.throat_quantity(1e-12, th, 0.78) + 2) < 1e-9


def test_desingularised_meridian_system():
    """With alpha' = -sin(alpha)/rho + G, G = 4 (-sin(alpha) psi_rho + cos(alpha) psi_z)/psi,
    u = rho cos(alpha), v = rho sin(alpha) and d/dsigma = rho d/ds:
    u' = rho - v P, v' = u P with P = 4 (u psi_z - v psi_rho)/psi, and
    d/dsigma (u^2 + v^2 - rho^2) = 0 on the system."""
    rho, z, a = sp.symbols("rho z alpha", positive=True)
    psi = sp.Function("psi")(rho, z)
    pr, pz = sp.diff(psi, rho), sp.diff(psi, z)
    G = 4 * (-sp.sin(a) * pr + sp.cos(a) * pz) / psi
    da = -sp.sin(a) / rho + G  # d alpha / ds
    ds = {rho: sp.cos(a), z: sp.sin(a), a: da}
    u, v = rho * sp.cos(a), rho * sp.sin(a)
    dsig = lambda f: rho * sum(sp.diff(f, k) * ds[k] for k in (rho, z, a))  # noqa: E731
    P = 4 * (u * pz - v * pr) / psi
    assert sp.simplify(dsig(u) - (rho - v * P)) == 0
    assert sp.simplify(dsig(v) - u * P) == 0
    U, V, R, PP = sp.symbols("U V R P")
    assert sp.expand(2 * U * (R - V * PP) + 2 * V * (U * PP) - 2 * R * U) == 0


def test_barriers_over_the_interval_of_separations():
    """Theorem F' (notes/12 section 7): both barriers hold for every d in [0.765625, 0.78125]
    at once (d an interval). The throat radius must shrink to 0.15: 0.155 fails at d = 0.766."""
    d = ["0.765625", "0.78125"]
    assert bg.certify_throat(d=d, s_star="0.15") < 0
    assert bg.certify_far(d=d, R_star="1.05") > 0
