"""The solutions and limits claimed in notes/03-05 really satisfy the field equations.

Units are G = c = 1, except that Newton's G is kept as `G_N` where it multiplies
matter, so the Friedmann equations can be read off directly.
"""

import pytest
import sympy as sp

from einstein import metrics, nonzero
from einstein.metrics import M, Lambda, Phi, a, epsilon, k, r, t, x, y, z

G_N, rho, p = sp.symbols("G_N rho p", positive=True)
ZERO = sp.zeros(4, 4)


def first_order(expr):
    """Coefficient of epsilon in the Taylor expansion of expr around epsilon = 0."""
    return sp.simplify(sp.diff(expr, epsilon).subs(epsilon, 0).doit())


def test_flat_space_in_spherical_coordinates():
    # Iterating a SymPy N-dim array yields sub-arrays, not components, so use nonzero().
    st = metrics.minkowski_spherical()
    assert len(nonzero(st.christoffel)) == 9, "curvilinear coordinates give nonzero Christoffels"
    assert nonzero(st.riemann) == [], "but flat space has zero curvature"


def test_schwarzschild_is_vacuum():
    assert metrics.schwarzschild().ricci == ZERO


def test_schwarzschild_kretschmann():
    # Finite at the horizon r = 2M, divergent at r = 0 (notes/02, notes/04).
    assert sp.simplify(metrics.schwarzschild().kretschmann - 48 * M**2 / r**6) == 0


def test_de_sitter_is_lambda_vacuum():
    st = metrics.de_sitter_static()
    assert (st.einstein + Lambda * st.metric).applyfunc(sp.simplify) == ZERO
    assert (st.ricci - Lambda * st.metric).applyfunc(sp.simplify) == ZERO


def test_flrw_friedmann_equations():
    st = metrics.flrw()
    G, g = st.einstein, st.metric
    H = a.diff(t) / a
    assert sp.simplify(G[0, 0] - 3 * (H**2 + k / a**2)) == 0
    spatial = -(2 * a.diff(t, 2) / a + H**2 + k / a**2)
    for i in (1, 2, 3):
        assert sp.simplify(G[i, i] - spatial * g[i, i]) == 0
    assert all(G[i, j] == 0 for i in range(4) for j in range(4) if i != j)


def test_flrw_acceleration_equation():
    """Perfect fluid at rest: T_tt = rho, T_ij = p g_ij. Solve for rho and p, then
    check that the acceleration equation follows from the other two."""
    st = metrics.flrw()
    G, g = st.einstein, st.metric
    rho_solved = (G[0, 0] + Lambda * g[0, 0]) / (8 * sp.pi * G_N)
    p_solved = (G[1, 1] + Lambda * g[1, 1]) / (8 * sp.pi * G_N * g[1, 1])
    acceleration = -4 * sp.pi * G_N / 3 * (rho_solved + 3 * p_solved) + Lambda / 3
    assert sp.simplify(acceleration - a.diff(t, 2) / a) == 0

    friedmann = 8 * sp.pi * G_N / 3 * rho_solved - k / a**2 + Lambda / 3
    assert sp.simplify(friedmann - (a.diff(t) / a) ** 2) == 0


def test_newtonian_limit():
    """To first order, R_00 = laplacian(Phi) and G_00 = 2 laplacian(Phi), so
    G_00 = 8 pi G rho is Poisson's equation laplacian(Phi) = 4 pi G rho (notes/03)."""
    st = metrics.weak_field()
    laplacian = sum(Phi.diff(v, 2) for v in (x, y, z))
    assert sp.simplify(first_order(st.ricci[0, 0]) - laplacian) == 0
    assert sp.simplify(first_order(st.einstein[0, 0]) - 2 * laplacian) == 0


@pytest.mark.parametrize("polarization", ["plus", "cross"])
def test_linearized_plane_wave_is_vacuum(polarization):
    st = metrics.plane_wave(polarization)
    assert all(first_order(st.ricci[i, j]) == 0 for i in range(4) for j in range(4))
