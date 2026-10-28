"""Identities that hold for every metric, checked on metrics with free functions."""

import pytest
import sympy as sp

from einstein import Spacetime, metrics


@pytest.mark.parametrize("make", [metrics.static_spherical, metrics.flrw])
def test_contracted_bianchi_identity(make):
    """nabla_a G^a_b = 0 identically, which forces nabla_a T^a_b = 0 (notes/01)."""
    st = make()
    assert st.divergence(st.einstein_mixed) == [0] * st.dim


def test_christoffel_symmetric_in_lower_indices():
    st = metrics.static_spherical()
    n = st.dim
    assert all(
        st.christoffel[a, b, c] == st.christoffel[a, c, b]
        for a in range(n) for b in range(n) for c in range(n)
    )


def test_two_sphere_has_positive_curvature():
    """Sign check for notes/00: with our conventions a sphere of radius rho has R = 2 / rho^2."""
    rho = sp.Symbol("rho", positive=True)
    theta, phi = sp.symbols("theta phi")
    sphere = Spacetime(sp.diag(rho**2, rho**2 * sp.sin(theta) ** 2), (theta, phi))
    assert sp.simplify(sphere.ricci_scalar - 2 / rho**2) == 0


def test_rejects_asymmetric_metric():
    t, x = sp.symbols("t x")
    with pytest.raises(ValueError):
        Spacetime(sp.Matrix([[-1, x], [0, 1]]), (t, x))
