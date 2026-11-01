"""Petrov type of the Brill-Lindquist development at the initial slice
(research/calc/petrov_bl.py, route map B/G): type I on an open dense set, for every d > 0."""

import sympy as sp

from einstein import Spacetime
from research.calc import petrov_bl as pb
from research.calc.petrov_bl import x, y, z


def test_conformally_flat_ricci_formula():
    """R_ij(psi^4 delta) from the general curvature code equals the closed formula used in
    petrov_bl, for a harmonic psi (a polynomial one keeps the general code fast; the formula is
    pointwise in psi and its first two derivatives)."""
    psi = 3 + x**2 - y**2 + x * z + sp.Rational(1, 2) * y * z
    assert sum(sp.diff(psi, c, 2) for c in (x, y, z)) == 0
    g = sp.diag(psi**4, psi**4, psi**4)
    R = Spacetime(g, (x, y, z), simplify=sp.cancel).ricci
    F = pb.ricci_conformally_flat(psi)
    for i in range(3):
        for j in range(3):
            assert sp.cancel(sp.sympify(R[i, j]) - F[i, j]) == 0


def test_type_I_at_an_exact_point_d078():
    """d = 39/50, point (4/5, 0, 0): three distinct rational eigenvalues, discriminant > 0
    exactly (rational arithmetic)."""
    lam, E = pb.equatorial_eigenvalues(sp.Rational(39, 50), sp.Rational(4, 5))
    assert all(v.is_Rational for v in lam) and len(set(lam)) == 3 and sum(lam) == 0
    D = pb.discriminant(E)
    assert D.is_Rational and D > 0


def test_equatorial_closed_forms_and_single_centre_limit():
    """The closed-form eigenvalue differences on z = 0 (r = distance to either puncture):
    lam_phi - lam_z = 3 d^2 / (2 r^4 (2r + 1)) vanishes only for d = 0 (Schwarzschild, type D)."""
    for d, rho in ((sp.Rational(39, 50), sp.Rational(4, 5)), (sp.Integer(3), sp.Integer(2))):
        lam, _ = pb.equatorial_eigenvalues(d, rho)
        r = sp.sqrt(rho**2 + d**2 / 4)
        assert sp.simplify(lam[1] - lam[2] - 3 * d**2 / (2 * r**4 * (2 * r + 1))) == 0
        assert sp.simplify(lam[0] - lam[1] + 3 * (2 * r - d) * (2 * r + d) / (r**3 * (2 * r + 1)**2)) == 0
    lam0, _ = pb.equatorial_eigenvalues(sp.Integer(0), sp.Integer(2))
    assert lam0[1] == lam0[2] and sp.simplify(lam0[0] + 2 * lam0[1]) == 0
