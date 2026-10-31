"""Petrov type of the Brill-Lindquist head-on spacetime near the initial slice (exact arithmetic).

Route B/G screening (research/route_map.md): every known family of exact radiative vacuum
solutions built from a shear-free geodesic null congruence (Kerr-Schild with one null vector,
Robinson-Trautman, Kundt, twisting type II) is algebraically special (Goldberg-Sachs). This module
checks that the vacuum development of equal-mass time-symmetric Brill-Lindquist data is NOT
algebraically special: it is of Petrov type I on an open set.

At a moment of time symmetry (K = 0) in vacuum, the magnetic part of the Weyl tensor vanishes and
the electric part is the Ricci tensor of the slice, E_ij = R_ij(gamma) (Gauss equation; R(gamma) = 0
by the Hamiltonian constraint). Then Q = E - iB = E is real, symmetric and traceless, and the Weyl
tensor is algebraically special exactly where E has a repeated eigenvalue, i.e. where the
discriminant
    Delta = 4 (tr E^2 / 2)^3 - 27 det(E)^2  (>= 0 for real symmetric E)
vanishes. Delta > 0 at a point of the slice is an open condition on the spacetime Weyl tensor, so
the development is of type I on a spacetime neighbourhood of that point.

For gamma = psi^4 delta (flat Laplacian Delta psi = 0) the Ricci tensor is (checked against the
general curvature code in tests/test_petrov_bl.py)
    R_ij = -2 psi^-1 d_i d_j psi + 6 psi^-2 d_i psi d_j psi - 2 psi^-2 delta_ij |d psi|^2 .

Run: python -m research.calc.petrov_bl
"""

import sympy as sp

x, y, z = sp.symbols("x y z", real=True)


def bl_psi(d, m=sp.Rational(1, 2)):
    """Brill-Lindquist conformal factor, bare masses m at (0, 0, +-d/2)."""
    rp = sp.sqrt(x**2 + y**2 + (z - d / 2) ** 2)
    rm = sp.sqrt(x**2 + y**2 + (z + d / 2) ** 2)
    return 1 + m / (2 * rp) + m / (2 * rm)


def ricci_conformally_flat(psi):
    """R_ij of psi^4 delta for harmonic psi (formula in the module docstring)."""
    X = (x, y, z)
    g = [sp.diff(psi, c) for c in X]
    grad2 = sum(a * a for a in g)
    return sp.Matrix(3, 3, lambda i, j: -2 * sp.diff(psi, X[i], X[j]) / psi
                     + 6 * g[i] * g[j] / psi**2 - (2 * grad2 / psi**2 if i == j else 0))


def discriminant(E):
    """4 (tr E^2/2)^3 - 27 det(E)^2 for a traceless symmetric 3x3 matrix (product of squared
    eigenvalue differences)."""
    p = (E * E).trace() / 2
    return 4 * p**3 - 27 * E.det()**2


def electric_weyl_at(d, point, m=sp.Rational(1, 2)):
    """E_ij = R_ij(gamma) of BL data at an exact point (exact arithmetic)."""
    psi = bl_psi(d, m)
    R = ricci_conformally_flat(psi)
    sub = dict(zip((x, y, z), point))
    return R.subs(sub).applyfunc(sp.nsimplify).applyfunc(sp.radsimp)


def equatorial_eigenvalues(d, rho, m=sp.Rational(1, 2)):
    """Eigenvalues of E at (rho, 0, 0) on the plane z = 0: there E is diagonal in (x, y, z) =
    (radial, azimuthal, axial)."""
    E = electric_weyl_at(d, (rho, 0, 0), m)
    assert all(sp.simplify(E[i, j]) == 0 for i in range(3) for j in range(3) if i != j)
    return [sp.radsimp(E[i, i]) for i in range(3)], E


def main():
    d = sp.Rational(39, 50)            # d = 0.78 (Theorem F); rho = 4/5 gives r = 89/100
    lam, E = equatorial_eigenvalues(d, sp.Rational(4, 5))
    D = sp.radsimp(discriminant(E))
    print("d = 39/50, point (4/5, 0, 0): eigenvalues of E =", [sp.nsimplify(v) for v in lam])
    print("  as decimals:", [float(v) for v in lam])
    print("  discriminant =", D, "=", float(D), "> 0:", bool(D > 0))
    # an off-axis, off-plane point with rational r_+ and r_- (d = 1): (x, z) = (4/5, 1/10) gives
    # r_+ = sqrt(0.64 + 0.16) irrational, so use symbolic radicals and exact sign evaluation
    E2 = electric_weyl_at(sp.Integer(1), (sp.Rational(3, 10), 0, sp.Rational(9, 10)))
    D2 = discriminant(E2)
    print("d = 1, point (3/10, 0, 9/10): discriminant =", sp.N(D2, 30))
    E3 = electric_weyl_at(sp.Integer(1), (0, 0, sp.Rational(3, 2)))
    print("on the axis, d = 1, (0, 0, 3/2): discriminant =", sp.simplify(discriminant(E3)),
          "(type D on the axis, by symmetry)")


if __name__ == "__main__":
    main()
