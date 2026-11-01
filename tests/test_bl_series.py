"""Route E (research/calc/bl_series.py): the geodesic-slicing development of Brill-Lindquist data as
a time series. Checked against the exact free fall of Schwarzschild, against the Hamiltonian
constraint order by order, and against the independent Cartesian formula for Ric(gamma_0)."""

import sympy as sp

from research.calc import bl_series as bs
from research.calc import petrov_bl as pb


def test_schwarzschild_free_fall_from_rest():
    """One centre (M = 1): the azimuthal metric coefficient at isotropic radius 1 evolves as
    (r(t)/R0)^2 with r'' = -M/r^2, r(0) = R0 = 9/4 (radial free fall from rest): all coefficients
    through t^8 agree to 1e-9 relative."""
    vals, _, ham = bs.series(1.0, 0.0, [(1.0, 0.0)], 8)
    exact = bs.schwarzschild_ratio(2.25, 1.0, 4)
    got = [v / vals["pp"][0] for v in vals["pp"]]
    for a, b in zip(got, exact):
        assert abs(a - b) <= 1e-9 * abs(b) + 1e-15
    assert max(abs(h) for h in ham) < 1e-12


def test_first_coefficient_is_minus_ricci():
    """BL data (d = 39/50): g_1 = -Ric(gamma_0) at a point, from the jet code in (rho, z, phi),
    equals minus the closed Cartesian formula of research/calc/petrov_bl.py (independent code):
    on the plane y = 0 the (rho, z) block is the (x, z) block and g_pp = rho^2 times the yy entry."""
    d, rho0, z0 = 0.78, 0.6, 0.15
    vals, _, ham = bs.series(rho0, z0, [(0.5, d / 2), (0.5, -d / 2)], 4)
    E = pb.electric_weyl_at(sp.Rational(39, 50), (sp.Rational(3, 5), 0, sp.Rational(3, 20)))
    assert abs(vals["rr"][1] + float(E[0, 0])) < 1e-10
    assert abs(vals["rz"][1] + float(E[0, 2])) < 1e-10
    assert abs(vals["zz"][1] + float(E[2, 2])) < 1e-10
    assert abs(vals["pp"][1] + rho0**2 * float(E[1, 1])) < 1e-10
    assert max(abs(h) for h in ham) < 1e-12


def test_magnetic_weyl_starts_at_t_cubed():
    """For conformally flat time-symmetric data, B = -t Cotton(gamma_0) + O(t^3) = O(t^3). BL
    (d = 0.78): the t coefficient vanishes, the t^3 one does not. Schwarzschild: B = 0 through t^5
    (spherical symmetry), a control for the code."""
    b = bs.magnetic_weyl(0.25, 0.1, [(0.5, 0.39), (0.5, -0.39)], 8)
    assert abs(b["rp"][0]) < 1e-12 and abs(b["zp"][0]) < 1e-12
    assert abs(b["zp"][1]) > 0.1
    c = bs.magnetic_weyl(1.0, 0.3, [(1.0, 0.0)], 8)
    assert max(abs(v) for v in c["rp"] + c["zp"]) < 1e-12


def test_momentum_constraint_holds_order_by_order():
    """The recursion uses only the evolution equations; the momentum constraint
    D_j K^j_i - D_i K must still vanish at every order the jets justify (BL, d = 0.78, and
    Schwarzschild)."""
    for centres, pt in (([(0.5, 0.39), (0.5, -0.39)], (0.25, 0.1)), ([(1.0, 0.0)], (1.0, 0.3))):
        m = bs.momentum_constraint(*pt, centres, 8)
        assert max(abs(v) for v in m["r"] + m["z"]) < 1e-10
