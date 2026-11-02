import sympy as sp
from research.calc.bl_gh_initial_data import symbolic_checks


def test_exact_BL_GH_initial_constraints_and_reduction():
    c=symbolic_checks()
    assert c["Hamiltonian_constraint_exact_off_punctures"]
    assert c["momentum_constraint_exact"]
    assert c["initial_harmonic_gauge_H_zero"]
    assert c["first_order_reduction_constraints_exact"]
    assert c["normal_derivative_GH_constraint_not_independently_verified"]
    assert c["Pi"]==sp.zeros(4)
    assert c["full_finite_slab_3plus1_certified"] is False



def test_BL_puncture_cutoff_GH_geometric_hyperbolicity():
    from research.calc.bl_gh_initial_data import initial_hyperbolicity_bounds
    from fractions import Fraction as Q
    b=initial_hyperbolicity_bounds()
    assert b["psi_max"]==Q(3)
    assert b["lapse_initial_lower"]==Q(1,9)
    assert b["GH_initial_symmetrizer_min"]==Q(1,81)
    assert b["metric_Neumann_multiplier"]==Q(243,500)
    assert b["GH_symmetrizer_tube_min"]==Q(500,40503)
    assert b["positive_time_full_GH_constraint_and_slab"] is False
    import pytest
    with pytest.raises(ValueError,match="Lorentzian"):
        initial_hyperbolicity_bounds(tube_H3_radius=Q(1,100))
