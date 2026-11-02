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
