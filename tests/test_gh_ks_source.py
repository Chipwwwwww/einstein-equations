import sympy as sp
from research.calc.gh_ks_source import verify_source


def test_exact_ks_gh_source():
    d=verify_source()
    assert d["exact_GH_constraint"]
    assert d["shell_r_3M_over_2_to_5M_over_2_symmetrizer_eigen_lower"]==sp.Rational(3,7)
    assert d["complete_H3_residual_bootstrap_or_vacuum_certificate"] is False
