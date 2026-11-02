import sympy as sp
import pytest
from research.calc.gh_50_principal import blocks,schwarzschild_ks_axis,exact_axis_report


def test_50_fields_at_horizon():
    d=exact_axis_report()
    assert d["field_count"]==50
    assert d["principal_symmetric_all_directions"]
    assert d["positive_symmetrizer_min_eigen"]==sp.Rational(1,2)
    assert d["shift_x"]==sp.Rational(1,2)
    assert d["lapse"]==1/sp.sqrt(2)
    assert d["full_einstein_source_and_constraints_verified"] is False


def test_ordinary_identity_fails_as_symmetrizer():
    N,beta,g=schwarzschild_ks_axis()
    S,As=blocks(N,beta,g)
    assert sp.simplify(As[0]-As[0].T)!=sp.zeros(50)
    assert sp.simplify(S*As[0]-As[0].T*S)==sp.zeros(50)
    assert As[0][:5,:5].eigenvals()=={
        0:2,sp.Rational(-1):1,sp.Rational(-1,2):2}


def test_different_radius_and_negative_controls():
    d=exact_axis_report(r=sp.Integer(3))
    assert d["positive_symmetrizer_min_eigen"]==sp.Rational(3,5)
    with pytest.raises(ValueError):
        blocks(sp.Integer(-1),[0,0,0],sp.eye(3))
    with pytest.raises(ValueError):
        blocks(sp.Integer(1),[0,0,0],sp.Matrix([[1,1,0],[0,1,0],[0,0,1]]))
