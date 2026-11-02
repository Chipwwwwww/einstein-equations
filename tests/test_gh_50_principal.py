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



def test_full_50_positive_gamma2_and_cross_term_correction():
    N,beta,ginv=schwarzschild_ks_axis()
    S,As=blocks(N,beta,ginv,gamma2=sp.Rational(1,4),Lambda2=sp.Integer(1))
    assert S[0,1]==-sp.Rational(1,4)
    assert S[1,0]==-sp.Rational(1,4)
    assert all(sp.simplify(S*A-A.T*S)==sp.zeros(50) for A in As)
    assert S[:5,:5].eigenvals()=={
        sp.Rational(3,4):1,sp.Rational(5,4):1,
        sp.Rational(1,2):1,sp.Integer(1):2}
    with pytest.raises(ValueError,match="Lambda2"):
        blocks(N,beta,ginv,gamma2=sp.Rational(3,2),Lambda2=sp.Integer(1))
