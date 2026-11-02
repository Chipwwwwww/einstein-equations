from fractions import Fraction as Q
import pytest
import sympy as sp
from research.calc.full_hs_bootstrap import verify


def test_closed_coupled_nonlinear_h3_bootstrap():
    c=verify()
    assert c["bootstrap_closed"]
    assert c["Sobolev"]=="H^3(T^3;R^2)"
    assert c["residual_H3_upper"]==Q(393,1250000)
    assert c["exponent_K"]==Q(75261,250000)
    assert c["closed_H3_error"]==Q(393,124623695)
    assert c["independent_actual_H3_error_upper"]<=c["closed_H3_error"]
    assert c["full_50_field_GH_vacuum_certified"] is False


def test_symmetric_noncommuting_principal_matrices():
    u0,u1=sp.symbols("u0 u1",real=True)
    A1=sp.Matrix([[u0,sp.Rational(1,20)],[sp.Rational(1,20),u1]])
    A2=sp.diag(u1,-u0)
    assert A1==A1.T and A2==A2.T
    assert (A1*A2-A2*A1)[0,1]!=0
    assert A1.det().subs({u0:0,u1:0})==-sp.Rational(1,400)


def test_fail_closed_on_bad_claims():
    with pytest.raises(ValueError,match="residual"):
        verify(residual_override=Q(0))
    with pytest.raises(ValueError,match="initial"):
        verify(require_exact_initial=False)
    with pytest.raises(ValueError,match="A-matrix"):
        verify(R=Q(1,2))
    with pytest.raises(ValueError,match="bootstrap"):
        verify(R=Q(1,10000000))
    with pytest.raises(ValueError):
        verify(a=-Q(1,100))


def test_iterated_slabs_are_sharper():
    c1=verify(nslabs=1)
    c2=verify(nslabs=2)
    c8=verify(nslabs=8)
    assert c8["closed_H3_error"]<c2["closed_H3_error"]<c1["closed_H3_error"]
    assert all(y<Q(1,1000) for y in c8["slab_errors"])
