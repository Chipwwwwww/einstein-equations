"""Independent positive/negative tests: full second GH jet, sharp C9 BL bounds.

These tests do not substitute for the mathematical proof in note31.
"""
from fractions import Fraction as Q
from sympy import (Rational as R, Integer, Matrix, diag, symbols, sqrt, diff,
                   simplify, zeros, factorial, binomial)
import sympy as sp
import pytest
from research.calc.bl_sharp_cutoff import (
    cutoff_coefficients, bernstein_derivative_bounds, derivative_envelopes,
    exact_ceil_log10, auxiliary_p)
from research.calc.bl_sharp_certificate import certificate
from research.calc.bl_full_second_jet import (
    physical_jets, initial_jets_from_p, full_50_field_state, metric,
    PAIRS, p_and_derivatives)

def test_C9_finite_polynomial_cutoff_and_exact_Bernstein():
    z=symbols("z");P=sum(sp.Rational(c.numerator,c.denominator)*z**k
                         for k,c in cutoff_coefficients(9).items())
    assert P.subs(z,0)==0 and P.subs(z,1)==1
    for i in range(1,10):
        assert diff(P,z,i).subs(z,0)==diff(P,z,i).subs(z,1)==0
    assert diff(P,z,10).subs(z,0)!=0  # C9, not C-infinity
    D=bernstein_derivative_bounds(cutoff_coefficients(9),8)
    assert len(D)==9 and all(d>0 for d in D)
    e=derivative_envelopes()
    assert e["H3_H7_log10_bounds"]=={3:22,4:27,5:32,6:37,7:42}
    assert e["k"]==8 and e["m"]==9

def test_exact_decimal_exponent_for_both_small_and_large_numbers():
    assert exact_ceil_log10(Q(1,1000))==-3
    assert exact_ceil_log10(Q(1,10000)*Q(3,2))==-3
    assert exact_ceil_log10(Q(74184649,10**9))==-1
    assert exact_ceil_log10(Q(100001))==6
    for k in range(-20,21):
        v=Q(10)**k
        assert exact_ceil_log10(v)==k
        assert exact_ceil_log10(v*Q(101,100))==k+1

def test_BL_auxiliary_exact_match_and_not_vacuum_outside():
    x=R(1,10);a=R(3831,10000)
    p=auxiliary_p(x,0,0)
    assert simplify(p-(1+1/(2*sqrt(x*x+a*a))))==0
    assert auxiliary_p(R(2),0,0)==1
    assert auxiliary_p(0,0,a)>1  # the other hole contributes; this puncture is smooth
    assert derivative_envelopes()["globally_vacuum"] is False

def test_complete_full_50_field_second_time_jet_has_missing_components():
    u0,u1,u2=physical_jets(R(1,10),0,0)
    assert len(u0)==len(u1)==len(u2)==50
    assert u1[11]==0
    assert u2[11]!=0 and u2[24]!=0 and u2[35]!=0 and u2[46]!=0
    assert u2[0]!=0 and u2[4]!=0

def test_all_Pi2_tensor_entries_offshell_nonharmonic_p():
    p=R(13,7);q=[R(2,5),-R(3,8),R(4,11)];L=-R(7,13)
    Hess=diag(L,0,0)
    u0,u1,u2=initial_jets_from_p(p,q,Hess,L)
    D=sum(a*a for a in q)
    for i in range(3):
        idx=10+PAIRS.index((0,i+1))
        assert simplify(u2[idx]-16*q[i]*(p*L+5*D)/p**13)==0
    for a,b in PAIRS:
        if not(a==0 and b!=0):
            assert u2[10+PAIRS.index((a,b))]==0

def test_Phi2_equals_gradient_metric2_at_nonsymmetric_point():
    x,y,z=symbols("x y z",real=True)
    a=R(3831,10000)
    p=1+1/(4*sqrt(x*x+y*y+(z-a)**2))+1/(4*sqrt(x*x+y*y+(z+a)**2))
    px,py,pz=[diff(p,v) for v in (x,y,z)]
    D=px*px+py*py+pz*pz
    point={x:R(1,10),y:0,z:0}
    u0,u1,u2=physical_jets(R(1,10),0,0)
    assert simplify(u2[20]-diff(-4*D/p**14,x).subs(point))==0
    assert simplify(u2[24]-diff((16*px*px-4*D)/p**6,x).subs(point))==0

def test_sharp_graph_certificate_and_fail_closed():
    c=certificate()
    assert all(c["flags"].values())
    assert c["reference_Hs_log10"]=={3:22,4:27,5:32,6:37,7:42}
    assert c["full_GH_source_graph_H3"]["F_l"]==538
    assert c["full_GH_source_graph_H3"]["F_h"]==714
    assert c["growth_K_log10"]==551
    assert c["residual_t2_coefficient_log10"]==2812
    assert c["time_T_over_M"]=="10^(-2832)"
    assert c["H3_error_upper"]=="10^(-5674)"
    with pytest.raises(ValueError,match="rejected"):
        certificate(time_exponent_override=1000)
    with pytest.raises(ValueError,match="positive integer"):
        certificate(time_exponent_override=-1)

def test_resolved_dynamics_and_exact_evaluator_domain_guards():
    a=Q(1,10); rmin=Q(39,100); rmax=Q(2,5)
    # r=sqrt(0.1^2+0.3831^2), so rmin<r<rmax exactly
    assert rmin*rmin < a*a+Q(3831,10000)**2 < rmax*rmax
    gxx2_lower=12*(a/(2*rmax**3))**2/(1+1/(2*rmin))**6
    assert gxx2_lower>Q(1,20)
    assert 10**10>480  # delta gxx>=T²/40>2*(6 E)
    t=R(1,10**2832)
    assert len(full_50_field_state(t,R(1,10),0,0))==50
    assert metric(t,R(1,10),0,0,precision=12).shape==(4,4)
    with pytest.raises(ValueError,match="time"):
        metric(t*2,R(1,10),0,0)
    with pytest.raises(ValueError,match="outside"):
        metric(t,R(1),0,0)
    with pytest.raises(ValueError,match="puncture"):
        metric(t,0,0,R(3831,10000))
