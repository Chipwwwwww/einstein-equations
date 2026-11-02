"""BL whole-domain analytic proof witness and generic full-GH source controls."""
from fractions import Fraction as F
import pytest
import sympy as sp
from research.calc.bl_local_cutoff import (
    theta_polynomials,cutoff_derivatives,derivative_envelopes,envelope_summary)
from research.calc.bl_torus_majorant import certificate,bound_graph
from research.calc.bl_gh_taylor_evaluator import (
    p_and_derivatives,physical_initial_and_first_jet,evaluate_metric)


def test_smooth_flat_cutoff_analytic_envelopes():
    b=envelope_summary()
    assert b['C_derivative_exponent']==36
    assert b['p_derivative_exponent']==95
    assert b['g_derivative_exponent']==789
    assert b['reference_H7_component_exponent']==791
    assert b['pmax_upper']=='31/6'
    p=theta_polynomials(3)
    assert p[0]=={0:1}
    assert p[1]=={2:1}
    assert p[2]=={4:1,3:-2}
    assert cutoff_derivatives(3)[1][0]==1


def test_explicit_binary_domain_and_closed_50_gh_energy():
    c=certificate()
    assert c['d']=='3831/5000'
    assert c['reference_H7_exp']==795
    assert c['time_T']=='10^(-29045)'
    assert c['metric_H3_error_upper']=='10^(-37098)'
    assert c['nonlinear_H3_radius']=='10^(-8050)'
    assert c['GH_residual_H3_upper']=='10^(-8060)'
    assert c['GH_speed_upper'].startswith('4 ')
    assert c['fallback_graph_speed_upper']=='10^(4019)'
    assert c['finite_speed_displacement_upper']=='4*10^(-29045)'
    assert all(c['proof_flags'].values())
    assert not c['local_exact_vacuum_outside_cutoff']
    assert not c['complete_binary_merger_or_extension']


def test_nonflat_BL_metric_and_true_nonzero_metric_time_jet():
    p,dp,ddp=p_and_derivatives(sp.Integer(0),sp.Integer(0),sp.Integer(0))
    assert sp.factor(p-sp.Rational(1)-sp.Rational(1,2)/sp.Rational(3831,10000))==0
    assert all(d==0 for d in dp) # reflection symmetry at center
    assert sp.factor(sum(ddp[i,i] for i in range(3)))==0
    g,pi1,g2=physical_initial_and_first_jet(sp.Integer(0),sp.Integer(0),sp.Integer(0))
    assert g[0,0]<0 and g[1,1]>25
    assert all(sp.simplify(g2[i,i])==0 for i in range(4)) # lapse choice cancels at symmetry center
    g_off,pi_off,g2_off=physical_initial_and_first_jet(sp.Rational(1,10),sp.Integer(0),sp.Integer(0))
    assert any(sp.simplify(g2_off[i,i])!=0 for i in range(4))
    T=sp.Rational(1,10**29045)
    approx=evaluate_metric(T,0,0,0)
    assert approx==g+T*T*g2/2
    assert approx==g # at midpoint this initial lapse cancels the second metric jet
    taylor_off=evaluate_metric(T,sp.Rational(1,10),0,0)
    assert taylor_off !=g_off
    with pytest.raises(ValueError,match='time'):
        evaluate_metric(sp.Rational(1,10**29044),0,0,0)
    with pytest.raises(ValueError,match='point'):
        evaluate_metric(T,0,0,sp.Rational(1,2))


def test_torus_proof_rejects_nongeometric_inputs():
    with pytest.raises(ValueError,match='BL punctures'):
        p_and_derivatives(0,0,sp.Rational(3831,10000))
    with pytest.raises(ValueError,match='unexpected'):
        from research.calc.bl_torus_majorant import E
        E(1).c(5000)


def test_exact_auxiliary_extension_equality_and_periodic_seam():
    from research.calc.bl_local_cutoff import bl_auxiliary_p
    x=sp.Rational(1,10)
    p=bl_auxiliary_p(x,sp.Integer(0),sp.Integer(0))
    assert sp.simplify(p-p_and_derivatives(x,0,0)[0])==0
    assert bl_auxiliary_p(2,0,0)==1
    assert bl_auxiliary_p(0,0,sp.Rational(3831,10000))>1


def test_BL_initial_normal_GH_gauge_constraint_direct_zero():
    p=sp.symbols('p',positive=True)
    dx,dy,dz=sp.symbols('dx dy dz',real=True)
    v=sp.Matrix([dx,dy,dz]);q=(v.T*v)[0]
    pi=sp.zeros(4)
    pi[0,0]=4*q/p**12
    for i in range(3):
        for j in range(3):
            pi[i+1,j+1]=4*(q*int(i==j)-4*v[i]*v[j])/p**4
    ginv=sp.diag(-p**4,p**-4,p**-4,p**-4)
    N=p**-2
    # d_t C_a = n^b d_t Pi_ba - (1/2)n_a g^{bc} d_t Pi_bc.
    # At t=0, C=0 as a spatial identity, g_t=Phi_t=0.
    dtC0=pi[0,0]/N + N*sp.trace(ginv*pi)/2
    dtCi=[pi[0,i]/N for i in range(1,4)]
    assert sp.factor(dtC0)==0
    assert all(sp.factor(v)==0 for v in dtCi)

def test_short_time_and_unverified_faster_slab_rejected():
    with pytest.raises(ValueError,match='rejected'):
        certificate(time_exponent_override=10000)
    with pytest.raises(ValueError,match='positive integer'):
        certificate(time_exponent_override=-1)


def test_all_fifty_field_evaluable_time_taylor():
    from research.calc.bl_gh_taylor_evaluator import evaluate_50_field_state
    T=sp.Rational(1,10**29045)
    u=evaluate_50_field_state(T,sp.Rational(1,10),0,0)
    assert len(u)==50
    assert any(z!=0 for z in u[10:20])
    p,dp,_=p_and_derivatives(sp.Rational(1,10),0,0)
    assert sp.simplify(u[20]-4*p**-5*dp[0])==0
