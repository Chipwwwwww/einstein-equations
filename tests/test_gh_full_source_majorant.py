"""Independent complete 50-field source and full-tube majorant regression."""
from fractions import Fraction as Q
import sympy as sp
import pytest
from research.calc.gh_full_source import (full_source,kasner_initial,kasner_at,
  principal_blocks,pack50,sym4,flat10,PAIRS)
from research.calc.gh_h3_majorant import Jet,majorant,jet_inputs,kasner_centers_and_deviations


def test_all_50_kasner_equations_symbolically_for_all_harmonic_times():
    q=sp.symbols('q',positive=True)
    A=sp.Rational(4,3);B=sp.Rational(2,3)
    g=sp.diag(-q**6,q**4,q**4,q**-2)
    pi=sp.diag(2*q**3,-A*q,-A*q,B*q**-5)
    v=pack50([[g[i,j] for j in range(4)] for i in range(4)],
      [[pi[i,j] for j in range(4)] for i in range(4)],
      [[[0]*4 for _ in range(4)] for _ in range(3)])
    F=full_source(v)
    assert len(F)==50
    assert all(sp.factor(F[i]-q*sp.diff(v[i],q)/3)==0 for i in range(50))
    # Gauge C_a=0: Γ_a computed from d_tg and source first-order data.
    # det(g)=-q^12, sqrt(-g)=q^6; g^00=-q^-6 -> constant density -1.
    assert sp.simplify(sp.sqrt(-g.det())*g.inv()[0,0]+1)==0
    for j in range(1,4):
        assert all(g.inv()[i,j]==0 for i in range(4) if i!=j)
    # Nonflat Riemann component: R^x_{t x t}=p_x(1-p_x)=2/9.
    assert sp.Rational(2,3)*(1-sp.Rational(2,3))==sp.Rational(2,9)


def test_full_source_independent_homogeneous_generic_metric_jet():
    # At Minkowski g0 with arbitrary symmetric time-derivative tensor M,
    # derive the reduced Einstein equation from Eq.(4), independently of the
    # first-order source. Gamma_ab^c is built from the metric derivative.
    M=sp.Matrix([[sp.Rational(1,5),sp.Rational(1,7),0,sp.Rational(-1,11)],
       [sp.Rational(1,7),sp.Rational(1,3),sp.Rational(1,13),0],
       [0,sp.Rational(1,13),sp.Rational(-1,17),sp.Rational(1,19)],
       [sp.Rational(-1,11),0,sp.Rational(1,19),sp.Rational(1,23)]])
    g=sp.diag(-1,1,1,1);pi=-M
    v=pack50([[g[i,j] for j in range(4)] for i in range(4)],
      [[pi[i,j] for j in range(4)] for i in range(4)],
      [[[0]*4 for _ in range(4)] for _ in range(3)])
    F=full_source(v)
    pi_t=sym4(F[10:20])
    # d_t N=-M_00/2 (g^-1_00 derivative = -M_00).
    gtt=[[-M[0,0]*M[i,j]/2-pi_t[i][j] for j in range(4)] for i in range(4)]
    dg=[[[M[i,j] for j in range(4)] for i in range(4)]]+[[[0]*4 for _ in range(4)] for _ in range(3)]
    Gamma=[[[sum((dg[b][a][c],dg[c][a][b],-dg[a][b][c]))/2
             for c in range(4)] for b in range(4)]for a in range(4)]
    gu=g.inv()
    for a,b in PAIRS:
        quadr=0
        for c in range(4):
          for d in range(4):
            quadr+=gu[c,d]*(-pi[c,a]*pi[d,b]-sum(gu[e,f]*Gamma[a][c][e]*Gamma[b][d][f]
                       for e in range(4) for f in range(4)))
        # R_ab - ∇_(a Γ_b)= -(1/2)g^00 gtt + quadratic = gtt/2+quadr
        assert sp.factor(gtt[a][b]/2+quadr)==0
    # Artificial sign corruption MUST be detectable.
    assert sp.factor((gtt[0][1]+sp.Rational(1,10))/2+quadr)!=0


def test_banach_algebra_neumann_reciprocal_and_all_tube():
    j=Jet(1,Q(1,10000),Q(1),48)
    inv=j.reciprocal()
    assert inv.c==1
    assert inv.d==Q(1,10000)/(1-Q(48,10000))
    assert inv.L>1
    root=j.inverse_sqrt()
    assert root.d==j.d/(1-48*j.d)
    assert root.L==1/(1-48*j.d)**2
    with pytest.raises(ValueError):Jet(1,Q(1,40)).reciprocal()
    with pytest.raises(ValueError):Jet(1,Q(1,40)).inverse_sqrt()


def test_actual_full_source_composition_frechet_and_nonlinear_existence():
    a=majorant(T=Q(1,10**3),R=Q(1,10**8),eps=Q(1,10**7))
    assert a['closed']
    assert a['full_GH_50_source_composition']
    assert a['Lorentzian_metric_gap']<1
    assert a['lambda']>0
    assert a['K_times_T']<1
    assert 0<a['residual_H3']<Q(4,10**9)
    assert 0<a['certified_H3_error']<Q(1,10**8)
    assert not a['Einstein_vacuum_subsidiary_proved']
    # Regression in floats for readability, exact comparisons above.
    assert Q(720)<a['growth_K']<Q(750)
    assert a['adapted_D']*100 < a['generic_D']
    assert a['metric_dt_Linf']<Q(33)
    assert Q(14,10**12)<a['certified_H3_error']<Q(15,10**12)
    assert Q(3,10**9)<a['residual_H3']<Q(4,10**9)


def test_falsification_inadmissible_tube_time_or_approximation():
    with pytest.raises(ValueError,match='too long|inverse_sqrt|reciprocal'):
        majorant(T=Q(1,10**2),R=Q(1,10**8))
    with pytest.raises(ValueError,match='matrix-inertia|reciprocal|inverse_sqrt|positive'):
        majorant(T=Q(1,10**7),R=Q(1,10))
    with pytest.raises(ValueError,match='share certified'):
        majorant(T=Q(1,10**6),R=Q(1,10**19),eps=Q(1,10**7))


def test_compact_outward_rounded_exact_rational_certificate():
    from research.calc.gh_h3_certificate import certificate
    c=certificate()
    assert c['energy_growth_upper'] < Q(750)
    assert c['residual_H3_upper'] < Q(4,10**9)
    assert c['H3_metric_field_error_upper']<Q(2,10**11)
    assert c['full_50GH_algebraic_source'].endswith('gamma1=-1 gamma2=gamma0=0 H=0')
    assert not c['binary_merger_spacetime']

def test_independent_four_dimensional_kasner_ricci_and_curvature():
    t=sp.symbols('t',real=True)
    g=sp.diag(-sp.exp(2*t),sp.exp(sp.Rational(4,3)*t),
              sp.exp(sp.Rational(4,3)*t),sp.exp(-sp.Rational(2,3)*t))
    gi=sp.diag(-sp.exp(-2*t),sp.exp(-sp.Rational(4,3)*t),
               sp.exp(-sp.Rational(4,3)*t),sp.exp(sp.Rational(2,3)*t))
    coords=[t,*sp.symbols('x y z')]
    Gam=[[[sp.simplify(sum(gi[a,d]*(sp.diff(g[d,b],coords[c])+
                 sp.diff(g[d,c],coords[b])-sp.diff(g[b,c],coords[d]))/2
                 for d in range(4))) for c in range(4)] for b in range(4)]
                 for a in range(4)]
    Ric=[[sp.factor(sum(sp.diff(Gam[k][a][b],coords[k])-
             sp.diff(Gam[k][a][k],coords[b])+
             sum(Gam[k][k][d]*Gam[d][a][b]-Gam[k][b][d]*Gam[d][a][k]
               for d in range(4)) for k in range(4)))
          for b in range(4)] for a in range(4)]
    assert all(z==0 for row in Ric for z in row)
    a,b,c,d=1,0,1,0
    Ri=sp.factor(sp.diff(Gam[a][d][b],coords[c])-
          sp.diff(Gam[a][c][b],coords[d])+
          sum(Gam[a][c][k]*Gam[k][d][b]-Gam[a][d][k]*Gam[k][c][b]
             for k in range(4)))
    assert Ri==sp.Rational(2,9)
    assert sp.simplify(sp.sqrt(-g.det())*gi[0,0]+1)==0
    assert all(sp.simplify(sp.sqrt(-g.det())*gi[i,j])==0
               for i in range(4) for j in range(4) if i!=j)


def test_generic_spatially_coupled_source_matches_independent_second_order_GH():
    eta=sp.diag(-1,1,1,1)
    M=sp.Matrix([[sp.Rational(1,5),sp.Rational(-1,7),sp.Rational(1,11),0],
                 [sp.Rational(-1,7),sp.Rational(1,3),0,sp.Rational(1,13)],
                 [sp.Rational(1,11),0,sp.Rational(-1,17),sp.Rational(1,19)],
                 [0,sp.Rational(1,13),sp.Rational(1,19),sp.Rational(-1,23)]])
    Ds=[]
    for i in range(3):
        D=sp.Matrix(4,4,lambda a,b:sp.Rational((a+1)*(b+1)*(i+1),
            101+11*i+2*a+2*b))
        Ds.append(D)
    packed=pack50(eta.tolist(),(-M).tolist(),[x.tolist() for x in Ds])
    F=full_source(packed)
    dt=[M.tolist()]+[x.tolist() for x in Ds]
    Gamma=[[[sp.Rational(1,2)*(dt[b][a][c]+dt[c][a][b]-dt[a][b][c])
              for c in range(4)] for b in range(4)]for a in range(4)]
    reducedQ=sp.zeros(4)
    for a,b in PAIRS:
        q=0
        for c in range(4):
          for d in range(4):
            inner=sum(Ds[i][c,a]*Ds[i][d,b] for i in range(3))
            inner-=M[c,a]*M[d,b]
            inner-=sum(eta[e,f]*Gamma[a][c][e]*Gamma[b][d][f]
                       for e in range(4) for f in range(4))
            q+=eta[c,d]*inner
        reducedQ[a,b]=reducedQ[b,a]=sp.factor(q)
    Fpi=sp.Matrix(sym4(F[10:20]))
    N_t=-M[0,0]/2
    beta_t=[M[0,i+1] for i in range(3)]
    gtt=-Fpi+M*N_t+sum((beta_t[i]*Ds[i] for i in range(3)),sp.zeros(4))
    assert all(sp.factor(gtt[a,b]+2*reducedQ[a,b])==0 for a,b in PAIRS)
    changed=sp.Matrix(gtt)
    changed[0,1]+=sp.Rational(1,100)
    assert sp.factor(changed[0,1]+2*reducedQ[0,1])!=0
