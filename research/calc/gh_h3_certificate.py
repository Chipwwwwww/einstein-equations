"""Small rational proof record for a complete undamped 50-field GH Kasner slab.

The large exact fractions of the underlying H3 automatic functional calculus
are checked and rounded OUTWARDS to a short human-auditable certificate.
The proven analytic GH constraint/subsidiary step is separate from arithmetic.
"""
from fractions import Fraction as Q
import json
from .gh_h3_majorant import majorant


def certificate():
    T=Q(1,1000);R=Q(1,10**8);eps=Q(1,10**7)
    m=majorant(T=T,R=R,eps=eps)
    # Strict rational outward enclosures of the computed full H3 graph:
    lam=Q(978,1000);Lam=Q(1022,1000)
    MF=Q(10);LF=Q(525)
    MA=[Q(21)]*3;LA=[Q(123),Q(123),Q(119)]
    Mc=[Q(66,1000),Q(66,1000),Q(58,1000)]
    JS=Q(10);JB=[Q(13)]*3
    gt=Q(33);gx=Q(65,1000)
    assert m['lambda']>=lam and m['Lambda']<=Lam
    assert m['source_H3_bound']<=MF
    assert m['source_Lipschitz_H3']<=LF
    for field,key in ((MA,'A_H3_bounds'),(LA,'A_Lipschitz_H3'),(Mc,'A_deviation_H3'),(JB,'J_B')):
        assert all(known<=upper for known,upper in zip(m[key],field)),key
    assert m['J_S']<=JS
    assert m['metric_dt_Linf']<=gt
    assert m['metric_dxi_Linf']<=gx
    assert m['Lorentzian_metric_gap']<Q(65,1000)
    assert m['V_H3']<Q(10)
    # All constants recomputed solely from short rational enclosures.
    D=JS*gt+sum(JB)*gx
    Vi=6*eps*T*T
    G=LF+48*LA[0]*Vi+144*sum(Mc)
    K=D/(2*lam)+Lam*G/lam
    rho=18*eps*T+48*sum(MA)*6*eps*T*T+LF*6*eps*T*T
    E=Lam*rho*T/(lam*(1-K*T))
    assert D < Q(333) and G < Q(553) and K < Q(750)
    assert rho<Q(4,10**9)
    assert K*T<Q(3,4)
    assert E<Q(2,10**11)<R
    assert m['residual_H3']<=rho
    assert m['certified_H3_error']<=E
    return {
      'classification':'CERTIFIED_FULL_50_GH_UNDAMPED_TORUS_KASNER_SLAB',
      'time_interval':[Q(0),T], 'R':R, 'amplitude':eps,
      'model':'Kasner p=(2/3,2/3,-1/3) harmonic t; approx perturbs three GH fields by eps*t^2*cos(x1)',
      'Sobolev':'H^3(T^3;R^50) normalized Fourier',
      'H3_product':48,'H3_commutator':144,'H3_C0':6,'H3_C1':12,
      'lambda_lower':lam,'Lambda_upper':Lam,
      'source_H3_upper':MF, 'source_H3_Jacobian_upper':LF,
      'A_H3_upper':MA, 'A_H3_Jacobian_upper':LA,
      'A_spatial_nonconstant_upper':Mc,
      'J_S_upper':JS, 'J_B_upper':JB,
      'metric_dt_Linf_upper':gt,'metric_dx_Linf_upper':gx,
      'adapted_divergence_upper':D,
      'G_upper':G,'energy_growth_upper':K,
      'residual_H3_upper':rho,
      'H3_metric_field_error_upper':E,
      'claimed_round_error_cap':Q(2,10**11),
      'exact_initial_data':'Kasner at harmonic time zero; constraint-exact',
      'subsidiary_propagation':'Lindblom et al. eqs 40-52 plus exact ADM constraints, torus no boundary',
      'full_50GH_algebraic_source':'Lindblom et al. eqs 35-39 gamma1=-1 gamma2=gamma0=0 H=0',
      'companion_known_Kasner_identity_used_for_residual_rewrite':True,
      'nonzero_residual':True,
      'binary_merger_spacetime':False,
      'physically_meaningful_black_hole_exterior_domain':False,
      'new_mathematics_relative_to_literature':False,
    }


def readable(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,list):return [readable(v) for v in x]
    if isinstance(x,dict):return {k:readable(v) for k,v in x.items()}
    return x


if __name__=='__main__':
    print(json.dumps(readable(certificate()),indent=2))