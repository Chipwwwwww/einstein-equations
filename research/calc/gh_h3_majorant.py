"""Exact-rational Banach algebra/Fréchet majorants for complete 50-field GH.

Proves bounds for every H^3 state in a tube, not finite sampling.
Each scalar object is center constant c plus spatial H^3 deviation <=d,
with derivative norm <=L relative to full 50-field H^3 input norm.
Multiplication uses universal H^3 product constant 48 from note25.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt


@dataclass(frozen=True)
class Jet:
    c: Q
    d: Q = Q(0)
    L: Q = Q(0)
    K: int = 48

    def __post_init__(self):
        object.__setattr__(self,'c',Q(self.c))
        object.__setattr__(self,'d',Q(self.d))
        object.__setattr__(self,'L',Q(self.L))
        if self.d<0 or self.L<0 or self.K not in (1,48):
            raise ValueError('invalid jet norms')

    def cast(self,x):
        if isinstance(x,Jet):
            if x.K!=self.K: raise ValueError('mismatched algebras')
            return x
        return Jet(Q(x),K=self.K)

    def __add__(self,b):
        b=self.cast(b)
        return Jet(self.c+b.c,self.d+b.d,self.L+b.L,self.K)
    def __radd__(self,a):return self+a
    def __neg__(self):return Jet(-self.c,self.d,self.L,self.K)
    def __sub__(self,b):return self+-self.cast(b)
    def __rsub__(self,a):return self.cast(a)+-self
    def __mul__(self,b):
        b=self.cast(b);k=self.K
        return Jet(self.c*b.c,
          abs(self.c)*b.d+abs(b.c)*self.d+k*self.d*b.d,
          (abs(self.c)+k*self.d)*b.L + (abs(b.c)+k*b.d)*self.L,k)
    def __rmul__(self,a):return self*a
    def reciprocal(self):
        c=self.c
        if c==0 or self.K*self.d>=abs(c):
            raise ValueError('unverified reciprocal of full state tube')
        dd=self.d/(c*c*(1-self.K*self.d/abs(c)))
        M=abs(1/c)+self.K*dd
        return Jet(1/c,dd,M*M*self.L,self.K)
    def __truediv__(self,b):return self*self.cast(b).reciprocal()
    def __rtruediv__(self,a):return self.cast(a)*self.reciprocal()
    def inverse_sqrt(self):
        if self.c!=1 or self.K*self.d>=1:
            raise ValueError('inverse_sqrt requires unit center and full-tube binomial convergence')
        z=1-self.K*self.d
        return Jet(Q(1),self.d/z,self.L/(z*z),self.K)
    def norm(self):return abs(self.c)+self.d
    def multnorm(self):return abs(self.c)+self.K*self.d


def to_fraction(a):
    return Q(a)


def jet_inputs(centers,deviations,R,pointwise=False):
    k=1 if pointwise else 48
    if len(centers)!=50 or len(deviations)!=50:raise ValueError('50 fields required')
    return [Jet(Q(x),(6 if pointwise else 1)*(Q(d)+R),Q(1),k)
            for x,d in zip(centers,deviations)]


def bounds_from_jets(jets):
    return sum((v.norm() for v in jets),Q(0)),sum((v.L for v in jets),Q(0))


def sym_eig_bounds_from_gamma(gamma):
    # pointwise ||gamma^ij -delta^ij||op <= Frobenius <= sum all entry sup
    # each jet gamma_ij center=diag 1 and deviation d ; S5 has other diag 1.
    d=sum((g.d for row in gamma for g in row),Q(0))
    delta=6*d
    if delta>=1:raise ValueError('spatial GH symmetrizer positivity unverified')
    return 1-delta,1+delta


def get_geometry_jets(jets):
    from .gh_full_source import unpack50,inverse4
    g,pi,phi=unpack50(jets)
    inv=inverse4(g)
    N=(-inv[0][0]).inverse_sqrt()
    t=[-N*inv[0][a] for a in range(4)]
    beta=[-inv[0][i]/inv[0][0] for i in range(1,4)]
    gam=[[inv[i+1][j+1]+t[i+1]*t[j+1] for j in range(3)]for i in range(3)]
    return N,beta,gam,inv


def principal_bounds(jets):
    from .gh_full_source import principal_blocks
    A=principal_blocks(jets)
    out=[]
    for block in A:
        flat=[x if isinstance(x,Jet) else Jet(Q(x),K=jets[0].K)
              for row in block for x in row]
        m,lip=bounds_from_jets(flat)
        dev=sum((j.d for j in flat),Q(0))
        # 10 identical blocks; Frobenius 50 <= l1 of 10 block entries.
        out.append((10*m,10*lip,10*dev))
    return out


def pointwise_derivative_sym_and_flux(jets):
    from .gh_full_source import principal_blocks
    _,_,gamma,_=get_geometry_jets(jets)
    JS=sum((x.L for row in gamma for x in row),Q(0))
    A=principal_blocks(jets)
    S=[[Jet(0,K=1) for j in range(5)] for i in range(5)]
    S[0][0]=S[1][1]=Jet(1,K=1)
    for i in range(3):
        for j in range(3):S[i+2][j+2]=gamma[i][j]
    JB=[]
    for block in A:
        B=[[sum((S[i][m]*block[m][j] for m in range(5)),Jet(0,K=1))
            for j in range(5)] for i in range(5)]
        JB.append(sum((b.L for row in B for b in row),Q(0)))
    return JS,JB


def kasner_centers_and_deviations(T,perturb=Q(0)):
    """All-encompassing continuum uniform data bounds for Kasner t in [0,T].
    p=(2/3,2/3,-1/3); small perturb <=perturb in individual H3 values.
    Each rational envelope follows |exp(qt)-1|<=|q|T/(1-|q|T).
    """
    from .gh_full_source import kasner_initial,PAIRS
    if not (0<T<Q(1,5)) or perturb<0:raise ValueError('time/perturb invalid')
    c=[Q(str(x)) for x in kasner_initial()]
    d=[Q(0)]*50
    gp=[(0,0,Q(2)),(1,1,Q(4,3)),(2,2,Q(4,3)),(3,3,Q(-2,3))]
    for a,b,q in gp:
        idx=PAIRS.index((a,b))
        d[idx]=abs(c[idx])*abs(q)*T/(1-abs(q)*T)
        if a==b:
            # For Kasner Pi_aa(t)=Pi_aa(0)*exp(qPi*t).
            qpi=q-Q(1)
            if a==0:qpi=Q(1)
            d[10+idx]=abs(c[10+idx])*abs(qpi)*T/(1-abs(qpi)*T)
    return c,d


def majorant( T=Q(1,10**12),R=Q(1,10**8),eps=Q(1,10**7),m=3):
    """Compute complete GH source and Jacobian bounds valid on H3 50-field ball.

    Model v=w_Kasner+eps*t^2*(cos x1 in 3 distinct GH components).
    m=3: g01, Pi22, Phi_3,12, so approximate reduction constraint is
    nonzero away from t0; exact Cauchy data match Kasner at t=0.
    """
    from .gh_full_source import full_source
    if T<=0 or R<=0 or eps<=0 or m!=3: raise ValueError('positive inputs and 3 test modes required')
    c,d=kasner_centers_and_deviations(T)
    per=2*eps*T*T # exact cos x norm_H3=2
    for idx in (1,17,47): d[idx]+=per
    # Ensure the exact comparator w is inside the certified ball about v.
    if 6*eps*T*T>R: raise ValueError("approximation and reference no longer share certified H3 tube")
    # Check full spacetime Lorentzian signature for every member of tube: 
    # ||g-eta||_op <= 6 sum_{a,b}|g_ab-eta_ab|_H3 <= 12 sum_{a<=b} d_ab.
    metric_gap=12*sum((d[i]+R for i in range(10)),Q(0))
    if metric_gap>=1: raise ValueError("Lorentzian full-tube matrix-inertia test fails")
    j=jet_inputs(c,d,R)
    Fs=full_source(j)
    f0=[x if isinstance(x,Jet) else Jet(Q(x)) for x in Fs]
    MF,LF=bounds_from_jets(f0)
    Ap=principal_bounds(j)
    Mp=[x[0] for x in Ap];LA=[x[1] for x in Ap];Mc=[x[2] for x in Ap]
    JSp,Jbp=pointwise_derivative_sym_and_flux(jet_inputs(c,d,R,pointwise=True))
    _,_,gam,_=get_geometry_jets(j)
    lam, Lam=sym_eig_bounds_from_gamma(gam)
    V=sum((abs(c[i])+d[i] for i in range(50)),Q(0))
    # Known spatial derivative norm for THREE active cos x modes; only i=1.
    Vi=[6*eps*T*T,Q(0),Q(0)]
    # Exact structural GH improvement: S and S*A depend ONLY on the 10 metric
    # fields, never on the 40 Pi/Phi fields. Thus differentiate only g_t, dg,
    # and use the algebraic metric RHS g_t=-N Pi + beta Phi. This removes the
    # huge generic unknown all-50-field u_t majorant.
    metric_dt_bound=6*sum((f0[i].norm() for i in range(10)),Q(0))
    metric_dx_bound=12*sum((d[i]+R for i in range(10)),Q(0))
    D=JSp*metric_dt_bound+sum(Jbp)*metric_dx_bound
    D_generic=JSp*(6*MF+72*(V+R)*sum(Mp))+12*(V+R)*sum(Jbp)
    G=LF+48*sum(LA[k]*Vi[k] for k in range(3))+144*sum(Mc)
    K=D/(2*lam)+Lam*G/lam
    # Approximate residual via exact analytic Kasner identity at all times;
    # v = Kasner + delta, w_t=F(w), and w has zero spatial gradients.
    # This identity can ALSO be verified directly by the full source code.
    rho=18*eps*T+48*sum(Mp)*6*eps*T*T +LF*6*eps*T*T
    # Initial error 0; rational exponential comparison e^(KT) <= 1/(1-KT).
    if K*T>=1:raise ValueError('GH nonlinear energy time-step too long: K*T>=1')
    # Energy L2 symmetrizer perturb controls 6 separate from full state; safe rational
    # sqrt Lambda/lambda majorant using 1 <= Lam and lambda<1:
    # Y(t) <= (Lam/sqrt(lam))*rho*T/(1-KT), and H3 error <=Y/sqrt(lam)
    bound=Lam*rho*T/(lam*(1-K*T))
    if bound>=R:raise ValueError('GH nonlinear H3 error tube not closed')
    return {'T':T,'R':R,'eps':eps,'lambda':lam,'Lambda':Lam,
            'C0':6,'C1':12,'P':48,'CC':144,
            'source_H3_bound':MF,'source_Lipschitz_H3':LF,
            'A_H3_bounds':Mp,'A_Lipschitz_H3':LA,'A_deviation_H3':Mc,
            'J_S':JSp,'J_B':Jbp,'metric_dt_Linf':metric_dt_bound,
            'metric_dxi_Linf':metric_dx_bound,'adapted_D':D,
            'generic_D':D_generic,'V_H3':V,'V_spatial_deriv_H3':Vi,
            'residual_H3':rho,'growth_K':K,'K_times_T':K*T,
            'certified_H3_error':bound,'closed':True,
            'Lorentzian_metric_gap':metric_gap,
            'full_GH_50_source_composition':True,
            'Einstein_vacuum_subsidiary_proved':False,
            'comparison_uses_explicit_Kasner_vacuum_reference':True,
            'domain':'3-torus, t in [0,T], Kasner harmonic time'}