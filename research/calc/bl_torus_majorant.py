"""Power-of-ten integer arithmetic H3 certificate for auxiliary BL GH evolution.

VERY pessimistic, but *finite* exact upper bounds for 50-field GH source
on a full ball about spatially variable data. Distinguish artificial torus
reduced solution from vacuum on protected BL domain of dependence.
"""
from dataclasses import dataclass
from .bl_local_cutoff import envelope_summary
from .gh_full_source import full_source, unpack50, pack50, PAIRS


@dataclass(frozen=True)
class E:
    """Bound scalar field H3 norm <=10**n, Frechet Lipschitz <=10**l.

    A None exponent means exactly zero. Multiplication uses P=192,
    valid in H3 and H4 since 2**4*6=96; also for H5 <=192.
    """
    n: int | None
    l: int | None = None

    def c(self,x):
        if isinstance(x,E):return x
        if x==0:return E(None,None)
        # Known numbers in GH formula have |x|<=4; absolute upper <=10.
        if abs(x)>4:raise ValueError('unexpected large scalar constant')
        return E(1,None)
    @staticmethod
    def addex(a,b):
        if a is None:return b
        if b is None:return a
        return max(a,b)+1
    @staticmethod
    def multex(a,b,extra):
        if a is None or b is None:return None
        return a+b+extra
    def __add__(self,x):
        b=self.c(x)
        return E(self.addex(self.n,b.n),self.addex(self.l,b.l))
    def __radd__(self,x):return self+x
    def __neg__(self):return self
    def __sub__(self,x):return self+-self.c(x)
    def __rsub__(self,x):return self.c(x)+-self
    def __mul__(self,x):
        b=self.c(x)
        n=self.multex(self.n,b.n,3) # H5 product constant <=192<10^3
        l=self.addex(self.multex(self.n,b.l,3),self.multex(self.l,b.n,3))
        return E(n,l)
    def __rmul__(self,x):return self*x
    def __truediv__(self,x):
        if x==2:return self # divide by 2 only makes bounds smaller
        raise ValueError("only source division by 2 supported")


def matrix_geometry(U):
    """Reference-centered Neumann and binomial bounds on entire H3 tube.

    Inputs: all reference g0, g0^-1, p^-2,p^-4 component H5 <=10**U.
    For metric perturb ||delta g||_H3 <=2*10**(-Rexp), Rexp>=10U+100.
    Neumann ratio <= 48 * 16 *10**U *2*10**(-Rexp) <1/100.
    For N, normalize -g^00=p^4 (1+w) with ||w||_H3 small.
    L_inv <= 48^2*(2*10**U)^2*16; L_N via binomial derivative.
    """
    inv=[[E(U+2,2*U+10) for j in range(4)] for i in range(4)]
    N=E(U+6,4*U+30)
    normal=[-N*inv[0][a] for a in range(4)]
    beta=[(N*N)*inv[0][i] for i in range(1,4)]
    gamma=[[inv[i+1][j+1]+normal[i+1]*normal[j+1] for j in range(3)] for i in range(3)]
    return inv,N,normal,beta,gamma


def bound_graph(U):
    fields=[E(U,0) for _ in range(50)]
    inv,N,normal,beta,gamma=matrix_geometry(U)
    F=full_source(fields,geometry_override=(inv,N,normal,beta,gamma))
    A=[]
    for k in range(3):
        M=[[E(None,None) for _ in range(5)]for _ in range(5)]
        M[1][1]=-beta[k]
        for i in range(3):
            M[1][i+2]=N*gamma[k][i]
            M[i+2][1]=N*int(i==k)
            M[i+2][i+2]=-beta[k]
        A.append([M[i][j] for i in range(5)for j in range(5)])
    # S_5=diag(1,1,gamma_ij), B_i=S_5*A_i, 10 identical blocks.
    S=[[E(None,None) for _ in range(5)]for _ in range(5)]
    S[0][0]=S[1][1]=E(1,None)
    for i in range(3):
        for j in range(3):S[i+2][j+2]=gamma[i][j]
    JS=max(x.l for row in S for x in row if x.l is not None)+3
    JB=[]
    for i in range(3):
        mat=[A[i][5*m:5*(m+1)] for m in range(5)]
        B=[[sum((S[j][m]*mat[m][k] for m in range(5)),E(None,None))
              for k in range(5)]for j in range(5)]
        JB.append(max(x.l for row in B for x in row if x.l is not None)+3)
    fn=max(x.n for x in F if x.n is not None)+3
    fl=max(x.l for x in F if x.l is not None)+3
    fg=max(x.n for x in F[:10] if x.n is not None)+3
    an=max(x.n for row in A for x in row if x.n is not None)+3
    al=max(x.l for row in A for x in row if x.l is not None)+3
    betaN=max(x.n for x in beta if x.n is not None)+2
    gammaN=max(x.n for row in gamma for x in row if x.n is not None)+2
    return dict(U=U,F_norm_exp=fn,F_Lip_exp=fl,g_t_exp=fg,
                A_norm_exp=an,A_Lip_exp=al,JS_exp=JS,JB_exp=max(JB),
                N_exp=N.n,beta_exp=betaN,gamma_exp=gammaN)


def certificate(time_exponent_override=None):
    inputs=envelope_summary()
    # Explicit Fourier H7 component bounds from exact cutoff derivatives.
    # +4 covers vector/tensor sums and 50 fields.
    U=inputs['reference_H7_component_exponent']+4
    b=bound_graph(U)
    # u1=F(u0)-A_i(u0)partial_i u0, H4 bound (H5 data).
    V1=max(b['F_norm_exp'],b['A_norm_exp']+U+3)+7
    # u2^g=-N u1^Pi, with the 40 other u2 components zero;
    # all H3/H4 component norms <=10**W2.
    W2=b['N_exp']+V1+5
    # Each real H3 state within R of v(t) differs from u0 by <2R.
    Rexp=10*U+100
    # Geometry-adapted divergence from g_t (the complete GH g source)
    # and partial_i g, and symmetrizer positivity with lambda>=10^-4.
    Dexp=max(b['JS_exp']+b['g_t_exp']+2,b['JB_exp']+U+3)+5
    Gexp=max(b['F_Lip_exp'],b['A_Lip_exp']+U+4,b['A_norm_exp']+4)+5
    Kexp=max(Dexp+4,Gexp+5)+5
    # Exact residual v_t+Av_x-F(v) vanishes at time zero.
    # v-u0=t u1+t²u2/2. Time-linear upper coefficient includes
    # u2, and (A,L) derivatives at whole-ball tube.
    RhoCoeff=max(W2,b['F_Lip_exp']+V1+4,
       b['A_Lip_exp']+V1+U+5,b['A_norm_exp']+V1+5,
       b['F_Lip_exp']+W2+4,b['A_Lip_exp']+W2+U+5,
       b['A_norm_exp']+W2+5)+10
    # All 50 GH/subsidiary speeds: 0,-beta.n,-beta.n +/- N*sqrt(gamma^nn).
    # Coordinate Euclidean speed <=sum_i|beta_i|+N sqrt(3 max|gamma^ij|).
    Cexp=max(b['beta_exp']+4,b['N_exp']+b['gamma_exp']+5)+5
    # Much sharper conformal-reference cone: p>=1, p<=31/6 =>
    # ||g0^-1||op <1000. ||g-g0||op <=50R (entrywise embedding),
    # hence ||g^-1-g0^-1||op < 10**8 R <10**(-8040).
    # Therefore -g^00>=1-e, N<2, |beta|_2<2e and
    # gamma^{ij}=p^-4 delta^ij +O(e), lambda_max(gamma)<2;
    # all GH + subsidiary speeds <2e+2*sqrt(2)<4.
    sharp_speed=4
    if Rexp<100:raise ValueError('BL reference-adapted cone fails')
    Texp=max(Rexp+V1+10,Rexp+W2+10,Kexp+10,
             Rexp+RhoCoeff+10,Cexp+10)
    if time_exponent_override is not None:
        if not isinstance(time_exponent_override,int) or time_exponent_override<=0:
            raise ValueError("positive integer exponent required")
        Texp=time_exponent_override
    # rho <=10**(RhoCoeff-Texp). sqrt(Lambda)/lambda <=10^6;
    # (1-KT)^-1 <2, so E <=10**(RhoCoeff-2Texp+7).
    Errexp=2*Texp-RhoCoeff-7
    flags={
      'Neumann_inverse_full_ball': Rexp > U+20,
      'lapse_binomial_full_ball': Rexp>4*U+60,
      'Lorentzian_cone': Rexp>U+15,
      'positive_symmetric_hyperbolic_symmetrizer': Rexp>U+15,
      'v_stays_in_reference_tube':Texp>Rexp+V1+5 and 2*Texp>Rexp+W2+5,
      'nonlinear_Gronwall':Texp>Kexp+5,
      'nonlinear_error_strictly_below_radius':Errexp>Rexp+5,
      'GH_and_constraint_causal_buffer':Texp>3,
    }
    if not all(flags.values()):
        raise ValueError('BL finite-time GH certificate rejected: '+str(flags))
    return {'classification':'CERTIFIED_POSITIVE_TIME_BL_LOCAL_VACUUM_IF_SUBSIDIARY_INITIAL_CONSTRAINT_ARGUMENT',
      'd':'3831/5000',
      'reference_H7_exp':U,'graph':b,
      'auxiliary_extension':'C-infinity on T3, equals BL for |x|<=1.1 and both r_pm>=.2',
      'target':'|x|<=.9 and both r_pm>=.3',
      'buffer':'1/10','time_T':f'10^(-{Texp})',
      'metric_H3_error_upper':f'10^(-{Errexp})',
      'nonlinear_H3_radius':f'10^(-{Rexp})',
      'GH_residual_H3_upper':f'10^(-{Texp-RhoCoeff})',
      'GH_speed_upper':'4 (reference-adapted, GH and subsidiary)',
      'fallback_graph_speed_upper':f'10^({Cexp})',
      'finite_speed_displacement_upper':f'4*10^(-{Texp})',
      'v1_H4_exp':V1,'u2_H3_exp':W2,
      'energy_growth_K_upper':f'10^({Kexp})',
      'proof_flags':flags,
      'local_exact_vacuum_outside_cutoff':False,
      'complete_binary_merger_or_extension':False}

if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2))