"""Exact-rational derivatives for a genuinely smooth torus extension of BL data.

No assertion of global vacuum outside physical protected region. C-infinity
cutoff has analytically bounded derivatives (not a mesh check).
"""
from fractions import Fraction as F
from math import comb, factorial

D=F(3831,10000) # separation d/2, d=0.7662
A=F(12,100); B=F(18,100)
PHYS_IN=F(20,100); TARGET_IN=F(30,100)
PHYS_OUT=F(11,10); OUT=F(3,2); TARGET_OUT=F(9,10)


def theta_polynomials(n):
    # P_n(r), theta^(n)(s)=exp(-r)P_n(r) with r=1/s>0
    p={0:1}
    result=[]
    for i in range(n+1):
        result.append(p.copy())
        q={}
        for m,c in p.items(): q[m+2]=q.get(m+2,0)+c
        for m,c in p.items():
            if m:q[m+1]=q.get(m+1,0)-m*c
        p={m:c for m,c in q.items() if c}
    return result


def cutoff_derivatives(k):
    """Integer rigorous sup |smoothstep^(n)| on real line, 0<=n<=k."""
    H=[]
    for p in theta_polynomials(k):
        H.append(sum(abs(c)*(m**m if m else 1) for m,c in p.items()))
    # max |exp(-r)r^m| <= m^m; denom>=exp(-2)>1/8 on [0,1].
    S=[1]
    for n in range(1,k+1):
        S.append(8*(H[n]+sum(comb(n,j)*2*H[j]*S[n-j] for j in range(1,n+1))))
    return H,S


def quadratic_composition_bound(k,first,second,S):
    """Explicit Faà di Bruno coefficient bound for any mixed spatial derivative.

    q(x) quadratic. Q(h)=first*(h1+h2+h3)+second/2*(h1²+h2²+h3²).
    For multi-index alpha, |D^alpha chi(q(x))| <= alpha! [h^alpha]
    sum_{j=1}^|alpha| S_j/j! Q(h)^j. We overbound by the
    univariate substitution h1=h2=h3=z so coefficient sums include all alpha.
    Thus the bound for any |alpha|=n is n!*[z^n]sum_j S_j/j!*(3first*z+3second/2*z²)^j.
    """
    f=3*first;g=3*second/2
    out=[F(1)]
    for n in range(1,k+1):
        z=F(0)
        for j in range(1,n+1):
            z+=F(S[j],factorial(j))*comb(j,n-j)*f**(2*j-n)*g**(n-j) if j<=n<=2*j else 0
        out.append(F(factorial(n))*z)
    return out


def derivative_envelopes(k=8):
    H,S=cutoff_derivatives(k)
    denin=B**2-A**2
    denout=OUT**2-PHYS_OUT**2
    inner=quadratic_composition_bound(k,2*B/denin,2/denin,S)
    outer=quadratic_composition_bound(k,2*OUT/denout,2/denout,S)
    inv=[F(6**n*factorial(n),1)*A**(-n-1) for n in range(k+1)]
    # 3-factor Leibniz: product derivatives of outer cutoff,
    # inner cutoff, inverse radius. B_n is uniform for all |alpha|<=n.
    C=max(inner+outer)
    V=max(inv)
    P=F(1)+max(F(1,2)*3**n*C*C*V for n in range(k+1))
    # p>=1; upper p<=1+1/(2A) since each cutoff<=1.
    pmax=1+1/(2*A)
    # Universal composition bound on all derivatives of p^4,p^-4,p^-2
    # to order k, from multivariate Faà di Bruno/Leibniz partition bound.
    # For n<=k, <=(k+4)!^3*(10 P)^k. Sufficient for powers +/-2,+/-4.
    G=F(factorial(k+4)**3)*(10*P)**k
    # H^m norm <=2^m*max_{|alpha|<=m} sup |D^alpha field|
    # (normalized torus Fourier norm). For Pi0=0, Phi_i=D_i g,
    # choose U for all H^{k-1} GH initial fields and inverse fields.
    U=F(2**(k-1))*G
    return dict(k=k,H=H,S=S,inner=inner,outer=outer,inv=inv,
                smoothstep_global_bound=C,reciprocal_bound=V,
                p_derivative_bound=P,p_pointwise_upper=pmax,
                g_derivative_bound=G,Hk_minus_1_component_bound=U)


def ceil_log10_fraction(x):
    if x<=0:raise ValueError("positive rational required")
    # exact decimal order from rational digit counts, no floats.
    n,d=x.numerator,x.denominator
    k=len(str(n))-len(str(d))
    while (n < d*10**k if k>=0 else n*10**(-k) < d):k-=1
    while (n > d*10**(k+1) if k+1>=0 else n*10**(-k-1) > d):k+=1
    return k+1


def smooth_step(z):
    """Explicit C-infinity flat transition; exact symbolic for rational z."""
    import sympy as sp
    z=sp.sympify(z)
    if z<=0:return sp.Integer(0)
    if z>=1:return sp.Integer(1)
    a=sp.exp(-1/z)
    b=sp.exp(-1/(1-z))
    return a/(a+b)


def bl_auxiliary_p(x,y,z):
    """Exact smooth positive auxiliary BL conformal factor on the 2pi-torus.

    Inputs are inside the canonical fundamental cube. Periodic extension is
    smooth because p_aux=1 for |x|>=1.5 < pi. Outside puncture cutoffs the
    same function equals the actual BL p on Omega_phys.
    """
    import sympy as sp
    x,y,z=map(sp.sympify,(x,y,z))
    rr=x*x+y*y+z*z
    if rr>=OUT**2:return sp.Integer(1)
    outer=smooth_step((sp.Rational(OUT**2)-rr)/sp.Rational(OUT**2-PHYS_OUT**2))
    sums=sp.Integer(0)
    for sign in (1,-1):
        r2=x*x+y*y+(z-sign*sp.Rational(D))**2
        if r2<=A**2: continue
        cc=smooth_step((r2-sp.Rational(A**2))/sp.Rational(B**2-A**2))
        sums+=cc/sp.sqrt(r2)
    return 1+outer*sums/4


def envelope_summary():
    b=derivative_envelopes()
    return {'C_derivative_exponent':ceil_log10_fraction(b['smoothstep_global_bound']),
      'p_derivative_exponent':ceil_log10_fraction(b['p_derivative_bound']),
      'g_derivative_exponent':ceil_log10_fraction(b['g_derivative_bound']),
      'reference_H7_component_exponent':ceil_log10_fraction(b['Hk_minus_1_component_bound']),
      'pmax_upper':str(b['p_pointwise_upper']),
      'physical_domain':'|x|<1.1 and r_+,r_->0.20',
      'target_domain':'|x|<=0.9 and r_+,r_->=0.30',
      'euclidean_buffer':str(min(PHYS_OUT-TARGET_OUT,TARGET_IN-PHYS_IN)),
      'exact_agreement_with_BL_inside_physical_domain':True,
      'global_auxiliary_Einstein_vacuum':False}

if __name__=='__main__':
    print(envelope_summary())