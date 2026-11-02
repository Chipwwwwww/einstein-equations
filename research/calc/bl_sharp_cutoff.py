"""Sharp exact-rational finite-Sobolev extension of equal-mass BL puncture data.

Finite smoothstep S_m has m continuous derivatives (m=9), not C-infinity.
Uniform multi-index derivative envelopes through order 8 use exact Bernstein
coefficients and exponential-generating-function Faa di Bruno bounds.
"""
from fractions import Fraction as Q
from math import factorial as fact, comb

A=Q(12,100)
B=Q(18,100)
OUT=Q(3,2)
PHYS_OUT=Q(11,10)
C=Q(3831,10000)


def cutoff_coefficients(m=9):
    if m<8:raise ValueError("at least eight continuous cutoff derivatives needed")
    k=Q(fact(2*m+1),fact(m)**2)
    return {m+j+1:k*Q((-1)**j*comb(m,j),m+j+1) for j in range(m+1)}


def bernstein_derivative_bounds(poly,k):
    """For Q(z)=sum a_j z^j=sum b_i binom(d,i) z^i(1-z)^(d-i),
    |Q(z)|<=max_i |b_i| on [0,1]. Exact rational conversion."""
    result=[Q(1)]
    for n in range(1,k+1):
        mono={j-n:c*Q(fact(j),fact(j-n)) for j,c in poly.items() if j>=n}
        d=max(mono)
        b=[sum((mono.get(j,Q(0))*Q(comb(i,j),comb(d,j))
                 for j in range(i+1)),Q(0)) for i in range(d+1)]
        result.append(max(abs(x) for x in b))
    return result


def quadratic_composition(k,first,second,S):
    f=3*first;g=3*second/2
    out=[Q(1)]
    for n in range(1,k+1):
        z=sum((Q(S[j],fact(j))*comb(j,n-j)*f**(2*j-n)*g**(n-j)
               for j in range(1,n+1) if j<=n<=2*j),Q(0))
        out.append(Q(fact(n))*z)
    return out


def convolution(a,b,k):
    return [sum((comb(n,j)*a[j]*b[n-j] for j in range(n+1)),Q(0))
            for n in range(k+1)]


def poly_convolution(a,b,k):
    return [sum((a[j]*b[n-j] for j in range(n+1)),Q(0))
            for n in range(k+1)]


def positive_power(p0,P,k,power):
    a=[Q(1)]+[Q(0)]*k
    b=[p0]+[P[j]/fact(j) for j in range(1,k+1)]
    for _ in range(power):a=poly_convolution(a,b,k)
    return [fact(n)*a[n] for n in range(k+1)]


def negative_power(P,k,power):
    # p>=1. This is the coefficientwise majorant (1-Q(z))^-power.
    q=[Q(0)]+[P[j]/fact(j) for j in range(1,k+1)]
    r=[Q(1)]+[Q(0)]*k
    for n in range(1,k+1):
        r[n]=sum((q[j]*r[n-j] for j in range(1,n+1)),Q(0))
    a=[Q(1)]+[Q(0)]*k
    for _ in range(power):a=poly_convolution(a,r,k)
    return [fact(n)*a[n] for n in range(k+1)]


def exact_ceil_log10(v):
    """Smallest integer exponent k with v<=10**k, for v>0."""
    v=Q(v)
    if v<=0:raise ValueError("positive rational required")
    k=len(str(v.numerator))-len(str(v.denominator))
    while v>Q(10)**k:k+=1
    while v<=Q(10)**(k-1):k-=1
    return k


def derivative_envelopes(k=8,m=9):
    S=bernstein_derivative_bounds(cutoff_coefficients(m),k)
    inner=quadratic_composition(k,2*B/(B*B-A*A),2/(B*B-A*A),S)
    outer=quadratic_composition(k,2*OUT/(OUT*OUT-PHYS_OUT*PHYS_OUT),
                                   2/(OUT*OUT-PHYS_OUT*PHYS_OUT),S)
    # Mixed derivatives of 1/r: <=6^n n!/r^(n+1), by polarization and
    # the Legendre polynomial bound for repeated directional derivatives.
    inv=[Q(6**n*fact(n))*A**(-n-1) for n in range(k+1)]
    triple=convolution(convolution(inner,outer,k),inv,k)
    P=[Q(1)+triple[0]/2]+[triple[n]/2 for n in range(1,k+1)]
    pmax=1+Q(1,2*A)
    powers={'+4':positive_power(pmax,P,k,4),'-4':negative_power(P,k,4),
            '-2':negative_power(P,k,2)}
    metric=[max(powers['+4'][j],powers['-4'][j]) for j in range(k+1)]
    U={s:exact_ceil_log10(50*2**s*max(metric[:s+2]))+1
       for s in range(3,min(8,k))}
    return {'k':k,'m':m,'step':S,'inner':inner,'outer':outer,'inverse_r':inv,
            'p_derivatives':P,'metric_derivatives':metric,'powers':powers,
            'H3_H7_log10_bounds':U,'pmax':pmax,
            'globally_vacuum':False,'auxiliary_regular_only_Cm':True}


def finite_step(z,m=9):
    import sympy as sp
    z=sp.sympify(z)
    if z.is_number and z<=0:return sp.Integer(0)
    if z.is_number and z>=1:return sp.Integer(1)
    s=sum((sp.Rational(a.numerator,a.denominator)*z**j
           for j,a in cutoff_coefficients(m).items()),sp.Integer(0))
    if z.is_number:return sp.factor(s)
    return sp.Piecewise((0,z<=0),(1,z>=1),(s,True))


def auxiliary_p(x,y,z,m=9):
    """Periodic C^9 field agreeing EXACTLY with BL on |x|<=1.1,r_+/- >=.2."""
    import sympy as sp
    x,y,z=map(sp.sympify,(x,y,z))
    rr=x*x+y*y+z*z
    if rr.is_number and rr>=OUT*OUT:return sp.Integer(1)
    co=finite_step((sp.Rational(OUT*OUT)-rr)/sp.Rational(OUT*OUT-PHYS_OUT*PHYS_OUT),m)
    a=sp.Rational(C.numerator,C.denominator)
    terms=[]
    for sign in (1,-1):
        q=x*x+y*y+(z-sign*a)**2
        if q.is_number and q<=A*A:
            terms.append(sp.Integer(0));continue
        ci=finite_step((q-sp.Rational(A*A))/sp.Rational(B*B-A*A),m)
        terms.append(ci/sp.sqrt(q))
    return 1+co*sum(terms)/4


if __name__=="__main__":
    import json
    x=derivative_envelopes()
    print(json.dumps({'cutoff_Cm':x['m'],'H3_to_H7_log10':x['H3_H7_log10_bounds'],
       'cutoff_derivative_log10':[exact_ceil_log10(v) for v in x['step']],
       'g_derivative_log10':[exact_ceil_log10(v) for v in x['metric_derivatives']],
       'global_auxiliary_vacuum':False},indent=2))
