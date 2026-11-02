"""Evaluable BL harmonic-gauge finite-time approximation on protected region.

Uses all 50 first-order GH equations for first time jet u1 and exact
second-time derivative of 10 metric components. Exact SymPy arithmetic
when inputs are rational. Not a numerical PDE solve or claimed merger.
"""
from fractions import Fraction
from sympy import Rational, Integer, sqrt, Matrix, zeros, simplify, N
from .gh_full_source import full_source, pack50, unpack50

C=Rational(3831,10000) # d/2 = .3831


def p_and_derivatives(x,y,z):
    p=Integer(1)
    grad=[Integer(0)]*3
    hess=zeros(3)
    for s in (+1,-1):
        dx=Matrix([x,y,z-s*C])
        r2=sum(q*q for q in dx)
        if r2==0:raise ValueError('BL punctures excluded')
        r=sqrt(r2)
        p+=1/(4*r)
        for i in range(3):
            grad[i]-=dx[i]/(4*r**3)
            for j in range(3):
                hess[i,j]+=(3*dx[i]*dx[j]/r**5-int(i==j)/r**3)/4
    return p,grad,hess


def physical_initial_and_first_jet(x,y,z,verify_against_source=False):
    """Exact full-GH initial jet in the BL vacuum region Delta p=0.

    Pi1_00=4 |grad p|^2/p^12, Pi1_0i=0,
    Pi1_ij=4 (|grad p|^2 delta_ij-4 p_i p_j)/p^4.
    No spherical/axial reduction is used to derive the 50-field equations.
    """
    p,dp,ddp=p_and_derivatives(x,y,z)
    g=Matrix.diag(-p**-4,p**4,p**4,p**4)
    D2=sum(d*d for d in dp)
    pi1=zeros(4)
    pi1[0,0]=4*D2/p**12
    for i in range(3):
        for j in range(3):
            pi1[i+1,j+1]=4*(D2*int(i==j)-4*dp[i]*dp[j])/p**4
    g2=-p**-2*pi1
    if verify_against_source:
        phi=[Matrix.diag(4*p**-5*dp[i],*[4*p**3*dp[i]]*3) for i in range(3)]
        initial=pack50(g.tolist(),zeros(4).tolist(),[m.tolist() for m in phi])
        Fg,Fpi,Fphi=unpack50(full_source(initial))
        lap=zeros(4)
        lap[0,0]=-20*p**-6*D2
        for j in range(1,4):lap[j,j]=12*p**2*D2
        derived=Matrix(Fpi)-p**-6*lap
        if any(simplify(q)!=0 for q in derived-pi1):
            raise ValueError('derived initial full-GH Pi jet disagrees with independent source')
        if any(simplify(q)!=0 for row in Fg for q in row):
            raise ValueError('metric first time jet nonzero')
        if any(simplify(q)!=0 for matrix in Fphi for row in matrix for q in row):
            raise ValueError('Phi first time jet nonzero')
    return g,pi1,g2


def evaluate_50_field_state(t,x,y,z):
    """Exact symbolic 50-field approximate state, defined on target at certified t."""
    metric=evaluate_metric(t,x,y,z)
    p,dp,_=p_and_derivatives(*map(lambda w: Rational(str(w)),(x,y,z)))
    _,pi1,_=physical_initial_and_first_jet(*map(lambda w: Rational(str(w)),(x,y,z)))
    phi=[Matrix.diag(4*p**-5*dp[i],*[4*p**3*dp[i]]*3) for i in range(3)]
    tt=t if isinstance(t,Rational) else Rational(str(t))
    return pack50(metric.tolist(),(tt*pi1).tolist(),[v.tolist() for v in phi])


def evaluate_metric(t,x,y,z,precision=None):
    t,x,y,z=map(lambda v: v if hasattr(v,'free_symbols') else Rational(str(v)),(t,x,y,z))
    from .bl_torus_majorant import certificate
    if not (t.is_number and t.is_real and t>=0):
        raise ValueError('time must be real nonnegative')
    Texp=int(certificate()['time_T'].split('-')[1].rstrip(')'))
    if t>Rational(1,10**Texp):
        raise ValueError('time outside the rigorously certified BL GH slab')
    # Exact rational domain test; never classify points using rounded floats.
    if not all(isinstance(v,Rational) for v in (t,x,y,z)):
        raise ValueError('evaluator requires rational time and coordinates')
    r1sq=x*x+y*y+(z-C)**2
    r2sq=x*x+y*y+(z+C)**2
    if not (x*x+y*y+z*z<=Rational(9,10)**2 and
            r1sq>=Rational(3,10)**2 and r2sq>=Rational(3,10)**2):
        raise ValueError('point not inside proven protected target; no BL time certificate here')
    g0,pi1,g2=physical_initial_and_first_jet(x,y,z)
    result=g0+t*t*g2/2
    return result.applyfunc(lambda v: N(v,precision) if precision else simplify(v))

if __name__=='__main__':
    g=evaluate_metric(Rational(1,10**29045),Rational(0),Rational(0),Rational(0),precision=12)
    print(g)