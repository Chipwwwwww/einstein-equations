"""All 50 GH second time-jet components for BL conformally flat K=0 data.

Signature (-,+,+,+), gamma1=-1 gamma2=gamma0=0, harmonic source H_a=0.
Uses p>0, p_i, p_ij, Lap p, grad Lap p; formulas hold also in
the nonvacuum artificial-cutoff region. Actual BL region has Lap p=0.
The evaluator rejects outside the rigorously protected space/time domain.
"""
from sympy import Matrix, zeros, diag, simplify, sqrt, Rational, Integer
C=Rational(3831,10000)
PAIRS=[(a,b) for a in range(4) for b in range(a,4)]

def pack(g,pi,phi):
    def ten(m):return [m[a,b] for a,b in PAIRS]
    return ten(g)+ten(pi)+sum([ten(m) for m in phi],[])

def initial_jets_from_p(p,q,H,L,Lgrad=None):
    """Return exact u0,u1,u2 where u2=DG(u0)[u1] in 50 GH fields.

    q_i=partial_i p; H_ij=partial_ij p; L=Delta p.
    Lgrad_i=partial_i Delta p. On the vacuum BL collar L=Lgrad=0.
    No symmetry of p is assumed.
    """
    p=Rational(p) if not hasattr(p,"free_symbols") else p
    q=Matrix(q);H=Matrix(H)
    Lgrad=Matrix([0,0,0]) if Lgrad is None else Matrix(Lgrad)
    S=sum(q[i]*q[i] for i in range(3))
    g0=diag(-p**-4,p**4,p**4,p**4)
    pi1=zeros(4)
    pi1[0,0]=4*S/p**12-4*L/p**11
    for i in range(3):
        for j in range(3):
            pi1[i+1,j+1]=4*(S*int(i==j)-4*q[i]*q[j])/p**4
            pi1[i+1,j+1]-=4*L*int(i==j)/p**3
    g2=-p**-2*pi1
    pi2=zeros(4)
    for i in range(3):
        pi2[0,i+1]=pi2[i+1,0]=16*q[i]*(p*L+5*S)/p**13
    phi0=[];phi2=[]
    for i in range(3):
        phi0.append(diag(4*q[i]/p**5,*([4*p**3*q[i]]*3)))
        mat=zeros(4)
        Hdot=sum(q[j]*H[i,j] for j in range(3))
        mat[0,0]=-8*Hdot/p**14+56*S*q[i]/p**15
        mat[0,0]+=4*Lgrad[i]/p**13-52*L*q[i]/p**14
        for j in range(3):
            for k in range(3):
                mat[j+1,k+1]=(16*(H[i,j]*q[k]+q[j]*H[i,k])
                               -8*Hdot*int(j==k))/p**6
                mat[j+1,k+1]-=6*q[i]*(16*q[j]*q[k]-4*S*int(j==k))/p**7
                mat[j+1,k+1]+=4*Lgrad[i]*int(j==k)/p**5-20*L*q[i]*int(j==k)/p**6
        phi2.append(mat)
    return pack(g0,zeros(4),phi0),pack(zeros(4),pi1,[zeros(4) for i in range(3)]),pack(g2,pi2,phi2)

def p_and_derivatives(x,y,z):
    p=Integer(1);q=zeros(3,1);H=zeros(3)
    for sign in (1,-1):
        xyz=Matrix([x,y,z-sign*C])
        r2=sum(v*v for v in xyz)
        if r2==0:raise ValueError("BL puncture excluded")
        r=sqrt(r2)
        p+=1/(4*r)
        for i in range(3):
            q[i]-=xyz[i]/(4*r**3)
            for j in range(3):
                H[i,j]+=(3*xyz[i]*xyz[j]/r**5-int(i==j)/r**3)/4
    return p,q,H

def physical_jets(x,y,z):
    p,q,H=p_and_derivatives(x,y,z)
    return initial_jets_from_p(p,q,H,Integer(0))

def auxiliary_jets_at_point(x0,y0,z0):
    """Definition of global torus u0,u1,u2 from the C9 auxiliary p_E.

    For off-vacuum points L and grad L are included. Symbolic derivatives
    may be expensive but are explicitly evaluable. Does not infer vacuum
    outside the physical agreement region.
    """
    import sympy as sp
    from .bl_sharp_cutoff import auxiliary_p
    x,y,z=sp.symbols("x y z",real=True)
    p=auxiliary_p(x,y,z)
    coords=(x,y,z)
    q=[sp.diff(p,v) for v in coords]
    H=[[sp.diff(p,v,w) for w in coords] for v in coords]
    L=sum(H[i][i] for i in range(3))
    dL=[sp.diff(L,v) for v in coords]
    args={x:sp.Rational(str(x0)),y:sp.Rational(str(y0)),z:sp.Rational(str(z0))}
    return initial_jets_from_p(p.subs(args),[qv.subs(args) for qv in q],
        [[h.subs(args) for h in row]for row in H],L.subs(args),
        [l.subs(args) for l in dL])

def _validated_input(t,x,y,z):
    # Do not round large rational times through string conversion.
    t,x,y,z=map(lambda q:q if isinstance(q,Rational) else Rational(str(q)),(t,x,y,z))
    from .bl_sharp_certificate import certificate
    e=int(certificate()["time_T_over_M"].split("-")[1].rstrip(")"))
    if t<0 or t>Rational(1,10**e):
        raise ValueError("time outside certified BL GH slab")
    if x*x+y*y+z*z>Rational(9,10)**2:
        raise ValueError("point outside certified target")
    if min(x*x+y*y+(z-C)**2,x*x+y*y+(z+C)**2)<Rational(3,10)**2:
        raise ValueError("point within uncertified puncture collar")
    return t,x,y,z

def full_50_field_state(t,x,y,z):
    """Complete second Taylor state on certified physical target, 50 fields."""
    t,x,y,z=_validated_input(t,x,y,z)
    u0,u1,u2=physical_jets(x,y,z)
    return [u0[i]+t*u1[i]+t*t*u2[i]/2 for i in range(50)]

def metric(t,x,y,z,precision=None):
    """Full GH approximate metric with strict input guard and error bound."""
    state=full_50_field_state(t,x,y,z)
    g=zeros(4)
    for i,(a,b) in enumerate(PAIRS):
        g[a,b]=g[b,a]=simplify(state[i]) if precision is None else state[i].evalf(precision)
    return g

if __name__=="__main__":
    from sympy import Rational as R
    g=metric(R(1,10**2832),R(1,10),0,0,precision=12)
    print(g)
