"""Complete 50-field Lindblom--Scheel--Kidder--Owen--Rinne GH RHS.

Signature (-,+,+,+), gamma1=-1,gamma2=gamma0=0 by default.
Source equations (35)-(39) of gr-qc/0512093 v3. H_a is prescribed
as a function of coordinates only, optionally with known partial derivatives.
All operations are scalar arithmetic, allowing SymPy or rigorous majorants.
Only pointwise algebra is implemented here; analytic GH constraints require
zero initial Einstein and subsidiary data and compatible boundary conditions.
"""
from fractions import Fraction
HALF = Fraction(1, 2)
import sympy as sp

PAIRS = [(a,b) for a in range(4) for b in range(a,4)]


def sym4(values):
    if len(values)!=10:
        raise ValueError('10 symmetric entries required')
    q=[[0 for j in range(4)] for i in range(4)]
    for (a,b),v in zip(PAIRS,values):
        q[a][b]=q[b][a]=v
    return q


def flat10(q):
    return [q[a][b] for a,b in PAIRS]


def unpack50(u):
    if len(u)!=50:
        raise ValueError('full GH system requires 50 entries')
    return sym4(u[:10]),sym4(u[10:20]),[sym4(u[20+10*i:30+10*i]) for i in range(3)]


def pack50(g,pi,phi):
    return flat10(g)+flat10(pi)+sum([flat10(m) for m in phi],[])


def inverse4(M):
    """Arithmetic-only inverse via exact fraction-free-like row reduction.

    No pivoting: for field validation a nonzero pivot must be certified on
    the entire state tube, not merely sampled at the central state.
    """
    a=[[M[i][j] for j in range(4)]+[int(i==j) for j in range(4)]
       for i in range(4)]
    for k in range(4):
        q=a[k][k].reciprocal() if hasattr(a[k][k], "reciprocal") else sp.S.One/a[k][k]
        a[k]=[q*x for x in a[k]]
        for i in range(4):
            if i==k:continue
            fac=a[i][k]
            a[i]=[a[i][j]-fac*a[k][j] for j in range(8)]
    return [row[4:] for row in a]


def full_source(u,H=None,dH=None,gamma0=0,geometry_override=None):
    """Return F in u_t + A^k(u)du/dx^k=F for gamma1=-1,gamma2=0.

    H_a(x) and dH[a][b]=partial_a H_b fixed known coordinate functions,
    not independent GH dynamical fields. For arbitrary algebraic H(x,g),
    dH must include its metric derivatives through its first-order fields.
    """
    if gamma0 != 0:
        raise NotImplementedError("Only gamma0=gamma2=0 is audited; no damping-source claim")
    g,pi,phi=unpack50(u)
    if geometry_override is None:
        gu=inverse4(g)
        N=(-gu[0][0]).inverse_sqrt() if hasattr(gu[0][0], "inverse_sqrt") else 1/sp.sqrt(-gu[0][0])
        normal=[-N*gu[0][a] for a in range(4)]
        beta=[-gu[0][i]/gu[0][0] for i in range(1,4)]
        gamma=[[gu[i+1][j+1]+normal[i+1]*normal[j+1] for j in range(3)] for i in range(3)]
    else:
        # Exact same algebraic source but inverse/lapse etc carry externally
        # proven reference-centered Banach-algebra bounds.
        gu,N,normal,beta,gamma=geometry_override
    dt=[[ -N*pi[a][b]+sum(beta[i]*phi[i][a][b] for i in range(3))
         for b in range(4)] for a in range(4)]
    dg=[dt]+phi
    Gamma=[[[ (dg[b][a][c]+dg[c][a][b]-dg[a][b][c])*HALF
              for c in range(4)] for b in range(4)] for a in range(4)]
    if H is None: H=[0]*4
    if dH is None: dH=[[0]*4 for _ in range(4)]
    if len(H)!=4 or len(dH)!=4 or any(len(row)!=4 for row in dH):
        raise ValueError("GH gauge input must have 4 covectors and a 4x4 derivative")
    Gamma_trace=[sum(gu[c][d]*Gamma[a][c][d] for c in range(4) for d in range(4)) for a in range(4)]
    Fpi=[[0]*4 for _ in range(4)]
    for a,b in PAIRS:
        quadr=0
        for c in range(4):
            for d in range(4):
                term=sum(gamma[i][j]*phi[i][c][a]*phi[j][d][b]
                         for i in range(3) for j in range(3))
                term-=pi[c][a]*pi[d][b]
                term-=sum(gu[e][f]*Gamma[a][c][e]*Gamma[b][d][f]
                          for e in range(4) for f in range(4))
                quadr+=gu[c][d]*term
        derivH=(dH[a][b]+dH[b][a])*HALF
        # covariant ∇_(a H_b) uses Γ^c_ab H_c
        derivH-=sum(gu[c][d]*Gamma[d][a][b]*H[c]
                    for c in range(4) for d in range(4))
        normalnormal=sum(normal[c]*normal[d]*pi[c][d] for c in range(4) for d in range(4))
        transport=sum(normal[c]*pi[c][i+1]*gamma[i][j]*phi[j][a][b]
                      for c in range(4) for i in range(3) for j in range(3))
        constraint=0
        if gamma0:
            constraint=sum((int(c==a)*(-N if b==0 else 0)+int(c==b)*(-N if a==0 else 0)
                 -g[a][b]*normal[c])*(H[c]+Gamma_trace[c])
                 for c in range(4))
            # t_b = -N delta_b0, and 2δ^c_(a t_b) above.
        val=2*N*quadr-2*N*derivH-N*normalnormal*pi[a][b]/2-N*transport+N*gamma0*constraint
        Fpi[a][b]=Fpi[b][a]=val
    Fphi=[]
    for i in range(3):
        mat=[[0]*4 for _ in range(4)]
        for a,b in PAIRS:
            q=sum(normal[c]*normal[d]*phi[i][c][d] for c in range(4) for d in range(4))
            term= N*q*pi[a][b]*HALF
            term+=N*sum(gamma[j][k]*normal[c]*phi[i][j+1][c]*phi[k][a][b]
                        for c in range(4) for j in range(3) for k in range(3))
            mat[a][b]=mat[b][a]=term
        Fphi.append(mat)
    return pack50(dt,Fpi,Fphi)


def principal_blocks(u):
    """Fifty by fifty A^k as 10 identical variable 5x5 blocks, gamma1=-1 gamma2=0."""
    g,_,_=unpack50(u);gu=inverse4(g)
    N=(-gu[0][0]).inverse_sqrt() if hasattr(gu[0][0], "inverse_sqrt") else 1/sp.sqrt(-gu[0][0])
    normal=[-N*gu[0][a] for a in range(4)]
    beta=[-gu[0][i]/gu[0][0] for i in range(1,4)]
    gamma=[[gu[i+1][j+1]+normal[i+1]*normal[j+1] for j in range(3)] for i in range(3)]
    # Each original 5 vector per ab=(psi,pi,phi1,phi2,phi3);
    # matrix blocks only, avoid allocating dense 50x50 matrices.
    blocks=[]
    for k in range(3):
        A=[[0]*5 for _ in range(5)]
        A[1][1]=-beta[k]
        for i in range(3):
            A[1][i+2]=N*gamma[k][i]
            A[i+2][1]=N*int(i==k)
            A[i+2][i+2]=-beta[k]
        blocks.append(A)
    return blocks


def kasner_initial():
    g=sym4([-1,0,0,0,1,0,0,1,0,1])
    p=sym4([2,0,0,0,sp.Rational(-4,3),0,0,sp.Rational(-4,3),0,sp.Rational(2,3)])
    return pack50(g,p,[sym4([0]*10) for _ in range(3)])


def kasner_at(t):
    from sympy import exp
    q=(sp.Rational(2,3),sp.Rational(2,3),sp.Rational(-1,3))
    g=sp.diag(-exp(2*t),*[exp(2*p*t) for p in q]);N=exp(t)
    Pi=sp.diag(2*exp(t),*[-2*p*exp((2*p-1)*t) for p in q])
    z=sp.zeros(4)
    return pack50([[g[i,j] for j in range(4)] for i in range(4)],
                  [[Pi[i,j] for j in range(4)] for i in range(4)],
                  [[[0]*4 for _ in range(4)] for _ in range(3)])