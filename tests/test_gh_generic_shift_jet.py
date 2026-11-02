"""Independently recover GH second-order reduced wave equation at generic off-shell jet.

All coordinates and first derivatives allowed; space-time metric has nonzero
shift and nonzero spatial offdiagonal entry. Spatial dPi,dPhi set zero at the
jet; first-order source supplies their time derivatives.
"""
from fractions import Fraction
import sympy as sp
from research.calc.gh_full_source import pack50,full_source,unpack50


def _make():
    h=sp.Matrix([[2,sp.Rational(1,10),0],[sp.Rational(1,10),3,sp.Rational(1,20)],
                 [0,sp.Rational(1,20),4]])
    b=sp.Matrix([sp.Rational(1,5),-sp.Rational(1,7),sp.Rational(1,11)])
    g=sp.zeros(4)
    g[0,0]=-1+(b.T*h*b)[0]
    for i in range(3):
        for j in range(3):g[i+1,j+1]=h[i,j]
        g[0,i+1]=g[i+1,0]=(h*b)[i]
    pi=sp.Matrix(4,4,lambda i,j:sp.Rational((i+1)*(j+1),31+i+j))
    ph=[sp.Matrix(4,4,lambda i,j:sp.Rational((k+1)*(i+1)*(j+1),97+k+i+j)) for k in range(3)]
    return g,pi,ph


def _rate_lapse_beta(g,dg):
    inv=g.inv(); dgi=-inv*dg*inv
    N=1/sp.sqrt(-inv[0,0]); Nd=sp.Rational(1,2)*N**3*dgi[0,0]
    beta=[-inv[0,i]/inv[0,0] for i in range(1,4)]
    betad=[-dgi[0,i]/inv[0,0]+inv[0,i]*dgi[0,0]/inv[0,0]**2 for i in range(1,4)]
    return Nd,betad


def test_generic_nonzero_shift_full_second_order_GH_identity():
    g,pi,ph=_make()
    u=pack50(g.tolist(),pi.tolist(),[x.tolist() for x in ph])
    src=full_source(u)
    fg,fpi,fphi=unpack50(src)
    inv=g.inv(); N=1/sp.sqrt(-inv[0,0]); beta=[-inv[0,i]/inv[0,0] for i in range(1,4)]
    assert N==1 and all(b!=0 for b in beta) and g[1,2]!=0
    dg0=-N*pi+sum((beta[i]*ph[i] for i in range(3)),sp.zeros(4))
    derivatives=[dg0]+ph
    n_t,b_t=_rate_lapse_beta(g,dg0)
    d2gtt=-n_t*pi - N*sp.Matrix(fpi)+sum((b_t[i]*ph[i]+beta[i]*sp.Matrix(fphi[i]) for i in range(3)),sp.zeros(4))
    d2gti=[]
    for i in range(3):
        ni,bi=_rate_lapse_beta(g,ph[i])
        d2gti.append(-ni*pi+sum((bi[j]*ph[j] for j in range(3)),sp.zeros(4)))
    christ=[[[ (derivatives[b][a,c]+derivatives[c][a,b]-derivatives[a][b,c])/2
              for c in range(4)]for b in range(4)]for a in range(4)]
    for a,b in [(i,j) for i in range(4) for j in range(i,4)]:
        q=sum(inv[c,d]*inv[e,f]*(derivatives[e][c,a]*derivatives[f][d,b]-
                  christ[a][c][e]*christ[b][d][f])
                 for c in range(4) for d in range(4) for e in range(4) for f in range(4))
        wave=inv[0,0]*d2gtt[a,b]+2*sum(inv[0,i+1]*d2gti[i][a,b] for i in range(3))
        assert sp.factor(wave-2*q)==0, (a,b,sp.factor(wave-2*q))