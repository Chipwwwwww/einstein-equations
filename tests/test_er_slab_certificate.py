"""Validated nonflat vacuum cylindrical-wave test, independent tensor checks."""
from fractions import Fraction as F
import sympy as sp
import pytest
from research.calc.er_slab_certificate import er_certificate


def test_rational_nonlinear_errored_geometry():
    result = er_certificate()
    assert result['radial_wave_residual_rho'] == F(1,50)
    assert result['gamma_r_constraint_residual'] == 0
    assert result['gamma_t_residual_bound'] == F(11,2500)
    assert result['null_derivative_error_bound'] == F(3,490)
    assert result['psi_error_bound'] == F(11,9800)
    assert result['independent_exact_psi_error'] == F(1,1600)
    assert result['independent_exact_gamma_error'] == F(21,20000)
    assert result['asymptotically_flat_binary_merger_certified'] is False
    assert result['metric_component_error_bound'] < F(1,100)
    assert er_certificate(slabs=4)['psi_error_bound'] < result['psi_error_bound']
    with pytest.raises(ValueError):
        er_certificate(rmin=F(0))
    with pytest.raises(ValueError):
        er_certificate(rmax=F(1))


def test_exact_cylindrical_reduced_equations_and_manufactured_residuals():
    t,r,a,eps = sp.symbols('t r a eps', real=True, positive=True)
    psi = (a+eps)*t**2 + a*r**2/2
    gamma = 2*(a+eps)**2*t**2*r**2 + a**2*r**4/4
    rho = sp.diff(psi,t,2)-sp.diff(psi,r,2)-sp.diff(psi,r)/r
    gr = sp.diff(gamma,r) -r*(sp.diff(psi,t)**2+sp.diff(psi,r)**2)
    gt = sp.diff(gamma,t) -2*r*sp.diff(psi,t)*sp.diff(psi,r)
    assert sp.factor(rho-2*eps)==0
    assert sp.factor(gr)==0
    assert sp.factor(gt-4*eps*(a+eps)*t*r**2)==0
    assert sp.factor(rho.subs(eps,0))==0
    assert sp.factor(gt.subs(eps,0))==0
    # The true data are identical up to first time derivatives at t=0.
    assert sp.simplify((psi-psi.subs(eps,0)).subs(t,0))==0
    assert sp.simplify(sp.diff(psi-psi.subs(eps,0),t).subs(t,0))==0
    assert sp.simplify((gamma-gamma.subs(eps,0)).subs(t,0))==0
    assert sp.simplify(sp.diff(gamma-gamma.subs(eps,0),t).subs(t,0))==0


def test_full_four_d_vacuum_and_nonflatness():
    t,r,z,ph = sp.symbols('t r z ph',real=True)
    a=sp.Rational(1,10)
    psi=a*(t**2+r**2/2)
    gamma=a*a*(2*t*t*r*r+r**4/4)
    chi=gamma-psi
    coords=(t,r,z,ph)
    g=sp.diag(-sp.exp(2*chi),sp.exp(2*chi),sp.exp(2*psi),r**2*sp.exp(-2*psi))
    inv=sp.diag(-sp.exp(-2*chi),sp.exp(-2*chi),sp.exp(-2*psi),sp.exp(2*psi)/r**2)
    G=[[[sp.simplify(sum(inv[i,k]*(sp.diff(g[k,j],coords[l])+
        sp.diff(g[k,l],coords[j])-sp.diff(g[j,l],coords[k]))/2
        for k in range(4))) for l in range(4)] for j in range(4)] for i in range(4)]
    Ric = [[sp.simplify(sum(sp.diff(G[k][i][j],coords[k])-
        sp.diff(G[k][i][k],coords[j])+
        sum(G[k][k][d]*G[d][i][j]-G[k][j][d]*G[d][i][k]
        for d in range(4)) for k in range(4))) for j in range(4)] for i in range(4)]
    assert all(z==0 for row in Ric for z in row)
    # Independent nonzero Riemann tensor component (upper z,t,z,t):
    # R^z_{t z t} = -psi_tt + psi_t*chi_t + psi_r*chi_r - psi_t^2.
    Rztzt=sp.simplify(-sp.diff(psi,t,2)+sp.diff(psi,t)*sp.diff(chi,t)+
                     sp.diff(psi,r)*sp.diff(chi,r)-sp.diff(psi,t)**2)
    assert Rztzt.subs({t:0,r:1}) == sp.Rational(-209,1000)


def test_independently_computed_riemann_curvature():
    t,r,z,ph=sp.symbols('t r z ph',real=True)
    a=sp.Rational(1,10)
    psi=a*(t**2+r**2/2)
    gamma=a*a*(2*t*t*r*r+r**4/4)
    chi=gamma-psi
    coords=[t,r,z,ph]
    g=sp.diag(-sp.exp(2*chi),sp.exp(2*chi),sp.exp(2*psi),r**2*sp.exp(-2*psi))
    inv=sp.diag(-sp.exp(-2*chi),sp.exp(-2*chi),sp.exp(-2*psi),sp.exp(2*psi)/r**2)
    Gamma = [[[sp.simplify(sum(inv[i,k]*(sp.diff(g[k,j],coords[l])+
        sp.diff(g[k,l],coords[j])-sp.diff(g[j,l],coords[k]))/2
        for k in range(4))) for l in range(4)] for j in range(4)] for i in range(4)]
    # R^z_{t z t} computed directly from Christoffels, no curvature ansatz.
    riem=sp.factor(sp.diff(Gamma[2][0][0],coords[2])-sp.diff(Gamma[2][2][0],coords[0])+
        sum(Gamma[2][2][k]*Gamma[k][0][0]-Gamma[2][0][k]*Gamma[k][2][0]
        for k in range(4)))
    assert sp.simplify(riem.subs({t:0,r:1})+sp.Rational(209,1000))==0
