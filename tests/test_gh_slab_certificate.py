"""Independent symbolic Einstein checks and exact-rational slab tests."""
from fractions import Fraction as F

import pytest
import sympy as sp

from research.calc.gh_slab_certificate import (
    NullError, NullSlab, gauge_wave_certificate, propagate,
)


def test_rational_slabs_match_direct_diamond_bound():
    c = gauge_wave_certificate()
    assert c["Happ_lower"] == F(1, 2)
    assert c["Happ_upper"] == F(2401, 1600)
    assert c["R_wave_bound"] == F(3, 25)
    assert c["Ric_tt_Ric_xx_bound"] == F(3, 50)
    assert c["log_error_bound"] == F(3, 800)
    assert c["null_derivative_bound"] == F(3, 100)
    assert c["metric_C0_error_bound"] == F(7203, 1275200)
    assert c["independently_known_actual_metric_C0_error"] == F(1, 1600)
    assert c["full_einstein_merger_certified"] is False
    for n in (1, 2, 3, 7, 40):
        d = gauge_wave_certificate(nslabs=n)
        assert d["log_error_bound"] == F(3, 800)
        assert d["null_derivative_bound"] == F(3, 100)
        assert d["metric_C0_error_bound"] == c["metric_C0_error_bound"]


def test_propagation_and_invalid_inputs():
    assert propagate(NullSlab(F(1, 4), F(1, 10)), NullError(F(0), F(0))) == NullError(F(1, 320), F(1, 40))
    with pytest.raises(ValueError):
        propagate(NullSlab(F(-1), F(1)), NullError(F(0), F(0)))
    with pytest.raises(ValueError):
        gauge_wave_certificate(a=F(2))


def test_symbolic_wave_residual_and_vacuum_solution():
    t, x, a, eps = sp.symbols("t x a eps", real=True)
    H = 1 + a * (x - t) + eps * t**2
    W = sp.diff(sp.log(H), t, 2) - sp.diff(sp.log(H), x, 2)
    rhs = 2 * eps * (1 + a * (x + t) - eps * t**2) / H**2
    assert sp.factor(W - rhs) == 0
    assert sp.factor(W.subs(eps, 0)) == 0
    assert sp.diff(H, t).subs(t, 0) == sp.diff(H.subs(eps, 0), t).subs(t, 0)


def test_full_four_dimensional_ricci_and_harmonic_gauge():
    t, x, y, z, eps = sp.symbols("t x y z eps", real=True)
    H = 1 + (x - t) / 2 + eps * t**2
    coords = (t, x, y, z)
    g = sp.diag(-H, H, 1, 1)
    inv = sp.diag(-1 / H, 1 / H, 1, 1)
    gamma = [[[sp.factor(sp.Rational(1, 2)*sum(
        inv[a, d] * (sp.diff(g[d,b], coords[c]) + sp.diff(g[d,c], coords[b]) - sp.diff(g[b,c], coords[d]))
        for d in range(4))) for c in range(4)] for b in range(4)] for a in range(4)]
    ric = [[sp.factor(sum(
        sp.diff(gamma[c][a][b], coords[c]) - sp.diff(gamma[c][a][c], coords[b]) +
        sum(gamma[c][c][d]*gamma[d][a][b] - gamma[c][b][d]*gamma[d][a][c] for d in range(4))
        for c in range(4))) for b in range(4)] for a in range(4)]
    residual = sp.diff(sp.log(H), t, 2) - sp.diff(sp.log(H), x, 2)
    assert sp.factor(ric[0][0] + residual / 2) == 0
    assert sp.factor(ric[1][1] - residual / 2) == 0
    assert all(ric[a][b] == 0 for a in range(4) for b in range(4) if (a,b) not in ((0,0),(1,1)))
    assert all(sp.simplify(sum(sp.diff(H*inv[mu,nu], coords[mu]) for mu in range(4))) == 0 for nu in range(4))
    assert all(sp.factor(expr.subs(eps, 0)) == 0 for row in ric for expr in row)


def test_strong_field_schwarzschild_kerr_schild_slice_control():
    # Horizon-penetrating spherical Kerr--Schild coordinates, NOT harmonic gauge.
    # Tests independent exact Ricci-flat geometry at r=2M, not a merger.
    t, r, th, ph, M = sp.symbols("t r th ph M", positive=True)
    f = 2*M/r
    g = sp.Matrix([[-1+f, f, 0, 0], [f, 1+f, 0, 0],
                   [0, 0, r**2, 0], [0, 0, 0, r**2*sp.sin(th)**2]])
    gi = sp.simplify(g.inv())
    assert sp.simplify(gi[0,0] + 1 + f) == 0
    assert sp.simplify(g[0,0].subs(r,2*M)) == 0
    assert sp.simplify(gi[0,0].subs(r,2*M) + 2) == 0
    assert sp.simplify(g[1,1].subs(r,2*M) - 2) == 0
    coords = [t,r,th,ph]
    Gamma = [[[sp.simplify(sum(gi[a,d] *
        (sp.diff(g[d,b],coords[c])+sp.diff(g[d,c],coords[b])-sp.diff(g[b,c],coords[d]))
        for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
    Ric = [[sp.simplify(sum(
        sp.diff(Gamma[c][a][b], coords[c])-sp.diff(Gamma[c][a][c],coords[b]) +
        sum(Gamma[c][c][d]*Gamma[d][a][b]-Gamma[c][b][d]*Gamma[d][a][c]
            for d in range(4)) for c in range(4)))
        for b in range(4)] for a in range(4)]
    assert all(sp.simplify(z) == 0 for row in Ric for z in row)
