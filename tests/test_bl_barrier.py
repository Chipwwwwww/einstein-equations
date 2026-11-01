"""Tests for the barrier certificate of notes/20 (research/calc/bl_barrier.py).

Fast tests check the ingredients: the general expansion formula (against an independent SymPy
level-set computation, and against the conformally flat formula of fold_proof), the source term,
the conventions (theta^(1) = h on a floating-point level surface), and that the true expansion from
the certified evaluator is close to h. Slow tests check that the copied integrator reproduces
fold_proof.shoot bit for bit, and the stored certificate record.
"""

import json
import math
import os

import pytest
import sympy as sp
from flint import arb, ctx

from research.calc import bl_barrier as bb
from research.calc import bl_metric as bm
from research.calc import fold_proof as fp

DATA = os.path.join(os.path.dirname(__file__), "..", "research", "calc", "data")


def _ellipsoid_case():
    x, y, z = sp.symbols("x y z", real=True)
    R2 = x**2 + y**2
    a = 1 + sp.Rational(3, 10) * z + sp.Rational(1, 5) * R2
    B = sp.Rational(2, 5) + sp.Rational(1, 7) * z
    C = sp.Rational(1, 3) + sp.Rational(1, 4) * R2
    E = sp.Rational(1, 6) + sp.Rational(1, 9) * z
    xp, e3 = [x, y, 0], [0, 0, 1]
    gam = sp.Matrix(3, 3, lambda i, j: a * (1 if i == j else 0) + B * xp[i] * xp[j] + C * e3[i] * e3[j]
                    + E * (xp[i] * e3[j] + e3[i] * xp[j]))
    Kt = sp.Matrix(3, 3, lambda i, j: (z / 3) * (1 if i == j else 0) + sp.Rational(1, 2) * xp[i] * xp[j]
                   + sp.Rational(1, 5) * e3[i] * e3[j] + sp.Rational(1, 4) * (xp[i] * e3[j] + e3[i] * xp[j]))
    return (x, y, z), gam, Kt


def test_theta_general_against_sympy_level_set():
    """A genuinely non-conformally-flat axisymmetric metric (with gamma_xz != 0) and K != 0: the
    expansion of an ellipsoid of revolution from theta_general equals tr_S K + div(nu) computed by
    SymPy from the level-set function."""
    (x, y, z), gam, Kt = _ellipsoid_case()
    A_, B_ = sp.Rational(7, 10), sp.Rational(5, 4)
    f = (x**2 + y**2) / A_**2 + z**2 / B_**2 - 1
    co = [x, y, z]
    gi = gam.inv()
    df = sp.Matrix([sp.diff(f, c) for c in co])
    nf = sp.sqrt((df.T * gi * df)[0])
    nu_up = gi * df / nf
    sq = sp.sqrt(gam.det())
    H = sum(sp.diff(sq * nu_up[i], co[i]) for i in range(3)) / sq
    trK = sum(gi[i, j] * Kt[i, j] for i in range(3) for j in range(3)) \
        - sum(Kt[i, j] * nu_up[i] * nu_up[j] for i in range(3) for j in range(3))
    th_sym = sp.lambdify((x, y, z), H + trK)
    keys = {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2), "xz": (0, 2), "xy": (0, 1), "yz": (1, 2)}
    for u in (0.3, 0.9, 1.3):
        a, b = float(A_), float(B_)
        rho, zz = a * math.sin(u), b * math.cos(u)
        rs, zs, rss, zss = a * math.cos(u), -b * math.sin(u), -a * math.sin(u), -b * math.cos(u)
        alpha = math.atan2(zs, rs)
        k1 = (zss * rs - rss * zs) / (rs**2 + zs**2) / math.hypot(rs, zs)
        k2 = math.sin(alpha) / rho
        sub = {x: rho, y: 0, z: zz}
        gv = {k: float(gam[i, j].subs(sub)) for k, (i, j) in keys.items() if k in ("xx", "yy", "zz", "xz")}
        Kv = {k: float(Kt[i, j].subs(sub)) for k, (i, j) in keys.items() if k in ("xx", "yy", "zz", "xz")}
        dg = {(l, k): float(sp.diff(gam[i, j], cl).subs(sub))
              for l, cl in (("x", x), ("y", y), ("z", z)) for k, (i, j) in keys.items()}
        mine = bb.theta_general(math.sin(alpha), math.cos(alpha), k1, k2, gv, dg, Kv, math.sqrt)
        assert abs(mine - th_sym(rho, 0, zz)) < 1e-12


def test_theta_general_conformally_flat_reduces_to_fold_proof_formula():
    """For gamma = psi^4 delta: H = psi^-2 (-alpha' - sin(alpha)/rho + 4 n.grad(psi)/psi)."""
    psi, px, pz = 2.1, -0.3, 0.45
    for al, k1, k2 in ((-0.4, 1.3, 0.7), (-1.2, -0.2, 2.5)):
        sa, ca = math.sin(al), math.cos(al)
        gam = {"xx": psi**4, "yy": psi**4, "zz": psi**4, "xz": 0.0}
        dg = {}
        for l, pl in (("x", px), ("z", pz)):
            for k in ("xx", "yy", "zz"):
                dg[(l, k)] = 4 * psi**3 * pl
            dg[(l, "xz")] = 0.0
        dg[("y", "xy")] = dg[("y", "yz")] = 0.0
        K = {"xx": 0.0, "yy": 0.0, "zz": 0.0, "xz": 0.0}
        ndl = (-sa * px + ca * pz) / psi
        want = (-k1 - k2 + 4 * ndl) / psi**2
        assert abs(bb.theta_general(sa, ca, k1, k2, gam, dg, K, math.sqrt) - want) < 1e-13


def test_psi2p_matches_fold_proof_t_coefficient():
    ctx.prec = 128
    for rho, z, al, d in ((0.31, 0.42, -0.7, 0.7662), (0.05, 0.6, -0.1, 0.7671)):
        a = [arb(v) for v in (rho, z, math.sin(al), math.cos(al), d)]
        ours = bb.psi2p(*a)
        theirs = fp.t_coefficient(*a)
        assert ours.overlaps(theirs) and abs(float((ours - theirs).mid())) < 1e-30


def _level_point_data(d, t, h, s):
    z_in, _, _, _ = bb.find_levels(d, t, h, 2e-6)
    sol = bb.level_shooter(d, t, h).shoot(z_in)
    rho, z, al = sol.sol(s)
    sa, ca = math.sin(al), math.cos(al)
    G = bb.make_source(arb(t), arb(h))
    Jb = lambda v: fp.J([arb(v)] + [arb(0)] * 5)  # noqa: E731
    k2 = sa / rho
    k1 = -k2 + G(Jb(rho), Jb(z), Jb(sa), Jb(ca), Jb(d)).c[0]
    return rho, z, sa, ca, k1, k2


def test_first_order_expansion_equals_level_on_level_surface():
    """Conventions: on a floating-point solution of the level-h shooting ODE, theta_general with the
    first-order data (g_0, d g_0, t g_1) from the jet recursion returns h (to ODE accuracy)."""
    ctx.prec = 128
    d, t, h = 0.7662, 1.2e-5, -0.75e-6
    for s in (0.05, 0.6):
        rho, z, sa, ca, k1, k2 = _level_point_data(d, t, h, s)
        Jt, g, _ = bm.series_jets((arb(rho), arb(0), arb(z)), arb(d), 2, arb(1))
        vals, d1, _ = bm.coefficients(Jt, g, arb(0), 2)
        K, ax = bb.KEYS, {"x": 0, "y": 1, "z": 2}
        gam0 = {k: vals[0][K[k]] for k in ("xx", "yy", "zz", "xz")}
        K0 = {k: arb(t) * vals[1][K[k]] for k in ("xx", "yy", "zz", "xz")}
        dg0 = {(l, k): d1[0][K[k]][ax[l]] for l in ("x", "z") for k in ("xx", "yy", "zz", "xz")}
        dg0[("y", "xy")], dg0[("y", "yz")] = d1[0][K["xy"]][1], d1[0][K["yz"]][1]
        th = bb.theta_general(arb(sa), arb(ca), k1, arb(k2), gam0, dg0, K0)
        assert abs(float(th.mid()) - h) < 1e-12


def test_true_expansion_close_to_level():
    """The true expansion (certified evaluator, order 4, tails included) on the level surface differs
    from h by far less than the margin |h| used by the certificate."""
    ctx.prec = 128
    d, t, h = 0.7662, 1.2e-5, -0.75e-6
    rho, z, sa, ca, k1, k2 = _level_point_data(d, t, h, 0.3)
    res = bm.evaluate_metric(t, rho, 0, z, d, order=4)
    assert res["status"] == "certified"
    g, dg = res["g"], res["dg"]
    K, ax = bb.KEYS, {"x": 0, "y": 1, "z": 2}
    gam = {k: g[K[k][0] + 1][K[k][1] + 1] for k in ("xx", "yy", "zz", "xz")}
    Kt = {k: dg[0][K[k][0] + 1][K[k][1] + 1] / 2 for k in ("xx", "yy", "zz", "xz")}
    dgt = {(l, k): dg[ax[l] + 1][K[k][0] + 1][K[k][1] + 1] for l in ("x", "z") for k in ("xx", "yy", "zz", "xz")}
    dgt[("y", "xy")], dgt[("y", "yz")] = dg[2][1][2], dg[2][2][3]
    th = bb.theta_general(arb(sa), arb(ca), k1, arb(k2), gam, dgt, Kt)
    assert abs(th - h).upper() < 1e-9


@pytest.mark.slow
def test_integrator_copy_reproduces_fold_proof():
    """bl_barrier.shoot with fold_proof's own source G gives fold_proof.shoot's result exactly."""
    ctx.prec = 256
    z, d = arb(0.6326595528058354), arb(0.766197428240025)
    F1, _ = fp.shoot(z, d)
    F2, _, _ = bb.shoot(z, d, fp.G)
    assert all(a.mid() == b.mid() and a.rad() == b.rad() for a, b in zip(F1.c, F2.c))


@pytest.mark.slow
def test_stored_certificate_record():
    """The stored record of the replay says certified, with the margins quoted in notes/20."""
    with open(os.path.join(DATA, "bl_barrier.json")) as f:
        rec = json.load(f)
    assert rec["certified"]
    assert rec["inner"]["theta_hi"] < 0 < rec["outer"]["theta_lo"]
    assert rec["inner"]["F_minus"][0] * rec["inner"]["F_plus"][1] < 0
    assert rec["geometry"]["nested"] and rec["geometry"]["Omega_in_U"]
