"""Route J (research/calc/ck_bound.py): the certified convergence time of the BL development near the
critical surface. These tests check the ingredients against independent code and exact solutions;
the certificate itself is `python -m research.calc.ck_bound` (about 2 minutes)."""

import math

import numpy as np
from flint import acb, arb, ctx
from scipy.integrate import solve_ivp

from research.calc import ck_bound as ck

R3 = range(3)


def _data(w, d, centres=None):
    psi, dpsi, d2 = ck.bl_jets(w, acb(d), centres)
    G = [[1 / psi**4 if i == j else acb(0) for j in R3] for i in R3]
    g0 = [[[4 * psi**3 * dpsi[k] if i == j else acb(0) for j in R3] for i in R3] for k in R3]
    Dg0 = [[[[4 * (3 * psi**2 * dpsi[m] * dpsi[k] + psi**3 * d2[m][k]) if i == j else acb(0)
              for j in R3] for i in R3] for k in R3] for m in R3]
    return psi, dpsi, d2, G, g0, Dg0


def test_reduced_ricci_is_the_ricci_tensor():
    """A + B in the variables (gamma^-1, g = d gamma, d g) equals the closed formula for
    Ric(psi^4 delta) (itself checked against the general curvature code in test_petrov_bl), at a
    real point and at a complex point (a holomorphic identity)."""
    with ck.working_precision(128):
        for w in [(acb("0.3"), acb("0.1"), acb("0.45")), (acb("0.2", "0.01"), acb(0, "-0.008"), acb("0.6", "0.005"))]:
            psi, dpsi, d2, G, g0, Dg0 = _data(w, "0.7662")
            A, B = ck.ricci_parts(G, g0, Dg0)
            R0 = ck.ricci0(psi, dpsi, d2)
            assert max(float(abs(A[i][j] + B[i][j] - R0[i][j]).upper()) for i in R3 for j in R3) < 1e-25


def _numpy_B(G, g):
    """Independent implementation of B with numpy (einsum), complex floats."""
    Gl = 0.5 * (np.einsum("ijl->lij", g) + np.einsum("jil->lij", g) - g)   # Gl[l,i,j]
    Gu = np.einsum("kl,lij->kij", G, Gl)
    dG = -np.einsum("ka,mab,bl->mkl", G, g, G)
    return (np.einsum("kkl,lij->ij", dG, Gl) - np.einsum("jkl,lik->ij", dG, Gl)
            + np.einsum("kkl,lij->ij", Gu, Gu) - np.einsum("kjl,lik->ij", Gu, Gu))


def _numpy_Q(G, K):
    return 2 * K @ G @ K - np.trace(G @ K) * K


def test_algebraic_lipschitz_bounds_hold_at_random_points():
    """For one covering cell, random complex points w in it and random u, v in the weighted ball:
    |B(u) - B(v)|_F <= Gam_op^2 JBG |du_gamma|_F + JBg |du_g|_F (and the same for Q), with B, Q
    from an independent numpy implementation. A wrong Jacobian or norm in ck_bound would show."""
    rng = np.random.default_rng(1)
    s0, radii = 0.012, (1.0, 20.0, 10.0)
    rho, z, d = 0.2, 0.3, 0.7662
    with ck.working_precision(128):
        b = ck.cell_bounds(arb(rho, 0.0025), arb(z, 0.0025), arb(d), s0, radii)
    J, gop2 = b["J"], b["gam_op"]**2
    for _ in range(20):
        w = np.array([rho, 0.0, z]) + rng.uniform(-0.0025, 0.0025, 3) + 1j * rng.uniform(-s0, s0, 3) / 2
        psi, g0 = 1.0 + 0j, np.zeros((3, 3, 3), complex)
        for zc in (d / 2, -d / 2):
            X = w - np.array([0, 0, zc])
            r = np.sqrt(np.sum(X * X))
            psi = psi + 0.25 / r
        dpsi = sum(-0.25 * (w - np.array([0, 0, zc])) / np.sqrt(np.sum((w - np.array([0, 0, zc]))**2))**3 for zc in (d / 2, -d / 2))
        for k in R3:
            g0[k] = 4 * psi**3 * dpsi[k] * np.eye(3)
        gam0 = psi**4 * np.eye(3)

        def rnd(shape, rad):
            x = rng.normal(size=shape) + 1j * rng.normal(size=shape)
            return x / np.linalg.norm(x) * rad * rng.uniform(0, 1)
        u = (rnd((3, 3), radii[0]), rnd((3, 3, 3), radii[1]), rnd((3, 3), radii[2]))
        v = (rnd((3, 3), radii[0]), rnd((3, 3, 3), radii[1]), rnd((3, 3), radii[2]))
        Gu, Gv = np.linalg.inv(gam0 + u[0]), np.linalg.inv(gam0 + v[0])
        dB = np.linalg.norm(_numpy_B(Gu, g0 + u[1]) - _numpy_B(Gv, g0 + v[1]))
        dQ = np.linalg.norm(_numpy_Q(Gu, u[2]) - _numpy_Q(Gv, v[2]))
        dgam, dgg, dK = (np.linalg.norm(u[i] - v[i]) for i in R3)
        # the two mean-value bounds used in cell_bounds
        assert dB <= gop2 * J["BG"] * dgam + J["Bg"] * dgg
        assert dQ <= gop2 * J["QG"] * dgam + J["QK"] * dK


def test_certify_refuses_when_a_hypothesis_fails():
    """No success-like output: a gamma-ball too large for invertibility (q >= 1), or a complex
    neighbourhood reaching a puncture (Re q > 0 fails), makes certify return ok = False."""
    box = [(0.0, 0.01, 0.6, 0.62)]
    r = ck.certify(0.012, (40.0, 20.0, 10.0), 0.7662, 0.7663, boxes=box, cell=0.01)
    assert r["ok"] is False
    near = [(0.0, 0.01, 0.37, 0.39)]                    # the upper puncture sits at z = 0.3831
    r = ck.certify(0.012, (1.0, 20.0, 10.0), 0.7662, 0.7663, boxes=near, cell=0.01)
    assert r["ok"] is False


def test_working_precision_is_restored():
    old = ctx.prec
    ck.certify(0.012, (1.0, 20.0, 10.0), 0.7662, 0.7663, boxes=[(0.25, 0.26, 0.0, 0.01)], cell=0.01)
    assert ctx.prec == old


def test_profile_boxes_contain_the_floating_point_critical_curve():
    """The stored certified profile (research/calc/data/ck_bound.json) contains the critical curve
    integrated in floating point from z0 = z*, d = d_c (independent integrator)."""
    import json
    boxes = json.load(open("research/calc/data/ck_bound.json"))["profile"]
    d, z0 = 0.766197428240, 0.632659552806

    def rhs(s, y):
        rho, zt, u, v = y
        psi, pr, pz = 1.0, 0.0, 0.0
        for Z in (zt, zt + d):
            ir = 1 / math.hypot(rho, Z)
            psi += 0.25 * ir
            pr -= 0.25 * rho * ir**3
            pz -= 0.25 * Z * ir**3
        P = 4 * (u * pz - v * pr) / psi
        return [u, v, rho - v * P, u * P]
    ev = lambda s, y: y[1] + d / 2  # noqa: E731
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 40), [1e-8, z0 - d / 2, 1e-8, 0.0], events=ev, rtol=1e-11, atol=1e-14,
                    dense_output=True, max_step=0.02)
    for t in np.linspace(5, sol.t[-1], 400):
        rho, zt = sol.sol(t)[:2]
        zz = zt + d / 2
        assert any(b[0] - 1e-6 <= rho <= b[1] + 1e-6 and b[2] - 1e-6 <= zz <= b[3] + 1e-6 for b in boxes)


def test_schwarzschild_benchmark():
    """One centre (M = 1), U = tube of radius 0.02 around the isotropic sphere r = 1. The exact data
    give |Ric(gamma_0)|_F = sqrt(6) psi^4 M / R^3 (coordinate components, R the areal radius): the
    certified sup K0 must enclose its maximum over U and stay within a modest factor of it (the
    complex neighbourhood and cells add only a little). The certified time must lie below the exact
    crunch time of geodesic slicing at the innermost point, pi sqrt(R^3 / 8M)."""
    boxes = []
    for k in range(64):
        th0, th1 = k * math.pi / 128, (k + 1) * math.pi / 128
        rs = [math.sin(th0), math.sin(th1)]
        zs = [math.cos(th0), math.cos(th1)]
        boxes.append((min(rs) - 1e-9, max(rs) + 1e-9, min(zs) - 1e-9, max(zs) + 1e-9))
    centres = [(acb(1) / 2, acb(0))]
    r = ck.certify(0.012, (1.0, 20.0, 10.0), 0.0, 0.0, boxes=boxes, cell=0.01, centres=centres)
    assert r["ok"]
    exact = []
    for rr in (0.98, 1.0, 1.02):
        psi = 1 + 0.5 / rr
        R = rr * psi**2
        exact.append(math.sqrt(6) * psi**4 / R**3)
    assert max(exact) <= r["K0"] <= 3 * max(exact)
    R_in = 0.98 * (1 + 0.5 / 0.98) ** 2
    assert r["T_min"] < math.pi * math.sqrt(R_in**3 / 8)
