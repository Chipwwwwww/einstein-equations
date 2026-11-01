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


def _A_matrix(G):
    """The linear map D -> A(G, D) (9 x 81) from ck_bound.ricci_parts itself."""
    cols = []
    for n in range(81):
        D = [0] * 81
        D[n] = 1
        Dg = [[[[D[((m * 3 + a) * 3 + b) * 3 + c] for c in R3] for b in R3] for a in R3] for m in R3]
        g = [[[0.0 for _ in R3] for _ in R3] for _ in R3]
        A, _ = ck.ricci_parts([[float(G[i][j]) for j in R3] for i in R3], g, Dg)
        cols.append([A[i][j] for i in R3 for j in R3])
    return np.array(cols).T


def test_exact_norms_of_the_principal_part():
    """||A(I, .)|| = 2 and |A(G, D)|_F <= (||T|| / 2) |G|_F |D|_F with ||T||^2 = 14 + 2 sqrt(33): the
    characteristic polynomials of the integer matrices 4 L L^T and T^T T, with L and T read off
    ck_bound.ricci_parts (A(E_kl, D)_ij = T_klij(D) / 2)."""
    from flint import fmpz_mat, fmpz_poly
    L = _A_matrix(np.eye(3))
    M = np.rint(4 * L @ L.T).astype(int)
    assert np.allclose(4 * L @ L.T, M)
    x = fmpz_poly([0, 1])
    assert fmpz_mat(M.tolist()).charpoly() == (x - 16) * (x - 10) ** 5 * (x - 6) ** 3
    rows = []
    for k in R3:
        for l in R3:
            E = np.zeros((3, 3))
            E[k, l] = 1
            rows.append(2 * _A_matrix(E))           # T_kl.. as a 9 x 81 block
    T = np.vstack(rows)
    TT = np.rint(T.T @ T).astype(int)
    assert np.allclose(T.T @ T, TT)
    assert fmpz_mat(TT.tolist()).charpoly() == x ** 54 * (x - 8) ** 15 * (x - 16) ** 6 * (x ** 2 - 28 * x + 64) ** 3
    assert abs(float(ck._half_norm_T().mid()) - math.sqrt(14 + 2 * math.sqrt(33)) / 2) < 1e-15


def test_sharp_principal_part_bounds_at_random_points():
    """|A(G, D)|_F <= (2 |c| + (||T||/2) |E|_F) |D|_F for G = c I + E, and <= (||T||/2) |G|_F |D|_F."""
    rng = np.random.default_rng(3)
    half = float(ck._half_norm_T().upper())
    for _ in range(50):
        c = complex(rng.normal(), rng.normal())
        E = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) * rng.uniform(0, 0.3)
        G = c * np.eye(3) + E
        D = rng.normal(size=81) + 1j * rng.normal(size=81)
        Dg = D.reshape(3, 3, 3, 3)
        A, _ = ck.ricci_parts([[G[i, j] for j in R3] for i in R3], [[[0j] * 3 for _ in R3] for _ in R3],
                              [[[[Dg[m, a, b, cc] for cc in R3] for b in R3] for a in R3] for m in R3])
        nA = np.linalg.norm(np.array([[A[i][j] for j in R3] for i in R3]))
        nD = np.linalg.norm(D)
        assert nA <= (2 * abs(c) + half * np.linalg.norm(E)) * nD * (1 + 1e-12)
        assert nA <= half * np.linalg.norm(G) * nD * (1 + 1e-12)


def test_outward_rounding_helpers_and_exact_cells():
    with ck.working_precision(128):
        x = arb(1) / 3
        assert arb(ck.up(x)) >= x and arb(ck.down(x)) <= x
        c = ck._cell(37, 0.005)
        assert c.lower() <= arb(37) * arb(0.005) and c.upper() >= arb(38) * arb(0.005)
        nxt = ck._cell(38, 0.005)
        assert c.upper() >= nxt.lower()                    # adjacent cells leave no gap
        cells = ck.covering_cells([(0.25, 0.26, 0.0, 0.01)], arb("0.032"), 0.005)
        ks = {k for _, k in cells}
        assert min(ks) * 0.005 <= -0.032 and (max(ks) + 1) * 0.005 >= 0.042


def test_baouendi_goulaouic_time_formula_is_the_corrected_one():
    """Baouendi-Goulaouic (12) reads a < R / (16 R C + 8 M) (s0 = 1). For the session-12 constants
    (C = 1.5745916, K0 = 194.45955, R_K = 10, s0 = 0.014) this is 5.115e-4, below the 5.557e-4 that the
    misread min(s0 / (16 C), R_K / (8 K0)) gave (notes/17 §1)."""
    a, B = ck.time_bg(1.5745916432444342, 19.445954549021442, 0.014)
    assert 5.11e-4 < a < 5.12e-4
    assert a < min(0.014 / (16 * 1.5745916432444342), 10 / (8 * 194.45954549021442))
    assert abs(B - 4 * a * 19.445954549021442) < 1e-12


def test_ck_prime_conditions_and_its_reduction_to_the_original():
    """Theorem CK' (notes/17 §2): with beta = 1/2 and no zero-order part, a_max = 1 / (16 C / s0 + 6 M)
    (the original structure with the invariant ball's radius optimised); every returned (a, B)
    satisfies 2 a k < 1 and 3 a M < (1 - beta)(1 - 2 a k), and B = a M / ((1 - beta)(1 - a k))."""
    C, M, s0 = 1.5, 20.0, 0.014
    a, B = ck.time_ck2({"g": C, "K": 0.0, "gamma": 0.0}, {"g": 0.0, "K": 0.0, "gamma": 0.0}, M, s0, 0.5, 0.999999)
    amax = 1 / (16 * C / s0 + 6 * M)
    assert 0.999 * amax < a <= amax
    for beta in (0.3, 0.6, 0.9):
        Crow, Lrow = {"gamma": 0.0, "g": 1.4, "K": 0.5}, {"gamma": 80.0, "g": 0.0, "K": 75.0}
        a, B = ck.time_ck2(Crow, Lrow, M, s0, beta, 0.95)
        with ck.working_precision(128):          # exact balls, not floats (B is an outward-rounded bound)
            ab, bb, Mb = arb(a), arb(beta), arb(M)
            k = max(4 * arb(Crow[A]) / (bb * arb(s0)) + arb(Lrow[A]) / (1 - bb) for A in Crow)
            assert 2 * ab * k < 1 and 3 * ab * Mb < (1 - bb) * (1 - 2 * ab * k)
            assert arb(B) >= ab * Mb / ((1 - bb) * (1 - ab * k)) or abs(B / float((ab * Mb / ((1 - bb) * (1 - ab * k))).mid()) - 1) < 1e-12


def test_lemma_constants_of_ck_prime_numerically():
    """Sampled diagnostic of the two integrals behind CK' (A = a(1 - s), 0 <= t < A):
    Lemma 1': int_0^t (1 - u/A)^(-beta) du <= A / (1 - beta);
    Lemma 2': (A/a)(1 - t/A)^beta * int_0^t 4 a^2 (A + u)^(beta - 1) (A - u)^(-1 - beta) du <= 4 a / beta."""
    from scipy.integrate import quad
    a = 1.0
    for beta in (0.3, 0.5, 0.8, 0.95):
        for s in (0.0, 0.3):
            A = a * (1 - s)
            for frac in (0.2, 0.9, 0.999):
                t = frac * A
                i1 = quad(lambda u: (1 - u / A) ** (-beta), 0, t)[0]
                assert i1 <= A / (1 - beta) * (1 + 1e-9)
                i2 = quad(lambda u: 4 * a * a * (A + u) ** (beta - 1) * (A - u) ** (-1 - beta), 0, t, limit=200)[0]
                assert (A / a) * (1 - t / A) ** beta * i2 <= 4 * a / beta * (1 + 1e-7)


def test_imaginary_quadrant_boxes_cover_the_hermitian_ball():
    """Every imaginary part with Im x, Im y >= 0 and |Im| < s0 lies in some box of imaginary_boxes."""
    rng = np.random.default_rng(5)
    s0 = 0.03
    with ck.working_precision(128):
        boxes = ck.imaginary_boxes(s0, 3)
        for _ in range(300):
            v = rng.normal(size=3)
            v = v / np.linalg.norm(v) * s0 * rng.uniform(0, 0.999)
            v[0], v[1] = abs(v[0]), abs(v[1])
            assert any(all(b[i].lower() <= v[i] <= b[i].upper() for i in range(3)) for b in boxes)


def test_zero_midpoint_squares_are_finite():
    """python-flint returns nan for x ** 2 when mid(x) = 0; frob_upper must not (it multiplies)."""
    assert math.isfinite(ck.frob_upper([acb(arb(0, 1e-5)), acb(0, arb(0, 1e-5)), acb(1)]))
