"""Route J: a certified lower bound on the convergence time of the Brill-Lindquist development near
the critical surface, by the abstract (Nirenberg-Nishida) Cauchy-Kovalevskaya theorem.

Mathematical statement (notes/15-route-j-convergence.md, Proposition T):

  Let U be the compact rotation- and reflection-invariant tube of radius W_TUBE around the critical
  minimal surface S_c (its meridian profile enclosed by a validated Lohner run, profile_boxes), and
  I = [d_lo, d_hi]. For every d in I the vacuum development of equal-mass time-symmetric BL data in
  geodesic slicing (lapse 1, shift 0), written as the first-order system below, has a unique real-
  analytic solution on U x {|t| < T_min}. Its Taylor series in t converges there, uniformly on U, with
  sup |u(t)|_w <= 4 T_min K0 / R_K < 1/2 (weighted norm below).

System. Unknowns gamma_ij, g_kij (= d_k gamma_ij), K_ij (repository sign, K_ij = (1/2) d_t gamma_ij):
    d_t gamma_ij = 2 K_ij,      d_t g_kij = 2 d_k K_ij,
    d_t K_ij     = -[A_ij + B_ij] + 2 K_ik G^kl K_lj - (G^kl K_kl) K_ij ,   G = gamma^-1,
where A + B is the Ricci tensor written in (gamma, g, dg): A is linear in dg with coefficients from G
(A_ij = (1/2) G^kl [d_k g_ijl + d_k g_jil - d_k g_lij - d_j g_ikl - d_j g_kil + d_j g_lik]) and B is
algebraic (Christoffel products, with d_m G^kl = -G^ka g_mab G^bl). The constraint g = d gamma
propagates exactly (d_t(g - d gamma) = 0), so solutions of the system solve the ADM evolution
equations; the ADM constraints then propagate (notes/15 §2).

Abstract CK theorem (Nishida 1977, J. Differential Geom. 12, 629-633; proof read in Baouendi-
Goulaouic, Journees EDP 1977, 1-8, numdam JEDP_1977____1_0): on a decreasing scale X_s (0 < s <= s0) of Banach
spaces, if for all 0 < s' < s <= s0 and u, v in the R-ball of X_s
    (H1) |F(u) - F(v)|_{s'} <= C |u - v|_s / (s - s'),     (H2) |F(0)|_s <= M / (s0 - s),
then du/dt = F(u), u(0) = 0 has a unique solution with sup |u(t)|_s < R for |t| < a (1 - s/s0), for any
a < s0 min(1/(16 C), R/(8 M)), and then sup |u(t)|_s <= 4 a M / s0 (constants re-derived in notes/15
§3 from the proof's inequalities (10)-(12)). With R = 1 in the weighted norm below and
M = s0 K0 / R_K, K0 = sup |Ric(gamma_0)|_F, this gives T_min = min(s0 / (16 C), R_K / (8 K0)).

Scale. X_s = bounded holomorphic maps on Omega_s = {w in C^3 : dist(w, U) < s} (Hermitian distance),
with the weighted sup norm |u|_w = max(|u_gamma|_F / R_gamma, |u_g|_F / R_g, |u_K|_F / R_K) (Frobenius
norms, radii = RADII) and ball radius R = 1. Cauchy's estimate gives |d_k u|_{s'} <= |u|_s / (s - s')
and |grad u|_F <= sqrt(3) |u|_s / (s - s').

Bounds (per covering cell, complex ball arithmetic, python-flint acb, with d an interval). Unweighted
block coefficients c_AB, in |F_A(u) - F_A(v)| <= sum_B c_AB |u_B - v_B|_s / (s - s'):
    c_gammaK = 2 s0,   c_gK = 2 sqrt(3),   c_KK = s0 J_QK,
    c_Kg     = 3 sqrt(3) Gam_F + s0 J_Bg,
    c_Kgamma = s0 (J_BG + J_QG) Gam_op^2 + 3 Gam_op^2 (s0 G1 + sqrt(3) R_g),
    Gam_op = sup |G|_op <= sup|psi|^-4 / (1 - q), q = R_gamma sup|psi|^-4 < 1, Gam_F = sqrt(3) Gam_op,
    G1 = sup |grad g0|_F, J_* = Frobenius norms of interval Jacobians (forward-mode duals over acb),
    K0 = sup |Ric(gamma_0)|_F.
The weighted Lipschitz constant is C = max_A sum_B (R_B / R_A) c_AB (lipschitz_C), M = s0 K0 / R_K.
Sharp variant (sharp=True, notes/16 §2): exact principal-part norms ||A(I, .)|| = 2 and
|A(G, D)|_F <= (||T|| / 2) |G|_F |D|_F, ||T||^2 = 14 + 2 sqrt(33), give
    c_Kg = sqrt(3) (2 sup|psi|^-4 + (||T||/2) sqrt(3) eps_op) + s0 J_Bg,   and ||T||/2 for 3 in c_Kgamma.
All scalar bounds are formed in ball arithmetic and rounded outward (up / down); cells are exact
intervals and the d-interval is an exact hull.
Reduction by symmetry: all these norms are invariant under real rotations about the axis and under
z -> -z (equal masses), so the sup over Omega_s0 is a sup over complex points whose real part lies in
the meridian half-plane y = 0: the covering uses cells in (rho, z) with imaginary parts in
[-s0, s0]^3 (a superset of the Hermitian s0-ball).

Run: python -m research.calc.ck_bound classic   (writes research/calc/data/ck_bound.json)
     python -m research.calc.ck_bound sharp     (writes research/calc/data/ck_bound_sharp.json)
"""

import json
import math
import os
import sys
import time
from contextlib import contextmanager

from flint import acb, arb, ctx

PREC = 128           # bits for all ball arithmetic here (set locally, restored afterwards)
W_TUBE = 0.02        # tube radius of U around the critical surface (units M = 1)
S0 = 0.012           # complex width of the largest domain of the analytic scale
RADII = (1.0, 20.0, 10.0)  # ball radii for (gamma, g, K) (weights of the norm); from a small search
CELL = 0.005         # covering cells in (rho, z)
D_LO = 0.7661974282399   # lower end of the d-interval (below d_c = 0.766197428240 +- 1.1e-13)
S0_SHARP = 0.014     # parameters of the sharp certificate (from a small search, notes/16 §2)
RADII_SHARP = (0.25, 22.0, 10.0)


@contextmanager
def working_precision(bits=PREC):
    old = ctx.prec
    ctx.prec = bits
    try:
        yield
    finally:
        ctx.prec = old


# ---------------------------------------------------------------------------
# The critical surface's meridian profile, enclosed by a validated Lohner run (Theorem C's box)
# ---------------------------------------------------------------------------
def profile_boxes(h=0.25, ds=0.004, L=14, sigma_max=60.0):
    """Boxes (rho_lo, rho_hi, z_lo, z_hi) whose union contains the profile of S_c in z >= 0: the
    sigma-system step boxes of the Lohner set started at z0 in the Theorem C box, d in its box, run
    until the whole set lies below the plane. The true critical curve is one of the curves of the
    set (Theorem C: z* and d_c lie in these boxes), and every point of it up to its first plane
    crossing lies in some step box."""
    from research.calc import gap3_fold as gf
    from research.calc import lohner
    slots = (0,)
    out = []
    with working_precision(128):
        f = gf.make_f_sigma(slots)
        S = gf.start(gf.Z_C_LO, gf.Z_C_HI, gf.D_C_LO, gf.D_C_HI, slots)
        blocks = [[0, 1, 2, 3]]
        half = arb(1) / 2
        sig = 0.0
        while sig < sigma_max:
            X = S.hull()
            hh = min(h, 2.0 ** math.floor(math.log2(ds / float(X[0].upper()))))
            S, B = lohner.step_param(f, S, hh, L, blocks)
            sig += hh
            zB = B[1] + B[4] * half
            pad = 1e-12
            out.append((float(B[0].lower()) - pad, float(B[0].upper()) + pad,
                        float(zB.lower()) - pad, float(zB.upper()) + pad))
            if gf.z_abs(S) < 0:
                return out
    raise AssertionError("profile run did not reach the plane")


def covering_cells(boxes, inflate, cell):
    """Grid cells (i, k), standing for [i cell, (i+1) cell] x [k cell, (k+1) cell], covering the union of
    the boxes inflated by `inflate` (componentwise, a superset of the Euclidean inflation). The index
    ranges are computed in ball arithmetic and rounded outward, so no needed cell is dropped by a
    floating-point floor()."""
    cells = set()
    a_inf, a_cell = arb(inflate), arb(cell)          # inflate may be an arb

    def lo(v):
        return math.floor(down((arb(v) - a_inf) / a_cell))

    def hi(v):
        return math.floor(up((arb(v) + a_inf) / a_cell))
    for r0, r1, z0, z1 in boxes:
        for i in range(lo(r0), hi(r1) + 1):
            for k in range(lo(z0), hi(z1) + 1):
                cells.add((i, k))
    return sorted(cells)


# ---------------------------------------------------------------------------
# Forward-mode duals over acb (gradient as a list), for interval Jacobians
# ---------------------------------------------------------------------------
class Du:
    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v, self.g = v, g

    def __add__(a, b):
        if isinstance(b, Du):
            return Du(a.v + b.v, [x + y for x, y in zip(a.g, b.g)])
        return Du(a.v + b, a.g)

    __radd__ = __add__

    def __neg__(a):
        return Du(-a.v, [-x for x in a.g])

    def __sub__(a, b):
        return a + (-b)

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        if isinstance(b, Du):
            return Du(a.v * b.v, [a.v * y + b.v * x for x, y in zip(a.g, b.g)])
        return Du(a.v * b, [x * b for x in a.g])

    __rmul__ = __mul__


def seed(values, offset, n):
    return [Du(v, [acb(1) if j == offset + i else acb(0) for j in range(n)])
            for i, v in enumerate(values)]


# ---------------------------------------------------------------------------
# The reduced Ricci tensor and the quadratic K terms (generic arithmetic: acb or Du)
# ---------------------------------------------------------------------------
R3 = range(3)


def ricci_parts(G, g, Dg=None):
    """A_ij (needs Dg[m][a][b][c] = d_m g_abc) and B_ij of the reduced Ricci tensor, with
    G = gamma^-1 (3x3 nested lists) and g[k][i][j] = d_k gamma_ij."""
    zero = G[0][0] * 0
    Gl = [[[(g[i][j][l] + g[j][i][l] - g[l][i][j]) * 0.5 for j in R3] for i in R3] for l in R3]
    Gu = [[[sum((G[k][l] * Gl[l][i][j] for l in R3), zero) for j in R3] for i in R3] for k in R3]
    dG = [[[-sum((G[k][a] * g[m][a][b] * G[b][l] for a in R3 for b in R3), zero) for l in R3]
           for k in R3] for m in R3]
    B = [[sum((dG[k][k][l] * Gl[l][i][j] - dG[j][k][l] * Gl[l][i][k]
               + Gu[k][k][l] * Gu[l][i][j] - Gu[k][j][l] * Gu[l][i][k]
               for k in R3 for l in R3), zero) for j in R3] for i in R3]
    if Dg is None:
        return None, B
    A = [[sum((G[k][l] * 0.5 * (Dg[k][i][j][l] + Dg[k][j][i][l] - Dg[k][l][i][j]
                                - Dg[j][i][k][l] - Dg[j][k][i][l] + Dg[j][l][i][k])
               for k in R3 for l in R3), zero) for j in R3] for i in R3]
    return A, B


def quad_K(G, K):
    zero = G[0][0] * 0
    tr = sum((G[k][l] * K[l][k] for k in R3 for l in R3), zero)
    return [[2 * sum((K[i][k] * G[k][l] * K[l][j] for k in R3 for l in R3), zero) - tr * K[i][j]
             for j in R3] for i in R3]


# ---------------------------------------------------------------------------
# Brill-Lindquist data on a complex ball
# ---------------------------------------------------------------------------
def bl_jets(w, d, centres=None):
    """psi, d_i psi, d_i d_j psi at the complex point(s) w = (x, y, z) (acb balls), bare masses 1/2
    at z = +-d/2 (or centres = [(m/2, z_c), ...] as acb, for psi = 1 + sum (m/2) / r). The
    holomorphic extension of 1/r is q^(-1/2), q = sum (w - p)^2, on the principal
    branch: valid (analytic and equal to the real function on real points) where Re q > 0, which is
    checked. Returns None if the check fails."""
    psi, dpsi = acb(1), [acb(0)] * 3
    d2 = [[acb(0)] * 3 for _ in R3]
    for m2, c in (centres or [(acb(1) / 4, d / 2), (acb(1) / 4, -d / 2)]):
        X = [w[0], w[1], w[2] - c]
        q = X[0] * X[0] + X[1] * X[1] + X[2] * X[2]
        if not q.real > 0:
            return None
        ir = q.rsqrt()
        ir3 = ir * ir * ir
        ir5 = ir3 * ir * ir
        m = m2                              # m/2 (bare mass m)
        psi = psi + m * ir
        dpsi = [dpsi[i] - m * X[i] * ir3 for i in R3]
        d2 = [[d2[i][j] + m * (3 * X[i] * X[j] * ir5 - (ir3 if i == j else 0)) for j in R3]
              for i in R3]
    return psi, dpsi, d2


def ricci0(psi, dpsi, d2):
    """Ric of psi^4 delta for harmonic psi (holomorphic identity; checked in tests)."""
    grad2 = sum((a * a for a in dpsi), acb(0))
    return [[-2 * d2[i][j] / psi + 6 * dpsi[i] * dpsi[j] / psi**2
             - (2 * grad2 / psi**2 if i == j else 0) for j in R3] for i in R3]


def up(x):
    """A float >= every point of the arb ball x (outward rounding of the upper endpoint)."""
    return math.nextafter(float(x.upper()), math.inf)


def down(x):
    """A float <= every point of the arb ball x."""
    return math.nextafter(float(x.lower()), -math.inf)


def frob_upper(entries):
    """A float >= sqrt(sum |e|^2) over complex balls e (outward rounded)."""
    s = sum((abs(e) * abs(e) for e in entries), arb(0))   # not ** 2: arb x ** 2 is nan for mid(x) = 0
    return up(arb(s.upper()).sqrt())        # the ball can straddle 0: use its upper end


# Exact norms of the principal part (notes/16 §2): A(I, D) has operator norm 2 (4 L L^T has
# characteristic polynomial (x - 16)(x - 10)^5 (x - 6)^3), and D -> T(D) (the six-term combination in
# A) has norm sqrt(14 + 2 sqrt(33)) (T^T T has characteristic polynomial
# x^54 (x - 8)^15 (x - 16)^6 (x^2 - 28 x + 64)^3). Checked exactly in tests/test_ck_bound.py.
NORM_L = 2
def _half_norm_T():
    return (arb(14) + 2 * arb(33).sqrt()).sqrt() / 2


def cell_bounds(rho, z, d, s0, radii, centres=None, sharp=False, im=None):
    """Unscaled block Lipschitz coefficients c_AB (F_A against u_B, in the form c / (s - s')) and
    K0 = sup |Ric(gamma_0)|_F over the complex box {Re x in rho, Re y = 0, Re z in z, imaginary
    parts in [-s0, s0]^3}, d in the real ball d, for the ball |u_gamma|_F < R_g0, |u_g|_F < R_g1,
    |u_K|_F < R_K with radii = (R_g0, R_g1, R_K). All scalar bounds are computed in ball arithmetic
    and rounded outward to floats (up / down). sharp=False reproduces the session-11 constants
    (bilinear constant 3 for A); sharp=True uses the exact principal-part norms (NORM_L, ||T|| / 2).
    Returns a dict, or None if a check fails (the principal-branch condition Re q > 0, |psi| > 0,
    or invertibility q < 1)."""
    Rgam, Rg, RK = (arb(r) for r in radii)
    S0 = arb(s0)
    if im is None:
        im = (arb(0, s0),) * 3
    w = (acb(rho, im[0]), acb(0, im[1]), acb(z, im[2]))
    jets = bl_jets(w, acb(d), centres)
    if jets is None:
        return None
    psi, dpsi, d2 = jets
    apsi = abs(psi)
    if not apsi > 0:
        return None
    ip4 = arb((arb(1) / arb(apsi.lower()) ** 4).upper())          # >= sup |psi|^-4, exact arb
    q = Rgam * ip4
    if not q < 1:
        return None
    gam_op = arb((ip4 / (1 - q)).upper())                          # >= sup |G|_op
    eps_op = arb((ip4 * q / (1 - q)).upper())                      # >= |G - psi^-4 I|_op
    eps = up(eps_op)                                               # entrywise radius (<= op norm)
    p4 = psi**4
    G = [[(1 / p4 if i == j else acb(0)) + acb(arb(0, eps), arb(0, eps)) for j in R3] for i in R3]
    g0 = [[[4 * psi**3 * dpsi[k] if i == j else acb(0) for j in R3] for i in R3] for k in R3]
    Dg0 = [[[[4 * (3 * psi**2 * dpsi[m] * dpsi[k] + psi**3 * d2[m][k]) if i == j else acb(0)
              for j in R3] for i in R3] for k in R3] for m in R3]
    G1 = frob_upper([Dg0[m][k][i][j] for m in R3 for k in R3 for i in R3 for j in R3])
    Ric = ricci0(psi, dpsi, d2)
    K0 = frob_upper([Ric[i][j] for i in R3 for j in R3])
    gball, Kball = acb(arb(0, radii[1]), arb(0, radii[1])), acb(arb(0, radii[2]), arb(0, radii[2]))
    flatG = [G[i][j] for i in R3 for j in R3]
    flatg = [g0[k][i][j] + gball for k in R3 for i in R3 for j in R3]
    sG, sg = seed(flatG, 0, 36), seed(flatg, 9, 36)
    _, B = ricci_parts([[sG[3 * i + j] for j in R3] for i in R3],
                       [[[sg[9 * k + 3 * i + j] for j in R3] for i in R3] for k in R3])
    JBG = frob_upper([B[i][j].g[c] for i in R3 for j in R3 for c in range(9)])
    JBg = frob_upper([B[i][j].g[c] for i in R3 for j in R3 for c in range(9, 36)])
    sG, sK = seed(flatG, 0, 18), seed([Kball] * 9, 9, 18)
    Q = quad_K([[sG[3 * i + j] for j in R3] for i in R3], [[sK[3 * i + j] for j in R3] for i in R3])
    JQG = frob_upper([Q[i][j].g[c] for i in R3 for j in R3 for c in range(9)])
    JQK = frob_upper([Q[i][j].g[c] for i in R3 for j in R3 for c in range(9, 18)])
    r3 = arb(3).sqrt()
    bil = _half_norm_T() if sharp else arb(3)        # |A(G, D)|_F <= bil |G|_F |D|_F
    if sharp:
        # A(G, D) = psi^-4 A(I, D) + A(E, D), |E|_F <= sqrt(3) eps_op; the sqrt(3) after the bracket
        # is the Cauchy estimate |grad f|_F <= sqrt(3) |f|_s / (s - s')
        aKg = r3 * (NORM_L * ip4 + bil * r3 * eps_op)
    else:
        aKg = 3 * r3 * (r3 * gam_op)                 # 3 sqrt(3) Gam_F, Gam_F = sqrt(3) Gam_op
    c = {"gamma_K": up(2 * S0),
         "g_K": up(2 * r3),
         "K_gamma": up(S0 * (arb(JBG) + arb(JQG)) * gam_op**2 + bil * gam_op**2 * (S0 * arb(G1) + r3 * Rg)),
         "K_g": up(aKg + S0 * arb(JBg)),
         "K_K": up(S0 * arb(JQK))}
    # the same bounds split into derivative-loss parts (form c / (s - s')) and zero-order parts
    # (pointwise Lipschitz at the output level), for Theorem CK' (notes/17 §2)
    loss = {"g_K": up(2 * r3),
            "K_gamma": up(bil * gam_op**2 * r3 * Rg),
            "K_g": up(aKg)}
    zero = {"gamma_K": 2.0,
            "K_gamma": up((arb(JBG) + arb(JQG)) * gam_op**2 + bil * gam_op**2 * arb(G1)),
            "K_g": JBg,
            "K_K": JQK}
    vals = list(c.values()) + list(loss.values()) + list(zero.values()) + [K0, G1, up(gam_op)]
    if not all(math.isfinite(v) for v in vals):        # a non-finite bound is a failed cell, never skipped
        return None
    return {"c": c, "loss": loss, "zero": zero, "K0": K0, "psi_abs_lower": down(apsi),
            "gam_op": up(gam_op), "G1": G1, "J": {"BG": JBG, "Bg": JBg, "QG": JQG, "QK": JQK}}


def lipschitz_C(c, radii):
    """C (rounded up) for the weighted norm |u| = max(|u_gamma|/R_g0, |u_g|/R_g1, |u_K|/R_K) (ball
    radius 1): C = max over output blocks A of sum_B (R_B / R_A) c_AB."""
    Rgam, Rg, RK = (arb(r) for r in radii)
    cc = {k: arb(v) for k, v in c.items()}
    return max(up(RK / Rgam * cc["gamma_K"]), up(RK / Rg * cc["g_K"]),
               up(Rgam / RK * cc["K_gamma"] + Rg / RK * cc["K_g"] + cc["K_K"]))


def _cell(i, cell):
    """The closed interval [i cell, (i + 1) cell] as an arb ball, with exact endpoints (adjacent cells
    share endpoints, so the union has no gaps)."""
    return (arb(i) * arb(cell)).union(arb(i + 1) * arb(cell))


def imaginary_boxes(s0, m):
    """Boxes covering the imaginary parts needed: by complex conjugation and the reflection y -> -y
    (both preserve every bound, which are norms of tensors built from real-analytic data), the sup over
    Omega_s0 is attained with Im x >= 0, Im y >= 0. The quadrant [0, s0]^2 x [-s0, s0] is split into
    m x m x 2m boxes; boxes with no point of norm < s0 (outside the Hermitian s0-ball) are dropped."""
    h = arb(s0) / m
    out = []
    for i in range(m):
        for j in range(m):
            for k in range(-m, m):
                bx = (arb(i) * h).union(arb(i + 1) * h)
                by = (arb(j) * h).union(arb(j + 1) * h)
                bz = (arb(k) * h).union(arb(k + 1) * h)
                mn = [arb(i) * h, arb(j) * h, arb(max(k, -k - 1, 0)) * h]   # min |coordinate|
                if (mn[0] * mn[0] + mn[1] * mn[1] + mn[2] * mn[2]) >= arb(s0) * arb(s0):
                    continue
                out.append((bx, by, bz))
    return out


def _cell_worker(args):
    i, k, cell, d_lo, d_hi, s0, radii, sharp, m = args
    with working_precision(PREC):
        d = arb(d_lo).union(arb(d_hi))
        res = None
        for im in (imaginary_boxes(s0, m) if m else [None]):
            b = cell_bounds(_cell(i, cell), _cell(k, cell), d, s0, radii, None, sharp, im)
            if b is None:
                return None
            if res is None:
                res = b
                continue
            for part in ("c", "loss", "zero"):
                for key in res[part]:
                    res[part][key] = max(res[part][key], b[part][key])
            res["K0"] = max(res["K0"], b["K0"])
            res["psi_abs_lower"] = min(res["psi_abs_lower"], b["psi_abs_lower"])
            res["gam_op"] = max(res["gam_op"], b["gam_op"])
        res.pop("J", None)
        return res


def time_bg(C, K0w, s0):
    """Baouendi-Goulaouic (12) (checked on the original page): a < R / (16 R C' + 8 M) with s0 = 1
    normalisation, C' = C / s0, M = sup |F(0)|_w = K0w, R = 1. Returns (a rounded down, B = 4 a M
    rounded up), B the bound sup |u(t)|_w of their proof."""
    a = down(1 / (16 * arb(C) / arb(s0) + 8 * arb(K0w)))
    return a, up(4 * arb(a) * arb(K0w))


def rows_CL(cl, cz, radii):
    """Weighted row constants: loss C_A = sum_B (R_B/R_A) c^loss_AB, zero order L_A likewise."""
    Rgam, Rg, RK = (arb(r) for r in radii)
    C = {"gamma": arb(0),
         "g": RK / Rg * arb(cl["g_K"]),
         "K": Rgam / RK * arb(cl["K_gamma"]) + Rg / RK * arb(cl["K_g"])}
    L = {"gamma": RK / Rgam * arb(cz["gamma_K"]),
         "g": arb(0),
         "K": Rgam / RK * arb(cz["K_gamma"]) + Rg / RK * arb(cz["K_g"]) + arb(cz["K_K"])}
    return {A: up(C[A]) for A in C}, {A: up(L[A]) for A in L}


def time_ck2(Crow, Lrow, K0w, s0, beta, theta):
    """Theorem CK' (notes/17 §2): with k = max_A (4 C_A / (beta s0) + L_A / (1 - beta)) and M = K0w,
    every a with 2 a k < 1 and 3 a M < (1 - beta)(1 - 2 a k) gives a holomorphic solution on
    |t| < a (1 - s/s0) with sup |u(t)|_w <= a M / ((1 - beta)(1 - a k)). Returns (a, B) for
    a = theta * a_max(beta) (rounded down), the conditions checked in ball arithmetic."""
    b, S0, M = arb(beta), arb(s0), arb(K0w)
    k = max((4 * arb(Crow[A]) / (b * S0) + arb(Lrow[A]) / (1 - b)) for A in Crow)
    k = arb(k.upper())
    amax = (1 - b) / (3 * M + 2 * (1 - b) * k)
    a = arb(down(arb(theta) * amax))
    assert 2 * a * k < 1 and 3 * a * M < (1 - b) * (1 - 2 * a * k), "CK' conditions"
    return float(a), up(a * M / ((1 - b) * (1 - a * k)))


def certify(s0, radii, d_lo, d_hi, w_tube=W_TUBE, cell=0.02, boxes=None, log=None, centres=None,
            sharp=False, im_split=0, workers=1):
    """Bounds over Omega_s0(U), U the profile boxes inflated by w_tube (rotated and reflected), d in
    [d_lo, d_hi] (exact hull): per-block coefficients maximised over cells (each separately), the
    combined C of note 15, the loss / zero-order split, and K0. Times:
      T_bg  = Baouendi-Goulaouic (12): 1 / (16 C / s0 + 8 K0 / R_K), with B = 4 T_bg K0 / R_K;
      schedule = Theorem CK' (notes/17 §2) pairs (beta, a, B) for several beta and fractions of a_max.
    T_min (the largest certified a) and sup_u_weighted_bound (its B) are kept for compatibility.
    im_split = 0: imaginary parts in the full cube [-s0, s0]^3 (session 11-12 covering); m >= 1: the
    quadrant covering of imaginary_boxes(s0, m). Returns ok = False with the reason if any check
    fails; there is no partial success."""
    t0 = time.time()
    with working_precision(PREC):
        if boxes is None:
            boxes = profile_boxes()
        cells = covering_cells(boxes, arb(w_tube) + arb(s0), cell)
        if workers > 1 and centres is None:
            import multiprocessing as mp
            import research.calc.ck_bound as _mod      # picklable by module path even when run as __main__
            with mp.get_context("fork").Pool(workers) as pool:
                results = pool.map(_mod._cell_worker, [(i, k, cell, d_lo, d_hi, s0, radii, sharp, im_split)
                                                       for (i, k) in cells], chunksize=8)
        else:
            results = []
            d = arb(d_lo).union(arb(d_hi))
            for (i, k) in cells:
                if im_split:
                    results.append(_cell_worker((i, k, cell, d_lo, d_hi, s0, radii, sharp, im_split))
                                   if centres is None else None)
                else:
                    results.append(cell_bounds(_cell(i, cell), _cell(k, cell), d, s0, radii, centres, sharp))
        agg = {"c": {}, "loss": {}, "zero": {}}
        K0, psi_low, gam_op, worst_cell = 0.0, math.inf, 0.0, None
        for (i, k), b in zip(cells, results):
            if b is None:
                return {"ok": False, "reason": f"check failed on cell {(i, k)}"}
            for part in agg:
                for key, v in b[part].items():
                    agg[part][key] = max(agg[part].get(key, 0.0), v)
            if b["K0"] > K0:
                K0, worst_cell = b["K0"], (i, k)
            psi_low = min(psi_low, b["psi_abs_lower"])
            gam_op = max(gam_op, b["gam_op"])
        cmax = agg["c"]
        C = lipschitz_C(cmax, radii)
        K0w = up(arb(K0) / arb(radii[2]))
        T_bg, B_bg = time_bg(C, K0w, s0)
        Crow, Lrow = rows_CL(agg["loss"], agg["zero"], radii)
        schedule = []
        for beta in (0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
            for theta in (0.25, 0.5, 0.7, 0.85, 0.95, 0.999):
                a, B = time_ck2(Crow, Lrow, K0w, s0, beta, theta)
                schedule.append({"beta": beta, "theta": theta, "a": a, "B": B})
        best = max(schedule, key=lambda e: e["a"])
        if log:
            log(f"  T_bg {T_bg:.4e}  CK' best a {best['a']:.4e} (beta {best['beta']})")
    return {"ok": True, "s0": s0, "radii": list(radii), "d_lo": d_lo, "d_hi": d_hi,
            "w_tube": w_tube, "cell": cell, "cells": len(cells), "profile_boxes": len(boxes),
            "sharp": sharp, "im_split": im_split, "c": cmax, "c_loss": agg["loss"], "c_zero": agg["zero"],
            "rows_C": Crow, "rows_L": Lrow, "C": C, "K0": K0, "K0w": K0w, "K0_cell": worst_cell,
            "psi_abs_lower": psi_low, "gam_op_upper": gam_op,
            "T_bg": T_bg, "B_bg": B_bg, "schedule": schedule,
            "T_min": best["a"], "sup_u_weighted_bound": best["B"], "theorem": "CK' (notes/17 §2)",
            "precision_bits": PREC, "rounding": "outward", "seconds": time.time() - t0}


# Parameters of the certificates written by the command line (notes/17 §3; chosen by the searches
# recorded there; any choice gives a valid certificate, the search only selects).
MODES = {
    # session 11 constants, Baouendi-Goulaouic (12) time corrected, plus the CK' schedule
    "classic": dict(s0=0.012, radii=(1.0, 20.0, 10.0), sharp=False, im_split=0, region="tube",
                    out="ck_bound.json"),
    # session 12 constants (exact principal-part norms), corrected time, plus the CK' schedule
    "sharp": dict(s0=0.014, radii=(0.25, 22.0, 10.0), sharp=True, im_split=0, region="tube",
                  out="ck_bound_sharp.json"),
    # session 13: CK', quadrant covering with sub-boxes, whole tube
    "ext": dict(s0=0.017, radii=(0.25, 22.0, 10.0), sharp=True, im_split=2, region="tube",
                out="ck_bound_ext.json"),
    # session 13: CK' on the waist region (profile boxes with z <= 0.15, near the symmetry plane)
    "waist": dict(s0=0.03, radii=(0.2, 18.0, 7.0), sharp=True, im_split=3, region="waist",
                  out="ck_bound_waist.json"),
}
WAIST_ZMAX = 0.15


if __name__ == "__main__":
    # python -m research.calc.ck_bound [classic|sharp|ext|waist]
    import platform
    import subprocess
    from research.calc import gap3_interval as gi
    mode = sys.argv[1] if len(sys.argv) > 1 else "ext"
    cfg = MODES[mode]
    boxes = profile_boxes()
    region = boxes if cfg["region"] == "tube" else [b for b in boxes if b[3] <= WAIST_ZMAX]
    res = certify(cfg["s0"], cfg["radii"], D_LO, gi.D_HOLE, cell=CELL, boxes=region, sharp=cfg["sharp"],
                  im_split=cfg["im_split"], workers=4, log=lambda m: print(m, flush=True))
    print(json.dumps({k: v for k, v in res.items() if k != "schedule"}, indent=1))
    if res["ok"]:
        res["profile"] = region
        res["region"] = cfg["region"] + (f" (profile boxes with z_hi <= {WAIST_ZMAX})" if cfg["region"] == "waist" else "")
        res["replay"] = "python -m research.calc.ck_bound " + mode
        try:
            res["source_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        except OSError:
            res["source_commit"] = None
        res["python"] = platform.python_version()
        with open(os.path.join(os.path.dirname(__file__), "data", cfg["out"]), "w") as f:
            json.dump(res, f, indent=1)
