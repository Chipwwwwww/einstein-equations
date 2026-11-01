"""An evaluable metric for the near-critical Brill-Lindquist development (category B, notes/16).

What is evaluated. The spacetime is the vacuum development of equal-mass, time-symmetric Brill-
Lindquist data,

    gamma_0 = psi^4 delta,   psi = 1 + 1/(4 r_+) + 1/(4 r_-),   r_pm = |x - (0, 0, +-d/2)|,   K = 0,

in units G = c = M_ADM = 1 (bare masses 1/2). The gauge is geodesic slicing: lapse 1, shift 0,
spatial coordinates comoving with the normal geodesics of the initial slice. So

    g = -dt^2 + gamma_ij(t, x) dx^i dx^j,    g_00 = -1, g_0i = 0 exactly,

with K_ij = (1/2) d_t gamma_ij (repository sign). The coordinates (t, x, y, z) are Cartesian on the
initial slice; they are regular on the symmetry axis, and they are the chart in which the
Cauchy-Kovalevskaya certificate (research/calc/ck_bound.py, notes/15) is proved.

How. gamma(t, x) = sum_k g_k(x) t^(2k) (even in t). The coefficients g_k(x) and their first and
second spatial derivatives at x are computed exactly in the chosen arithmetic (python-flint arb
balls, or exact rationals at points with rational distances to the punctures) by a truncated-jet
recursion of the evolution equations

    (2k + 1) kappa_k = [-Ric(gamma)]_k + [2 kappa gamma^-1 kappa - tr(gamma^-1 kappa) kappa]_(k-1),
    g_(k+1) = kappa_k / (k + 1),       K = t sum_k kappa_k t^(2k),

with jets in (X, Y, Z, s = t^2) of weight |alpha| + 2 deg_s <= W. There is no spatial approximation:
truncated-polynomial arithmetic gives the Taylor coefficients at the point exactly.

Error control. The infinite series converges on the certified domain (Proposition T, notes/15):
for x in U (the closed tube of radius 0.02 about the critical surface S_c), d in I and |t| < T_min,
with Taylor-coefficient bounds |u_n|_w <= B T_s^-n on the complex s-neighbourhood of U,
T_s = T_min (1 - s/s0), B = 4 T_min K0 / R_K. tail_bounds() turns these into explicit bounds for the
metric, its first and second derivatives (spatial ones by Cauchy estimates, with the spatial loss
s), and evaluate_metric() adds them to the ball enclosures of the truncated sums. Points, times or
separations outside the certified domain are rejected, or (allow_uncertified=True) returned with
status "uncertified" and no tail bound.

Separately reported: (1) coefficient arithmetic error (ball radii of the truncated sum); (2)
spatial approximation error (none); (3) time-truncation error (tail bound); (4) parameter and
imported-certificate uncertainty (d may be a ball; the certificate constants are imported from
research/calc/data/ck_bound_sharp.json or ck_bound.json, with their replay command).

Run: python -m research.calc.bl_metric    (the demonstration of notes/16 §4, writes
research/calc/data/bl_metric_demo.json)
"""

import json
import math
import os
import time
from contextlib import contextmanager

from flint import arb, arb_mat, ctx, fmpq

DATA = os.path.join(os.path.dirname(__file__), "data")
CERTS = ("ck_bound_waist.json", "ck_bound_ext.json", "ck_bound_sharp.json", "ck_bound.json")


@contextmanager
def working_precision(bits):
    old = ctx.prec
    ctx.prec = bits
    try:
        yield
    finally:
        ctx.prec = old


# ---------------------------------------------------------------------------
# Truncated jets in (X, Y, Z, s), generic coefficients (arb, fmpq, float)
# ---------------------------------------------------------------------------
class Jets:
    """Truncated polynomials sum c_m X^a Y^b Z^c s^k with weight a + b + c + sw k <= W, stored as
    dicts {monomial index: coefficient}; monomials are sorted by weight."""

    def __init__(self, W, sw=2):
        self.W, self.sw = W, sw
        mons = []
        for k in range(W // sw + 1):
            for deg in range(W - sw * k + 1):
                for a in range(deg, -1, -1):
                    for b in range(deg - a, -1, -1):
                        mons.append((a, b, deg - a - b, k))
        mons.sort(key=lambda m: (m[0] + m[1] + m[2] + sw * m[3], m))
        self.mons = mons
        self.index = {m: i for i, m in enumerate(mons)}
        self.wt = [m[0] + m[1] + m[2] + sw * m[3] for m in mons]
        self.cut = [sum(1 for w in self.wt if w <= r) for r in range(W + 1)]   # prefix by weight
        self.add = []
        for i, m in enumerate(mons):
            room = W - self.wt[i]
            self.add.append([self.index[(m[0] + n[0], m[1] + n[1], m[2] + n[2], m[3] + n[3])]
                             for n in mons[: self.cut[room]]])
        self.dmap = []
        for ax in range(3):
            dm = {}
            for i, m in enumerate(mons):
                if m[ax]:
                    n = list(m)
                    n[ax] -= 1
                    dm[i] = (self.index[tuple(n)], m[ax])
            self.dmap.append(dm)
        self.smap = {i: self.index[(m[0], m[1], m[2], m[3] + 1)] for i, m in enumerate(mons)
                     if self.wt[i] + sw <= W}

    # arithmetic -----------------------------------------------------------
    def const(self, c):
        return {0: c}

    def var(self, ax, c=1):
        m = [0, 0, 0, 0]
        m[ax] = 1
        return {self.index[tuple(m)]: c}

    @staticmethod
    def add_(a, b):
        out = dict(a)
        for k, v in b.items():
            out[k] = out[k] + v if k in out else v
        return out

    @staticmethod
    def sub(a, b):
        out = dict(a)
        for k, v in b.items():
            out[k] = out[k] - v if k in out else -v
        return out

    @staticmethod
    def scale(a, c):
        return {k: v * c for k, v in a.items()}

    def mul(self, a, b):
        if not a or not b:
            return {}
        if len(a) > len(b):
            a, b = b, a
        bi = sorted(b.items())
        out = {}
        add = self.add
        for i, ai in a.items():
            row = add[i]
            L = len(row)
            for j, bj in bi:
                if j >= L:
                    break
                k = row[j]
                if k in out:
                    out[k] = out[k] + ai * bj
                else:
                    out[k] = ai * bj
        return out

    def d(self, a, ax):
        """d/dX (ax 0), d/dY (1), d/dZ (2). The result is exact to weight W - 1 (it has no weight-W
        terms)."""
        dm = self.dmap[ax]
        out = {}
        for i, v in a.items():
            if i in dm:
                j, e = dm[i]
                out[j] = v * e
        return out

    def times_s(self, a):
        return {self.smap[i]: v for i, v in a.items() if i in self.smap}

    def s_slice(self, a, k):
        """The coefficient of s^k, as a jet in s^0 (spatial jet)."""
        idx = self.index
        return {idx[(m[0], m[1], m[2], 0)]: v for i, v in a.items()
                for m in (self.mons[i],) if m[3] == k}

    def put_s(self, a, k, spatial, max_weight):
        """Set the s^k slice of a (in place) from a spatial jet, keeping weights <= max_weight."""
        idx = self.index
        for i, v in spatial.items():
            m = self.mons[i]
            if m[3] == 0 and self.wt[i] + self.sw * k <= max_weight:
                a[idx[(m[0], m[1], m[2], k)]] = v

    def drop_above(self, a, w):
        return {i: v for i, v in a.items() if self.wt[i] <= w}

    def inv(self, a, one):
        """1/a, for a jet with an invertible constant term."""
        a0 = a[0]
        ia0 = one / a0
        x = self.scale({k: v for k, v in a.items() if k != 0}, -ia0)
        acc, term = {0: one}, {0: one}
        for _ in range(self.W):
            term = self.mul(term, x)
            if not term:
                break
            acc = self.add_(acc, term)
        return self.scale(acc, ia0)

    def coef(self, a, m, zero):
        return a.get(self.index[m], zero)


# ---------------------------------------------------------------------------
# Brill-Lindquist data as jets (Cartesian)
# ---------------------------------------------------------------------------
def binom_half(n):
    """binomial(-1/2, n) as an exact rational."""
    c = fmpq(1)
    for j in range(n):
        c = c * fmpq(-1 - 2 * j, 2 * (j + 1))
    return c


def psi_jet(J, x0, centres, one):
    """psi = 1 + sum c / |x - p| as a spatial jet at x0 = (x, y, z) (field elements), centres =
    [(c, (px, py, pz), r0), ...] with r0 = |x0 - p| supplied in the field (arb sqrt, or exact for
    rational-distance points). 1/|x - p| = r0^-1 (1 + e)^-1/2 with e = (2 (x0 - p).X + |X|^2) / r0^2."""
    psi = {0: one}
    X = [J.var(i, one) for i in range(3)]
    XX = J.add_(J.add_(J.mul(X[0], X[0]), J.mul(X[1], X[1])), J.mul(X[2], X[2]))
    for c, p, r0 in centres:
        a = [x0[i] - p[i] for i in range(3)]
        ir2 = one / (r0 * r0)
        lin = {}
        for i in range(3):
            lin = J.add_(lin, J.scale(X[i], 2 * a[i]))
        e = J.scale(J.add_(lin, XX), ir2)
        acc, term = {0: one}, {0: one}
        for n in range(1, J.W + 1):
            term = J.mul(term, e)
            acc = J.add_(acc, J.scale(term, one * binom_half(n)))
        psi = J.add_(psi, J.scale(acc, c / r0))
    return psi


def bl_centres(d, one, x0=None, exact_r=None):
    """The two punctures (bare masses 1/2, coefficient 1/4) at z = +-d/2. r0 = |x0 - p| is computed
    by arb sqrt unless exact_r supplies it."""
    pts = [(0 * one, 0 * one, d / 2), (0 * one, 0 * one, -d / 2)]
    out = []
    for n, p in enumerate(pts):
        if exact_r is not None:
            r0 = exact_r[n]
        else:
            q = sum(((x0[i] - p[i]) * (x0[i] - p[i]) for i in range(3)), 0 * one)   # not ** 2 (nan at mid 0)
            r0 = q.sqrt()
        out.append((one / 4, p, r0))
    return out


PAIRS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def sym_get(S, i, j):
    return S[(i, j) if i <= j else (j, i)]


def inverse3(J, g, one):
    """Inverse of a symmetric 3x3 matrix of jets (adjugate over determinant)."""
    m = lambda a, b: J.mul(a, b)  # noqa: E731
    G = lambda i, j: sym_get(g, i, j)  # noqa: E731
    cof = {}
    for (i, j) in PAIRS:
        r = [k for k in range(3) if k != j]
        c = [k for k in range(3) if k != i]
        minor = J.sub(m(G(r[0], c[0]), G(r[1], c[1])), m(G(r[0], c[1]), G(r[1], c[0])))
        cof[(i, j)] = minor if (i + j) % 2 == 0 else J.scale(minor, -1)
    det = {}
    for j in range(3):
        det = J.add_(det, m(G(0, j), cof[(0, j)] if j >= 0 else None))
    idet = J.inv(det, one)
    return {pq: J.mul(cof[pq], idet) for pq in PAIRS}


def ricci3(J, g, one):
    """Ric_bd = d_a Gam^a_bd - d_d Gam^a_ab + Gam^a_ae Gam^e_bd - Gam^a_de Gam^e_ab (repository
    convention), for a symmetric 3x3 matrix of jets in Cartesian coordinates. Exact to weight W - 2."""
    Gi = inverse3(J, g, one)
    dg = {(pq, k): J.d(g[pq], k) for pq in PAIRS for k in range(3)}
    D = lambda i, j, k: dg[((i, j) if i <= j else (j, i), k)]  # noqa: E731  d_k g_ij
    low = {}   # Gam_{e, bc} = (1/2)(d_c g_eb + d_b g_ec - d_e g_bc)
    for e in range(3):
        for (b, c) in PAIRS:
            low[(e, b, c)] = J.scale(J.sub(J.add_(D(e, b, c), D(e, c, b)), D(b, c, e)), one / 2)
    Gam = {}
    for a in range(3):
        for (b, c) in PAIRS:
            acc = {}
            for e in range(3):
                acc = J.add_(acc, J.mul(sym_get(Gi, a, e), low[(e, b, c)]))
            Gam[(a, b, c)] = acc
    gm = lambda a, b, c: Gam[(a, b, c) if b <= c else (a, c, b)]  # noqa: E731
    tr = [J.add_(J.add_(gm(0, 0, e), gm(1, 1, e)), gm(2, 2, e)) for e in range(3)]  # Gam^a_ae
    out = {}
    for (b, dd) in PAIRS:
        acc = {}
        for a in range(3):
            acc = J.add_(acc, J.d(gm(a, b, dd), a))
        acc = J.sub(acc, J.d(tr[b], dd))
        for e in range(3):
            acc = J.add_(acc, J.mul(tr[e], gm(e, b, dd)))
            for a in range(3):
                acc = J.sub(acc, J.mul(gm(a, dd, e), gm(e, a, b)))
        out[(b, dd)] = acc
    return out


def matmul_sym(J, A, B, C):
    """A B C for symmetric matrices given as dicts over PAIRS (the result is symmetric when A = C)."""
    AB = {(i, j): {} for i in range(3) for j in range(3)}
    for i in range(3):
        for j in range(3):
            acc = {}
            for k in range(3):
                acc = J.add_(acc, J.mul(sym_get(A, i, k), sym_get(B, k, j)))
            AB[(i, j)] = acc
    out = {}
    for (i, j) in PAIRS:
        acc = {}
        for k in range(3):
            acc = J.add_(acc, J.mul(AB[(i, k)], sym_get(C, k, j)))
        out[(i, j)] = acc
    return out


def trace_prod(J, A, B):
    acc = {}
    for i in range(3):
        for j in range(3):
            acc = J.add_(acc, J.mul(sym_get(A, i, j), sym_get(B, j, i)))
    return acc


_JET_CACHE = {}
_JET_CACHE_MAX = 64


def series_jets_cached(x0, d, N, one, prec):
    """series_jets for arb inputs, memoised on (x0, d, N, precision): the jets do not depend on t, so a
    scan over times reuses them."""
    key = (tuple(str(v) for v in x0), str(d), N, prec)
    if key not in _JET_CACHE:
        if len(_JET_CACHE) >= _JET_CACHE_MAX:
            _JET_CACHE.pop(next(iter(_JET_CACHE)))
        _JET_CACHE[key] = series_jets(x0, d, N, one)
    return _JET_CACHE[key]


def series_jets(x0, d, N, one, exact_r=None, centres=None, W=None):
    """The geodesic-slicing development as jets at x0: returns (J, g, kap) where g[(i,j)] is the
    (X, Y, Z, s)-jet of gamma_ij = sum_k g_k s^k and kap that of kappa (K = t kappa). With
    W = 2 N (default), g_k is exact with its spatial derivatives up to order 2N - 2k: values for
    k <= N, first and second derivatives for k <= N - 1."""
    W = 2 * N if W is None else W
    J = Jets(W)
    cs = centres if centres is not None else bl_centres(d, one, x0, exact_r)
    psi = psi_jet(J, x0, cs, one)
    p2 = J.mul(psi, psi)
    p4 = J.mul(p2, p2)
    g = {pq: (dict(p4) if pq[0] == pq[1] else {}) for pq in PAIRS}
    kap = {pq: {} for pq in PAIRS}
    for k in range(J.W // 2):
        Ric = ricci3(J, g, one)
        Gi = inverse3(J, g, one)
        quad_a = matmul_sym(J, kap, Gi, kap)
        trK = trace_prod(J, Gi, kap)
        for pq in PAIRS:
            quad = J.sub(J.scale(quad_a[pq], 2), J.mul(trK, kap[pq]))
            rhs = J.add_(J.scale(Ric[pq], -1), J.times_s(quad))
            sl = J.s_slice(rhs, k)                        # exact to spatial weight W - 2 - 2k
            sl = J.drop_above(sl, J.W - 2 - 2 * k)
            kk = J.scale(sl, one / (2 * k + 1))
            J.put_s(kap[pq], k, kk, J.W - 2)
            J.put_s(g[pq], k + 1, J.scale(kk, one / (k + 1)), J.W)
    return J, g, kap


def coefficients(J, g, zero, N):
    """From the jets: lists over k of {pair: value}, {pair: [d_x, d_y, d_z]}, {pair: 3x3 second
    derivatives}. Values for k <= N, derivatives for k <= N - 1."""
    vals, d1, d2 = [], [], []
    for k in range(N + 1):
        vals.append({pq: J.coef(g[pq], (0, 0, 0, k), zero) for pq in PAIRS})
        if k <= N - 1:
            e = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            d1.append({pq: [J.coef(g[pq], e[a] + (k,), zero) for a in range(3)] for pq in PAIRS})
            sec = {}
            for pq in PAIRS:
                M = [[None] * 3 for _ in range(3)]
                for a in range(3):
                    for b in range(3):
                        m = [0, 0, 0]
                        m[a] += 1
                        m[b] += 1
                        c = J.coef(g[pq], tuple(m) + (k,), zero)
                        M[a][b] = c * 2 if a == b else c
                sec[pq] = M
            d2.append(sec)
    return vals, d1, d2


def constraint_residuals(J, g, kap, one):
    """Hamiltonian R + K^2 - K_ij K^ij and momentum D_j K^j_i - D_i K constraints of the jet solution
    at the point, order by order (coefficients of s^k, resp. t s^k), for the orders the jets
    justify. They vanish by constraint propagation; nonzero values would expose a code error."""
    Ric = ricci3(J, g, one)
    Gi = inverse3(J, g, one)
    R = trace_prod(J, Gi, Ric)
    trK = trace_prod(J, Gi, kap)
    KK = trace_prod(J, matmul_sym(J, Gi, kap, Gi), kap)
    ham = J.add_(R, J.times_s(J.sub(J.mul(trK, trK), KK)))
    zero = 0 * one
    hres = [J.coef(ham, (0, 0, 0, k), zero) for k in range(J.W // 2)]
    # momentum: D_j K^j_i - D_i K = G^jk (D_k kap_ji - D_i kap_jk)
    dg = {(pq, k): J.d(g[pq], k) for pq in PAIRS for k in range(3)}
    D = lambda i, j, k: dg[((i, j) if i <= j else (j, i), k)]  # noqa: E731
    Gam = {}
    for a in range(3):
        for b in range(3):
            for c in range(3):
                acc = {}
                for e in range(3):
                    t_ = J.sub(J.add_(D(e, b, c), D(e, c, b)), D(b, c, e))
                    acc = J.add_(acc, J.mul(sym_get(Gi, a, e), t_))
                Gam[(a, b, c)] = J.scale(acc, one / 2)

    def DK(k, l, j):  # D_k kap_lj
        acc = J.d(sym_get(kap, l, j), k)
        for m in range(3):
            acc = J.sub(acc, J.add_(J.mul(Gam[(m, k, l)], sym_get(kap, m, j)),
                                    J.mul(Gam[(m, k, j)], sym_get(kap, l, m))))
        return acc
    mom = []
    for i in range(3):
        acc = {}
        for j in range(3):
            for k in range(3):
                acc = J.add_(acc, J.mul(sym_get(Gi, j, k), J.sub(DK(k, j, i), DK(i, j, k))))
        mom.append([J.coef(acc, (0, 0, 0, k), zero) for k in range((J.W - 3) // 2 + 1)])
    return hres, mom


class JV:
    """A jet with arithmetic operators, so that generic code (ck_bound.ricci_parts, quad_K) runs on
    jets."""
    __slots__ = ("J", "a")

    def __init__(self, J, a):
        self.J, self.a = J, a

    def _lift(self, b):
        return b.a if isinstance(b, JV) else ({0: b} if b != 0 else {})

    def __add__(self, b):
        return JV(self.J, self.J.add_(self.a, self._lift(b)))

    __radd__ = __add__

    def __sub__(self, b):
        return JV(self.J, self.J.sub(self.a, self._lift(b)))

    def __rsub__(self, b):
        return JV(self.J, self.J.sub(self._lift(b), self.a))

    def __neg__(self):
        return JV(self.J, self.J.scale(self.a, -1))

    def __mul__(self, b):
        if isinstance(b, JV):
            return JV(self.J, self.J.mul(self.a, b.a))
        return JV(self.J, self.J.scale(self.a, b) if b != 0 else {})

    __rmul__ = __mul__


def first_order_coefficients(x0, d, W, one):
    """Independent cross-check: Taylor coefficients in t of the first-order system that the
    Cauchy-Kovalevskaya certificate is proved for (ck_bound: d_t gamma = 2K, d_t g = 2 dK,
    d_t K = -(A + B)(G, g, dg) + Q(G, K), with ck_bound's own ricci_parts and quad_K), computed by
    u_(n+1) = [F(u)]_n / (n + 1) on jets in (X, Y, Z, t) (t of weight 1), with no parity assumption.
    Returns lists over n of {pair: gamma_n(x0)} and {pair: K_n(x0)}, exact for n <= W - 1."""
    from research.calc import ck_bound as ck
    J = Jets(W, sw=1)
    psi = psi_jet(J, x0, bl_centres(d, one, x0), one)
    p4 = J.mul(J.mul(psi, psi), J.mul(psi, psi))
    R = range(3)
    gam = [[JV(J, dict(p4) if i == j else {}) for j in R] for i in R]
    gg = [[[JV(J, J.d(gam[i][j].a, k)) for j in R] for i in R] for k in R]
    K = [[JV(J, {}) for j in R] for i in R]

    def slice_t(v, n):
        return {i: c for i, c in v.items() if J.mons[i][3] == n}

    def shift_up(v, n):   # spatial jet of the t^n slice -> t^(n+1) slice
        out = {}
        for i, c in v.items():
            m = J.mons[i]
            key = (m[0], m[1], m[2], n + 1)
            if key in J.index:
                out[J.index[key]] = c
        return out
    for n in range(W):
        Gs = inverse3(J, {pq: gam[pq[0]][pq[1]].a for pq in PAIRS}, one)
        G = [[JV(J, sym_get(Gs, i, j)) for j in R] for i in R]
        Dg = [[[[JV(J, J.d(gg[k][i][j].a, m)) for j in R] for i in R] for k in R] for m in R]
        A, B = ck.ricci_parts(G, gg, Dg)
        Q = ck.quad_K(G, K)
        FK = [[-(A[i][j] + B[i][j]) + Q[i][j] for j in R] for i in R]
        inv = one / (n + 1)
        for i in R:
            for j in R:
                gam[i][j] = gam[i][j] + JV(J, J.scale(shift_up(slice_t(J.scale(K[i][j].a, 2), n), n), inv))
                newK = J.scale(shift_up(slice_t(FK[i][j].a, n), n), inv)
                K[i][j] = K[i][j] + JV(J, newK)
                for k in R:
                    gg[k][i][j] = gg[k][i][j] + JV(J, J.scale(shift_up(slice_t(J.scale(J.d(K[i][j].a, k), 2), n), n), inv))
    zero = 0 * one
    gvals = [{pq: J.coef(gam[pq[0]][pq[1]].a, (0, 0, 0, n), zero) for pq in PAIRS} for n in range(W)]
    kvals = [{pq: J.coef(K[pq[0]][pq[1]].a, (0, 0, 0, n), zero) for pq in PAIRS} for n in range(W)]
    return gvals, kvals


# ---------------------------------------------------------------------------
# The certified domain and the tail bounds
# ---------------------------------------------------------------------------
def load_certificates():
    """Every valid certificate in research/calc/data (each certifies its own compact region U_c, its
    d-interval and its times). The evaluator combines all that apply to a point (notes/17 §4)."""
    out = []
    for fn in CERTS:
        path = os.path.join(DATA, fn)
        if os.path.exists(path):
            with open(path) as f:
                c = json.load(f)
            if c.get("ok"):
                c["file"] = fn
                out.append(c)
    if not out:
        raise FileNotFoundError("no valid ck_bound certificate in research/calc/data")
    return out


def load_certificate(name=None):
    """Constants of one certificate (imported, with their replay command); by default the first of
    CERTS that exists."""
    for fn in ([name] if name else CERTS):
        path = os.path.join(DATA, fn)
        if os.path.exists(path):
            with open(path) as f:
                c = json.load(f)
            if not c.get("ok"):
                continue
            c["file"] = fn
            return c
    raise FileNotFoundError("no valid ck_bound certificate in research/calc/data")


def in_domain(cert, x, y, z):
    """True if every point of the (ball) point (x, y, z) lies in the certified spatial domain U of the
    certificate: the points whose meridian coordinates (rho, |z|) lie in a profile box inflated by
    w_tube componentwise. This U is compact, contains the closed tube of radius w_tube about S_c
    (each profile box contains a point of S_c's meridian curve), and its complex s0-neighbourhood is
    exactly what the certificate's cells cover (notes/16 §1). Also returns the margin: the largest
    componentwise inflation still containing the point, minus w_tube (<= 0 when inside)."""
    rr = x * x + y * y
    r_hi = arb(rr.upper()).sqrt().upper()                 # rho <= r_hi
    r_lo = arb(rr.lower()).sqrt().lower() if rr.lower() > 0 else arb(0)
    az = abs(z)
    z_lo, z_hi = az.lower(), az.upper()
    best = math.inf
    for b in cert["profile"]:
        need = max((arb(b[0]) - r_lo).upper(), (r_hi - arb(b[1])).upper(),
                   (arb(b[2]) - z_lo).upper(), (z_hi - arb(b[3])).upper(), arb(0))
        best = min(best, math.nextafter(float(need), math.inf))
    return best <= cert["w_tube"], best


def _geom_even(q, n0):
    """sum_{n >= n0, n = n0 mod 2} q^n = q^n0 / (1 - q^2)."""
    return q ** n0 / (1 - q * q)


def _dgeom_odd(q, j):
    """sum_{i >= j} (2i + 1) q^(2i) = d/dq [q^(2j+1) / (1 - q^2)]."""
    return ((2 * j + 1) * q ** (2 * j) * (1 - q * q) + 2 * q ** (2 * j + 2)) / (1 - q * q) ** 2


def _dgeom_even(q, N):
    """sum_{j >= N} 2j q^(2j-1) = d/dq [q^(2N) / (1 - q^2)]."""
    return (2 * N * q ** (2 * N - 1) * (1 - q * q) + 2 * q ** (2 * N + 1)) / (1 - q * q) ** 2


def _tails_one(T, B, cert, ta, N):
    """Tail bounds for one certified pair (a = T, B): sup_{Omega_s(U)} |u(t)|_w <= B for complex
    |t| < T (1 - s/s0). See tail_bounds."""
    Rgam, Rg, RK = (arb(r) for r in cert["radii"])
    s0 = arb(cert["s0"])
    q = ta / T
    out = {
        # gamma: even n >= 2N + 2 (gamma block)
        "value": Rgam * B * _geom_even(q, 2 * N + 2),
        # d_x gamma = the g block (no Cauchy loss): even n >= 2N
        "d_space": Rg * B * _geom_even(q, 2 * N),
        # d_t gamma = 2 K: odd n >= 2N + 1
        "d_time": 2 * RK * B * _geom_even(q, 2 * N + 1),
        # d_t^2 gamma = 2 sum n K_n t^(n-1), odd n >= 2N + 1 (Cauchy in t on the K block)
        "d2_time": 2 * RK * B / T * _dgeom_odd(q, N),
        # d_t d_x gamma = d_t g = sum n g_n t^(n-1), even n >= 2N (Cauchy in t on the g block)
        "d2_time_space": Rg * B / T * _dgeom_even(q, N),
    }
    best_ss, best_ts = None, None
    qf = float(q.upper())
    for f in (1 / 64, 1 / 32, 1 / 16, 1 / 8, 3 / 16, 1 / 4, 3 / 8, 1 / 2, 5 / 8, 3 / 4, (1 - qf) / 2):
        if not 0 < f < 1:
            continue
        s = s0 * arb(f)
        Ts = T * (1 - arb(f))
        if not ta < Ts:
            continue
        qs = ta / Ts
        # d_x d_x gamma: Cauchy in x (radius s, per direction) on the g block, even n >= 2N
        ss = arb(3).sqrt() * Rg * B / s * _geom_even(qs, 2 * N)
        # d_t d_x gamma = 2 d_x K: Cauchy in x on the K block, odd n >= 2N - 1
        ts = 2 * arb(3).sqrt() * RK * B / s * _geom_even(qs, 2 * N - 1)
        if best_ss is None or ss.upper() < best_ss.upper():
            best_ss = ss
        if best_ts is None or ts.upper() < best_ts.upper():
            best_ts = ts
    out["d2_space"] = best_ss if best_ss is not None else arb("inf")
    if best_ts is not None and best_ts.upper() < out["d2_time_space"].upper():
        out["d2_time_space"] = best_ts
    return out


def certified_pairs(cert):
    """The certified (a, B) pairs: every pair of the CK' schedule (notes/17), or the single
    (T_min, B) of older certificates."""
    if "schedule" in cert:
        return [(e["a"], e["B"]) for e in cert["schedule"]]
    return [(cert["T_min"], cert["sup_u_weighted_bound"])]


def tail_bounds(cert, t_abs, N):
    """Upper bounds (Frobenius norms of the coordinate components) of the series tails left out by
    evaluate_metric at order N (metric values from g_0..g_N, derivatives from g_0..g_(N-1)), for
    |t| <= t_abs. Each certified pair (a, B) gives |u_n|_w <= B T_s^-n on Omega_s(U), T_s = a (1 - s/s0)
    (notes/16 §3, notes/17 §3); spatial derivatives use Cauchy estimates in x with spatial loss s
    (optimised on a grid), time derivatives Cauchy estimates in t, and every quantity takes the
    minimum over all valid pairs and methods (each bound is valid by itself). Returns None if
    t_abs is not below any certified a."""
    ta = arb(t_abs)
    best = None
    for a, B in certified_pairs(cert):
        T = arb(a)
        if not ta < T:
            continue
        one = _tails_one(T, arb(B), cert, ta, N)
        if best is None:
            best = one
        else:
            best = {k: (one[k] if one[k].upper() < best[k].upper() else best[k]) for k in best}
    if best is None:
        return None
    return {k: math.nextafter(float(v.upper()), math.inf) for k, v in best.items()}


# ---------------------------------------------------------------------------
# The evaluator
# ---------------------------------------------------------------------------
class DomainError(ValueError):
    pass


def _ball(v):
    return v if isinstance(v, arb) else arb(v)


def evaluate_metric(t, x, y, z, d, tolerance=None, precision=128, order=None, max_order=6,
                    allow_uncertified=False, certificate=None, tolerance_on="metric"):
    """The metric g_{mu nu} at (t, x, y, z) (comoving geodesic-slicing Cartesian coordinates, units
    M_ADM = 1) of the BL development with separation d, with certified enclosures.

    Arguments may be floats, strings (parsed exactly by arb) or arb balls (a ball means "every value
    in it": the enclosures then hold uniformly, e.g. over a spatial patch or a time interval).
    order N: metric values from g_0..g_N (through t^(2N)), derivatives from g_0..g_(N-1). With
    tolerance and no order, N increases from 1 until every radius of the quantities named by
    tolerance_on ("metric": g and g_inv; "first": also dg; "all": also ddg) is <= tolerance; if
    max_order is reached first, DomainError explains which bound is too large. Outside the
    certified domain (d not in the certificate's interval, x not in U, |t| >= T_min) a DomainError
    is raised, unless allow_uncertified=True: then the truncated series is returned with status
    "uncertified" and no tail bound.

    Returns a dict:
      status          "certified" or "uncertified" (with "reasons")
      g, g_inv        4x4 lists of arb balls (g_00 = -1 and g_0i = 0 are exact)
      dg              dg[l][m][n] = d_l g_mn (l = 0 is d_t)
      ddg             ddg[k][l][m][n] = d_k d_l g_mn
      error_budget    coefficient arithmetic (max ball radius before the tail is added), spatial
                      approximation (0), time-truncation bounds per quantity, parameter (radius of d)
      order, precision, domain, certificate, seconds
    """
    t0 = time.time()
    cert = certificate or load_certificates()
    if order is None and tolerance is None:
        order = 3
    orders = [order] if order is not None else list(range(1, max_order + 1))
    last_err = None
    for N in orders:
        res = _evaluate(t, x, y, z, d, N, precision, allow_uncertified, cert)
        res["seconds"] = time.time() - t0
        if tolerance is None:
            return res
        worst = res["max_radius_by"][tolerance_on]
        if worst <= tolerance:
            return res
        last_err = (N, worst, res["error_budget"])
    N, worst, budget = last_err
    raise DomainError(f"tolerance {tolerance:g} not reached up to order {N}: largest returned radius "
                      f"{worst:.3g}; error budget {budget}")


def _evaluate(t, x, y, z, d, N, precision, allow_uncertified, certs):
    if isinstance(certs, dict):
        certs = [certs]
    with working_precision(precision):
        t, x, y, z, d = (_ball(v) for v in (t, x, y, z, d))
        t_abs = math.nextafter(float(abs(t).upper()), math.inf)
        reasons, tails, used, dist = [], None, [], math.inf
        for cert in certs:
            why = []
            if not (d.lower() >= arb(cert["d_lo"]) and d.upper() <= arb(cert["d_hi"])):
                why.append(f"d outside [{cert['d_lo']}, {cert['d_hi']}]")
            inside, dc = in_domain(cert, x, y, z)
            dist = min(dist, dc)
            if not inside:
                why.append(f"point not inside U (needs inflation {dc:.4g} > {cert['w_tube']})")
            tc = tail_bounds(cert, t_abs, N) if not why else None
            if not why and tc is None:
                why.append(f"|t| = {t_abs:.4g} not below the certified time {cert['T_min']:.6g}")
            if why:
                reasons.append(f"{cert.get('file', 'certificate')}: " + "; ".join(why))
                continue
            used.append(cert)
            tails = tc if tails is None else {k: min(tails[k], tc[k]) for k in tails}
        if used:
            reasons = []
            cert = used[0]
        else:
            cert = certs[0]
        if reasons and not allow_uncertified:
            raise DomainError(" | ".join(reasons))
        one = arb(1)
        J, g, kap = series_jets_cached((x, y, z), d, N, one, precision)
        zero = arb(0)
        vals, d1, d2 = coefficients(J, g, zero, N)
        if any(v.rad() > 0 for v in (x, y, z)):
            # a box of points: mean-value form around the centre (Taylor coefficients at the centre
            # plus derivative enclosures over the box times the half-widths), intersected with the
            # direct enclosure. Both are valid; the mean-value form avoids most of the dependency
            # overestimation of the recursion on a box.
            xc = [arb(v.mid()) for v in (x, y, z)]
            h = [arb(0, float(arb(v.rad()).upper())) for v in (x, y, z)]
            Jc, gc, _ = series_jets_cached(tuple(xc), d, N, one, precision)
            vc, d1c, _ = coefficients(Jc, gc, zero, N)
            for k in range(N):
                for pq in PAIRS:
                    mv = vc[k][pq] + sum((d1[k][pq][a] * h[a] for a in range(3)), zero)
                    vals[k][pq] = vals[k][pq].intersection(mv) if vals[k][pq].overlaps(mv) else mv
                    for a in range(3):
                        mv1 = d1c[k][pq][a] + sum((d2[k][pq][a][b] * h[b] for b in range(3)), zero)
                        d1[k][pq][a] = d1[k][pq][a].intersection(mv1) if d1[k][pq][a].overlaps(mv1) else mv1
        # time sums ----------------------------------------------------------
        tp = [one]
        for _ in range(2 * N + 2):
            tp.append(tp[-1] * t)
        gam, dt, dtt, dx, dtx, dxx = {}, {}, {}, {}, {}, {}
        for pq in PAIRS:
            gam[pq] = sum((vals[k][pq] * tp[2 * k] for k in range(N + 1)), zero)
            dt[pq] = sum((2 * k * vals[k][pq] * tp[2 * k - 1] for k in range(1, N + 1)), zero)
            dtt[pq] = sum((2 * k * (2 * k - 1) * vals[k][pq] * tp[2 * k - 2] for k in range(1, N + 1)), zero)
            dx[pq] = [sum((d1[k][pq][a] * tp[2 * k] for k in range(N)), zero) for a in range(3)]
            dtx[pq] = [sum((2 * k * d1[k][pq][a] * tp[2 * k - 1] for k in range(1, N)), zero) for a in range(3)]
            dxx[pq] = [[sum((d2[k][pq][a][b] * tp[2 * k] for k in range(N)), zero) for b in range(3)]
                       for a in range(3)]
        coef_rad = 0.0
        for pq in PAIRS:
            for v in [gam[pq], dt[pq], dtt[pq]] + dx[pq] + dtx[pq] + [e for r in dxx[pq] for e in r]:
                coef_rad = max(coef_rad, float(v.rad()))
        if tails is not None:
            add = lambda v, e: v + arb(0, e)  # noqa: E731
            for pq in PAIRS:
                gam[pq] = add(gam[pq], tails["value"])
                dt[pq] = add(dt[pq], tails["d_time"])
                dtt[pq] = add(dtt[pq], tails["d2_time"])
                dx[pq] = [add(v, tails["d_space"]) for v in dx[pq]]
                dtx[pq] = [add(v, tails["d2_time_space"]) for v in dtx[pq]]
                dxx[pq] = [[add(v, tails["d2_space"]) for v in r] for r in dxx[pq]]
        # 4-metric ----------------------------------------------------------
        Z = arb(0)
        g4 = [[arb(-1), Z, Z, Z]] + [[Z] + [sym_get(gam, i, j) for j in range(3)] for i in range(3)]
        sp = arb_mat([[sym_get(gam, i, j) for j in range(3)] for i in range(3)])
        lor = _positive_definite(sp)
        try:
            spi = sp.inv()
            g4i = [[arb(-1), Z, Z, Z]] + [[Z] + [spi[i, j] for j in range(3)] for i in range(3)]
        except ZeroDivisionError:
            spi, g4i = None, None
        if tails is not None and spi is not None:
            # a second certified enclosure of d_t^2 gamma from the field equation; keep the intersection
            Gi = [[spi[i, j] for j in range(3)] for i in range(3)]
            eq = ddt_from_equation(gam, dx, dxx, dt, Gi)
            for pq in PAIRS:
                if eq[pq].overlaps(dtt[pq]):
                    dtt[pq] = dtt[pq].intersection(eq[pq])
                else:      # disjoint certified enclosures would be a bug, never a result
                    raise AssertionError(f"inconsistent enclosures of d_t^2 gamma_{pq}")
        dg = [[[Z] * 4 for _ in range(4)] for _ in range(4)]
        ddg = [[[[Z] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
        for i in range(3):
            for j in range(3):
                pq = (i, j) if i <= j else (j, i)
                dg[0][i + 1][j + 1] = dt[pq]
                ddg[0][0][i + 1][j + 1] = dtt[pq]
                for a in range(3):
                    dg[a + 1][i + 1][j + 1] = dx[pq][a]
                    ddg[0][a + 1][i + 1][j + 1] = ddg[a + 1][0][i + 1][j + 1] = dtx[pq][a]
                    for b in range(3):
                        ddg[a + 1][b + 1][i + 1][j + 1] = dxx[pq][a][b]
        rad = lambda vs: max(float(v.rad()) for v in vs)  # noqa: E731
        r_metric = rad([e for r in g4 for e in r] + ([e for r in g4i for e in r] if g4i else []))
        r_first = max(r_metric, rad([e for a in dg for r in a for e in r]))
        r_all = max(r_first, rad([e for a in ddg for b in a for r in b for e in r]))
        max_rad = r_all
        return {
            "status": "certified" if not reasons else "uncertified",
            "reasons": reasons,
            "g": g4, "g_inv": g4i, "dg": dg, "ddg": ddg,
            "lorentzian": lor if not reasons else None,
            "order": N, "t_power": 2 * N, "precision": precision,
            "error_budget": {
                "coefficient_arithmetic": coef_rad,
                "spatial_approximation": 0.0,
                "time_truncation": tails,
                "parameter_d_radius": float(d.rad()),
                "certificate": [c.get("file", "?") for c in used],
            },
            "max_radius": max_rad,
            "max_radius_by": {"metric": r_metric, "first": r_first, "all": r_all},
            "domain": {"inflation_needed": dist, "w_tube": cert["w_tube"], "t_abs": t_abs,
                       "certified_times": {c.get("file", "?"): c["T_min"] for c in used},
                       "d_interval": [cert["d_lo"], cert["d_hi"]]},
            "coordinates": "t = proper time along the normal geodesics of the initial slice; "
                           "(x, y, z) Cartesian comoving labels; punctures at (0, 0, +-d/2); G = c = M_ADM = 1",
        }


def ricci_pointwise(gam, dx, dxx, Gi):
    """Ric_ij at a point from the values gamma_ij, d_a gamma_ij, d_a d_b gamma_ij (arb balls) and an
    enclosure Gi of gamma^-1: the repository convention R_bd = d_a Gam^a_bd - d_d Gam^a_ab +
    Gam^a_ae Gam^e_bd - Gam^a_de Gam^e_ab, with d_a G = -G (d_a gamma) G."""
    R = range(3)
    D = lambda i, j, k: sym_get(dx, i, j)[k]  # noqa: E731
    DD = lambda i, j, k, l: sym_get(dxx, i, j)[k][l]  # noqa: E731
    low = [[[(D(e, b, c) + D(e, c, b) - D(b, c, e)) / 2 for c in R] for b in R] for e in R]
    dlow = [[[[(DD(e, b, c, a) + DD(e, c, b, a) - DD(b, c, e, a)) / 2 for a in R] for c in R]
             for b in R] for e in R]
    dG = [[[-sum((Gi[k][m] * D(m, n, a) * Gi[n][l] for m in R for n in R), arb(0)) for a in R]
           for l in R] for k in R]
    Gam = [[[sum((Gi[a][e] * low[e][b][c] for e in R), arb(0)) for c in R] for b in R] for a in R]
    dGam = [[[[sum((dG[a][e][m] * low[e][b][c] + Gi[a][e] * dlow[e][b][c][m] for e in R), arb(0))
               for m in R] for c in R] for b in R] for a in R]
    out = {}
    for (b, dd) in PAIRS:
        acc = sum((dGam[a][b][dd][a] for a in R), arb(0)) - sum((dGam[a][a][b][dd] for a in R), arb(0))
        for a in R:
            for e in R:
                acc += Gam[a][a][e] * Gam[e][b][dd] - Gam[a][dd][e] * Gam[e][a][b]
        out[(b, dd)] = acc
    return out


def ddt_from_equation(gam, dx, dxx, dt, Gi):
    """d_t^2 gamma = 2 d_t K = 2 (-Ric(gamma) + 2 K G K - tr(G K) K), K = (1/2) d_t gamma, evaluated on
    enclosures: valid wherever the inputs enclose the CK solution, which satisfies these equations
    exactly."""
    R = range(3)
    Ric = ricci_pointwise(gam, dx, dxx, Gi)
    K = [[sym_get(dt, i, j) / 2 for j in R] for i in R]
    KG = [[sum((K[i][k] * Gi[k][j] for k in R), arb(0)) for j in R] for i in R]
    trGK = sum((KG[i][i] for i in R), arb(0))
    out = {}
    for (i, j) in PAIRS:
        KGK = sum((KG[i][k] * K[k][j] for k in R), arb(0))
        out[(i, j)] = 2 * (-Ric[(i, j)] + 2 * KGK - trGK * K[i][j])
    return out


def _positive_definite(A):
    """True if every symmetric matrix in the ball matrix A is positive definite (Sylvester)."""
    m1 = A[0, 0]
    m2 = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    m3 = A.det()
    return bool(m1 > 0 and m2 > 0 and m3 > 0)


# ---------------------------------------------------------------------------
# Demonstration (notes/16 §4)
# ---------------------------------------------------------------------------
ORDER_DEMO = 5      # metric through t^10, derivatives through t^8


def demo(out_path=os.path.join(DATA, "bl_metric_demo.json")):
    cert = load_certificate("ck_bound_ext.json")       # times are fractions of the whole-tube certificate
    d = "0.7662"
    T = cert["T_min"]
    pts = {"axis, near the top of S_c": ("0", "0", "0.62"),
           "plane z = 0, near the waist of S_c": ("0.26", "0", "0"),
           "generic point near S_c": ("0.17", "0.11", "0.5")}
    rows = []
    for name, (x, y, z) in pts.items():
        for frac in ("0.1", "0.25", "0.5"):
            with working_precision(128):
                t = arb(frac) * arb(T)
            t1 = time.time()
            r = evaluate_metric(t, x, y, z, d, order=ORDER_DEMO)
            rows.append({"point": name, "x": [x, y, z], "t_over_T": frac, "t": str(t),
                         "status": r["status"], "lorentzian": r["lorentzian"],
                         "gamma": {f"{i}{j}": str(r["g"][i + 1][j + 1]) for i, j in PAIRS},
                         "g_inv": {f"{i}{j}": str(r["g_inv"][i + 1][j + 1]) for i, j in PAIRS},
                         "d_t gamma": {f"{i}{j}": str(r["dg"][0][i + 1][j + 1]) for i, j in PAIRS},
                         "d_x gamma_xx": [str(r["dg"][a + 1][1][1]) for a in range(3)],
                         "d_t^2 gamma_xx": str(r["ddg"][0][0][1][1]),
                         "max_radius_by": r["max_radius_by"], "budget": r["error_budget"],
                         "seconds": time.time() - t1})
            print(name, frac, r["status"], {k: f"{v:.2e}" for k, v in r["max_radius_by"].items()},
                  f"{time.time() - t1:.1f}s", flush=True)
    # spatial patches (boxes around the generic point) over a time interval
    patch = []
    for hw in ("0.001", "0.0001"):
        t1 = time.time()
        box = evaluate_metric(arb(0, 0.5 * T), arb("0.17", hw), arb("0.11", hw), arb("0.5", hw), d,
                              order=ORDER_DEMO)
        patch.append({"box": f"x in 0.17 +- {hw}, y in 0.11 +- {hw}, z in 0.5 +- {hw}, |t| <= T_min/2",
                      "half_width": hw, "status": box["status"], "lorentzian": box["lorentzian"],
                      "gamma": {f"{i}{j}": str(box["g"][i + 1][j + 1]) for i, j in PAIRS},
                      "max_radius_by": box["max_radius_by"], "budget": box["error_budget"],
                      "seconds": time.time() - t1})
        print("patch", hw, box["status"], box["max_radius_by"], f"{time.time() - t1:.1f}s", flush=True)
    # exact low-order coefficients at a point with rational distances to both punctures
    t1 = time.time()
    dq = fmpq(3831, 5000)
    p = dq / 2
    k = fmpq(15, 8)
    xq = (p * (k * k - 1) / (2 * k), fmpq(0), fmpq(0))
    rq = p * (k * k + 1) / (2 * k)
    J, g, _ = series_jets(xq, dq, 2, fmpq(1), exact_r=[rq, rq])
    vals, _, _ = coefficients(J, g, fmpq(0), 2)
    exact = {"point": [str(v) for v in xq], "r_plus = r_minus": str(rq), "d": str(dq),
             "coefficients": [{f"{i}{j}": str(vals[n][(i, j)]) for (i, j) in PAIRS} for n in range(3)],
             "decimal": [{f"{i}{j}": float(vals[n][(i, j)]) for (i, j) in PAIRS} for n in range(3)],
             "seconds": time.time() - t1}
    out = {"d": d, "certificate": cert["file"], "T_min": T, "order": ORDER_DEMO, "precision": 128,
           "evaluations": rows, "patch": patch, "exact_low_order": exact,
           "replay": "python -m research.calc.bl_metric"}
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    return out


def compare(out_path=os.path.join(DATA, "bl_metric_compare.json"), orders=(5, 7)):
    """Before/after comparison of notes/17 §4, at fixed fractions of the session-12 T0 = 5.557e-4.
    'before_corrected': the session-12 constants (ck_bound_sharp.json) with the corrected
    Baouendi-Goulaouic time (one certified pair); 'after': every certificate, CK' schedules, combined.
    The session-12 numbers as originally claimed (invalid time constant) are in
    bl_metric_baseline_session12.json (produced at the start of session 13 by the session-12 code)."""
    import resource
    T0 = 5.556996340950146e-4
    sharp = load_certificate("ck_bound_sharp.json")
    corrected = {k: v for k, v in sharp.items() if k != "schedule"}
    corrected["T_min"], corrected["sup_u_weighted_bound"] = sharp["T_bg"], sharp["B_bg"]
    corrected["file"] = "ck_bound_sharp.json (BG (12) time only)"
    configs = {"before_corrected": [corrected], "after": load_certificates()}
    pts = {"generic (0.17, 0.11, 0.5)": ("0.17", "0.11", "0.5"), "axis (0, 0, 0.62)": ("0", "0", "0.62"),
           "waist (0.26, 0, 0.05)": ("0.26", "0", "0.05"),
           "box (0.17, 0.11, 0.5) +- 1e-4": (arb("0.17", "0.0001"), arb("0.11", "0.0001"), arb("0.5", "0.0001"))}
    rows = []
    for cname, certs in configs.items():
        for N in orders:
            if cname != "after" and N != orders[0]:
                continue
            for pname, (x, y, z) in pts.items():
                for f in ("0.1", "0.25", "0.5", "1.0", "2.0"):
                    with working_precision(128):
                        t = arb(f) * arb(T0)
                    t1 = time.time()
                    try:
                        r = evaluate_metric(t, x, y, z, "0.7662", order=N, certificate=certs)
                    except DomainError as e:
                        rows.append({"config": cname, "order": N, "point": pname, "t/T0": f, "status": "refused",
                                     "reason": str(e)[:200]})
                        continue
                    rows.append({"config": cname, "order": N, "point": pname, "t/T0": f, "status": r["status"],
                                 "metric": r["max_radius_by"]["metric"], "first": r["max_radius_by"]["first"],
                                 "second": r["max_radius_by"]["all"], "tails": r["error_budget"]["time_truncation"],
                                 "arith": r["error_budget"]["coefficient_arithmetic"],
                                 "certs": r["error_budget"]["certificate"], "seconds": round(time.time() - t1, 2)})
                    print(cname, N, pname, f, {k: f"{rows[-1][k]:.2e}" for k in ("metric", "first", "second")},
                          rows[-1]["seconds"], flush=True)
    out = {"T0_reference": T0, "rows": rows, "maxrss_MB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
           "replay": "python -m research.calc.bl_metric compare"}
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "compare":
        compare()
    else:
        demo()
