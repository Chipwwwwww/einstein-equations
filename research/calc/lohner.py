"""Lohner's method for validated integration of autonomous ODEs, in Arb ball arithmetic.

Plain interval Taylor steps (research/calc/fold_proof.py) suffer from the wrapping effect: a
box carried by a rotating or hyperbolic flow is re-enclosed in a larger axis-parallel box at
every step. On the Brill-Lindquist shooting curves of gap 3, which wind around the individual
horizon and pass close to the axis, this blows up even point enclosures. Lohner's method keeps
the enclosure as a parallelepiped m + A r (A a moving, nearly orthogonal frame from a QR
factorisation, r a small box) and propagates it with the Jacobian of the step map:

    phi_h(x) in phi_h(m) + D phi_h(X) (x - m),   x in X = m + A r.

Per step:
  * an a-priori enclosure B of all trajectories from the hull of X over [0, h] (Picard);
  * phi_h(m) = sum_{k<=p} c_k(m) h^k + c_{p+1}(B) h^{p+1}  (Lagrange remainder);
  * D phi_h(X) = sum_{k<=p} Dc_k(X) h^k + Dc_{p+1}(B) Phi_box h^{p+1}, where c_k are the Taylor
    coefficients of the flow, Dc_k their gradients (forward-mode dual numbers on the Taylor
    series) and Phi_box encloses D phi_t over [0, h] (Picard for Phi' = Df Phi);
  * C = D phi_h(X) A,  A' = Q from a floating-point QR of mid(C) (inverted rigorously),
    r' = (A'^-1 C) r + A'^-1 (phi_h(m) - m').
Checked against flows with exact solutions in tests/test_lohner.py.
"""

import math

import numpy as np
from flint import arb, arb_mat, arb_series, ctx

ctx.cap = max(ctx.cap, 64)  # series truncation (python-flint default 10 is too short)


class D:
    """A value with a gradient (forward-mode dual number); components arb or arb_series."""

    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v, self.g = v, g

    def __add__(a, b):
        if isinstance(b, D):
            return D(a.v + b.v, [x + y for x, y in zip(a.g, b.g)])
        return D(a.v + b, a.g)

    __radd__ = __add__

    def __neg__(a):
        return D(-a.v, [-x for x in a.g])

    def __sub__(a, b):
        return a + (-b)

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        if isinstance(b, D):
            return D(a.v * b.v, [a.v * y + b.v * x for x, y in zip(a.g, b.g)])
        return D(a.v * b, [x * b for x in a.g])

    __rmul__ = __mul__

    def inv(a):
        i = 1 / a.v
        return D(i, [-x * i * i for x in a.g])

    def __truediv__(a, b):
        return a * b.inv() if isinstance(b, D) else D(a.v / b, [x / b for x in a.g])

    def __rtruediv__(a, b):
        return a.inv() * b

    def sqrt(a):
        s = a.v.sqrt()
        return D(s, [x / (2 * s) for x in a.g])


class Dl:
    """Tangent (dual) number: value + one directional derivative, for writing variational
    equations as f evaluated on Dl inputs. Components may be arb, arb_series or D (whose own
    gradients are then taken with respect to the whole state)."""

    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v, self.g = v, g

    def __add__(a, b):
        if isinstance(b, Dl):
            return Dl(a.v + b.v, a.g + b.g)
        return Dl(a.v + b, a.g)

    __radd__ = __add__

    def __neg__(a):
        return Dl(-a.v, -a.g)

    def __sub__(a, b):
        return a + (-b)

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        if isinstance(b, Dl):
            return Dl(a.v * b.v, a.v * b.g + a.g * b.v)
        return Dl(a.v * b, a.g * b)

    __rmul__ = __mul__

    def inv(a):
        i = 1 / a.v
        return Dl(i, -(a.g * i * i))

    def __truediv__(a, b):
        return a * b.inv() if isinstance(b, Dl) else Dl(a.v / b, a.g / b)

    def __rtruediv__(a, b):
        return a.inv() * b

    def sqrt(a):
        s = a.v.sqrt()
        return Dl(s, a.g / (2 * s))


def _integral(x):
    return D(x.v.integral(), [g.integral() for g in x.g]) if isinstance(x, D) else x.integral()


def taylor(f, y0, L, grad=False):
    """Taylor series (length L) of the solution of y' = f(y) through y0 (arb balls). With
    grad=True the components are dual numbers carrying d/dy0 (Jacobian of each coefficient)."""
    n = len(y0)
    if grad:
        const = [D(arb_series([y], prec=L), [arb_series([1 if i == j else 0], prec=L)
                                             for j in range(n)]) for i, y in enumerate(y0)]
    else:
        const = [arb_series([y], prec=L) for y in y0]
    ys = const
    for _ in range(L):
        F = f(ys)
        ys = [c + _integral(Fi) for c, Fi in zip(const, F)]
    return ys


def _coef(s, k):
    c = s.coeffs()
    return c[k] if k < len(c) else arb(0)


def _picard_box(f, X, h):
    """B containing every trajectory from the box X over [0, h]: B >= X + [0, h] f(B)."""
    H = arb(h / 2, h / 2)
    fx = f(X)
    B = [x + fi * H for x, fi in zip(X, fx)]
    B = [arb(b.mid(), b.rad() * 1.5 + arb(1e-30)) for b in B]
    for _ in range(30):
        Bn = [x + fi * H for x, fi in zip(X, f(B))]
        if all(b.contains(bn) for b, bn in zip(B, Bn)):
            return B
        # widen only the components that are not yet contained: widening all of them can
        # diverge when some widths depend quadratically on others
        B = [b if b.contains(bn) else
             arb(b.union(bn).mid(), b.union(bn).rad() * 1.5 + arb(1e-30)) for b, bn in zip(B, Bn)]
    raise RuntimeError("a-priori enclosure failed")


def _jacobian(f, B):
    """Interval Jacobian Df over the box B."""
    n = len(B)
    seeded = [D(b, [arb(1 if i == j else 0) for j in range(n)]) for i, b in enumerate(B)]
    F = f(seeded)
    return arb_mat([[F[i].g[j] for j in range(n)] for i in range(n)])


def _phi_box(Df, h):
    """Enclosure of D phi_t, t in [0, h]: Phi = I + [0, h] Df Phi (Picard)."""
    n = Df.nrows()
    H = arb(h / 2, h / 2)
    eye = arb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
    P = arb_mat([[arb(1 if i == j else 0, 0.1) for j in range(n)] for i in range(n)])
    for _ in range(40):
        Pn = eye + (Df * P) * H
        if all(P[i, j].contains(Pn[i, j]) for i in range(n) for j in range(n)):
            return P
        P = arb_mat([[arb(P[i, j].union(Pn[i, j]).mid(), P[i, j].union(Pn[i, j]).rad() * 1.5 + arb(1e-30))
                      for j in range(n)] for i in range(n)])
    raise RuntimeError("variational a-priori enclosure failed")


class Set:
    """The parallelepiped m + A r."""

    def __init__(self, m, A, r):
        self.m, self.A, self.r = m, A, r

    @classmethod
    def from_box(cls, X):
        n = len(X)
        m = [arb(x.mid()) for x in X]
        A = arb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
        return cls(m, A, [x - mi for x, mi in zip(X, m)])

    def hull(self):
        n = len(self.m)
        return [self.m[i] + sum((self.A[i, j] * self.r[j] for j in range(n)), arb(0)) for i in range(n)]


def step(f, S, h, L):
    """One Lohner step of length h. Returns (new Set, a-priori enclosure B over the step)."""
    n, p = len(S.m), L - 2
    X = S.hull()
    B = _picard_box(f, X, h)
    hb = arb(h)
    # centre trajectory
    Tm = taylor(f, S.m, L)
    TB = taylor(f, B, L)
    xm = []
    for i in range(n):
        acc = arb(0)
        for k in range(p, -1, -1):
            acc = acc * hb + _coef(Tm[i], k)
        xm.append(acc + _coef(TB[i], p + 1) * hb ** (p + 1))
    # Jacobian of the step map over X
    TX = taylor(f, X, L, grad=True)
    TBg = taylor(f, B, L, grad=True)
    Pbox = _phi_box(_jacobian(f, B), h)
    J = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            acc = arb(0)
            for k in range(p, -1, -1):
                acc = acc * hb + _coef(TX[i].g[j], k)
            J[i, j] = acc
    Rem = arb_mat([[_coef(TBg[i].g[j], p + 1) for j in range(n)] for i in range(n)]) * Pbox
    J = J + Rem * hb ** (p + 1)
    C = J * S.A
    # new frame: QR of mid(C), columns ordered by the width they carry
    Cm = np.array([[float(C[i, j].mid()) for j in range(n)] for i in range(n)])
    widths = np.array([float(S.r[j].rad()) * np.linalg.norm(Cm[:, j]) for j in range(n)])
    order = np.argsort(-widths)
    Q, _ = np.linalg.qr(Cm[:, order])
    A_new = arb_mat([[float(Q[i, j]) for j in range(n)] for i in range(n)])
    Ainv = A_new.inv()
    m_new = [arb(x.mid()) for x in xm]
    M = Ainv * C
    r_new = []
    for i in range(n):
        acc = arb(0)
        for j in range(n):
            acc += M[i, j] * S.r[j]
        for j in range(n):
            acc += Ainv[i, j] * (xm[j] - m_new[j])
        r_new.append(acc)
    return Set(m_new, A_new, r_new), B


# ---------------------------------------------------------------------------
# Parameters kept out of the frame. With constant parameters p treated as extra state
# variables (p' = 0), the QR frames mix parameter and state directions, and the parameter box
# gets re-wrapped at every step: on the thin-neck system of gap 3 the set blows up after
# sigma ~ 0.7 (research/idea_log.md, session 7). Here the set is
#     x = m + A r + P q,   p = p0 + q,   q in a FIXED box,
# so the parameter dependence is carried linearly by the point matrix P and never re-boxed:
#     x_new in phi(m, p0) + Jx (A r + P q) + Jp q
#           = xm + (Jx A) r + (Jx P + Jp) q,
#     P' = mid(Jx P + Jp),   r' = A'^-1 [ (Jx A) r + (Jx P + Jp - P') q + (xm - m') ].
# f acts on the augmented state (x, p) and returns zeros for the parameter components.
# ---------------------------------------------------------------------------
class PSet:
    """The set m + A r + P q of states, for parameters p0 + q; q is fixed for the whole run."""

    def __init__(self, m, A, r, P, p0, q):
        self.m, self.A, self.r, self.P, self.p0, self.q = m, A, r, P, p0, q

    @classmethod
    def from_box(cls, X, params):
        n, k = len(X), len(params)
        m = [arb(x.mid()) for x in X]
        p0 = [arb(p.mid()) for p in params]
        A = arb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
        return cls(m, A, [x - mi for x, mi in zip(X, m)], arb_mat(n, k), p0,
                   [p - pi for p, pi in zip(params, p0)])

    def hull(self):
        n, k = len(self.m), len(self.p0)
        return [self.m[i] + sum((self.A[i, j] * self.r[j] for j in range(n)), arb(0))
                + sum((self.P[i, l] * self.q[l] for l in range(k)), arb(0)) for i in range(n)]

    def params(self):
        return [a + b for a, b in zip(self.p0, self.q)]


def _qr_frame(Cm, widths, blocks):
    """Orthogonal frame from QR of the point matrix Cm, columns ordered by the width they carry.
    With blocks (lists of state indices), the frame is block diagonal: one QR per diagonal
    block. For a triangular system (variational equations, jets) this keeps the widths of later
    blocks out of the earlier ones: with one full QR, the large widths of the jets leak into the
    base state through the rotated frame and the set blows up (research/calc/gap3_fold.py)."""
    n = Cm.shape[0]
    Q = np.zeros((n, n))
    for blk in blocks or [list(range(n))]:
        idx = np.array(blk)
        order = idx[np.argsort(-widths[idx])]
        q, _ = np.linalg.qr(Cm[np.ix_(idx, order)])
        Q[np.ix_(idx, idx)] = q
    return Q


def step_param(f, S, h, L, blocks=None):
    """One Lohner step of length h for a PSet; f is the augmented right-hand side.
    Returns (new PSet, a-priori enclosure B of the augmented state over the step).
    blocks: optional partition of the state indices for a block-diagonal frame (_qr_frame)."""
    n, k, p = len(S.m), len(S.p0), L - 2
    N = n + k
    Y = S.hull() + S.params()
    B = _picard_box(f, Y, h)
    hb = arb(h)
    Tm = taylor(f, S.m + S.p0, L)
    TB = taylor(f, B, L)
    xm = []
    for i in range(n):
        acc = arb(0)
        for kk in range(p, -1, -1):
            acc = acc * hb + _coef(Tm[i], kk)
        xm.append(acc + _coef(TB[i], p + 1) * hb ** (p + 1))
    TX = taylor(f, Y, L, grad=True)
    TBg = taylor(f, B, L, grad=True)
    Pbox = _phi_box(_jacobian(f, B), h)
    J = arb_mat(n, N)
    for i in range(n):
        for j in range(N):
            acc = arb(0)
            for kk in range(p, -1, -1):
                acc = acc * hb + _coef(TX[i].g[j], kk)
            J[i, j] = acc
    Rem = arb_mat([[_coef(TBg[i].g[j], p + 1) for j in range(N)] for i in range(n)]) * Pbox
    J = J + Rem * hb ** (p + 1)
    Jx = arb_mat([[J[i, j] for j in range(n)] for i in range(n)])
    Jp = arb_mat([[J[i, n + l] for l in range(k)] for i in range(n)])
    C1 = Jx * S.A
    C2 = Jx * S.P + Jp
    Cm = np.array([[float(C1[i, j].mid()) for j in range(n)] for i in range(n)])
    widths = np.array([float(S.r[j].rad()) * np.linalg.norm(Cm[:, j]) for j in range(n)])
    Q = _qr_frame(Cm, widths, blocks)
    A_new = arb_mat([[float(Q[i, j]) for j in range(n)] for i in range(n)])
    Ainv = A_new.inv()
    P_new = arb_mat([[arb(C2[i, l].mid()) for l in range(k)] for i in range(n)])
    m_new = [arb(x.mid()) for x in xm]
    M = Ainv * C1
    E = Ainv * (C2 - P_new)
    r_new = []
    for i in range(n):
        acc = arb(0)
        for j in range(n):
            acc += M[i, j] * S.r[j]
        for l in range(k):
            acc += E[i, l] * S.q[l]
        for j in range(n):
            acc += Ainv[i, j] * (xm[j] - m_new[j])
        r_new.append(acc)
    return PSet(m_new, A_new, r_new, P_new, S.p0, S.q), B
