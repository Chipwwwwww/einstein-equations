"""Computer-assisted verification of the fold hypothesis (F) (notes/09, H4).

RESULT (session 3, 2026-10-29). With F(z0, d, t) the shooting function described below:

  (1) Krawczyk: F = dF/dz0 = 0 (t = 0) has exactly one solution in the box
      X = FOLD_CENTRE +- 1e-12, with z0_c = 0.632659552806 +- 3e-13 and
      d_c = 0.766197428240 +- 1.1e-13 (ADM mass M = 1), i.e. 1.532394856480 +- 2e-13 in
      units of the bare mass.
  (2) At the fold point: d^2F/dz0^2 in [240.87, 242.70] > 0, dF/dd in [1.2315, 1.2364] > 0,
      dF/dt in [-0.51816, -0.51655] < 0 (lapse N = 1), so the onset slope dF/dd / (-dF/dt)
      lies in [2.3768, 2.3936].
  (3) Along the curve: the variation field dX/dz0 has positive normal component (so the
      Jacobi field phi_0 > 0), sin(alpha) < 0 and rho > 0 after the axis, and z > 0 where the
      height becomes the independent variable: the critical surface is an embedded sphere
      enclosing both punctures.

Rigour rests on Arb's ball arithmetic and on this implementation, which is checked against
ODEs with exact solutions at several orders, jets against SymPy derivatives, F and its jets
against the independent floating-point code research/calc/bl_fold.py, and the whole
evaluation against a second, finer step schedule (all enclosures overlap). Run main() or the
slow test tests/test_fold_proof.py::test_fold_certificate (~5 min).

History, kept because it matters for trust:
- Session 2: a first "verified" run was INVALID. python-flint caps power series at length 10
  by default (ctx.cap), silently dropping Taylor terms and the remainder. The cap is now set
  explicitly and the known-solution tests catch any recurrence.
- Session 2 also blamed a "structural barrier" at the axis for the failure that followed.
  That diagnosis was WRONG. The blow-up came from low-order Taylor steps once the step size
  grew. The fix was the piecewise step schedule MU_SCHEDULE: order 8 near the axis with
  h <= s/64, order 18 further out with h <= min(s/32, 2^-9).

Shooting function. Start the meridian curve of an axisymmetric surface at the axis point
(0, z0) with horizontal tangent, solve theta_+ = H_S + t p = 0 (first order in t, lapse 1)
until the curve meets z = 0, and record F = alpha + pi/2 there. A reflection-symmetric common
MOTS is a zero of F.

Method. The ODE alpha' = -sin(alpha)/rho + G has a singular point on the axis.
  * Axis segment [0, s1], s1 = 2^-26: with alpha = s A, rho = s + s^3 Q, z = z0 + s^2 Z,
    A(s) = s^-2 Int_0^s sigma Phi,  Q(s) = -s^-3 Int_0^s sigma^2 (A^2/2) sinc^2(sigma A/2),
    Z(s) = s^-2 Int_0^s sigma A sinc(sigma A),  Phi = G + A (1 - sinc(s A) / (1 + s^2 Q)).
    These are weighted averages, so A, Q, Z lie in (1/2) hull Phi, -(1/6) hull(..),
    (1/2) hull(..). A box X with T(X) in int(X) encloses them, by continuation from s = 0.
  * mu = rho sin(alpha), for which (rho sin alpha)' = rho cos(alpha) G: the singular term
    cancels exactly. cos(alpha) = +sqrt(1 - (mu/rho)^2) is justified because the enclosure
    keeps |mu/rho| < 1.
  * then (rho, z, alpha); finally z itself as the independent variable (tau in [0, 1],
    z = z_switch (1 - tau)), so the end condition z = 0 is at a fixed point.
  * Each step uses Moore's method: a Taylor polynomial of the solution, a Picard a-priori
    enclosure B (y0 + [0,h] f(B) inside B), and a Lagrange remainder from the Taylor
    coefficient of order p+1 evaluated on B.
  * Jets in (z0, d, t) with monomials 1, z, d, t, z^2, zd (a quotient ring) carry F_z, F_d,
    F_t, F_zz, F_zd.
  * The Krawczyk test runs on Psi = (F, F_z): Psi at the centre comes from a point run, and
    the enclosure of D Psi over the box from a run with ball-valued parameters.

Run: python -m research.calc.fold_proof   (~10 min: Newton refinement, then the certificate)
"""

import math

from flint import arb, arb_series, ctx

ctx.prec = 256
ctx.cap = 64  # python-flint truncates power series at this length (default 10!)

# ---------------------------------------------------------------------------
# Jets in the parameters (z0, d, t): monomials 1, z, d, t, z^2, z d.
# Products are truncated to these monomials (a quotient ring, since the excluded
# monomials d^2, t^2, zt, dt, z^3, z^2 d generate an ideal). Coefficient of z^2 is
# F_zz / 2; of z d is F_zd.
# ---------------------------------------------------------------------------
_TABLE = [(0, k, k) for k in range(6)] + [(k, 0, k) for k in range(1, 6)] + [(1, 1, 4), (1, 2, 5), (2, 1, 5)]


class J:
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = c

    def __add__(a, b):
        if isinstance(b, J):
            return J([x + y for x, y in zip(a.c, b.c)])
        c = list(a.c)
        c[0] = c[0] + b
        return J(c)

    __radd__ = __add__

    def __neg__(a):
        return J([-x for x in a.c])

    def __sub__(a, b):
        return a + (-b)

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        if isinstance(b, J):
            out = [None] * 6
            for i, j, k in _TABLE:
                prod = a.c[i] * b.c[j]
                out[k] = prod if out[k] is None else out[k] + prod
            return J(out)
        return J([x * b for x in a.c])

    __rmul__ = __mul__

    def __truediv__(a, b):
        return a * inv(b) if isinstance(b, J) else J([x / b for x in a.c])

    def __rtruediv__(a, b):
        return inv(a) * b


def _compose(x, f0, f1, f2):
    """f(x) for a jet x, given f, f', f'' at the constant part (delta^3 = 0 in the ring)."""
    zero = x.c[0] * 0
    dlt = J([zero] + list(x.c[1:]))
    return J([f0] + [zero] * 5) + dlt * f1 + (dlt * dlt) * (f2 * 0.5)


def _sc(b):
    if isinstance(b, arb_series):
        return b.sin_cos()
    m = b.abs_upper()
    if m < 0.01:
        # Arb's ball sin/cos are crude when the radius is comparable to the midpoint (as on
        # the axis segment); use alternating-series enclosures, valid for |x| <= 1:
        # sin x = x - x^3/6 + R, |R| <= |x|^5/120;  cos x = 1 - x^2/2 + R, |R| <= x^4/24.
        x2 = b * b
        return (b - b * x2 / 6 + arb(0, m**5 / 120), 1 - x2 / 2 + arb(0, m**4 / 24))
    return b.sin(), b.cos()


def sin_cos(x):
    s, c = _sc(x.c[0])
    return _compose(x, s, c, -s), _compose(x, c, -s, -c)


def sqrt(x):
    r = x.c[0].sqrt()
    return _compose(x, r, 1 / (2 * r), -1 / (4 * r * r * r))


def inv(x):
    i = 1 / x.c[0]
    return _compose(x, i, -i * i, 2 * i * i * i)


def asin(x):  # arb base only
    x0 = x.c[0]
    q = 1 / (1 - x0 * x0).sqrt()
    return _compose(x, x0.asin(), q, x0 * q * q * q)


def sinc_small(x):
    """sinc of a jet whose constant part is a small ball (possibly containing 0), arb base.
    Alternating-series bounds valid for |x| <= 1:
      sinc = 1 - x^2/6 + R0, |R0| <= x^4/120;  sinc' = -x/3 + R1, |R1| <= |x|^3/30;
      sinc'' = -1/3 + R2, |R2| <= x^2/10."""
    x0 = x.c[0]
    m = x0.abs_upper()
    assert m < 1
    f0 = 1 - x0 * x0 / 6 + arb(0, m**4 / 120)
    f1 = -x0 / 3 + arb(0, m**3 / 30)
    f2 = arb(-1) / 3 + arb(0, m**2 / 10)
    return _compose(x, f0, f1, f2)


def ball(mid, rad=0):
    return arb(mid, rad) if rad else arb(mid)


def contained(new, box, interior=False):
    test = (lambda b, n: b.contains_interior(n)) if interior else (lambda b, n: b.contains(n))
    return all(test(bj, nj) for nb, bb in zip(new, box) for nj, bj in zip(nb.c, bb.c))


def inflate(v, factor=2, eps=1e-60):
    return J([arb(x.mid(), x.rad() * factor + arb(eps)) for x in v.c])


def hull_union(a, b):
    return J([x.union(y) for x, y in zip(a.c, b.c)])


# ---------------------------------------------------------------------------
# Brill-Lindquist geometry: bare masses 1/2 at z = +-d/2, lapse N = 1.
# ---------------------------------------------------------------------------
HALF_M = arb(1) / 4  # m/2 with m = 1/2


def psi_parts(rho, z, d):
    """psi, (d psi/d rho)/rho, d psi/dz for jets (or base objects)."""
    zp, zm = z - d * 0.5, z + d * 0.5
    i1 = inv(sqrt(rho * rho + zp * zp))
    i2 = inv(sqrt(rho * rho + zm * zm))
    c1, c2 = i1 * i1 * i1, i2 * i2 * i2
    return 1 + HALF_M * (i1 + i2), -HALF_M * (c1 + c2), -HALF_M * (zp * c1 + zm * c2)


def t_coefficient(rho, z, sa, ca, d):
    """psi^2 p = (Kdot_tt + Kdot_yy) / psi^2 at base (constant-part) values, N = 1,
    Kdot_ij = -R_ij, R_ij = -2 psi_ij/psi + 6 psi_i psi_j/psi^2 - 2 |grad psi|^2 delta_ij/psi^2."""
    psi, px, pz = 1, 0, 0
    pxx = pyy = pzz = pxz = 0
    for zc in (d * 0.5, -(d * 0.5)):
        Zr = z - zc
        r2 = rho * rho + Zr * Zr
        r = r2.sqrt()
        r3, r5 = r2 * r, r2 * r2 * r
        psi = psi + HALF_M / r
        px = px - HALF_M * rho / r3
        pz = pz - HALF_M * Zr / r3
        pxx = pxx + HALF_M * (3 * rho * rho - r2) / r5
        pyy = pyy - HALF_M / r3
        pzz = pzz + HALF_M * (3 * Zr * Zr - r2) / r5
        pxz = pxz + HALF_M * 3 * rho * Zr / r5
    g2 = px * px + pz * pz
    Rxx = -2 * pxx / psi + 6 * px * px / psi**2 - 2 * g2 / psi**2
    Ryy = -2 * pyy / psi - 2 * g2 / psi**2
    Rzz = -2 * pzz / psi + 6 * pz * pz / psi**2 - 2 * g2 / psi**2
    Rxz = -2 * pxz / psi + 6 * px * pz / psi**2
    ktt = -(ca * ca * Rxx + 2 * ca * sa * Rxz + sa * sa * Rzz)
    return (ktt - Ryy) / psi**2


def G(rho, z, sa, ca, d):
    """theta_+ ODE source: 4 n.grad(psi)/psi + t psi^2 p (the t part as a t-monomial)."""
    psi, prr, pz = psi_parts(rho, z, d)
    g = 4 * (-sa * rho * prr + ca * pz) / psi
    tc = t_coefficient(rho.c[0], z.c[0], sa.c[0], ca.c[0], d.c[0])
    zero = tc * 0
    return g + J([zero, zero, zero, tc, zero, zero])


# ---------------------------------------------------------------------------
# Validated Taylor integration (Moore's method with Picard a-priori enclosure).
# ---------------------------------------------------------------------------
def _const_series(v, L):
    return J([arb_series([x], prec=L) for x in v.c])


def taylor(rhs, y, L):
    """Taylor coefficients (as series of length L) of the solution through state y."""
    const = [_const_series(v, L) for v in y]
    ys = const
    for _ in range(L):
        f = rhs(ys)
        ys = [J([cc + fc.integral() for cc, fc in zip(c.c, fj.c)]) for c, fj in zip(const, f)]
    return ys


def _coeffs(series, L):
    c = series.coeffs()
    return c + [arb(0)] * (L - len(c))


def _eval(series, x, L, upto):
    c = _coeffs(series, L)
    acc = arb(0)
    for k in range(upto, -1, -1):
        acc = acc * x + c[k]
    return acc


def step(rhs, y, h, L):
    """One validated step of length h (exact float). Returns (y(h) enclosure, enclosure of
    y over [0, h])."""
    p = L - 2
    P = taylor(rhs, y, L)
    H = arb(h / 2, h / 2)
    guess = [J([_eval(s, H, L, p) for s in v.c]) for v in P]
    B = [inflate(g, 1.25, 1e-40) for g in guess]
    for _ in range(20):
        f = rhs(B)
        Bn = [yv + fv * H for yv, fv in zip(y, f)]
        if contained(Bn, B):
            break
        # enlarge only the components that are not yet contained
        B = [J([bj if bj.contains(nj) else arb(bj.union(nj).mid(), bj.union(nj).rad() * 1.25 + arb(1e-40))
                for bj, nj in zip(b.c, bn.c)]) for b, bn in zip(B, Bn)]
    else:
        raise RuntimeError("a-priori enclosure failed; reduce the step")
    PB = taylor(rhs, B, L)
    hb = arb(h)
    out = []
    for v, vb in zip(P, PB):
        comps = []
        for s, sb in zip(v.c, vb.c):
            comps.append(_eval(s, hb, L, p) + _coeffs(sb, L)[p + 1] * hb ** (p + 1))
        out.append(J(comps))
    return out, B


def integrate(rhs, y, t0, t1, h, L, check=None, geometric=False, ratio=None):
    """Integrate from t0 to t1 (state component 0 is the independent variable) with dyadic
    steps of at most h (last step shortened). With geometric=True the step is also at most
    the largest power of 2 not exceeding t/GEOMETRIC_RATIO, as needed near the axis. `check(B)` is called
    on every step enclosure."""
    t = t0
    n = 0
    while t < t1:
        hh = min(h, t1 - t)
        if geometric:
            hh = min(hh, 2.0 ** math.floor(math.log2(t / (ratio or GEOMETRIC_RATIO))))
        while True:
            try:
                y_new, B = step(rhs, y, hh, L)
                break
            except (RuntimeError, ValueError):
                # enclosure failed (or a ball reached a singular value): halve the step
                hh /= 2
                if hh < 1e-12:
                    raise
        y = y_new
        if check is not None:
            check(B)
        t += hh
        n += 1
    return y, n


# ---------------------------------------------------------------------------
# The three regular segments and the axis enclosure.
# ---------------------------------------------------------------------------
def make_rhs_mu(d):
    """State (s, rho, z, mu), mu = rho sin(alpha):  rho' = cos a, z' = sin a,
    mu' = rho cos a G."""
    def rhs(y):
        s, rho, z, mu = y
        sa = mu / rho
        ca = sqrt(1 - sa * sa)
        one = s * 0 + 1
        return [one, ca, sa, rho * ca * G(rho, z, sa, ca, d)]
    return rhs


def make_rhs_alpha(d):
    """State (s, rho, z, alpha)."""
    def rhs(y):
        s, rho, z, a = y
        sa, ca = sin_cos(a)
        one = s * 0 + 1
        return [one, ca, sa, -sa / rho + G(rho, z, sa, ca, d)]
    return rhs


def make_rhs_tau(d, zsw):
    """State (tau, rho, alpha) with z = zsw (1 - tau) on tau in [0, 1]."""
    def rhs(y):
        tau, rho, a = y
        sa, ca = sin_cos(a)
        z = zsw * (1 - tau)
        one = tau * 0 + 1
        ds = -zsw / sa
        return [one, ds * ca, ds * (-sa / rho + G(rho, z, sa, ca, d))]
    return rhs


def axis_enclosure(z0, d, s1):
    """Rigorous boxes for A, Q, Z on [0, s1] (see the module docstring)."""
    sig = arb(s1 / 2, s1 / 2)
    zero_j = J([arb(0)] * 6)
    rho_a, z_a = zero_j + 0, z0
    sa0, ca0 = zero_j + 0, zero_j + 1
    A = G(rho_a, z_a, sa0, ca0, d) * 0.5
    A = inflate(A, 1, 1e-6)
    Q = J([arb(0, 1)] + [arb(0, 1)] * 5)
    Z = inflate(A * 0.5, 1, 1e-6)
    for _ in range(40):
        x = A * sig
        sa, ca = sin_cos(x)
        rho = sig + Q * sig**3
        z = z0 + Z * sig**2
        Phi = G(rho, z, sa, ca, d) + A * (1 - sinc_small(x) / (1 + Q * sig**2))
        An = Phi * 0.5
        sh = sinc_small(x * 0.5)
        Qn = (A * A) * (sh * sh) * (arb(-1) / 6)
        Zn = A * sinc_small(x) * 0.5
        if contained([An, Qn, Zn], [A, Q, Z], interior=True):
            # Verified: T(X) inside the interior of X. The solution starts at the centre and,
            # by continuation, stays in T(X); so it also lies in T(T(X)), ... : tighten.
            A, Q, Z = An, Qn, Zn
            for _ in range(8):
                x = A * sig
                sa, ca = sin_cos(x)
                Phi = G(sig + Q * sig**3, z0 + Z * sig**2, sa, ca, d) \
                    + A * (1 - sinc_small(x) / (1 + Q * sig**2))
                sh = sinc_small(x * 0.5)
                A, Q, Z = Phi * 0.5, (A * A) * (sh * sh) * (arb(-1) / 6), A * sinc_small(x) * 0.5
            return A, Q, Z
        A, Q, Z = (inflate(hull_union(o, n), 1.5, 1e-40) for o, n in ((A, An), (Q, Qn), (Z, Zn)))
    raise RuntimeError("axis enclosure did not verify")


class Checks:
    """Collects rigorous bounds along the curve: normal component of dX/dz0 (Jacobi field),
    sin(alpha) and rho."""

    def __init__(self):
        self.axis_A_max = math.nan
        self.z_switch_min = math.nan
        self.jacobi_min = math.inf
        self.sin_max = -math.inf
        self.rho_min = math.inf

    def _update(self, rho, sa, ca, rho_z, z_z):
        n_dot = -sa * rho_z + ca * z_z
        self.jacobi_min = min(self.jacobi_min, float(n_dot.lower()))
        self.sin_max = max(self.sin_max, float(sa.upper()))
        self.rho_min = min(self.rho_min, float(rho.lower()))

    def mu(self, B):
        s, rho, z, mu = B
        sa = mu.c[0] / rho.c[0]
        self._update(rho.c[0], sa, (1 - sa * sa).sqrt(), rho.c[1], z.c[1])

    def alpha(self, B):
        s, rho, z, a = B
        sa, ca = a.c[0].sin(), a.c[0].cos()
        self._update(rho.c[0], sa, ca, rho.c[1], z.c[1])

    def tau(self, zsw):
        def f(B):
            tau, rho, a = B
            sa, ca = a.c[0].sin(), a.c[0].cos()
            self._update(rho.c[0], sa, ca, rho.c[1], zsw.c[1] * (1 - tau.c[0]))
        return f


S1 = 2.0**-26
GEOMETRIC_RATIO = 32  # near the axis, steps of at most s / 32 (dyadic), order L_axis
S_A = 0.25
S_SW = 0.375


# Step schedule. Near the axis the Taylor coefficients of the near-singular field carry
# radii growing like (C/s)^k, so steps must stay small relative to s; further out, low order
# with large steps lets the per-step remainder dominate. Each entry is
# (end of piece, Taylor length L = order + 2, step <= s / ratio, maximal step).
MU_SCHEDULE = [(2.0**-9, 10, 64, 1.0), (S_A, 20, 32, 2.0**-9)]
ALPHA_SCHEDULE = (20, 2.0**-9)   # (L, step) on [S_A, S_SW]
TAU_SCHEDULE = (20, 2.0**-7)     # (L, step) in tau on [0, 1]


def shoot(z0, d, checks=None, mu_schedule=None, alpha_schedule=None, tau_schedule=None):
    """Validated shooting. z0, d: arb balls. Returns the jet of F = alpha_end + pi/2 and the
    number of steps per segment."""
    mu_schedule = mu_schedule or MU_SCHEDULE
    La, ha = alpha_schedule or ALPHA_SCHEDULE
    Lt, ht = tau_schedule or TAU_SCHEDULE
    zj = J([z0, arb(1), arb(0), arb(0), arb(0), arb(0)])
    dj = J([d, arb(0), arb(1), arb(0), arb(0), arb(0)])
    A, Q, Z = axis_enclosure(zj, dj, S1)
    s1 = arb(S1)
    rho = s1 + Q * s1**3
    z = zj + Z * s1**2
    a = A * s1
    sa, _ = sin_cos(a)
    y = [J([s1] + [arb(0)] * 5), rho, z, rho * sa]
    if checks is not None:
        # Axis segment, 0 < s <= s1: alpha = s A with A < 0 gives sin(alpha) < 0, and
        # rho = s (1 + s^2 Q) > 0; the Jacobi field there is
        # -sin(alpha) rho_z + cos(alpha) z_z with rho_z = s^3 Q_z, z_z = 1 + s^2 Z_z.
        checks.axis_A_max = float(A.c[0].upper())
        assert A.c[0] < 0 and (Q.c[0] * S1**2).abs_upper() < 1
        sig = arb(S1 / 2, S1 / 2)
        aa = (A * sig).c[0]
        n_dot = -aa.sin() * Q.c[1] * sig**3 + aa.cos() * (1 + Z.c[1] * sig**2)
        checks.jacobi_min = min(checks.jacobi_min, float(n_dot.lower()))
    t, n1 = S1, 0
    rhs = make_rhs_mu(dj)
    for s_end, L, ratio, hmax in mu_schedule:
        y, n = integrate(rhs, y, t, s_end, hmax, L, checks.mu if checks else None,
                         geometric=True, ratio=ratio)
        t, n1 = s_end, n1 + n
    s, rho, z, mu = y
    y = [s, rho, z, asin(mu / rho)]
    y, n2 = integrate(make_rhs_alpha(dj), y, S_A, S_SW, ha, La, checks.alpha if checks else None)
    s, rho, z, a = y
    zsw = z
    if checks is not None:
        # z decreases along the whole curve (sin(alpha) < 0), so it first meets z = 0 in the
        # tau segment provided z > 0 at the switch.
        checks.z_switch_min = float(zsw.c[0].lower())
    y = [J([arb(0)] + [arb(0)] * 5), rho, a]
    y, n3 = integrate(make_rhs_tau(dj, zsw), y, 0.0, 1.0, ht, Lt,
                      checks.tau(zsw) if checks else None)
    F = y[2] + arb.pi() / 2
    return F, (n1, n2, n3)


NAMES = ["F", "F_z", "F_d", "F_t", "F_zz/2", "F_zd"]


def newton_refine(z, d, iters=3):
    """Floating-point Newton on Psi = (F, F_z) using validated point evaluations."""
    for _ in range(iters):
        F, _ = shoot(arb(z), arb(d))
        f, fz, fd, fzz, fzd = (float(F.c[k].mid()) for k in (0, 1, 2, 4, 5))
        fzz *= 2
        det = fz * fzd - fd * fzz
        dz = -(fzd * f - fd * fz) / det
        dd = -(-fzz * f + fz * fz) / det
        z, d = z + dz, d + dd
    return z, d


def krawczyk(z, d, r):
    """Krawczyk test for Psi(z0, d) = (F, F_z) = 0 on X = [z - r, z + r] x [d - r, d + r].
    Returns (verified, K, point jet, box jet, checks over X)."""
    Fp, _ = shoot(arb(z), arb(d))
    chk = Checks()
    Fb, steps = shoot(arb(z, r), arb(d, r), checks=chk)
    psi = [Fp.c[0], Fp.c[1]]
    D = [[Fb.c[1], Fb.c[2]], [2 * Fb.c[4], Fb.c[5]]]  # dPsi over X
    m = [[float(D[i][j].mid()) for j in range(2)] for i in range(2)]
    det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
    C = [[arb(m[1][1] / det), arb(-m[0][1] / det)], [arb(-m[1][0] / det), arb(m[0][0] / det)]]
    rad = arb(0, r)
    xh = [arb(z), arb(d)]
    K = []
    for i in range(2):
        val = xh[i] - (C[i][0] * psi[0] + C[i][1] * psi[1])
        for j in range(2):
            coef = (1 if i == j else 0) - (C[i][0] * D[0][j] + C[i][1] * D[1][j])
            val = val + coef * rad
        K.append(val)
    X = [arb(z, r), arb(d, r)]
    ok = all(X[i].contains_interior(K[i]) for i in range(2))
    return ok, K, Fp, Fb, chk, steps


# Newton-refined centre (from main(); recorded so the certificate can be re-run quickly).
FOLD_CENTRE = (0.6326595528058354, 0.766197428240025)
BOX_RADIUS = 1e-12


def certify_fold(z=FOLD_CENTRE[0], d=FOLD_CENTRE[1], r=BOX_RADIUS):
    """Run the Krawczyk test on X = (z, d) +- r, then re-run the box evaluation on the
    (much smaller) Krawczyk image K to get tight enclosures at the fold point.
    Returns a dict of rigorous results."""
    ok, K, Fp, Fb, chk, steps = krawczyk(z, d, r)
    rK = max(float(K[0].rad()), float(K[1].rad())) * 1.01
    zK, dK = float(K[0].mid()), float(K[1].mid())
    chkK = Checks()
    FK, _ = shoot(arb(zK, rK), arb(dK, rK), checks=chkK)
    # K is contained in the ball (zK, dK) +- rK, so these enclosures hold at the fold point.
    slope = FK.c[2] / (-FK.c[3])
    return {"verified": ok, "K": K, "point": Fp, "box": Fb, "checks_X": chk, "steps": steps,
            "at_fold": FK, "checks_K": chkK, "K_ball": (zK, dK, rK), "slope": slope}


def radius_growth_diagnostic(z=0.62, d=0.7662, s_end=0.05, L=8):
    """Integrate the near-axis segment from an EXACT state at s1 (zero radii) and print how
    the ball radius of mu grows: the overestimation that blocks the proof."""
    zj = J([arb(z), arb(1), arb(0), arb(0), arb(0), arb(0)])
    dj = J([arb(d), arb(0), arb(1), arb(0), arb(0), arb(0)])
    A, Q, Z = axis_enclosure(zj, dj, S1)
    s1 = arb(S1)
    rho = s1 + Q * s1**3
    sa, _ = sin_cos(A * s1)
    state = [J([s1] + [arb(0)] * 5), rho, zj + Z * s1**2, rho * sa]
    y = [J([arb(c.mid()) for c in v.c]) for v in state]
    rhs = make_rhs_mu(dj)
    t, k, marks = S1, 0, []
    while t < s_end:
        h = min(2.0**-6, s_end - t, 2.0 ** math.floor(math.log2(t / GEOMETRIC_RATIO)))
        y, _ = step(rhs, y, h, L)
        t += h
        k += 1
        if k % 80 == 0:
            marks.append((t, float(y[3].c[0].rad()), abs(float(y[3].c[0].mid()))))
    for t, rad, mu in marks:
        print(f"   s = {t:.2e}: radius of mu = {rad:.2e} (|mu| = {mu:.2e}, relative {rad / mu:.1e})")


def main():
    z, d = 0.6326595577095465, 0.7661974282424348  # floating-point fold from bl_fold
    print("Newton refinement with validated point evaluations ...", flush=True)
    z, d = newton_refine(z, d)
    print(f"   centre z0 = {z!r}, d = {d!r}", flush=True)
    res = certify_fold(z, d, BOX_RADIUS)
    K = res["K"]
    print(f"\nKrawczyk on X = centre +- {BOX_RADIUS:g} (steps per segment {res['steps']}): "
          f"verified = {res['verified']}")
    print(f"   the unique fold point lies in K: z0_c = {K[0].str(15, radius=True)}, "
          f"d_c = {K[1].str(15, radius=True)}")
    print("   jet at the centre of X (point evaluation):")
    for nm, c in zip(NAMES, res["point"].c):
        print(f"      {nm:7s} {c.str(14, radius=True)}")
    print("   enclosures over all of X:")
    for nm, c in zip(NAMES, res["box"].c):
        print(f"      {nm:7s} [{float(c.lower()):+.10e}, {float(c.upper()):+.10e}]")
    chk = res["checks_X"]
    print(f"   over X: Jacobi field >= {chk.jacobi_min:.6f}; sin(alpha) <= {chk.sin_max:.3e} "
          f"after the axis segment (A <= {chk.axis_A_max:.6f} on it); rho >= {chk.rho_min:.3e}; "
          f"z at the switch >= {chk.z_switch_min:.4f}")
    print("\nenclosures at the fold point (box run on a ball containing K):")
    for nm, c in zip(NAMES, res["at_fold"].c):
        print(f"      {nm:7s} [{float(c.lower()):+.12e}, {float(c.upper()):+.12e}]")
    sl = res["slope"]
    print(f"   onset slope F_d / (-F_t) in [{float(sl.lower()):.8f}, {float(sl.upper()):.8f}]")
    print(f"   d_c in bare-mass units: {(2 * K[1]).str(15, radius=True)}")


if __name__ == "__main__":
    main()
