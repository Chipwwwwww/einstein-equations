"""An explicit apparent-horizon merger in the near-critical Brill-Lindquist development (notes/20).

RESULT. See notes/20 (Theorem M) and research/calc/data/bl_barrier.json. In the geodesic-slicing
development of equal-mass BL data at an explicit separation d (default d = 0.7662, just above
d_c = 0.766197428240 +- 1.1e-13), the slice Sigma_t at an explicit time t (default t = 1.2e-5 M)
contains two explicit closed surfaces S_in inside S_out, both enclosing both punctures, with

    theta_+[S_in] < 0  (outer trapped),        theta_+[S_out] > 0  (outer untrapped),

where theta_+ = tr_S K + H is the TRUE outgoing null expansion of the slice data (gamma(t), K(t)), not a
truncation. The region Omega between them lies in the domain U of the Cauchy-Kovalevskaya certificate
(notes/15-17). Andersson-Metzger (CMP 290 (2009), Thms 3.1 and 4.1) then give a smooth, embedded, stable
MOTS in Omega homologous to S_out: a common apparent horizon of Sigma_t. Theorem F' (notes/12) says
Sigma_0 has none. So the common apparent horizon is born in (0, t].

Method (no Lyapunov-Schmidt reduction, no spectral gap):

1. Level surfaces of the first-order expansion. theta_+^(1) := H[gamma_0] + t tr_S[gamma_0](kappa_0),
   kappa_0 = -Ric(gamma_0) = g_1, is the expansion to first order in t (notes/15 §7.1). For a level h the
   profile curve of an axisymmetric surface with theta_+^(1) = h solves the shooting ODE of fold_proof.py
   with source  G = 4 n.grad(psi)/psi + psi^2 (t p - h),  p = tr_S kappa_0.  Validated shooting
   (fold_proof's Moore/Picard steps, axis enclosure, three segments) at z0 - delta and z0 + delta gives
   mismatches F of opposite signs; by continuity some z0* in between closes the curve at a right angle on
   the plane z = 0, so S_h (rotation + reflection of that curve) is a smooth closed surface on which
   theta_+^(1) = h EXACTLY. A run with z0 the whole ball [z0 - delta, z0 + delta] encloses every such curve.

2. The true expansion on S_h is theta_+ = h + R with R = Theta(p0 + dp) - Theta(p0), where Theta is the
   expansion as a function of the pointwise data p = (gamma_ij, d_l gamma_ij, K_ij) (theta_general below,
   valid for any axisymmetric metric with gamma_xy = gamma_yz = 0 on the meridian plane y = 0), p0 is the
   first-order data (g_0, d g_0, t g_1) and dp = p_true - p0. Componentwise
       d gamma  in  t^2 g_1 + [-V, V],      V  = tail of gamma beyond g_0 + t^2 g_1  (order N = 1),
       d dgamma in  t^2 d g_1 + [-D, D],    D  = tail of d_x gamma beyond d g_0 + t^2 d g_1 (N = 2),
       d K      in  2 t^3 g_2 + [-E, E],    E  = (1/2) tail of d_t gamma beyond 2 t g_1 + 4 t^3 g_2 (N = 2),
   with g_1, d g_1, g_2 enclosed on each box by the exact jet recursion of bl_metric.series_jets and V, D,
   E the certified CK tail bounds (bl_metric.tail_bounds, minimum over the certificates whose domain
   contains the box). By the mean-value theorem R = f'(sigma) for f(sigma) = Theta(p0 + sigma dp), some
   sigma in [0, 1]; an arb_series evaluation in sigma with sigma = [0, 1] and dp the balls above encloses
   f' on the whole box. Every box of the profile (the axis segment [0, s1] from the axis enclosure, then
   sub-boxes of every Taylor step) is checked.

3. Geometry: rho > 0 and sin(alpha) < 0 after the axis (the profile is a graph over z, so S_h is
   embedded), x.n > 0 (S_h is a radial graph about the origin), and on a fine polar subdivision
   max r(S_in) < min r(S_out) per angular bin (so S_in lies inside S_out) and every polar box between them
   lies in the certified domain U.

Integrators. The Moore/Picard integrator below is fold_proof.py's (same step, same schedules, same axis
enclosure), copied only to take the source term G as a parameter and to record each step's Taylor
polynomial; tests/test_bl_barrier.py checks that it reproduces fold_proof.shoot bit for bit. It gives
the sign change of F (point runs, enclosures ~5e-14) and the enclosure of the bracket's curves used for
the theta_+ bounds (ball run). The bracket must be tiny: its width grows by wrapping by a factor ~3e7
along the curve (delta = 1e-8 blows up; delta = 1e-10 reaches the plane with rho-widths ~6e-3). So the
centre is refined by validated secant steps, and the thin enclosures needed for the nesting test come
from Lohner sets over the same bracket, with sub-step enclosures (research/calc/bl_foliate.py), which
also re-check theta_+ on their own boxes.

Run: python -m research.calc.bl_barrier   (writes research/calc/data/bl_barrier.json; ~40 min, 2 cores)
"""

import json
import math
import os
import time

from flint import arb, arb_series, ctx

from research.calc import bl_metric as bm
from research.calc import fold_proof as fp
from research.calc.fold_proof import (J, asin, contained, hull_union, inflate, inv, sin_cos, sinc_small,
                                      sqrt, taylor, _coeffs, _eval)

DATA = os.path.join(os.path.dirname(__file__), "data")
HALF_M = arb(1) / 4
DC_ENCLOSURE = (0.766197428240 - 1.1e-13, 0.766197428240 + 1.1e-13)   # Theorem C


# ---------------------------------------------------------------------------
# The expansion of a surface of revolution in an axisymmetric metric
# ---------------------------------------------------------------------------
def theta_general(sa, ca, k1, k2, gam, dg, K, sqrt_=None):
    """theta_+ = tr_S K + H of the surface of revolution whose profile passes through the meridian
    point (rho, 0, z) with flat tangent T = (cos a, 0, sin a), flat normal n = (-sin a, 0, cos a)
    (outward), k1 = alpha' (flat arclength) and k2 = sin(alpha)/rho, for a metric with
    gamma_xy = gamma_yz = 0 at the point (true for every metric invariant under rotations about the z
    axis and the reflection y -> -y). Repository conventions: K_ij = (1/2) d_t gamma_ij for unit lapse,
    theta_+ = tr_S K + H, H = div nu with nu the outward unit normal.

    gam: {"xx","yy","zz","xz"};  K: same keys;  dg: {(l, ij)} for l in "x","z" and ij in
    "xx","yy","zz","xz", plus ("y","xy") and ("y","yz") (the only other components that enter).

        H = -(1/|n|) [ (k1 + n.Gamma(T,T)) / gamma_TT + (k2 + n.Gamma_yy) / gamma_yy ],
    from A_ab = -nu_i (d_a d_b X^i + Gamma^i_jk d_a X^j d_b X^k) on the parametrisation
    X(s, phi) = (rho cos phi, rho sin phi, z), where n.Gamma(V,W) = n^l Gamma_{l,jk} V^j W^k,
    n^l = gamma^{lk} n_k. Checked against a SymPy level-set computation in tests/test_bl_barrier.py.
    Works for any number type with +, -, *, / (arb, arb_series, float) given a sqrt."""
    sq = sqrt_ or (lambda v: v.sqrt())
    gxx, gyy, gzz, gxz = gam["xx"], gam["yy"], gam["zz"], gam["xz"]
    det = gxx * gzz - gxz * gxz
    nx, nz = -sa, ca
    nux = (gzz * nx - gxz * nz) / det
    nuz = (-gxz * nx + gxx * nz) / det
    n2 = nx * nux + nz * nuz
    gTT = ca * ca * gxx + 2 * ca * sa * gxz + sa * sa * gzz
    T = {"x": ca, "z": sa}

    def dgv(l, a, b):
        key = a + b if a + b in ("xx", "zz", "xz") else b + a
        return dg[(l, key)]

    GTT = {}
    for l in ("x", "z"):
        acc = 0
        for j in ("x", "z"):
            for k in ("x", "z"):
                acc = acc + T[j] * T[k] * (dgv(j, l, k) - dgv(l, j, k) * 0.5)
        GTT[l] = acc
    Gyy = {"x": dg[("y", "xy")] - dg[("x", "yy")] * 0.5, "z": dg[("y", "yz")] - dg[("z", "yy")] * 0.5}
    nGTT = nux * GTT["x"] + nuz * GTT["z"]
    nGyy = nux * Gyy["x"] + nuz * Gyy["z"]
    H = -((k1 + nGTT) / gTT + (k2 + nGyy) / gyy) / sq(n2)
    KTT = ca * ca * K["xx"] + 2 * ca * sa * K["xz"] + sa * sa * K["zz"]
    return H + KTT / gTT + K["yy"] / gyy


# ---------------------------------------------------------------------------
# The source term of the level-h first-order expansion equation
# ---------------------------------------------------------------------------
def _sqrt(x):
    return sqrt(x) if isinstance(x, J) else x.sqrt()


def _inv(x):
    return inv(x) if isinstance(x, J) else 1 / x


def psi2p(rho, z, sa, ca, d):
    """psi^2 p = psi^2 tr_S kappa_0 = (kappa_TT + kappa_yy) / psi^2 with kappa_0 = -Ric(psi^4 delta),
    Ric_ij = -2 psi_ij/psi + 6 psi_i psi_j/psi^2 - 2 |grad psi|^2 delta_ij/psi^2 (lapse 1), at (rho, 0, z).
    Same formula as fold_proof.t_coefficient, written with generic operations so that it acts on jets
    (fold_proof evaluates it only on constant parts, because there t is a jet variable)."""
    psi, px, pz = 1, 0, 0
    pxx = pyy = pzz = pxz = 0
    for zc in (d * 0.5, -(d * 0.5)):
        Zr = z - zc
        ir = _inv(_sqrt(rho * rho + Zr * Zr))
        ir3 = ir * ir * ir
        ir5 = ir3 * ir * ir
        psi = psi + HALF_M * ir
        px = px - HALF_M * rho * ir3
        pz = pz - HALF_M * Zr * ir3
        pxx = pxx + HALF_M * (3 * rho * rho * ir5 - ir3)
        pyy = pyy - HALF_M * ir3
        pzz = pzz + HALF_M * (3 * Zr * Zr * ir5 - ir3)
        pxz = pxz + HALF_M * 3 * rho * Zr * ir5
    ip = _inv(psi)
    ip2 = ip * ip
    g2 = px * px + pz * pz
    Rxx = -2 * pxx * ip + (6 * px * px - 2 * g2) * ip2
    Ryy = -2 * pyy * ip - 2 * g2 * ip2
    Rzz = -2 * pzz * ip + (6 * pz * pz - 2 * g2) * ip2
    Rxz = -2 * pxz * ip + 6 * px * pz * ip2
    ktt = -(ca * ca * Rxx + 2 * ca * sa * Rxz + sa * sa * Rzz)
    return (ktt - Ryy) * ip2


def make_source(t, h):
    """G for theta_+^(1) = h at time t: 4 n.grad(psi)/psi + t psi^2 p - h psi^2 (t, h: arb)."""
    def G(rho, z, sa, ca, d):
        psi, prr, pz = fp.psi_parts(rho, z, d)
        g = 4 * (-sa * rho * prr + ca * pz) / psi
        return g + psi2p(rho, z, sa, ca, d) * t - psi * psi * h
    return G


# ---------------------------------------------------------------------------
# fold_proof's validated integrator, with the source as a parameter and step recording
# ---------------------------------------------------------------------------
def step_rec(rhs, y, h, L):
    """fold_proof.step, also returning the Taylor polynomials P (at the step start) and the order-(p+1)
    coefficients over the a-priori enclosure B, so that sub-intervals can be enclosed later."""
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
    rec = {"h": h, "L": L, "P": [_coeffs(v.c[0], L) for v in P], "R": [_coeffs(vb.c[0], L)[p + 1] for vb in PB],
           "B": [b.c[0] for b in B], "Bz": [b.c[1] for b in B]}
    return out, B, rec


def integrate_rec(rhs, y, t0, t1, h, L, kind, steps, geometric=False, ratio=None):
    """fold_proof.integrate with step_rec; appends (kind, t_start, record) to steps."""
    t = t0
    while t < t1:
        hh = min(h, t1 - t)
        if geometric:
            hh = min(hh, 2.0 ** math.floor(math.log2(t / (ratio or fp.GEOMETRIC_RATIO))))
        while True:
            try:
                y_new, B, rec = step_rec(rhs, y, hh, L)
                break
            except (RuntimeError, ValueError):
                hh /= 2
                if hh < 1e-12:
                    raise
        steps.append((kind, t, rec))
        y = y_new
        t += hh
    return y


def rhs_mu(d, G):
    def rhs(y):
        s, rho, z, mu = y
        sa = mu / rho
        ca = sqrt(1 - sa * sa)
        return [s * 0 + 1, ca, sa, rho * ca * G(rho, z, sa, ca, d)]
    return rhs


def rhs_alpha(d, G):
    def rhs(y):
        s, rho, z, a = y
        sa, ca = sin_cos(a)
        return [s * 0 + 1, ca, sa, -sa / rho + G(rho, z, sa, ca, d)]
    return rhs


def rhs_tau(d, zsw, G):
    def rhs(y):
        tau, rho, a = y
        sa, ca = sin_cos(a)
        z = zsw * (1 - tau)
        ds = -zsw / sa
        return [tau * 0 + 1, ds * ca, ds * (-sa / rho + G(rho, z, sa, ca, d))]
    return rhs


def axis_enclosure(z0, d, s1, G):
    """fold_proof.axis_enclosure with the source G as a parameter."""
    sig = arb(s1 / 2, s1 / 2)
    zero_j = J([arb(0)] * 6)
    A = G(zero_j + 0, z0, zero_j + 0, zero_j + 1, d) * 0.5
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


def shoot(z0, d, G, record=False, jets=True):
    """fold_proof.shoot with source G (same segments and schedules). Returns (F jet, steps, axis) where
    F = alpha_end + pi/2, steps = [(kind, start, record)] if record, axis = (A, Q, Z, z0 jet, d jet).
    jets=False seeds no parameter derivatives (all jet slots stay exactly 0): only the values are
    integrated, and wide jet enclosures cannot spoil the Picard containment tests."""
    one = arb(1) if jets else arb(0)
    zj = J([z0, one, arb(0), arb(0), arb(0), arb(0)])
    dj = J([d, arb(0), one, arb(0), arb(0), arb(0)])
    S1 = fp.S1
    A, Q, Z = axis_enclosure(zj, dj, S1, G)
    s1 = arb(S1)
    rho = s1 + Q * s1**3
    z = zj + Z * s1**2
    a = A * s1
    sa, _ = sin_cos(a)
    y = [J([s1] + [arb(0)] * 5), rho, z, rho * sa]
    steps = []
    t = S1
    rhs = rhs_mu(dj, G)
    for s_end, L, ratio, hmax in fp.MU_SCHEDULE:
        y = integrate_rec(rhs, y, t, s_end, hmax, L, "mu", steps, geometric=True, ratio=ratio)
        t = s_end
    s, rho, z, mu = y
    y = [s, rho, z, asin(mu / rho)]
    La, ha = fp.ALPHA_SCHEDULE
    y = integrate_rec(rhs_alpha(dj, G), y, fp.S_A, fp.S_SW, ha, La, "alpha", steps)
    s, rho, z, a = y
    zsw = z
    Lt, ht = fp.TAU_SCHEDULE
    y = [J([arb(0)] + [arb(0)] * 5), rho, a]
    y = integrate_rec(rhs_tau(dj, zsw, G), y, 0.0, 1.0, ht, Lt, "tau:", steps)
    F = y[2] + arb.pi() / 2
    return F, (steps if record else None), (A, Q, Z, zj, dj, zsw)


# ---------------------------------------------------------------------------
# Boxes of the profile curve
# ---------------------------------------------------------------------------
def _horner(coeffs, u, p):
    acc = arb(0)
    for k in range(p, -1, -1):
        acc = acc * u + coeffs[k]
    return acc


def sub_states(rec, m):
    """Enclosures of the (constant-part) state on m equal sub-intervals of one recorded step."""
    h, L = rec["h"], rec["L"]
    p = L - 2
    out = []
    for i in range(m):
        u = arb((i + 0.5) * h / m, 0.5 * h / m)
        up = u ** (p + 1)
        out.append([_horner(c, u, p) + r * up for c, r in zip(rec["P"], rec["R"])])
    return out


def geometry(kind, state, zsw0):
    """(rho, z, sin a, cos a, k2 = sin(a)/rho) from a constant-part state of a segment."""
    if kind == "mu":
        s, rho, z, mu = state
        sa = mu / rho
        ca = (1 - sa * sa).sqrt()
        k2 = mu / (rho * rho)
    elif kind == "alpha":
        s, rho, z, a = state
        sa, ca = a.sin(), a.cos()
        k2 = sa / rho
    else:
        tau, rho, a = state
        z = zsw0 * (1 - tau)
        sa, ca = a.sin(), a.cos()
        k2 = sa / rho
    return rho, z, sa, ca, k2


def axis_boxes(axis, n=4):
    """Boxes of the axis segment 0 <= s <= s1 from the axis enclosure: rho = s + Q s^3, z = z0 + Z s^2,
    alpha = A s, sin(alpha)/rho = A sinc(A s) / (1 + Q s^2)."""
    A, Q, Z, zj, dj, _ = axis
    out = []
    for i in range(n):
        sig = arb((i + 0.5) * fp.S1 / n, 0.5 * fp.S1 / n)
        a = A.c[0] * sig
        x = J([a] + [arb(0)] * 5)
        sa, ca = (v.c[0] for v in sin_cos(x))
        k2 = A.c[0] * sinc_small(x).c[0] / (1 + Q.c[0] * sig**2)
        rho = sig + Q.c[0] * sig**3
        z = zj.c[0] + Z.c[0] * sig**2
        out.append((rho.union(arb(0)), z, sa, ca, k2))
    return out


# ---------------------------------------------------------------------------
# The remainder R = theta_+(true) - h on one box
# ---------------------------------------------------------------------------
KEYS = {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2), "xz": (0, 2), "xy": (0, 1), "yz": (1, 2)}


def tails_for_box(certs, t, x, z):
    """Certified tails (V, D, E) at |t| for the box: minimum over the certificates whose domain U
    contains it (and whose d-interval contains d, checked by the caller). None if none applies."""
    best = None
    for c, tb1, tb2 in certs:
        inside, _ = bm.in_domain(c, x, arb(0), z)
        if not inside:
            continue
        v = (tb1["value"], tb2["d_space"], tb2["d_time"] / 2)
        best = v if best is None else tuple(min(a, b) for a, b in zip(best, v))
    return best


def remainder(box, d, t, h, G, certs):
    """Enclosure of R on the box (rho, z, sa, ca, k2): theta_+(true) in h + R. Also returns the
    enclosure of Theta(p0) - h (should contain 0) as a consistency check."""
    rho, z, sa, ca, k2 = box
    Jb = lambda v: J([v] + [arb(0)] * 5)  # noqa: E731
    Gt = G(Jb(rho), Jb(z), Jb(sa), Jb(ca), Jb(d)).c[0]
    k1 = -k2 + Gt                                    # the ODE: alpha' = -sin(a)/rho + G
    tails = tails_for_box(certs, t, rho, z)
    if tails is None:
        raise bm.DomainError(f"box rho={rho}, z={z} not in any certified domain")
    V, D, E = (arb(0, v) for v in tails)
    Jt, g, _ = bm.series_jets((rho, arb(0), z), d, 2, arb(1))
    vals, d1, _ = bm.coefficients(Jt, g, arb(0), 2)
    t2, t3 = t * t, t * t * t
    sig = arb_series([arb(0.5, 0.5), 1], prec=2)
    gam, dg, K = {}, {}, {}
    gam0, dg0, K0 = {}, {}, {}
    for k in ("xx", "yy", "zz", "xz"):
        pq = KEYS[k]
        gam0[k], K0[k] = vals[0][pq], t * vals[1][pq]
        gam[k] = gam0[k] + sig * (t2 * vals[1][pq] + V)
        K[k] = K0[k] + sig * (2 * t3 * vals[2][pq] + E)
    comps = [(l, k) for l in ("x", "z") for k in ("xx", "yy", "zz", "xz")] + [("y", "xy"), ("y", "yz")]
    ax = {"x": 0, "y": 1, "z": 2}
    for l, k in comps:
        pq = KEYS[k]
        dg0[(l, k)] = d1[0][pq][ax[l]]
        dg[(l, k)] = dg0[(l, k)] + sig * (t2 * d1[1][pq][ax[l]] + D)
    R = theta_general(sa, ca, k1, k2, gam, dg, K).coeffs()[1]
    check = theta_general(sa, ca, k1, k2, gam0, dg0, K0) - h
    return R, check


# ---------------------------------------------------------------------------
# One barrier surface
# ---------------------------------------------------------------------------
def _down(v):
    return math.nextafter(float(v), -math.inf)


def _up(v):
    return math.nextafter(float(v), math.inf)


def _hull(lo, hi):
    return arb(lo).union(arb(hi))


def parameter_endpoints(value, radius=None):
    """Outward binary64 endpoints enclosing the SAME Arb parameter ball.

    In particular, float('6.5e-6') is not the exact decimal 6.5e-6.
    Use this bridge when one certificate hands a parameter to another.
    """
    ball = arb(value) if radius is None else arb(value, radius)
    return _down(ball.lower()), _up(ball.upper())


def polar(rho, z):
    r = (rho * rho + z * z).sqrt()
    phi = arb.atan2(rho, z)
    return phi, r


def certify_barrier(d, t, h, z0, delta, side, m_check=4, m_geo=32, prec=128, log=print):
    """Certify one barrier: S_h exists (F changes sign on [z0 - delta, z0 + delta]) and on every box
    theta_+(true) < 0 (side = "in") or > 0 (side = "out"). Returns a dict with the bounds and the
    fine polar boxes of the profile (for nesting and the domain check)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        d, t, h = (v if isinstance(v, arb) else arb(v) for v in (d, t, h))
        G = make_source(t, h)
        T0 = time.time()
        # validated secant refinement of the centre: the floating-point root is only good to ~3e-10,
        # while the bracket must stay ~1e-10 wide for the ball run (wider brackets blow up)
        zc = z0
        for it in range(5):
            Fm, _, _ = shoot(arb(zc - delta), d, G, jets=False)
            Fp, _, _ = shoot(arb(zc + delta), d, G, jets=False)
            fm, fpl = Fm.c[0], Fp.c[0]
            sign_change = (fm < 0 and fpl > 0) or (fm > 0 and fpl < 0)
            log(f"  [{side}] z0 = {zc!r}: F(z0 - delta) = {fm.str(6, radius=True)}, F(z0 + delta) = "
                f"{fpl.str(6, radius=True)}: sign change {sign_change} ({time.time() - T0:.0f} s)")
            if sign_change:
                break
            a, b = float(fm.mid()), float(fpl.mid())
            zc = zc - (a + b) / 2 * (2 * delta) / (b - a)
        else:
            raise AssertionError("no certified sign change of the shooting function")
        z0 = zc
        Fb, steps, axis = shoot(arb(z0, delta), d, G, record=True, jets=False)
        log(f"  [{side}] ball run: {len(steps)} steps ({time.time() - T0:.0f} s)")
        certs = []
        for c in bm.load_certificates():
            if not (d.lower() >= arb(c["d_lo"]) and d.upper() <= arb(c["d_hi"])):
                continue
            t_abs = _up(abs(t).upper())
            tb1, tb2 = bm.tail_bounds(c, t_abs, 1), bm.tail_bounds(c, t_abs, 2)
            if tb1 is None or tb2 is None:
                continue
            certs.append((c, tb1, tb2))
        if not certs:
            raise bm.DomainError("no certificate covers this d and t")
        zsw0 = axis[5].c[0]
        boxes = axis_boxes(axis)
        geo = []
        for kind, _, rec in steps:
            for st in sub_states(rec, m_check):
                boxes.append(geometry(kind, st, zsw0))
            for st in sub_states(rec, m_geo):
                geo.append(geometry(kind, st, zsw0))
        theta_hi, theta_lo, Rmax, check_ok = -math.inf, math.inf, 0.0, True
        sin_max_after_axis, rho_min, xn_min = -math.inf, math.inf, math.inf
        for i, b in enumerate(boxes):
            R, chk = remainder(b, d, t, h, G, certs)
            th = h + R
            theta_hi = max(theta_hi, float(th.upper()))
            theta_lo = min(theta_lo, float(th.lower()))
            Rmax = max(Rmax, float(abs(R).upper()))
            check_ok = check_ok and chk.contains(0)
        for j, (rho, z, sa, ca, k2) in enumerate(geo):
            sin_max_after_axis = max(sin_max_after_axis, float(sa.upper()))
            rho_min = min(rho_min, float(rho.lower()))
            xn_min = min(xn_min, float((-rho * sa + z * ca).lower()))
        ok = (theta_hi < 0) if side == "in" else (theta_lo > 0)
        log(f"  [{side}] {len(boxes)} boxes: theta_+ in [{theta_lo:.6e}, {theta_hi:.6e}], |R| <= {Rmax:.3e}, "
            f"first-order check contains 0: {check_ok}; after the axis segment sin(alpha) <= {sin_max_after_axis:.3e}, "
            f"rho >= {rho_min:.3e}; x.n >= {xn_min:.4f}; certified: {ok} ({time.time() - T0:.0f} s)")
        A = axis[0].c[0]
        return {
            "side": side, "h": float(h.mid()), "z0": z0, "delta": delta,
            "F_minus": [float(fm.lower()), float(fm.upper())], "F_plus": [float(fpl.lower()), float(fpl.upper())],
            "steps": len(steps), "boxes": len(boxes), "theta_lo": theta_lo, "theta_hi": theta_hi,
            "R_max": Rmax, "first_order_check": check_ok, "axis_A_max": float(A.upper()),
            "sin_alpha_max": sin_max_after_axis, "rho_min": rho_min, "x_dot_n_min": xn_min,
            "certified": bool(ok and check_ok and sin_max_after_axis < 0 and rho_min > 0 and xn_min > 0
                              and A < 0),
            "_geo": [(_down(rho.lower()), _up(rho.upper()), _down(z.lower()), _up(z.upper()))
                     for (rho, z, *_r) in geo],
            "_z0ball": (_down(axis[3].c[0].lower()), _up(axis[3].c[0].upper())),
        }
    finally:
        ctx.prec = old


NBINS = 4000


def radial_profile(g, nbins=NBINS):
    """Per angular bin of [0, pi/2]: outward-rounded (min r, max r) over the boxes of a profile
    (g["_geo"]: float boxes (rho_lo, rho_hi, z_lo, z_hi); g["_z0ball"]: the axis point)."""
    # Bins have exact dyadic width. Cover pi/2 outward and select every
    # intersected bin using interval division, never a rounded float ratio.
    width = _up(arb.pi().upper() / (2 * nbins))
    lo, hi = [math.inf] * nbins, [-math.inf] * nbins
    for (r0, r1, z0, z1) in g["_geo"]:
        phi, r = polar(_hull(r0, r1), _hull(z0, z1))
        ratio = phi / arb(width)
        a = max(0, int(ratio.lower().floor().fmpz()))
        b = min(nbins - 1, int(ratio.upper().floor().fmpz()))
        for k in range(a, b + 1):
            lo[k] = min(lo[k], _down(r.lower()))
            hi[k] = max(hi[k], _up(r.upper()))
    # the axis segment 0 <= s <= 2^-26 is not in _geo: phi < 2^-25 (bin 0), polar radius within
    # 2^-50 of the z0 ball
    z0lo, z0hi = g["_z0ball"]
    lo[0] = min(lo[0], z0lo - 1e-12)
    hi[0] = max(hi[0], z0hi + 1e-12)
    if not all(math.isfinite(v) for v in lo + hi):
        raise AssertionError("a bin is not covered by the profile boxes")
    return lo, hi, width


def radial_gap(inner, outer, nbins=NBINS):
    """min over bins of (min r of outer) - (max r of inner); > 0 means inner lies inside outer (both
    radial graphs)."""
    _, ih, _ = radial_profile(inner, nbins)
    ol, _, _ = radial_profile(outer, nbins)
    return min(_down((arb(o) - arb(i)).lower()) for o, i in zip(ol, ih))


def nesting_and_domain(inner, outer, nbins=NBINS, prec=128):
    """S_in inside S_out and Omega (between them) inside U. Both profiles are radial graphs over
    phi in [0, pi/2] (x.n > 0). Per angular bin: max r over boxes of S_in meeting the bin < min r over
    boxes of S_out meeting it; the polar box [bins] x [min r_in, max r_out] lies in U."""
    old = ctx.prec
    ctx.prec = prec
    try:
        r_in_lo, r_in_hi, width = radial_profile(inner, nbins)
        r_out_lo, r_out_hi, _ = radial_profile(outer, nbins)
        gap = min(_down((arb(o) - arb(i)).lower()) for o, i in zip(r_out_lo, r_in_hi))
        certs = [c for c in bm.load_certificates() if c.get("region") == "tube"]
        worst, dom_ok = math.inf, True
        coarse = 200
        per = nbins // coarse
        for j in range(coarse):
            ks = range(j * per, (j + 1) * per)
            rl = min(r_in_lo[k] for k in ks)
            rh = max(r_out_hi[k] for k in ks)
            phi = _hull(arb(j * per) * arb(width), arb((j + 1) * per) * arb(width))
            r = _hull(rl, rh)
            rho, z = r * phi.sin(), r * phi.cos()
            ins = [bm.in_domain(c, rho, arb(0), z) for c in certs]
            if not any(i for i, _ in ins):
                dom_ok = False
            need = min(m for _, m in ins)            # inflation of U's profile boxes needed here
            worst = need if worst == math.inf else max(worst, need)
        return {"nested": gap > 0, "min_radial_gap": gap, "Omega_in_U": dom_ok,
                "largest_inflation_needed": worst, "w_tube": certs[0]["w_tube"], "bins": nbins}
    finally:
        ctx.prec = old


# ---------------------------------------------------------------------------
# Non-existence in a shell before the merger: a foliation by leaves with theta_+ > 0
# ---------------------------------------------------------------------------
def certify_leaves(d, t_hi, h, z0, r, prec=128, m_check=4, m_geo=32, log=print):
    """One ball of the foliation. Leaves: the level-h shooting curves from (0, z0'), z0' in
    [z0 - r, z0 + r], at every time t in [0, t_hi], reflected in z = 0 and rotated. Certifies on the
    whole ball: F > 0 (the curve meets the plane at the angle F > 0 from the vertical, so the leaf has an
    outward, i.e. convex, equatorial ridge); theta_+(true) = h + R > 0 on every box; the normal
    component of dX/dz0 is > 0 (leaves strictly nested); rho > 0, sin(alpha) < 0 after the axis and
    x.n > 0 (embedded radial graphs)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        d, h = (v if isinstance(v, arb) else arb(v) for v in (d, h))
        tb = arb(t_hi / 2, t_hi / 2)            # every t in [0, t_hi]
        G = make_source(tb, h)
        T0 = time.time()
        F, steps, axis = shoot(arb(z0, r) if r else arb(z0), d, G, record=True)
        f = F.c[0]
        certs = _certs_for(d, t_hi)
        A, Q, Z, zj, dj, zsw = axis
        zsw0 = zsw.c[0]
        boxes = axis_boxes(axis)
        geo = []
        jac_min = math.inf
        # axis segment Jacobi field: -sin(A s) rho_z + cos(A s) z_z, rho_z = s^3 Q_z, z_z = 1 + s^2 Z_z
        sig = arb(fp.S1 / 2, fp.S1 / 2)
        aa = A.c[0] * sig
        jac_min = min(jac_min, float((-aa.sin() * Q.c[1] * sig**3 + aa.cos() * (1 + Z.c[1] * sig**2)).lower()))
        for kind, _, rec in steps:
            for st in sub_states(rec, m_check):
                boxes.append(geometry(kind, st, zsw0))
            for st in sub_states(rec, m_geo):
                geo.append(geometry(kind, st, zsw0))
            Bc, Bz = rec["B"], rec["Bz"]
            if kind == "mu":
                sa = Bc[3] / Bc[1]
                ca = (1 - sa * sa).sqrt()
                rho_z, z_z = Bz[1], Bz[2]
            elif kind == "alpha":
                sa, ca = Bc[3].sin(), Bc[3].cos()
                rho_z, z_z = Bz[1], Bz[2]
            else:
                sa, ca = Bc[2].sin(), Bc[2].cos()
                rho_z, z_z = Bz[1], zsw.c[1] * (1 - Bc[0])
            jac_min = min(jac_min, float((-sa * rho_z + ca * z_z).lower()))
        theta_lo, Rmax = math.inf, 0.0
        for b in boxes:
            R, _ = remainder(b, d, tb, h, G, certs)
            theta_lo = min(theta_lo, float((h + R).lower()))
            Rmax = max(Rmax, float(abs(R).upper()))
        sin_max = max(float(g[2].upper()) for g in geo)
        rho_min = min(float(g[0].lower()) for g in geo)
        xn_min = min(float((-g[0] * g[2] + g[1] * g[3]).lower()) for g in geo)
        ok = bool(f > 0 and theta_lo > 0 and jac_min > 0 and sin_max < 0 and rho_min > 0 and xn_min > 0
                  and A.c[0] < 0)
        log(f"  [leaves z0 = {z0} +- {r:g}] F in [{float(f.lower()):.4e}, {float(f.upper()):.4e}]; "
            f"theta_+ >= {theta_lo:.4e} (|R| <= {Rmax:.2e}); Jacobi >= {jac_min:.4f}; x.n >= {xn_min:.4f}; "
            f"sin(alpha) <= {sin_max:.2e}; certified {ok} ({time.time() - T0:.0f} s)")
        return {"z0": z0, "r": r, "F": [_down(f.lower()), _up(f.upper())], "theta_lo": theta_lo, "R_max": Rmax,
                "jacobi_min": jac_min, "x_dot_n_min": xn_min, "sin_alpha_max": sin_max, "rho_min": rho_min,
                "certified": ok,
                "_geo": [(_down(g[0].lower()), _up(g[0].upper()), _down(g[1].lower()), _up(g[1].upper()))
                         for g in geo],
                "_z0ball": (_down(zj.c[0].lower()), _up(zj.c[0].upper()))}
    finally:
        ctx.prec = old


def _certs_for(d, t_abs):
    out = []
    for c in bm.load_certificates():
        if not (d.lower() >= arb(c["d_lo"]) and d.upper() <= arb(c["d_hi"])):
            continue
        tb1, tb2 = bm.tail_bounds(c, t_abs, 1), bm.tail_bounds(c, t_abs, 2)
        if tb1 is not None and tb2 is not None:
            out.append((c, tb1, tb2))
    if not out:
        raise bm.DomainError("no certificate covers this d and t")
    return out


# ---------------------------------------------------------------------------
# Floating-point location of the level surfaces (not part of the proof)
# ---------------------------------------------------------------------------
Z0_FOLD = 0.6326595528058354
_FLOAT = {}


def level_shooter(d, t, h):
    """bl_fold.BL with the level term: floating-point shooting for theta_+^(1) = h (locates z0 only)."""
    key = (d, t, h)
    if key in _FLOAT:
        return _FLOAT[key]
    from scipy.integrate import solve_ivp
    from research.calc.bl_fold import BL

    class BLh(BL):
        def _rhs(self, s, y):
            rho, z, al = y
            extra = self.psi(rho, z) ** 2 * (self.t * self.p(rho, z, al) - h)
            return [math.cos(al), math.sin(al), -math.sin(al) / rho + 4 * self.n_dlogpsi(rho, z, al) + extra]

        def shoot(self, z0):
            k0 = (4 * self.n_dlogpsi(0.0, z0, 0.0) + self.psi(0.0, z0) ** 2 * (self.t * self.p(0.0, z0, 0.0) - h)) / 2
            ds = 1e-7
            y0 = [ds, z0 + k0 * ds**2 / 2, k0 * ds]

            def hit_axis(s, y):
                return y[0] - 1e-10
            hit_axis.terminal, hit_axis.direction = True, -1

            def hit_plane(s, y):
                return y[1]
            hit_plane.terminal, hit_plane.direction = True, -1
            return solve_ivp(self._rhs, (ds, 20.0), y0, method="DOP853", rtol=self.rtol,
                             atol=self.rtol * 1e-2, events=[hit_axis, hit_plane], dense_output=True)

    _FLOAT[key] = BLh(d, t)
    return _FLOAT[key]


def find_levels(d, t, h_in, h_out):
    """Floating-point z0 of the outer branch of theta^(1) = h_out and the inner branch of h_in."""
    from scipy.optimize import brentq, minimize_scalar
    fin = level_shooter(d, t, h_in).mismatch
    res = minimize_scalar(fin, bounds=(Z0_FOLD - 1e-3, Z0_FOLD + 1e-3), method="bounded",
                          options={"xatol": 1e-12})
    if res.fun >= 0:
        raise AssertionError("no inner level surface: h_in too negative or t too small")
    z_in = brentq(fin, Z0_FOLD - 3e-3, res.x, xtol=1e-15)
    z_out = brentq(level_shooter(d, t, h_out).mismatch, res.x, Z0_FOLD + 5e-3, xtol=1e-15)
    return z_in, z_out, res.x, res.fun


def _barrier_job(d, t, h, z0, delta, side):
    return certify_barrier(d, t, h, z0, delta, side, log=lambda m: print(m, flush=True))


def main(d="0.7662", t="1.2e-5", h_in=-0.75e-6, h_out=2e-6, z_in=None, z_out=None, delta=1e-10,
         out_path=os.path.join(DATA, "bl_barrier.json")):
    print(f"Barrier certificate: d = {d}, t = {t}, levels h_in = {h_in}, h_out = {h_out}", flush=True)
    T0 = time.time()
    if z_in is None or z_out is None:
        z_in, z_out, zmin, fmin = find_levels(float(d), float(t), h_in, h_out)
        print(f"  floating-point roots: z_in = {z_in!r}, z_out = {z_out!r} "
              f"(minimum {fmin:.3e} of the h_in mismatch at {zmin:.10f})", flush=True)
    import multiprocessing as mp
    with mp.get_context("fork").Pool(2) as pool:
        jobs = [pool.apply_async(_barrier_job, (d, t, h, z, delta, side))
                for h, z, side in ((h_in, z_in, "in"), (h_out, z_out, "out"))]
        inner, outer = (j.get() for j in jobs)
    # nesting and Omega in U need enclosures much thinner than the Moore ball run's (whose widths grow
    # by wrapping to ~6e-3 at the plane): Lohner sets over the same brackets (bl_foliate)
    from research.calc import bl_foliate as bf
    with mp.get_context("fork").Pool(2) as pool:
        lj = [pool.apply_async(bf._barrier_geometry_job, ((d, t, h, g["z0"], delta),))
              for h, g in ((h_in, inner), (h_out, outer))]
        lin, lout = (j.get() for j in lj)
    for g, lg, name in ((inner, lin, "in"), (outer, lout, "out")):
        print(f"  [{name}] Lohner enclosure: {len(lg['_geo'])} boxes, theta_+ in [{lg['theta_lo']:.6e}, "
              f"{lg['theta_hi']:.6e}], Jacobi >= {lg['jacobi_min']:.4f}, x.n >= {lg['x_dot_n_min']:.4f}", flush=True)
        g["lohner"] = {k: v for k, v in lg.items() if not k.startswith("_")}
        g["_geo"], g["_z0ball"] = lg["_geo"], lg["_z0ball"]
    lohner_ok = lin["theta_hi"] < 0 < lout["theta_lo"] and min(lin["x_dot_n_min"], lout["x_dot_n_min"]) > 0 \
        and min(lin["rho_min"], lout["rho_min"]) > 0
    geo = nesting_and_domain(inner, outer)
    geo["lohner_recheck"] = bool(lohner_ok)
    print(f"  nesting: S_in inside S_out: {geo['nested']} (radial gap >= {geo.get('min_radial_gap', float('nan')):.3e});"
          f" Omega inside U: {geo.get('Omega_in_U')} (inflation needed <= "
          f"{geo.get('largest_inflation_needed', float('nan')):.4f} of {geo.get('w_tube')})", flush=True)
    certified = bool(inner["certified"] and outer["certified"] and geo["nested"] and geo.get("Omega_in_U")
                     and geo["lohner_recheck"])
    print(f"\nRESULT: common MOTS in Sigma_t certified: {certified} ({time.time() - T0:.0f} s)")
    rec = {
        "theorem": "notes/20, Theorem M", "d": d, "t": t, "h_in": h_in, "h_out": h_out, "delta": delta,
        "z_in": z_in, "z_out": z_out,
        "inner": {k: v for k, v in inner.items() if not k.startswith("_")},
        "outer": {k: v for k, v in outer.items() if not k.startswith("_")},
        "geometry": geo, "certified": certified, "precision_bits": 128, "seconds": time.time() - T0,
        "certificates": [c["file"] for c in bm.load_certificates()],
        "replay": "python -m research.calc.bl_barrier",
    }
    if out_path:
        with open(out_path, "w") as f:
            json.dump(rec, f, indent=1)
    return rec


if __name__ == "__main__":
    main()
