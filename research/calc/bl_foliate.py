"""Theorem N (notes/20): no MOTS in an explicit shell about S_c before the merger, by a foliation.

Leaves. For z0 in [Z_A, Z_B] and t in [0, t0], the leaf S_t(z0) is the surface of revolution (and
reflection in z = 0) of the curve leaving the axis point (0, z0) horizontally and solving the
level-h first-order expansion equation theta_+^(1) = h (research/calc/bl_barrier.py, notes/20 §2.2),
up to its first crossing of the plane z = 0. Notes/20 Lemma 3 excludes every MOTS from the shell
swept by the leaves at time t, provided that for every leaf:

  (a) F = alpha_end + pi/2 > 0, i.e. u = rho cos(alpha) > 0 at the crossing: the reflected leaf has
      an outward (convex) equatorial ridge, and is smooth elsewhere;
  (b) theta_+(true) = h + R > 0 on the smooth part (R from bl_barrier.remainder, with the certified CK
      tails at |t| <= t0);
  (c) the leaves are strictly nested and embedded radial graphs: the normal component of dX/dz0 is
      > 0 (Jacobi field), rho > 0, sin(alpha) < 0, x.n > 0.

Method. gap3_fold.py's Lohner sets (sigma-system in (rho, zt = z - d/2, u = rho cos a,
v = rho sin a), d sigma = ds / rho, then z as the clock on the last segment), with jets in z0 (slot 1,
for the Jacobi field), block frames, and z0 as a dummy parameter. Two changes: the source of the
level-h equation,

    rho G = 4 (u psi_z - v psi_rho) / psi + [ t ((u^2 k_xx + 2 u v k_xz + v^2 k_zz) / rho + rho k_yy)
            - h rho psi^4 ] / psi^2,          k = kappa_0 = -Ric(psi^4 delta),

and the time t as a third set parameter (t' = 0), so that one run covers every t in [0, t0].

Status (session 16): NOT COMPLETED. On test pieces, (b) and (c) are certified for all t in [0, t0], but
(a) fails on z0-pieces near the fold: the Lohner remainder A r inflates like (piece width)^2 times
steps (9.6e-6 on a 6e-5 piece, against a true second-order term ~3e-8), so u_end > 0 is not certified
there. A single leaf (z0 a point, t the whole interval) passes. The fix is a convexity argument in
z0 with second jets, F >= F(z_c) - F_z(z_c)^2 / (2 inf F_zz) (notes/20 §4). barrier_geometry (used by
bl_barrier for the nesting test) and step_param_sub are complete and tested.

Run: python -m research.calc.bl_foliate   (would write research/calc/data/bl_foliate.json)
"""

import json
import math
import os
import time

from flint import arb, arb_mat, ctx

from research.calc import bl_barrier as bb
from research.calc import fold_proof as fp
from research.calc import gap3_fold as gf
from research.calc import lohner
from research.calc.fold_proof import J
from research.calc.gap3_interval import QUARTER, _ball

ctx.prec = 128
DATA = os.path.join(os.path.dirname(__file__), "data")
SLOTS = (0, 1)          # base state and its z0-derivative (Jacobi field)
NPAR = 3                # parameters (d, z0, t); z0 is a dummy


def _terms(rho, zt, d):
    """psi and its first and second Cartesian derivatives at (rho, 0, z), zt = z - d/2."""
    psi, pr, pz = 1, 0, 0
    pxx = pyy = pzz = pxz = 0
    for Z in (zt, zt + d):
        ir = 1 / (rho * rho + Z * Z).sqrt()
        ir2 = ir * ir
        ir3 = ir2 * ir
        ir5 = ir3 * ir2
        psi = psi + QUARTER * ir
        pr = pr - QUARTER * rho * ir3
        pz = pz - QUARTER * Z * ir3
        pxx = pxx + QUARTER * (3 * rho * rho * ir5 - ir3)
        pyy = pyy - QUARTER * ir3
        pzz = pzz + QUARTER * (3 * Z * Z * ir5 - ir3)
        pxz = pxz + QUARTER * 3 * rho * Z * ir5
    return psi, pr, pz, pxx, pyy, pzz, pxz


def f_state(rho, zt, u, v, d, t, h):
    """The level-h sigma-system: [rho', zt', u', v'] = [u, v, rho - v P, u P], P = rho G."""
    psi, pr, pz, pxx, pyy, pzz, pxz = _terms(rho, zt, d)
    ip = 1 / psi
    g2 = pr * pr + pz * pz
    kxx = 2 * pxx * ip + (2 * g2 - 6 * pr * pr) * ip * ip
    kyy = 2 * pyy * ip + 2 * g2 * ip * ip
    kzz = 2 * pzz * ip + (2 * g2 - 6 * pz * pz) * ip * ip
    kxz = 2 * pxz * ip - 6 * pr * pz * ip * ip
    p2 = psi * psi
    src = t * ((u * u * kxx + 2 * u * v * kxz + v * v * kzz) / rho + rho * kyy) - h * rho * p2 * p2
    P = 4 * (u * pz - v * pr) * ip + src * ip * ip
    return [u, v, rho - v * P, u * P]


def _unpack(y, width):
    k = len(SLOTS)
    d, z0, t = y[width * k], y[width * k + 1], y[width * k + 2]
    zero = d * 0
    c = [[zero] * 6 for _ in range(width)]
    for a, s in enumerate(SLOTS):
        for i in range(width):
            c[i][s] = y[width * a + i]
    return [J(ci) for ci in c], J([d, zero, zero + 1, zero, zero, zero]), J([t] + [zero] * 5), zero


def make_f_sigma(h):
    def f(y):
        X, dj, t, zero = _unpack(y, 4)
        F = f_state(*X, dj, t, h)
        return [F[i].c[s] for s in SLOTS for i in range(4)] + [zero] * NPAR
    return f


def make_f_tau(h):
    def f(y):
        X, dj, t, zero = _unpack(y, 5)
        F = f_state(*X[:4], dj, t, h)
        c = -X[4] / X[3]
        G = [F[0] * c, -X[4], F[2] * c, F[3] * c]
        return [(G[i].c[s] if i < 4 else zero) for s in SLOTS for i in range(5)] + [zero] * NPAR
    return f


def start(z_lo, z_hi, d, t_lo, t_hi, h):
    """PSet (parameters d, z0, t) at arclength s1 = 2^-26 from (0, z0), z0 in [z_lo, z_hi], t in [0, t0];
    the axis enclosure uses bl_barrier's source with t the whole interval. Returns (set, axis data)."""
    z0 = _ball(z_lo, z_hi)
    tb = _ball(t_lo, t_hi)
    zm = arb(z0.mid())
    one, zero = arb(1), arb(0)
    G = bb.make_source(tb, h)
    zj = J([z0, one, zero, zero, zero, zero])
    dj = J([d, zero, one, zero, zero, zero])
    A, Q, Z = bb.axis_enclosure(zj, dj, fp.S1, G)
    s1 = arb(fp.S1)
    rho = s1 + Q * s1**3
    zt = J([zm, one, zero, zero, zero, zero]) + Z * s1**2 - dj * (one / 2)
    sa, ca = fp.sin_cos(A * s1)
    X = [rho, zt, rho * ca, rho * sa]
    comps = [X[i].c[s] for s in SLOTS for i in range(4)]
    S = lohner.PSet.from_box(comps, [d, z0, tb])
    S.P[1, 0] = -one / 2
    S.P[1, 1] = one
    return S, (A, Q, Z, zj, dj, None)


def to_tau(S):
    """gap3_fold.to_tau for the 3-parameter set (z_sw = zt + d/2 appended to each block)."""
    k, n = len(SLOTS), len(S.m)
    n2 = n + k
    half = arb(1) / 2
    new = [5 * (j // 4) + j % 4 for j in range(n)]
    m, r = [None] * n2, [arb(0)] * n2
    A, P = arb_mat(n2, n2), arb_mat(n2, NPAR)
    for i in range(n):
        m[new[i]], r[new[i]] = S.m[i], S.r[i]
        for j in range(n):
            A[new[i], new[j]] = S.A[i, j]
        for l in range(NPAR):
            P[new[i], l] = S.P[i, l]
    for a, s in enumerate(SLOTS):
        zt, sw = 4 * a + 1, 5 * a + 4
        m[sw] = S.m[zt] + (half * S.p0[0] if s == 0 else 0)
        for j in range(n):
            A[sw, new[j]] = S.A[zt, j]
        A[sw, sw] = 1
        P[sw, 0] = S.P[zt, 0] + (half if s == 0 else 0)
        for l in range(1, NPAR):
            P[sw, l] = S.P[zt, l]
    return lohner.PSet(m, A, r, P, S.p0, S.q)


def z_abs(S):
    n = len(S.m)
    half = arb(1) / 2
    acc = S.m[1] + half * S.p0[0] + sum((S.A[1, j] * S.r[j] for j in range(n)), arb(0))
    acc += (S.P[1, 0] + half) * S.q[0]
    for l in range(1, NPAR):
        acc += S.P[1, l] * S.q[l]
    return acc


class Collector:
    """Checks (b) and (c) on every step box, and the polar boxes of the leaves."""

    def __init__(self, d, t_lo, t_hi, h, keep_geo):
        self.d, self.h = d, arb(h)
        self.tb = _ball(t_lo, t_hi)
        self.G = bb.make_source(self.tb, self.h)
        self.certs = bb._certs_for(d, t_hi)
        self.theta_lo, self.theta_hi, self.R_max = math.inf, -math.inf, 0.0
        self.jac_min = self.xn_min = self.rho_min = math.inf
        self.v_max = -math.inf
        self.geo = [] if keep_geo else None

    def box(self, rho, z, u, v, rho_z, z_z):
        sa, ca = v / rho, u / rho
        if u > 0 and abs(sa) < 1:
            # cos(alpha) >= 0 here, so cos = sqrt(1 - sin^2): much tighter than u / rho when the step
            # hull is long in rho (near the axis); both are valid, keep the intersection
            c2 = (1 - sa * sa).sqrt()
            ca = ca.intersection(c2) if ca.overlaps(c2) else c2
        k2 = v / (rho * rho)
        R, _ = bb.remainder((rho, z, sa, ca, k2), self.d, self.tb, self.h, self.G, self.certs)
        self.theta_lo = min(self.theta_lo, float((self.h + R).lower()))
        self.theta_hi = max(self.theta_hi, float((self.h + R).upper()))
        self.R_max = max(self.R_max, float(abs(R).upper()))
        self.jac_min = min(self.jac_min, float((-v * rho_z + u * z_z).lower()))
        self.xn_min = min(self.xn_min, float((-rho * v + z * u).lower()))
        self.rho_min = min(self.rho_min, float(rho.lower()))
        self.v_max = max(self.v_max, float(v.upper()))
        if self.geo is not None:
            self.geo.append((bb._down(rho.lower()), bb._up(rho.upper()), bb._down(z.lower()), bb._up(z.upper())))

    def axis(self, axis):
        A, Q, Z, zj, dj, _ = axis
        for b in bb.axis_boxes(axis):
            R, _ = bb.remainder(b, self.d, self.tb, self.h, self.G, self.certs)
            self.theta_lo = min(self.theta_lo, float((self.h + R).lower()))
            self.theta_hi = max(self.theta_hi, float((self.h + R).upper()))
            self.R_max = max(self.R_max, float(abs(R).upper()))
        sig = arb(fp.S1 / 2, fp.S1 / 2)
        aa = A.c[0] * sig
        self.jac_min = min(self.jac_min,
                           float((-aa.sin() * Q.c[1] * sig**3 + aa.cos() * (1 + Z.c[1] * sig**2)).lower()))
        if not A.c[0] < 0:
            self.v_max = math.inf


def step_param_sub(f, S, h, L, blocks, m_sub):
    """lohner.step_param (same computation), also returning enclosures of the state on m_sub equal
    sub-intervals [u_i, u_i+1] of [0, h]: the centre trajectory's Taylor polynomial with Lagrange
    remainder over the a-priori box B, plus the Lohner mean-value term at the sub-time (see below)."""
    from flint import arb_mat
    import numpy as np
    n, k, p = len(S.m), len(S.p0), L - 2
    N = n + k
    Y = S.hull() + S.params()
    B = lohner._picard_box(f, Y, h)
    hb = arb(h)
    Tm = lohner.taylor(f, S.m + S.p0, L)
    TB = lohner.taylor(f, B, L)
    xm = []
    for i in range(n):
        acc = arb(0)
        for kk in range(p, -1, -1):
            acc = acc * hb + lohner._coef(Tm[i], kk)
        xm.append(acc + lohner._coef(TB[i], p + 1) * hb ** (p + 1))
    TX = lohner.taylor(f, Y, L, grad=True)
    TBg = lohner.taylor(f, B, L, grad=True)
    Pbox = lohner._phi_box(lohner._jacobian(f, B), h)
    Jm = arb_mat(n, N)
    for i in range(n):
        for j in range(N):
            acc = arb(0)
            for kk in range(p, -1, -1):
                acc = acc * hb + lohner._coef(TX[i].g[j], kk)
            Jm[i, j] = acc
    Rem = arb_mat([[lohner._coef(TBg[i].g[j], p + 1) for j in range(N)] for i in range(n)]) * Pbox
    Jm = Jm + Rem * hb ** (p + 1)
    Jx = arb_mat([[Jm[i, j] for j in range(n)] for i in range(n)])
    Jp = arb_mat([[Jm[i, n + l] for l in range(k)] for i in range(n)])
    C1 = Jx * S.A
    C2 = Jx * S.P + Jp
    Cm = np.array([[float(C1[i, j].mid()) for j in range(n)] for i in range(n)])
    widths = np.array([float(S.r[j].rad()) * np.linalg.norm(Cm[:, j]) for j in range(n)])
    Q = lohner._qr_frame(Cm, widths, blocks)
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
    S_new = lohner.PSet(m_new, A_new, r_new, P_new, S.p0, S.q)
    # sub-step enclosures, in Lohner form at the sub-time u:
    #   x(u) in phi_u(m) + Jx(u) (A r) + (Jx(u) P + Jp(u)) q,
    # Jx, Jp the gradient Taylor series of the flow over the set's hull (the same series as for the step
    # map), evaluated at u, with the Lagrange term through Pbox. The correlations of the set are kept, so
    # quantities like v ~ rho^2 near the axis stay relatively tight.
    RemB = arb_mat([[lohner._coef(TBg[i].g[j], p + 1) for j in range(N)] for i in range(n)]) * Pbox
    subs = []
    for isub in range(m_sub):
        u = arb((isub + 0.5) * h / m_sub, 0.5 * h / m_sub)
        up = u ** (p + 1)
        Ju = arb_mat(n, N)
        for i in range(n):
            for j in range(N):
                acc = arb(0)
                for kk in range(p, -1, -1):
                    acc = acc * u + lohner._coef(TX[i].g[j], kk)
                Ju[i, j] = acc + RemB[i, j] * up
        Jxu = arb_mat([[Ju[i, j] for j in range(n)] for i in range(n)])
        Jpu = arb_mat([[Ju[i, n + l] for l in range(k)] for i in range(n)])
        C1u = Jxu * S.A
        C2u = Jxu * S.P + Jpu
        st = []
        for i in range(n):
            acc = arb(0)
            for kk in range(p, -1, -1):
                acc = acc * u + lohner._coef(Tm[i], kk)
            acc += lohner._coef(TB[i], p + 1) * up
            acc += sum((C1u[i, j] * S.r[j] for j in range(n)), arb(0))
            acc += sum((C2u[i, l] * S.q[l] for l in range(k)), arb(0))
            st.append(acc)
        subs.append((u, st))
    return S_new, B, subs


def _step_sub(f, S, hh, L, blocks, m_sub, what):
    while True:
        try:
            return step_param_sub(f, S, hh, L, blocks, m_sub) + (hh,)
        except (RuntimeError, ValueError, ZeroDivisionError):
            hh /= 2
            if hh < 1e-6:
                raise AssertionError(f"{what}: Lohner step fails")


def run(z_lo, z_hi, d, t0, h, ds=0.004, L=14, hmax=0.25, sigma_max=60.0, z_switch=0.1, h_tau=2.0**-5,
        keep_geo=False, t_lo=0.0, h_axis=2.0**-5, rho_axis=1e-3, m_sub=8):
    """One piece: z0 in [z_lo, z_hi], t in [t_lo, t0]. Returns a dict with the certified bounds (or raises
    AssertionError). Steps in sigma are at most h_axis while rho < rho_axis (tight hulls near the axis)."""
    d = arb(d) if not isinstance(d, arb) else d
    hh_ = arb(h)
    blocks = [list(range(4 * a, 4 * a + 4)) for a in range(len(SLOTS))]
    f = make_f_sigma(hh_)
    S, axis = start(z_lo, z_hi, d, t_lo, t0, h)
    col = Collector(d, t_lo, t0, h, keep_geo)
    col.axis(axis)
    half = arb(1) / 2
    sig, n = 0.0, 0
    while True:
        if sig > sigma_max:
            raise AssertionError("no switch before sigma_max")
        X = S.hull()
        rho_hi = float(X[0].upper())
        hh = min(hmax if rho_hi > rho_axis else h_axis, 2.0 ** math.floor(math.log2(ds / rho_hi)))
        S, B, subs, hh = _step_sub(f, S, hh, L, blocks, m_sub, f"sigma = {sig}")
        if not B[1] + B[8] * half > 0:
            raise AssertionError(f"a step box reaches the plane before the switch (sigma {sig})")
        for _, X in subs:
            col.box(X[0], X[1] + S.p0[0] * half + S.q[0] * half, X[2], X[3], X[4], X[5])
        sig, n = sig + hh, n + 1
        if z_abs(S) < z_switch:
            break
    zsw = z_abs(S)
    if not zsw > 0:
        raise AssertionError("z_sw > 0 not certified")
    S = to_tau(S)
    blocks = [list(range(5 * a, 5 * a + 5)) for a in range(len(SLOTS))]
    f = make_f_tau(hh_)
    tau = 0.0
    while tau < 1:
        S, B, subs, hh = _step_sub(f, S, min(h_tau, 1 - tau), L, blocks, m_sub, f"tau = {tau}")
        # z = z_sw (1 - tau) on the step; its z0-derivative is z_sw,z (1 - tau) = zt,z
        for u, X in subs:
            col.box(X[0], X[4] * (1 - (arb(tau) + u)), X[2], X[3], X[5], X[6])
        tau, n = tau + hh, n + 1
    Y = S.hull()
    u_end, v_end = Y[2], Y[3]
    ok = bool(u_end > 0 and v_end < 0 and col.theta_lo > 0 and col.jac_min > 0 and col.xn_min > 0
              and col.rho_min > 0 and col.v_max < 0)
    out = {"z0": [z_lo, z_hi], "t": [t_lo, t0], "u_end": [float(u_end.lower()), float(u_end.upper())],
           "u_end_lo": float(u_end.lower()), "theta_lo": col.theta_lo, "theta_hi": col.theta_hi, "R_max": col.R_max,
           "jacobi_min": col.jac_min, "x_dot_n_min": col.xn_min, "rho_min": col.rho_min,
           "sin_alpha_max_sign": col.v_max, "steps": n, "certified": ok}
    if keep_geo:
        out["_geo"] = col.geo
        out["_z0ball"] = (bb._down(_ball(z_lo, z_hi).lower()), bb._up(_ball(z_lo, z_hi).upper()))
    return out


def _piece_job(args):
    lo, hi, d, t0, h = args
    t1 = time.time()
    try:
        r = run(lo, hi, d, t0, h)
    except AssertionError as e:
        r = {"z0": [lo, hi], "certified": False, "error": str(e)}
    r["seconds"] = time.time() - t1
    return r


def plan_pieces(z_a, z_b, z_fold, near=1.2e-4, fine=2.0**-14, coarse=1.4e-4):
    """Pieces tiling [z_a, z_b]: width <= fine within `near` of the fold (where F is smallest and
    second-order terms in z0 must stay below it), <= coarse elsewhere."""
    cuts = sorted({z_a, z_b, max(z_a, min(z_b, z_fold - near)), max(z_a, min(z_b, z_fold + near))})
    out = []
    for lo, hi in zip(cuts, cuts[1:]):
        w = fine if (lo >= z_fold - near - 1e-15 and hi <= z_fold + near + 1e-15) else coarse
        n = max(1, math.ceil((hi - lo) / w))
        out += [(lo + (hi - lo) * i / n, lo + (hi - lo) * (i + 1) / n if i < n - 1 else hi) for i in range(n)]
    return out


def certify(z_a, z_b, d="0.7662", t0=4.0e-6, h=1e-7, todo=None, min_width=2.0**-18, processes=3,
            log=print):
    """All pieces (default: plan_pieces is not used; one piece per 2^-12), bisecting failures. Returns
    the list of certified pieces, which tile [z_a, z_b]."""
    from multiprocessing import get_context
    if todo is None:
        n = max(1, math.ceil((z_b - z_a) / 2.0**-12))
        w = (z_b - z_a) / n
        todo = [(z_a + i * w, z_a + (i + 1) * w if i < n - 1 else z_b) for i in range(n)]
    done = []
    with get_context("fork").Pool(processes) as pool:
        while todo:
            nxt = []
            for r in pool.imap(_piece_job, [(lo, hi, d, t0, h) for lo, hi in todo], chunksize=1):
                lo, hi = r["z0"]
                log(f"  [{lo!r}, {hi!r}] certified {r['certified']}  "
                    + (f"u_end >= {r['u_end_lo']:.3e}, theta >= {r['theta_lo']:.3e}, |R| <= {r['R_max']:.2e}, "
                       f"Jacobi >= {r['jacobi_min']:.4f}, x.n >= {r['x_dot_n_min']:.4f}" if "u_end_lo" in r
                       else r.get("error", "")) + f"  ({r['seconds']:.0f} s)")
                if r["certified"]:
                    done.append(r)
                elif hi - lo > min_width:
                    m = (lo + hi) / 2
                    nxt += [(lo, m), (m, hi)]
                else:
                    raise AssertionError(f"foliation fails on [{lo}, {hi}]")
            todo = nxt
    done.sort(key=lambda r: r["z0"][0])
    for a, b in zip(done, done[1:]):
        assert a["z0"][1] == b["z0"][0], "pieces do not tile [z_a, z_b]"
    return done


def barrier_geometry(d, t, h, z0, delta, m_sub=32):
    """Fine Lohner enclosure of every barrier curve with z0' in [z0 - delta, z0 + delta] at the single time
    t (level h): no wrapping blow-up, so the boxes are short (arclength <= ds) and thin. Also re-checks
    the barrier inequality on its own boxes (theta_lo, theta_hi), the Jacobi field, x.n and rho."""
    r = run(z0 - delta, z0 + delta, d, float(t), float(h), keep_geo=True, t_lo=float(t), m_sub=m_sub)
    return r


def _barrier_geometry_job(args):
    return barrier_geometry(*args)


def _extreme_job(args):
    z0, d, t0, h, m_sub = args
    return run(z0, z0, d, t0, h, keep_geo=True, m_sub=m_sub)


def _profile_job(args):
    rec, side = args
    g = rec["inner"] if side == "in" else rec["outer"]
    h = rec["h_in"] if side == "in" else rec["h_out"]
    return barrier_geometry(rec["d"], float(rec["t"]), h, g["z0"], rec["delta"])


def main(d="0.7662", t0=5.0e-6, h=5e-8, margin=2e-4, m_extreme=32,
         barrier_record=os.path.join(DATA, "bl_barrier_near.json"), out_path=os.path.join(DATA, "bl_foliate.json")):
    """Theorem N: the foliation of [z_a, z_b] = [z_in - margin, z_out + margin] for t in [0, t0], and
    W_t contains Omega (the region of the Theorem M record barrier_record): S_t(z_a) inside S_in and S_out
    inside S_t(z_b), by radial comparison of Lohner sub-step boxes. All runs go through one pool."""
    from multiprocessing import get_context
    with open(barrier_record) as f:
        brec = json.load(f)
    assert brec["d"] == d and brec["certified"], "needs a certified Theorem M record at the same d"
    z_in, z_out = brec["inner"]["z0"], brec["outer"]["z0"]
    z_a, z_b = z_in - margin, z_out + margin
    # the fold of the h-level family at t0 (floating point, only to place the fine pieces)
    from scipy.optimize import minimize_scalar
    z_fold = minimize_scalar(bb.level_shooter(float(d), t0, h).mismatch, bounds=(bb.Z0_FOLD - 1e-3, bb.Z0_FOLD + 1e-3),
                             method="bounded", options={"xatol": 1e-10}).x
    todo = plan_pieces(z_a, z_b, z_fold)
    print(f"Foliation certificate: d = {d}, t in [0, {t0}], level h = {h}, z0 in [{z_a!r}, {z_b!r}], "
          f"{len(todo)} pieces (fold of the level family at {z_fold:.9f})", flush=True)
    T0 = time.time()
    log = lambda m: print(m, flush=True)  # noqa: E731
    with get_context("fork").Pool(4) as pool:
        side = [pool.apply_async(_extreme_job, ((z_a, d, t0, h, m_extreme),)),
                pool.apply_async(_extreme_job, ((z_b, d, t0, h, m_extreme),)),
                pool.apply_async(_profile_job, ((brec, "in"),)),
                pool.apply_async(_profile_job, ((brec, "out"),))]
        pieces, nxt = [], []
        for r in pool.imap_unordered(_piece_job, [(lo, hi, d, t0, h) for lo, hi in todo]):
            lo, hi = r["z0"]
            log(f"  [{lo!r}, {hi!r}] certified {r['certified']}  "
                + (f"u_end >= {r['u_end_lo']:.3e}, theta >= {r['theta_lo']:.3e}, |R| <= {r['R_max']:.2e}, "
                   f"Jacobi >= {r['jacobi_min']:.3e}, x.n >= {r['x_dot_n_min']:.3e}" if "u_end_lo" in r
                   else r.get("error", "")) + f"  ({r['seconds']:.0f} s)")
            (pieces if r["certified"] else nxt).append(r)
        la, lb, pin, pout = (j.get() for j in side)
    if nxt:
        log(f"  {len(nxt)} pieces failed; bisecting")
        pieces += certify(z_a, z_b, d, t0, h, todo=[tuple(r["z0"]) for r in nxt], log=log)
    pieces.sort(key=lambda r: r["z0"][0])
    tiled = pieces[0]["z0"][0] == z_a and pieces[-1]["z0"][1] == z_b and all(
        a["z0"][1] == b["z0"][0] for a, b in zip(pieces, pieces[1:]))
    gap_a = bb.radial_gap(la, pin)      # S_t(z_a) inside S_in
    gap_b = bb.radial_gap(pout, lb)     # S_out inside S_t(z_b)
    leaves_ok = la["certified"] and lb["certified"]
    contains = bool(leaves_ok and gap_a > 0 and gap_b > 0)
    print(f"  pieces tile [z_a, z_b]: {tiled}; extreme leaves certified: {leaves_ok}; W_t contains Omega: "
          f"{contains} (radial gaps {gap_a:.3e}, {gap_b:.3e}); barrier re-check: S_in theta <= {pin['theta_hi']:.3e}, "
          f"S_out theta >= {pout['theta_lo']:.3e}", flush=True)
    rec = {"theorem": "notes/20, Theorem N", "d": d, "t0": t0, "h": h, "z_a": z_a, "z_b": z_b, "z_fold": z_fold,
           "pieces": pieces, "pieces_tile": bool(tiled),
           "extreme_leaves": [{k: v for k, v in x.items() if not k.startswith("_")} for x in (la, lb)],
           "gap_leaf_a_to_S_in": gap_a, "gap_S_out_to_leaf_b": gap_b, "W_contains_Omega": contains,
           "certified": bool(all(p["certified"] for p in pieces) and tiled and contains),
           "barrier_record": os.path.basename(barrier_record), "seconds": time.time() - T0,
           "replay": "python -m research.calc.bl_foliate"}
    print(f"RESULT: Theorem N certified: {rec['certified']} ({rec['seconds']:.0f} s)", flush=True)
    if out_path:
        with open(out_path, "w") as f:
            json.dump(rec, f, indent=1)
    return rec


if __name__ == "__main__":
    main()
