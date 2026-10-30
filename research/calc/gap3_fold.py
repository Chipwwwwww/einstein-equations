"""Gap 3 near the fold: the local lemma that takes Theorem F' down to d_c.

Let F(z0, d) = alpha + pi/2 at the first crossing of the plane z = 0 by the curve leaving the
axis point (0, z0) (the shooting function of research/calc/fold_proof.py at t = 0; the profile
of a common minimal surface has F = 0). Theorem C (notes/09) certifies the fold: a unique point
(z*, d*) with F = F_z = 0 in a box of radius 1e-12, z* = 0.632659552806 +- 3e-13,
d* = d_c = 0.766197428240 +- 1.1e-13.

Local lemma. Let Z = [z_lo, z_hi] contain z*, and d_1 > d_c. Suppose that for every z0 in Z and
d in [D_C_LO, d_1] (D_C_LO <= d*) the curve reaches the plane with v < 0, and
    (i)  F_zz > 0,    (ii) F_d > 0.
Then F(z0, d*) = (1/2) F_zz(xi, d*) (z0 - z*)^2 >= 0 on Z (Taylor at the fold point, where
F = F_z = 0), and F(z0, d) = F(z0, d*) + int_{d*}^{d} F_d > 0 for every d in (d*, d_1]. So for
those d no curve from Z meets the plane at a right angle, and none is the profile of a common
minimal surface. gap3_interval's slab sweep covers the z0 outside Z (it leaves out exactly the
window HOLE, whose z0-image fold_window() contains).

Enclosures. The sigma-system in the relative height zt = z - d/2 (gap3_interval.f_rel) is
carried as jets in (z0, d) (fold_proof.J: slots 1, z0, d, t, z0^2 [holds x_zz / 2], z0 d), and
the state with the jet slots in use goes through Lohner's method (lohner.step_param). Three
choices make the enclosures tight enough (session 9; each was needed):
  * Block-diagonal frame, one block per jet slot. The system is triangular (the base state
    does not see the jets); a full QR frame lets the large jet widths leak into the base state
    and the set blows up at sigma ~ 15.
  * z0 is a second (dummy) parameter of the PSet next to d. Lohner's method then carries the
    first-order dependence on z0 in the column P[:, z0], exactly correlated across the jet
    blocks, and r keeps only the second-order remainder. (With z0 in r, the blocks lose that
    correlation and F_zz comes out as +-300.)
  * Near the plane the independent variable becomes tau in [0, 1] with z = z_sw (1 - tau), as in
    fold_proof's last segment: z_sw (the height at the switch, a constant of each curve, with
    its jets) joins the state, and dx/dtau = -z_sw f(x) / v. Every curve reaches z = 0 exactly
    at tau = 1, where the jets are the derivatives of the crossing state:
        F_z = alpha_x . x_z,   F_zz = alpha_xx(x_z, x_z) + alpha_x . x_zz,   F_d = alpha_x . x_d,
    alpha = atan2(v, u), evaluated in mean-value form on the final Lohner set. (Evaluating the
    crossing formulas with the time corrections s', s'' on the sigma-steps that straddle the
    plane gives a union over a window of steps along the flow: F_d +-0.8, F_zz +-400.)
The reparametrisation needs v < 0 on the tau segment: 1/v appears in the right-hand side, so
the Lohner steps fail if v can vanish, and z_sw > 0 is checked at the switch. Before it, every
step box has z > 0, and z = z_sw (1 - tau) > 0 for tau < 1: tau = 1 is the first crossing.

Run: python -m research.calc.gap3_fold [processes]   (the lemma on fold_window(), ~2 CPU-hours)
"""

import math
import sys
import time

from flint import arb

from research.calc import fold_proof as fp
from research.calc import lohner
from research.calc.fold_proof import J
from research.calc.gap3_interval import QUARTER, _ball

# Theorem C (notes/09): the fold point, M = 1
D_C_LO, D_C_HI = 0.76619742823990, 0.76619742824010
Z_C = 0.632659552806
Z_C_LO, Z_C_HI = 0.6326595528055, 0.6326595528065


def _psi_terms(rho, zt, d):
    psi, pr, pz = 1, 0, 0
    for Z in (zt, zt + d):
        ir = 1 / (rho * rho + Z * Z).sqrt()
        ir3 = ir * ir * ir
        psi = psi + QUARTER * ir
        pr = pr - QUARTER * rho * ir3
        pz = pz - QUARTER * Z * ir3
    return psi, pr, pz


def f_state(rho, zt, u, v, d):
    """gap3_interval.f_rel without the parameter component; works on arb, D and J."""
    psi, pr, pz = _psi_terms(rho, zt, d)
    P = 4 * (u * pz - v * pr) / psi
    return [u, v, rho - v * P, u * P]


def _jets(y, slots, width):
    """Unpack the state into J's: `width` components per slot (slot 0 first), then (d, z0)."""
    k = len(slots)
    d = y[width * k]
    zero = d * 0
    c = [[zero] * 6 for _ in range(width)]
    for a, s in enumerate(slots):
        for i in range(width):
            c[i][s] = y[width * a + i]
    return [J(ci) for ci in c], J([d, zero, zero + 1, zero, zero, zero]), zero


def make_f_sigma(slots):
    """sigma-system on (rho, zt, u, v) per jet slot in `slots`, then the parameters d and z0
    (z0 is a dummy: f does not depend on it, see start). Unused slots are set to zero. That is
    exact as long as no slot in use is fed by one that is not: products feed slot 4 (z0^2)
    from slot 1 and slot 5 (z0 d) from slots 1, 2, and nothing else mixes, so slot 4 needs
    slot 1 and any other choice containing 0 is fine."""
    def f(y):
        X, dj, zero = _jets(y, slots, 4)
        F = f_state(*X, dj)
        return [F[i].c[s] for s in slots for i in range(4)] + [zero, zero]
    return f


def make_f_tau(slots):
    """tau-system on (rho, zt, u, v, z_sw) per jet slot: dx/dtau = -z_sw f(x) / v, except
    dzt/dtau = -z_sw exactly (z = z_sw (1 - tau)); z_sw is constant."""
    def f(y):
        X, dj, zero = _jets(y, slots, 5)
        F = f_state(*X[:4], dj)
        c = -X[4] / X[3]
        G = [F[0] * c, -X[4], F[2] * c, F[3] * c]
        return [(G[i].c[s] if i < 4 else zero) for s in slots for i in range(5)] + [zero, zero]
    return f


def start(z_lo, z_hi, d_lo, d_hi, slots):
    """PSet with parameters (d, z0) of the jet states at arclength s1 = 2^-26 from (0, z0),
    z0 in [z_lo, z_hi], d in [d_lo, d_hi], from fold_proof.axis_enclosure. The state depends on
    z0 and d at first order only through zt = z0 - d/2 + Z s1^2; those terms go into the
    parameter columns P (entries 1 and -1/2) and the rest into r."""
    z0, d = _ball(z_lo, z_hi), _ball(d_lo, d_hi)
    zm, dm = arb(z0.mid()), arb(d.mid())
    one, zero = arb(1), arb(0)
    zj = J([z0, one, zero, zero, zero, zero])
    A, Q, Z = fp.axis_enclosure(zj, J([d, zero, one, zero, zero, zero]), fp.S1)
    s1 = arb(fp.S1)
    rho = s1 + Q * s1**3
    zt = J([zm, one, zero, zero, zero, zero]) + Z * s1**2 \
        - J([dm, zero, one, zero, zero, zero]) * (one / 2)
    sa, ca = fp.sin_cos(A * s1)
    X = [rho, zt, rho * ca, rho * sa]
    comps = [X[i].c[s] for s in slots for i in range(4)]
    S = lohner.PSet.from_box(comps, [d, z0])
    S.P[1, 0] = -one / 2
    S.P[1, 1] = one
    return S


def z_abs(S):
    """z = zt + d/2 over a sigma-phase set, from its parts (keeps the correlations)."""
    n = len(S.m)
    half = arb(1) / 2
    acc = S.m[1] + half * S.p0[0] + sum((S.A[1, j] * S.r[j] for j in range(n)), arb(0))
    return acc + (S.P[1, 0] + half) * S.q[0] + S.P[1, 1] * S.q[1]


def to_tau(S, slots):
    """The sigma-phase set with z_sw = z = zt + d/2 appended to every jet block (as the 5th
    component; its slot-0 part is zt + d/2, slot 2 gets + 1/2, the others copy zt). z_sw is an
    affine function of the set's coordinates, so the new set is m' + A' r' + P' q with one new
    frame column per block and r' = 0 there. Rows and columns are renumbered alike (old index
    4 a + i -> 5 a + i, the new one 5 a + 4), so the block-diagonal frames stay aligned."""
    from flint import arb_mat
    k, n = len(slots), len(S.m)
    n2 = n + k
    half = arb(1) / 2
    new = [5 * (j // 4) + j % 4 for j in range(n)]
    m, r = [None] * n2, [arb(0)] * n2
    A, P = arb_mat(n2, n2), arb_mat(n2, 2)
    for i in range(n):
        m[new[i]], r[new[i]] = S.m[i], S.r[i]
        for j in range(n):
            A[new[i], new[j]] = S.A[i, j]
        for l in range(2):
            P[new[i], l] = S.P[i, l]
    for a, s in enumerate(slots):
        zt, sw = 4 * a + 1, 5 * a + 4
        m[sw] = S.m[zt] + (half * S.p0[0] if s == 0 else 0) + (half if s == 2 else 0)
        for j in range(n):
            A[sw, new[j]] = S.A[zt, j]
        A[sw, sw] = 1
        P[sw, 0] = S.P[zt, 0] + (half if s == 0 else 0)
        P[sw, 1] = S.P[zt, 1]
    return lohner.PSet(m, A, r, P, S.p0, S.q)


def _alpha_parts(u, v):
    r2 = u * u + v * v
    r4 = r2 * r2
    return -v / r2, u / r2, 2 * u * v / r4, (v * v - u * u) / r4, -2 * u * v / r4


def crossing(y, slots):
    """F_z, F_zz, F_d (those the slots allow, else None) from a tau = 1 state y."""
    jet = {s: [y[5 * a + i] for i in range(4)] for a, s in enumerate(slots)}
    u, v = jet[0][2], jet[0][3]
    au, av, auu, auv, avv = _alpha_parts(u, v)
    Fz = Fzz = Fd = None
    if 1 in jet:
        uz, vz = jet[1][2], jet[1][3]
        Fz = au * uz + av * vz
        if 4 in jet:
            uzz, vzz = 2 * jet[4][2], 2 * jet[4][3]
            Fzz = auu * uz * uz + 2 * auv * uz * vz + avv * vz * vz + au * uzz + av * vzz
    if 2 in jet:
        Fd = au * jet[2][2] + av * jet[2][3]
    return Fz, Fzz, Fd


def crossing_enclosure(S, slots):
    """crossing() over the tau = 1 set S = m + A r + P q in mean-value form: with c the point
    m + p0 and Y the box hull(S) united with c (it contains the segment from c to every point),
        G(y) in G(c) + grad G(Y) . (A r + P q, q)."""
    n, k = len(S.m), len(S.p0)
    c = S.m + S.p0
    Y = [a.union(b) for a, b in zip(S.hull() + S.params(), c)]
    seeded = [lohner.D(b, [arb(1 if i == j else 0) for j in range(n + k)]) for i, b in enumerate(Y)]
    out = []
    for gc, gY in zip(crossing(c, slots), crossing(seeded, slots)):
        if gc is None:
            out.append(None)
            continue
        g = gY.g
        acc = gc + sum((g[n + l] * S.q[l] for l in range(k)), arb(0))
        for j in range(n):
            acc += sum((g[i] * S.A[i, j] for i in range(n)), arb(0)) * S.r[j]
        for l in range(k):
            acc += sum((g[i] * S.P[i, l] for i in range(n)), arb(0)) * S.q[l]
        out.append(acc)
    return out


def _step(f, S, hh, L, blocks, what):
    while True:
        try:
            return lohner.step_param(f, S, hh, L, blocks) + (hh,)
        except (RuntimeError, ValueError, ZeroDivisionError):
            hh /= 2
            if hh < 1e-6:
                raise AssertionError(f"{what}: Lohner step fails")


def run(z_lo, z_hi, d_lo, d_hi, slots, h=0.25, ds=0.004, L=14, sigma_max=60.0, z_switch=0.1,
        h_tau=2.0**-5, log=None):
    """Curves from z0 in [z_lo, z_hi], d in [d_lo, d_hi] to the plane z = 0: sigma-steps until the
    set lies below z_switch (every step box with z > 0), then tau-steps to tau = 1. Returns
    (F_z, F_zz, F_d, steps) at the crossing (None where the slots do not allow), or raises
    AssertionError."""
    blocks = [list(range(4 * a, 4 * a + 4)) for a in range(len(slots))]
    f = make_f_sigma(slots)
    S = start(z_lo, z_hi, d_lo, d_hi, slots)
    half = arb(1) / 2
    sig, n = 0.0, 0
    while True:
        if sig > sigma_max:
            raise AssertionError("no switch before sigma_max")
        X = S.hull()
        hh = min(h, 2.0 ** math.floor(math.log2(ds / float(X[0].upper()))))
        S, B, hh = _step(f, S, hh, L, blocks, f"sigma = {sig}")
        if not B[1] + B[4 * len(slots)] * half > 0:
            raise AssertionError(f"a step box reaches the plane before the switch (sigma {sig})")
        sig, n = sig + hh, n + 1
        z = z_abs(S)
        if log and n % 25 == 0:
            log(f"  sigma {sig:.3f}  rho {S.hull()[0].str(5)}  z {z.str(5)}")
        if z < z_switch:
            break
    if not z > 0:
        raise AssertionError("z_sw > 0 not certified")
    S = to_tau(S, slots)
    blocks = [list(range(5 * a, 5 * a + 5)) for a in range(len(slots))]
    f = make_f_tau(slots)
    tau = 0.0
    while tau < 1:
        S, B, hh = _step(f, S, min(h_tau, 1 - tau), L, blocks, f"tau = {tau}")
        tau, n = tau + hh, n + 1
    Fz, Fzz, Fd = crossing_enclosure(S, slots)
    return Fz, Fzz, Fd, n


# ---------------------------------------------------------------------------
# Driver: the local lemma on the whole fold window that the slab sweep leaves out
# (gap3_interval.HOLE, image gap3_interval.hole_z0() = [0.629036.., 0.636869..] in z0).
# ---------------------------------------------------------------------------
PIECE = 2.0**-12
SLOTS = (0, 1, 2, 4)  # F_z, F_zz, F_d in one run


def fold_window():
    """Dyadic z0-window [k PIECE, l PIECE] containing hole_z0() and the fold point."""
    from research.calc.gap3_interval import hole_z0
    lo, hi = hole_z0()
    lo, hi = min(lo, Z_C_LO), max(hi, Z_C_HI)
    return math.floor(lo / PIECE) * PIECE, math.ceil(hi / PIECE) * PIECE


def check_piece(piece):
    """(piece, F_zz lower bound, F_d lower bound, steps or error, seconds) on z0 in piece,
    d in [D_C_LO, D1]; the bounds are None if the run fails."""
    from research.calc.gap3_interval import D1
    lo, hi = piece
    t0 = time.time()
    try:
        _, Fzz, Fd, n = run(lo, hi, D_C_LO, D1, SLOTS)
    except AssertionError as e:
        return piece, None, None, str(e), time.time() - t0
    return piece, float(Fzz.lower()), float(Fd.lower()), n, time.time() - t0


def certify_window(processes=4, min_width=2.0**-16, log=print):
    """Conditions (i) and (ii) of the local lemma on fold_window(), bisecting failed pieces.
    Returns (pieces, worst F_zz lower bound, worst F_d lower bound)."""
    from multiprocessing import Pool
    a, b = fold_window()
    todo = [(a + i * PIECE, a + (i + 1) * PIECE) for i in range(round((b - a) / PIECE))]
    done, wzz, wd = [], math.inf, math.inf
    with Pool(processes) as pool:
        while todo:
            nxt = []
            for piece, fzz, fd, n, sec in pool.imap(check_piece, todo, chunksize=1):
                ok = fzz is not None and fzz > 0 and fd > 0
                if log:
                    log(f"  {piece[0]!r} {piece[1]!r}  F_zz >= {fzz}  F_d >= {fd}  steps {n}  "
                        f"{sec:.0f} s  {'ok' if ok else 'SPLIT'}")
                if ok:
                    done.append(piece)
                    wzz, wd = min(wzz, fzz), min(wd, fd)
                elif piece[1] - piece[0] > min_width:
                    m = (piece[0] + piece[1]) / 2
                    nxt += [(piece[0], m), (m, piece[1])]
                else:
                    raise AssertionError(f"fold lemma fails on {piece}")
            todo = nxt
    return sorted(done), wzz, wd


if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t0 = time.time()
    pieces, wzz, wd = certify_window(procs, log=lambda m: print(m, flush=True))
    a, b = fold_window()
    print(f"Fold lemma: z0 in [{a!r}, {b!r}], d in [{D_C_LO!r}, D1]: F_zz >= {wzz:.4g}, "
          f"F_d >= {wd:.4g} on {len(pieces)} pieces ({time.time() - t0:.0f} s)")
