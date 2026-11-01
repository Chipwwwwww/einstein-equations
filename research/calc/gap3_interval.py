"""Gap 3 over an interval of separations: no common minimal surface of equal-mass Brill-Lindquist
data for every d in (d_c, D2] (notes/12, Theorem F'), certified slab by slab upward from d_c.

Same proof as Theorem F (research/calc/gap3_sweep.py, notes/12), with the separation d as a
parameter everywhere:

  * barriers: certify_throat / certify_far of research/calc/bl_global.py with d an interval;
    over [0.765625, 0.78125] the throat radius must shrink to s_* = 0.15 (0.155 fails at
    d = 0.766), the far radius R_* = 1.05 still works;
  * thin necks, radius <= LAM0_INT = 3 * 2^-8 (approach lemma, event N): research/calc/gap3_neck.py,
    certify_region_rel, in the coordinates (lam, zeta = z_w - d/2) for all d at once;
  * medium necks (event W): the curve is integrated to its waist; the waist box
    (lam, zeta, direction) is recorded and the thin-neck sign pattern is certified on it
    afterwards (research/calc/gap3_neck_lohner.py, certify_rel, on dyadic cells, for all
    d in [D_LO, D_HI] at once);
  * the fold window FOLD_WINDOW (z0 near the fold point) for d <= D_HOLE: left to the local fold
    lemma, research/calc/gap3_fold.py;
  * the sweep: Lohner sets with d kept out of the moving frame (lohner.PSet). The initial sets
    are band-aligned, z0 = Z0_REF + delta + SLOPE (d - D_REF), so that a (delta, d)-box follows
    the top of the individual horizon as d changes (z_top(d) has slope ~0.543). Box sizes come
    from the calibrated table BANDS.

Events (whole enclosure):
  E1  z < 0, every step that allowed z = 0 having u != 0 (steps are halved until it does);
  E2  outside B(0, R_*);
  E3  inside B(p_+-, s_*), with z - d/2 and z + d/2 taken from the set itself (not from hulls);
  N   approach lemma, as in gap3_sweep, in the variable zeta = z - d/2;
  W   u < 0 on the whole set after some step and u > 0 on it after a later one: every curve has
      its waist (u = 0) in one of the step boxes in between, with lam = rho_w and zeta_w in their
      unions, direction sign(v) over them. Excluded once the neck lemma is certified there.

Run: python -m research.calc.gap3_interval 4            (sweep, resumes from gap3_interval.log)
     python -m research.calc.gap3_interval progress LOG (largest d_k with [D0, d_k] covered)
     python -m research.calc.gap3_interval verify LOG [d_k]  (cover, holes, band certificate)
"""

import math
import time

from flint import arb, ctx

from research.calc import fold_proof as fp
from research.calc import lohner
from research.calc.fold_proof import J
from research.calc.gap3_neck import D_HI, D_LO, LAM0_INT, ZE_HI, ZE_LO
from research.calc.gap3_sweep import QUARTER

ctx.prec = 128  # fold_proof sets 256; enclosures here are ~1e-12 at best, 128 bits is ample

S_STAR, R_STAR = arb("0.15"), arb("1.05")
RHO_N = arb("0.02")
Z0_REF, D_REF, SLOPE = 0.5896559, 0.78, 0.543  # z_top(0.78) and its slope in d (floating point)
GRID = 2.0**-10  # starting cells of the medium-neck certificate in (lam, zeta)
GRID_MIN = 2.0**-15  # cells that fail are quartered down to this size
LAM_MAX = 0.05  # largest waist radius handed to the medium-neck certificate


def _ball(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)


def k_bound_rel(rho_e=RHO_N, zlo=ZE_LO, zhi=ZE_HI, d_lo=D_LO, d_hi=D_HI, pieces=64):
    """Upper bound for 4 (|psi_z| + |psi_rho|) / psi on rho in [0, rho_e], zeta = z - d/2 in
    [zlo, zhi], d in [d_lo, d_hi]."""
    rho = arb(rho_e) / 2 + arb(0, 1) * (arb(rho_e) / 2)
    d = _ball(d_lo, d_hi)
    best = arb(0)
    for k in range(pieces):
        a, b = zlo + (zhi - zlo) * k / pieces, zlo + (zhi - zlo) * (k + 1) / pieces
        ze = _ball(a, b)
        psi, pr, pz = arb(1), arb(0), arb(0)
        for Z in (ze, ze + d):
            ir = 1 / (rho * rho + Z * Z).sqrt()
            psi += QUARTER * ir
            pr -= QUARTER * rho * ir**3
            pz -= QUARTER * Z * ir**3
        best = max(best, (4 * (abs(pz) + abs(pr)) / psi).upper())
    return best


K_REL = k_bound_rel()


def f_rel(y):
    """sigma-system in the variables (rho, zt, u, v), zt = z - d/2 (height above the upper
    puncture), with the parameter d appended (d' = 0). Then d enters only through the lower
    puncture at zt = -d, which keeps a box in d from decorrelating z and d inside psi."""
    rho, zt, u, v, d = y
    psi, pr, pz = 1, 0, 0
    for Z in (zt, zt + d):
        ir = 1 / (rho * rho + Z * Z).sqrt()
        ir3 = ir * ir * ir
        psi = psi + QUARTER * ir
        pr = pr - QUARTER * rho * ir3
        pz = pz - QUARTER * Z * ir3
    P = 4 * (u * pz - v * pr) / psi
    return [u, v, rho - v * P, u * P, d * 0]


def start_pset(delta_lo, delta_hi, d_lo, d_hi):
    """Lohner set (lohner.PSet, parameter d) of the states at arclength s1 = 2^-26 from the axis
    points (0, z0), z0 = Z0_REF + delta + SLOPE (d - D_REF), for all delta, d in the box.
    The state is (s1 + Q s1^3, z0 + Z s1^2, rho cos(A s1), rho sin(A s1)) with A, Q, Z enclosed
    over the box (fold_proof.axis_enclosure); its only O(1) dependence on d is z0's, which goes
    into P (slope SLOPE), everything else into r."""
    s1 = arb(fp.S1)

    def state(delta, d):
        z0 = arb(Z0_REF) + delta + arb(SLOPE) * (d - arb(D_REF))
        zero = [arb(0)] * 5
        A, Q, Z = fp.axis_enclosure(J([z0] + zero), J([d] + zero), fp.S1)
        A, Q, Z = A.c[0], Q.c[0], Z.c[0]
        rho = s1 + Q * s1**3
        a = A * s1
        return [rho, Z * s1**2, rho * a.cos(), rho * a.sin()]  # z minus z0 in slot 1
    db, ddb = _ball(delta_lo, delta_hi), _ball(d_lo, d_hi)
    dm, ddm = arb(db.mid()), arb(ddb.mid())
    centre, box = state(dm, ddm), state(db, ddb)
    # zt = z - d/2 = Z0_REF + delta + SLOPE (d - D_REF) - d/2 + Z s1^2: linear in d with slope
    # SLOPE - 1/2, which goes into P; the rest into r
    half = arb(1) / 2
    zt0c = arb(Z0_REF) + dm + arb(SLOPE) * (ddm - arb(D_REF)) - ddm * half
    m = [arb(centre[0].mid()), arb((zt0c + centre[1]).mid()), arb(centre[2].mid()),
         arb(centre[3].mid())]
    r = [box[0] - m[0],
         (arb(Z0_REF) + db + arb(SLOPE) * (ddm - arb(D_REF)) - ddm * half) + box[1] - m[1],
         box[2] - m[2], box[3] - m[3]]
    S = lohner.PSet.from_box(m, [ddb])
    S.r = r
    S.P[1, 0] = arb(SLOPE) - half
    return S


def zt_plus(S, c=0):
    """zt + c d over the set, from its parts (keeps the correlation of zt with d). With
    zt = z - d/2: c = 0 is the height above the upper puncture, c = 1/2 the height z, c = 1 the
    height above the lower puncture."""
    n = len(S.m)
    c = arb(c)
    acc = S.m[1] + c * S.p0[0]
    acc += sum((S.A[1, j] * S.r[j] for j in range(n)), arb(0))
    acc += (S.P[1, 0] + c) * S.q[0]
    return acc


def neck_event(X, ze):
    """Approach lemma (gap3_sweep.neck_event) with zeta = z - d/2 in place of z."""
    rho, _, u, v = X[:4]
    if not (u < 0 and rho < RHO_N):
        return False
    re, K = rho.upper(), K_REL
    c = 1 - K * re - K * K * re * re / 4
    if not c > 0:
        return False
    dz = (2 * re / arb(1).exp() + K * re * re / 2) / c.sqrt()
    return bool(ze - dz > ZE_LO and ze + dz < ZE_HI and abs(v) + K * re * re / 2 < LAM0_INT)


def run_box(box, h=0.25, ds=0.004, L=14, sigma_max=80.0, max_steps=60000):
    """box = (delta_lo, delta_hi, d_lo, d_hi). Returns (box, event, sigma, steps, seconds, waist)
    with waist = (lam_lo, lam_hi, zeta_lo, zeta_hi, s) for event W, else None."""
    t0 = time.time()
    try:
        S = start_pset(*box)
    except (RuntimeError, ValueError, ZeroDivisionError) as e:
        return box, f"error_{type(e).__name__}", 0.0, 0, time.time() - t0, None
    sig, n = 0.0, 0
    half = arb(1) / 2
    pre = None  # (zt at the last step end with u < 0 on the whole set, B_rho union, B_v union, H)
    while sig < sigma_max and n < max_steps:
        X = S.hull()
        hh = min(h, 2.0 ** math.floor(math.log2(ds / float(X[0].upper()))))
        if sig > 14 and float(X[0].upper()) < 0.06:
            hh = min(hh, 2.0**-5)  # at a neck: long steps inflate the set
        if X[2] < 0 and float(X[2].upper()) > -0.3 * float(X[0].lower()):
            hh = min(hh, 2.0**-6)  # near a waist: short steps keep the waist box small
        while True:
            try:
                S_new, B = lohner.step_param(f_rel, S, hh, L)
            except (RuntimeError, ValueError, ZeroDivisionError):
                S_new = None
            # z = zt + d/2 on the step box: accept if z > 0 there, or u has one sign there
            if S_new is not None and (B[1] + B[4] * half > 0 or B[2] > 0 or B[2] < 0):
                break
            hh /= 2
            if hh < 1e-7:
                ev = "stuck" if S_new is None else "plane-undecided"
                return box, ev, sig, n, time.time() - t0, None
        Xn = S_new.hull()
        # Waist window: from the last step end with u < 0 on the whole set to the first with
        # u > 0 on it, every curve passes u = 0 (a waist) at a state inside one of the step
        # boxes B in between: rho_w in their rho-union, z_w - d/2 in zt(start) + [0, H] v-union.
        # During a step, rho = rho(start) + [0, hh] u and zt = zt(start) + [0, hh] v, with u, v
        # in the step box B; the waist lies in one of the window's steps.
        step_zt = zt_plus(S) + arb(hh / 2, hh / 2) * B[3]
        step_rho = X[0].union(Xn[0]) + arb(0, 1) * (hh * abs(B[2]).upper())
        if X[2] < 0:
            pre = (step_zt, step_rho, B[3])
        elif pre is not None:
            pre = (pre[0].union(step_zt), pre[1].union(step_rho), pre[2].union(B[3]))
        if pre is not None and Xn[2] > 0:
            ze, brho, bv = pre
            if (bv > 0 or bv < 0) and brho.upper() < LAM_MAX:
                s = 1 if bv > 0 else -1
                waist = (float(brho.lower()), float(brho.upper()), float(ze.lower()),
                         float(ze.upper()), s)
                return box, "W", sig + hh, n + 1, time.time() - t0, waist
            pre = None
        S = S_new
        sig += hh
        n += 1
        X = Xn
        rho_e = X[0]
        z_e, zt_up, zt_lo = zt_plus(S, half), zt_plus(S, 0), zt_plus(S, 1)
        if z_e < 0:
            return box, "E1", sig, n, time.time() - t0, None
        if rho_e * rho_e + z_e * z_e > R_STAR * R_STAR:
            return box, "E2", sig, n, time.time() - t0, None
        if any(rho_e * rho_e + w * w < S_STAR * S_STAR for w in (zt_up, zt_lo)):
            return box, "E3", sig, n, time.time() - t0, None
        if neck_event(X, zt_up):
            return box, "N", sig, n, time.time() - t0, None
    return box, "timeout", sig, n, time.time() - t0, None


def grid_pieces(waist):
    """The dyadic GRID pieces (i, k, s) in (lam, zeta) that cover a waist box."""
    lam_lo, lam_hi, ze_lo, ze_hi, s = waist
    i0, i1 = math.floor(lam_lo / GRID), math.floor(lam_hi / GRID)
    k0, k1 = math.floor(ze_lo / GRID), math.floor(ze_hi / GRID)
    return [(i, k, s) for i in range(i0, i1 + 1) for k in range(k0, k1 + 1)]


def certify_piece(piece):
    """The thin-neck sign pattern on the grid piece (i, k, s), for all d in [D_LO, D_HI]."""
    from research.calc.gap3_neck_lohner import certify_rel
    i, k, s = piece
    if i < 0:
        return piece, None
    res = certify_rel(_ball(i * GRID, (i + 1) * GRID), _ball(k * GRID, (k + 1) * GRID),
                      _ball(D_LO, D_HI), s, sigma_max=3.0)
    return piece, (min(res) if res else None)


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------
D0 = 0.7661974282399  # d_c - 1.0e-13 <= d_c (Theorem C: d_c = 0.766197428240 +- 1.1e-13)
D1 = 0.7662353515625  # = 0.765625 + 5 * 2^-13: the fold lemma's first part ends here
D_HOLE = 0.767578125  # below this the fold window is left to the fold lemma (gap3_fold.py)
D2 = D_HI  # 0.78125, contains 0.78
DELTA_MIN, DELTA_MAX = -0.0625, 0.46875  # covers z0 in [Z_LO, Z_HI] = [136/256, 269/256]
TOP = (-0.00390625, 0.0078125)  # neighbourhood of the individual horizon's top: necks
# The fold window in z0, certified by research/calc/gap3_fold.py: there F > 0 for every
# d in (d_c, D_HOLE] (F_zz > 0 and F_d > 0 up to D1, F_d > 0 from D1 to D_HOLE).
FOLD_WINDOW = (0.62890625, 0.636962890625)

# Partition of [DELTA_MIN, DELTA_MAX] x d into boxes: bands (lo, hi, delta-width, d-width).
# Calibrated on test boxes at both ends of the d-range (session 9): wider boxes get stuck,
# because Lohner's first-order enclosure explodes where the curves turn. The closer to the top
# of the individual horizon (delta = 0), the narrower; the single-d sweep at d = 0.78 needed
# similar widths (research/calc/data/gap3_partition_d078.txt).
BANDS = [
    (DELTA_MIN, -0.046875, 2.0**-10, 2.0**-9),
    (-0.046875, -0.015625, 2.0**-12, 2.0**-10),
    (-0.015625, -0.0078125, 2.0**-13, 2.0**-10),
    (-0.0078125, TOP[0], 2.0**-14, 2.0**-10),
    (TOP[0], TOP[1], 2.0**-14, 2.0**-10),
    (TOP[1], 0.01171875, 2.0**-15, 2.0**-10),
    (0.01171875, 0.013671875, 2.0**-14, 2.0**-10),
    (0.013671875, 0.0625, 2.0**-13, 2.0**-10),
    (0.0625, 0.0703125, 2.0**-12, 2.0**-10),
    (0.0703125, 0.15625, 2.0**-11, 2.0**-10),
    (0.15625, DELTA_MAX, 2.0**-10, 2.0**-9),
]


def slabs(w, d1=D0, d2=D2):
    """The d-slabs [D2 - (k+1) w, D2 - k w] (a grid anchored at D2) with top <= d2, the lowest
    one cut at d1."""
    out, k = [], 0
    while D2 - k * w > d1:
        hi, lo = D2 - k * w, max(d1, D2 - (k + 1) * w)
        if hi <= d2:
            out.append((lo, hi))
        k += 1
    return out[::-1]


def hole(dl, dh, w):
    """For a slab below D_HOLE: the delta-range on the w-grid whose z0-image over d in [dl, dh]
    lies inside FOLD_WINDOW (rounded inward, with interval arithmetic). Else None."""
    if dh > D_HOLE:
        return None
    a, b = FOLD_WINDOW
    lo = (arb(a) - arb(Z0_REF) - arb(SLOPE) * (arb(dl) - arb(D_REF))) / w
    hi = (arb(b) - arb(Z0_REF) - arb(SLOPE) * (arb(dh) - arb(D_REF))) / w
    return math.ceil(float(lo.upper()) + 1e-9) * w, math.floor(float(hi.lower()) - 1e-9) * w


def hole_inside(dl, dh, h):
    """Interval check: the z0-image of h x [dl, dh] lies in FOLD_WINDOW."""
    lo = arb(Z0_REF) + arb(h[0]) + arb(SLOPE) * (arb(dl) - arb(D_REF))
    hi = arb(Z0_REF) + arb(h[1]) + arb(SLOPE) * (arb(dh) - arb(D_REF))
    return bool(lo >= FOLD_WINDOW[0] and hi <= FOLD_WINDOW[1])


def initial_boxes(d1=D0, d2=D2):
    """All boxes of the partition (BANDS) in [DELTA_MIN, DELTA_MAX] x [d1, d2], minus the fold
    holes, ordered by d (so the sweep certifies (d_c, d_k] for growing d_k)."""
    boxes = []
    for lo, hi, w, dw in BANDS:
        for dl, dh in slabs(dw, d1, d2):
            h = hole(dl, dh, w)
            x = lo
            while x < hi:
                if h and h[0] <= x < h[1]:
                    x = h[1]
                    continue
                boxes.append((x, x + w, dl, dh))
                x += w
    return sorted(boxes, key=lambda b: (b[2], b[0]))


EXCLUDED = ("E1", "E2", "E3", "N", "W")


def split(box):
    a, b, c, d = box
    if b - a >= 2.0**-17:
        m = (a + b) / 2
        return [(a, m, c, d), (m, b, c, d)]
    if d - c >= 2.0**-16:
        m = (c + d) / 2
        return [(a, b, c, m), (a, b, m, d)]
    return []


def sweep(boxes=None, processes=4, logfile="gap3_interval.log", resume=False, chunk=64):
    """Run all boxes, lowest d first, bisecting failures (split) as soon as they occur, so that
    the certified range (d_c, d_k] grows steadily. Appends one line per run to logfile:
    delta_lo delta_hi d_lo d_hi event sigma steps seconds [lam_lo lam_hi ze_lo ze_hi s].
    resume: reuse logfile. A box logged as excluded is skipped; one logged as not excluded is
    replaced by its children (exclusions stay valid when the neck windows or LAM_MAX grow)."""
    import heapq
    import os
    from multiprocessing import Pool
    boxes = list(boxes if boxes is not None else initial_boxes())
    old = read_log(logfile) if resume and os.path.exists(logfile) else {}
    todo = []

    def expand(b):
        r = old.get(b)
        if r is None:
            todo.append(b)
        elif r[0] not in EXCLUDED:
            kids = split(b)
            for k in kids:
                expand(k)
            if not kids:
                todo.append(b)
    for b in boxes:
        expand(b)
    heap = [(b[2], b[0], b) for b in todo]
    heapq.heapify(heap)
    results, failed = [], []
    with Pool(processes) as pool, open(logfile, "a") as out:
        while heap:
            batch = [heapq.heappop(heap)[2] for _ in range(min(chunk, len(heap)))]
            for r in pool.imap(run_box, batch, chunksize=1):
                box, ev, sig, n, sec, waist = r
                results.append(r)
                extra = " ".join(repr(x) for x in waist) if waist else ""
                out.write(" ".join(repr(x) for x in box) + f" {ev} {sig:.4f} {n} {sec:.1f} {extra}\n")
                out.flush()
                if ev not in EXCLUDED:
                    kids = split(box)
                    for k in kids:
                        heapq.heappush(heap, (k[2], k[0], k))
                    if not kids:
                        failed.append(r)
    return results, failed


def certify_band(waists, processes=4, log=print):
    """Certify the medium-neck sign pattern on every GRID cell (lam, zeta, s) that meets a waist
    box, for all d in [D_LO, D_HI]; quarter cells that fail (down to GRID_MIN). Returns the set of
    certified cells as (lam_lo, lam_hi, ze_lo, ze_hi, s) and the list of uncertified ones."""
    from multiprocessing import Pool
    cells = sorted({p for w in waists for p in grid_pieces(w)})
    todo = [(i * GRID, (i + 1) * GRID, k * GRID, (k + 1) * GRID, s) for i, k, s in cells]
    ok, bad = [], []
    with Pool(processes) as pool:
        while todo:
            res = pool.map(_certify_cell, todo, chunksize=1)
            nxt = []
            for c, m in zip(todo, res):
                if m is not None:
                    ok.append(c)
                elif c[1] - c[0] > GRID_MIN:
                    l0, l1, z0, z1, s = c
                    lm, zm = (l0 + l1) / 2, (z0 + z1) / 2
                    nxt += [(a, b, e, f, s) for a, b in ((l0, lm), (lm, l1)) for e, f in ((z0, zm), (zm, z1))]
                else:
                    bad.append(c)
            if log:
                log(f"  band round: {len(todo)} cells, {len(nxt)} new after splitting, "
                    f"{len(ok)} certified, {len(bad)} failed")
            todo = nxt
    return ok, bad


def _certify_cell(c):
    from research.calc.gap3_neck_lohner import certify_rel
    l0, l1, z0, z1, s = c
    if l0 < 0:
        return None
    res = certify_rel(_ball(l0, l1), _ball(z0, z1), _ball(D_LO, D_HI), s, sigma_max=3.0)
    return min(res) if res else None


def covered_by(waist, cells):
    """Is the waist box inside the union of certified cells with the same direction? cells:
    the certified dyadic cells (lam_lo, lam_hi, ze_lo, ze_hi, s) of certify_band, all on the
    grids GRID / 2^j. Every GRID cell meeting the box must be certified or, recursively, every
    one of its quarters that meets the box."""
    lam_lo, lam_hi, ze_lo, ze_hi, s = waist
    cert = cells if isinstance(cells, set) else {(c[0], c[2], c[1] - c[0], c[4]) for c in cells}

    def ok(l0, z0, size):
        if (l0, z0, size, s) in cert:
            return True
        if size <= GRID_MIN:
            return False
        h = size / 2
        return all(ok(a, b, h) for a in (l0, l0 + h) for b in (z0, z0 + h)
                   if a <= lam_hi and a + h >= lam_lo and b <= ze_hi and b + h >= ze_lo)
    return all(ok(i * GRID, k * GRID, GRID) for i, k, _ in grid_pieces(waist))


def read_log(path):
    """Parse a sweep log into {box: (event, waist)} (the last run of each box wins)."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            box = tuple(float(x) for x in t[:4])
            waist = None
            if len(t) >= 13:
                waist = (float(t[8]), float(t[9]), float(t[10]), float(t[11]), int(t[12]))
            out[box] = (t[4], waist)
    return out


def check_cover(runs, d1=D0, d2=D2):
    """Every initial box in [d1, d2] is excluded, or split (split()) into children that are,
    recursively. Returns the list of uncovered boxes (empty: [DELTA_MIN, DELTA_MAX] x [d1, d2]
    minus the fold holes is covered by excluded boxes)."""
    bad = []

    def ok(box):
        r = runs.get(box)
        if r is None:
            bad.append(box)
            return False
        if r[0] in EXCLUDED:
            return True
        kids = split(box)
        if not kids:
            bad.append(box)
            return False
        return all([ok(k) for k in kids])
    for b in initial_boxes(d1, d2):
        ok(b)
    return bad


def holes_inside(d1=D0, d2=D2):
    """Every fold hole of the partition maps into FOLD_WINDOW (interval check)."""
    return all(hole_inside(dl, dh, h) for _, _, w, dw in BANDS for dl, dh in slabs(dw, d1, d2)
               for h in [hole(dl, dh, w)] if h)


def z0_window_covered(d1=D0, d2=D2):
    """For every d in [d1, d2], z0 = Z0_REF + delta + SLOPE (d - D_REF) with delta in
    [DELTA_MIN, DELTA_MAX] covers [d/2 + s_*, R_*] (interval arithmetic; the map is monotone
    in d, so the endpoints suffice)."""
    d = _ball(d1, d2)
    lo = arb(Z0_REF) + arb(DELTA_MIN) + arb(SLOPE) * (d - arb(D_REF))
    hi = arb(Z0_REF) + arb(DELTA_MAX) + arb(SLOPE) * (d - arb(D_REF))
    return bool(lo < d / 2 + S_STAR and hi > R_STAR)


def certified_up_to(runs):
    """The largest slab top d_k such that the partition on [D0, d_k] is covered (check_cover)."""
    tops = sorted(set.intersection(*[{dh for _, dh in slabs(dw)} for _, _, _, dw in BANDS]))
    best = None
    for t in tops:
        if check_cover(runs, D0, t):
            break
        best = t
    return best


def verify(logfile, d2=D2, processes=4, log=print):
    """Post-processing of a sweep log for d in [D0, d2]: cover check, holes inside the fold
    window, z0-window check, and the medium-neck band certificate for every waist event in
    that range. Returns a summary dict."""
    runs = {b: r for b, r in read_log(logfile).items() if b[3] <= d2}
    uncovered = check_cover(runs, D0, d2)
    waists = [w for ev, w in runs.values() if ev == "W"]
    cells, bad_cells = certify_band(waists, processes=processes, log=log)
    cert = {(c[0], c[2], c[1] - c[0], c[4]) for c in cells}
    uncertified = [w for w in waists if not covered_by(w, cert)]
    counts = {}
    for ev, _ in runs.values():
        counts[ev] = counts.get(ev, 0) + 1
    return {"d_range": (D0, d2), "runs": len(runs), "events": counts, "uncovered": uncovered,
            "holes_inside_fold_window": holes_inside(D0, d2),
            "z0_window": z0_window_covered(D0, d2), "waists": len(waists),
            "band_cells": len(cells), "band_failed": len(bad_cells),
            "uncertified_waists": uncertified}


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2 and sys.argv[1] == "verify":
        d2 = float(sys.argv[3]) if len(sys.argv) > 3 else D2
        summary = verify(sys.argv[2], d2=d2, processes=int(sys.argv[4]) if len(sys.argv) > 4 else 4)
        for k, v in summary.items():
            print(f"{k}: {v if not isinstance(v, list) or len(v) < 20 else str(len(v)) + ' items'}")
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[1] == "progress":
        print("covered up to d =", certified_up_to(read_log(sys.argv[2])))
        sys.exit(0)
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t0 = time.time()
    res, bad = sweep(processes=procs, resume=True)
    print(f"sweep: {len(res)} runs in {time.time() - t0:.0f} s; {len(bad)} boxes failed at the "
          f"minimal size")
