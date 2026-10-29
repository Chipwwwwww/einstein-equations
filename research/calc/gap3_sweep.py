"""Gap 3, validated sweep: exclude a common minimal surface of Brill-Lindquist data at d = 0.78,
shooting curve by shooting curve, with Lohner enclosures (research/calc/lohner.py).

Setting (research/calc/bl_global.py). An outermost common minimal surface S* would be a stable,
embedded, G-invariant sphere whose profile runs from an axis point (0, z0) to the plane z = 0,
meets that plane once and at a right angle, stays in B(0, R_*) and outside the throat balls
B(p+-, s_*) (Lemmas G1, G2), with z0 in [d/2 + s_*, R_*] = [0.545, 1.05].

Meridian ODE (H = 0, metric psi^4 delta) in sigma-time, ds = rho dsigma,
(u, v) = rho (cos alpha, sin alpha):

    rho' = u,  z' = v,  u' = rho - v P,  v' = u P,   P = 4 (u psi_z - v psi_rho) / psi.

Regular up to the axis (tests/test_bl_global.py); u^2 + v^2 = rho^2 is invariant and the axis
points are fixed points, so a curve never reaches the axis at finite sigma. Start: the box of
states at arclength 2^-26 from (0, z0) for all z0 in a box (fold_proof.axis_enclosure).

Events, checked for the whole z0 box (every trajectory in the enclosure):
  E1  z < 0 at a step end, and every step enclosure B that allowed z = 0 had u != 0: the first
      plane crossing is not at a right angle;
  E2  outside B(0, R_*) at a step end;
  E3  inside a throat ball at a step end;
  N   (approach lemma + thin-neck lemma, research/calc/gap3_neck.py) at a step end
      u < 0, rho <= rho_e <= RHO_N, and with K the bound below,
          z -+ dz inside [ZW_LO, ZW_HI],  dz = (2 rho_e / e + K rho_e^2 / 2) / sqrt(c),
          |v| + K rho_e^2 / 2 < LAM0,   c = 1 - K rho_e - K^2 rho_e^2 / 4 > 0.

Approach lemma. Let K >= 4 (|psi_z| + |psi_rho|) / psi on [0, RHO_N] x [ZW_LO, ZW_HI]; then
|P| <= K rho there (|u|, |v| <= rho). Follow a curve from a state with u < 0, rho <= rho_e. While
u < 0, rho decreases and dv/drho = P, so |dv/drho| <= K rho. Either u < 0 forever, and then
rho -> 0 (else u' >= rho (1 - K rho) stays positive): the curve closes up on the axis at a
point with z > 0, and is not the profile of S*; or u reaches 0 at a waist with
rho_w = |v_w| <= |v_e| + K rho_e^2 / 2 < LAM0. Between, |v| <= rho_w + K x / 2 with
x = rho^2 - rho_w^2, so rho^2 - v^2 >= c x and
    |dz| <= int (rho_w + K x / 2) / sqrt(c x) drho <= (rho_w arccosh(rho_e / rho_w) + K rho_e^2 / 2) / sqrt(c),
and rho_w arccosh(rho_e / rho_w) <= rho_w log(2 rho_e / rho_w) <= 2 rho_e / e. So the curve stays
in the region where K holds (continuity) and the waist height z_w is in [ZW_LO, ZW_HI]. The
thin-neck lemma then gives a Jacobi field with two zeros around the waist: the surface is
unstable, so it is not S*.

Run: python -m research.calc.gap3_sweep [workers]   (the full sweep, several hours on 4 cores;
results appended to the log file named in sweep()).
"""

import math
import time

from flint import arb

from research.calc import fold_proof as fp
from research.calc import lohner
from research.calc.fold_proof import J
from research.calc.gap3_neck import LAM0, ZW_HI, ZW_LO

D = arb("0.78")
S_STAR, R_STAR = arb("0.155"), arb("1.05")
QUARTER = arb(1) / 4  # m/2 with m = 1/2
RHO_N = arb("0.02")  # largest section radius for event N


def f_bl(y, d=D):
    """sigma-system right-hand side for (rho, z, u, v); works on arb, arb_series and duals."""
    rho, z, u, v = y
    psi, pr, pz = 1, 0, 0
    for c in (d / 2, -d / 2):
        Z = z - c
        ir = 1 / (rho * rho + Z * Z).sqrt()
        ir3 = ir * ir * ir
        psi = psi + QUARTER * ir
        pr = pr - QUARTER * rho * ir3
        pz = pz - QUARTER * Z * ir3
    P = 4 * (u * pz - v * pr) / psi
    return [u, v, rho - v * P, u * P]


def k_bound(rho_e=RHO_N, zlo=ZW_LO, zhi=ZW_HI, pieces=64, d=D):
    """Upper bound for 4 (|psi_z| + |psi_rho|) / psi on [0, rho_e] x [zlo, zhi]."""
    rho = arb(rho_e) / 2 + arb(0, 1) * (arb(rho_e) / 2)
    zlo, zhi = arb(zlo), arb(zhi)
    best = arb(0)
    for k in range(pieces):
        a, b = zlo + (zhi - zlo) * k / pieces, zlo + (zhi - zlo) * (k + 1) / pieces
        z = (a + b) / 2 + arb(0, 1) * ((b - a) / 2)
        psi, pr, pz = arb(1), arb(0), arb(0)
        for c in (d / 2, -d / 2):
            Z = z - c
            ir = 1 / (rho * rho + Z * Z).sqrt()
            psi += QUARTER * ir
            pr -= QUARTER * rho * ir**3
            pz -= QUARTER * Z * ir**3
        best = max(best, (4 * (abs(pz) + abs(pr)) / psi).upper())
    return best


K_N = k_bound()


def neck_event(X):
    """Event N for the hull X = (rho, z, u, v) (see the module docstring)."""
    rho, z, u, v = X
    if not (u < 0 and rho < RHO_N):
        return False
    re, K = rho.upper(), K_N
    c = 1 - K * re - K * K * re * re / 4
    if not c > 0:
        return False
    dz = (2 * re / arb(1).exp() + K * re * re / 2) / c.sqrt()
    return bool(z - dz > ZW_LO and z + dz < ZW_HI and abs(v) + K * re * re / 2 < LAM0)


def start_box(z0, d=D):
    """Box of states (rho, z, u, v) at arclength 2^-26 from the axis point (0, z0), for all z0
    in the ball z0."""
    zj = J([z0, arb(0), arb(0), arb(0), arb(0), arb(0)])
    dj = J([d, arb(0), arb(0), arb(0), arb(0), arb(0)])
    A, Q, Z = fp.axis_enclosure(zj, dj, fp.S1)
    s1 = arb(fp.S1)
    rho = (s1 + Q * s1**3).c[0]
    z = (zj + Z * s1**2).c[0]
    a = (A * s1).c[0]
    return [rho, z, rho * a.cos(), rho * a.sin()]


def _split(S):
    """Cover the Lohner set S by two sets, halving its widest coordinate box r_j."""
    j = max(range(len(S.r)), key=lambda k: S.r[k].rad())
    rj = S.r[j]
    lo, mid, hi = arb(rj.lower()), arb(rj.mid()), arb(rj.upper())
    out = []
    for piece in (lo.union(mid), mid.union(hi)):
        r = list(S.r)
        r[j] = piece
        out.append(lohner.Set(S.m, S.A, r))
    return out


def run_lohner(z0, d=D, h=0.25, ds=0.004, L=14, sigma_max=80.0, max_steps=40000, log=None,
               split_r=None, max_pieces=4096):
    """Integrate the z0 box (an arb ball) until every trajectory in it is excluded.
    Returns (event, sigma, number of steps); event is E1, E2, E3, N, a combination such as
    "E2+N" when the box was split, or a failure (stuck, plane-undecided, timeout,
    too-many-pieces). With split_r, a set whose widest coordinate r_j exceeds split_r is cut
    into two halves that are followed separately (this is z0-bisection that keeps the common
    part of the integration). Event N is certified for d = D only; for another d
    (falsification runs, tests/test_gap3.py) it is switched off, and E2, E3 use the radii
    of D, which are then just test values."""
    def f(y):
        return f_bl(y, d)
    neck = d is D
    stack = [(lohner.Set.from_box(start_box(z0, d)), 0.0)]
    events, n, sig_max, pieces = set(), 0, 0.0, 1
    while stack:
        S, sig = stack.pop()
        ev = None
        while ev is None:
            if sig >= sigma_max or n >= max_steps:
                return "timeout", sig, n
            if split_r is not None and max(float(r.rad()) for r in S.r) > split_r:
                pieces += 1
                if pieces > max_pieces:
                    return "too-many-pieces", sig, n
                S, other = _split(S)
                stack.append((other, sig))
            X = S.hull()
            hh = min(h, 2.0 ** math.floor(math.log2(ds / float(X[0].upper()))))
            while True:
                try:
                    S_new, B = lohner.step(f, S, hh, L)
                    break
                except (RuntimeError, ValueError, ZeroDivisionError):
                    hh /= 2
                    if hh < 1e-6:
                        return "stuck", sig, n
            if not B[1] > 0 and not (B[2] > 0 or B[2] < 0):
                return "plane-undecided", sig, n
            S = S_new
            sig += hh
            n += 1
            X = S.hull()
            rho_e, z_e = X[0], X[1]
            if log is not None:
                log(sig, S)
            if z_e < 0:
                ev = "E1"
            elif rho_e * rho_e + z_e * z_e > R_STAR * R_STAR:
                ev = "E2"
            elif any(rho_e * rho_e + (z_e - c) ** 2 < S_STAR * S_STAR for c in (d / 2, -d / 2)):
                ev = "E3"
            elif neck and neck_event(X):
                ev = "N"
        events.add(ev)
        sig_max = max(sig_max, sig)
    return "+".join(sorted(events)), sig_max, n


EXCLUDED = ("E1", "E2", "E3", "N")


def excluded(event):
    return all(e in EXCLUDED for e in event.split("+"))


Z_LO, Z_HI = 139 / 256, 269 / 256  # contains [d/2 + s_*, R_*] = [0.545, 1.05]


def run_box(box, d=D, split_r=None):
    """Worker: box = (lo, hi), dyadic floats. Returns (lo, hi, event, sigma, steps, seconds)."""
    lo, hi = box
    t0 = time.time()
    z0 = (arb(lo) + arb(hi)) / 2 + arb(0, 1) * ((arb(hi) - arb(lo)) / 2)
    try:
        ev, sig, n = run_lohner(z0, d, split_r=split_r)
    except (RuntimeError, ValueError, ZeroDivisionError) as e:
        ev, sig, n = f"error {type(e).__name__}", 0.0, 0
    return lo, hi, ev, sig, n, time.time() - t0


def sweep(lo=Z_LO, hi=Z_HI, width=2.0**-8, min_width=2.0**-20, processes=4,
          logfile="gap3_sweep.log", boxes=None, split_r=None):
    """Cover [lo, hi] by dyadic boxes of the given width (or start from the list boxes);
    run each, bisect every box that is not excluded. Each result is appended to logfile as
    it arrives. Returns the list of results and the list of boxes that failed at min_width
    (empty means the whole range is excluded)."""
    from functools import partial
    from multiprocessing import Pool
    if boxes is None:
        n = round((hi - lo) / width)
        boxes = [(lo + k * width, lo + (k + 1) * width) for k in range(n)]
    todo, results, failed = list(boxes), [], []
    work = partial(run_box, split_r=split_r)
    with Pool(processes) as pool, open(logfile, "a") as out:
        while todo:
            nxt = []
            for r in pool.imap_unordered(work, todo, chunksize=1):
                results.append(r)
                b_lo, b_hi, ev = r[:3]
                out.write(f"{b_lo!r} {b_hi!r} {ev.replace(' ', '_')} {r[3]:.4f} {r[4]} {r[5]:.1f}\n")
                out.flush()
                if not excluded(ev):
                    if b_hi - b_lo <= min_width:
                        failed.append(r)
                    else:
                        m = (b_lo + b_hi) / 2
                        nxt += [(b_lo, m), (m, b_hi)]
            todo = sorted(nxt)
    return results, failed


def covered(results, lo=Z_LO, hi=Z_HI):
    """The gaps in [lo, hi] left by the excluded boxes among results (empty: all covered)."""
    boxes = sorted((r[0], r[1]) for r in results if excluded(r[2]))
    gaps, cur = [], lo
    for a, b in boxes:
        if a > cur:
            gaps.append((cur, a))
        cur = max(cur, b)
    if cur < hi:
        gaps.append((cur, hi))
    return gaps


PARTITION = "research/calc/data/gap3_partition_d078.txt"


def read_partition(path=PARTITION):
    with open(path) as f:
        return [tuple(float(x) for x in line.split()) for line in f if not line.startswith("#")]


def main():
    """python -m research.calc.gap3_sweep [workers] [partition | uniform]: run the certificate
    over the stored partition (default), or from uniform 2^-8 boxes, bisecting failures."""
    import sys
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    mode = sys.argv[2] if len(sys.argv) > 2 else PARTITION
    t0 = time.time()
    boxes = None if mode == "uniform" else read_partition(mode)
    res, bad = sweep(processes=procs, boxes=boxes)
    gaps = covered(res)
    counts = {}
    for r in res:
        if excluded(r[2]):
            for e in r[2].split("+"):
                counts[e] = counts.get(e, 0) + 1
    print(f"{len(res)} box runs in {time.time() - t0:.0f} s; excluded boxes by event: {counts}")
    print(f"{len(bad)} boxes not excluded at the minimal width; gaps in [Z_LO, Z_HI]: {gaps}")
    if not bad and not gaps:
        print(f"Every z0 in [{Z_LO}, {Z_HI}] is excluded: no common minimal surface at d = 0.78.")


if __name__ == "__main__":
    main()
