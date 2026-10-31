"""B-to-A structural experiments on the series representation of the BL development (notes/16 §5).

Each experiment fixes a candidate class in advance and decides it by an exact or ball-arithmetic
computation:

  A  coefficient structure: is the t^4 coefficient a universal covariant expression in the data's
     curvature, gamma_2 = a Delta Ric + b Ric.Ric + c |Ric|^2 gamma_0 (with Ric.Ric = Ric gamma_0^-1 Ric)?
     The coefficients (a, b, c) = (-1/12, 5/6, -1/6) are derived by hand (notes/16 §5.1); the
     identity is checked exactly (rational arithmetic) at rational-distance points of BL data and
     of Schwarzschild, and in ball arithmetic at generic points.
  C1 finite ansatz "the geodesic slices stay conformally flat" (gamma(t) = Phi^4 x flat): decided by
     the s^1 coefficient of the Cotton tensor of gamma(t). A ball excluding 0 proves it nonzero.
  C2 finite ansatz "silent flow": the magnetic Weyl tensor relative to the comoving geodesic
     observers vanishes (B_ij = 0), which reduces the field equations to ODEs along each worldline
     (it contains Schwarzschild free fall and every 'local free fall' ansatz). Decided by the t^3
     coefficient of B.
  D  Pade resummation in s = t^2 of the coefficients (ball midpoints; exploratory only): pole locations
     against the certified radius and the Schwarzschild crunch.

Run: python -m research.calc.bl_structure   (writes research/calc/data/bl_structure.json)
"""

import json
import os
import time

import numpy as np
from flint import arb, fmpq

from research.calc import bl_metric as bm

PAIRS = bm.PAIRS
R3 = range(3)


def christoffel(J, g, Gi, one):
    """Gam^a_bc (dict over (a, b, c), symmetric in b, c) of a symmetric matrix of jets."""
    dg = {(pq, k): J.d(g[pq], k) for pq in PAIRS for k in R3}
    D = lambda i, j, k: dg[((i, j) if i <= j else (j, i), k)]  # noqa: E731
    G = {}
    for a in R3:
        for b in R3:
            for c in R3:
                acc = {}
                for e in R3:
                    acc = J.add_(acc, J.mul(bm.sym_get(Gi, a, e), J.sub(J.add_(D(e, b, c), D(e, c, b)), D(b, c, e))))
                G[(a, b, c)] = J.scale(acc, one / 2)
    return G


def cov_sym(J, S, Gam):
    """nabla_k S_ij (dict over (i, j, k), i, j full) of a symmetric 2-tensor of jets."""
    out = {}
    for i in R3:
        for j in R3:
            for k in R3:
                acc = J.d(bm.sym_get(S, i, j), k)
                for m in R3:
                    acc = J.sub(acc, J.add_(J.mul(Gam[(m, k, i)], bm.sym_get(S, m, j)),
                                            J.mul(Gam[(m, k, j)], bm.sym_get(S, i, m))))
                out[(i, j, k)] = acc
    return out


def cov_3(J, T, Gam):
    """nabla_l T_ijk of a 3-tensor of jets given as a dict over (i, j, k)."""
    out = {}
    for i in R3:
        for j in R3:
            for k in R3:
                for l in R3:
                    acc = J.d(T[(i, j, k)], l)
                    for m in R3:
                        acc = J.sub(acc, J.add_(J.add_(J.mul(Gam[(m, l, i)], T[(m, j, k)]),
                                                       J.mul(Gam[(m, l, j)], T[(i, m, k)])),
                                                J.mul(Gam[(m, l, k)], T[(i, j, m)])))
                    out[(i, j, k, l)] = acc
    return out


def t4_identity(x0, d, one, exact_r=None, centres=None):
    """gamma_2 (the t^4 coefficient) from the recursion, and the covariant candidate
    -(1/12) Delta Ric + (5/6) Ric.Ric - (1/6) |Ric|^2 gamma_0 of the data, at x0. Returns (gamma_2,
    candidate, pieces) as dicts over PAIRS of field elements (exact for fmpq input)."""
    J, g, _ = bm.series_jets(x0, d, 2, one, exact_r=exact_r, centres=centres, W=6)
    zero = 0 * one
    g2 = {pq: J.coef(g[pq], (0, 0, 0, 2), zero) for pq in PAIRS}
    g0 = {pq: J.s_slice(g[pq], 0) for pq in PAIRS}
    Gi = bm.inverse3(J, g0, one)
    Ric = bm.ricci3(J, g0, one)
    Gam = christoffel(J, g0, Gi, one)
    dR = cov_sym(J, Ric, Gam)
    ddR = cov_3(J, dR, Gam)
    at = lambda a: J.coef(a, (0, 0, 0, 0), zero)  # noqa: E731
    Gv = {pq: at(Gi[pq]) for pq in PAIRS}
    Rv = {pq: at(Ric[pq]) for pq in PAIRS}
    g0v = {pq: at(g0[pq]) for pq in PAIRS}
    lap = {pq: sum((bm.sym_get(Gv, l, k) * at(ddR[(pq[0], pq[1], k, l)]) for k in R3 for l in R3), zero)
           for pq in PAIRS}
    RR = {pq: sum((bm.sym_get(Rv, pq[0], k) * bm.sym_get(Gv, k, l) * bm.sym_get(Rv, l, pq[1])
                   for k in R3 for l in R3), zero) for pq in PAIRS}
    R2 = sum((bm.sym_get(Rv, i, j) * bm.sym_get(Gv, j, k) * bm.sym_get(Rv, k, l) * bm.sym_get(Gv, l, i)
              for i in R3 for j in R3 for k in R3 for l in R3), zero)
    cand = {pq: -lap[pq] / 12 + RR[pq] * 5 / 6 - R2 * g0v[pq] / 6 for pq in PAIRS}
    return g2, cand, {"lap": lap, "RR": RR, "R2": R2, "g0": g0v}


def cotton_s1(x0, d, one, centres=None):
    """The s^1 coefficient (t^2) of the Cotton tensor C_ijk = nabla_k S_ij - nabla_j S_ik,
    S = Ric - (R/4) gamma, of the slice metric gamma(t), at x0 (all 27 components)."""
    J, g, _ = bm.series_jets(x0, d, 3, one, centres=centres, W=7)
    Gi = bm.inverse3(J, g, one)
    Ric = bm.ricci3(J, g, one)
    R = bm.trace_prod(J, Gi, Ric)
    S = {pq: J.sub(Ric[pq], J.scale(J.mul(R, g[pq]), one / 4)) for pq in PAIRS}
    Gam = christoffel(J, g, Gi, one)
    dS = cov_sym(J, S, Gam)
    zero = 0 * one
    out = {}
    for i in R3:
        for j in R3:
            for k in R3:
                C = J.sub(dS[(i, j, k)], dS[(i, k, j)])
                out[(i, j, k)] = [J.coef(C, (0, 0, 0, n), zero) for n in range(2)]
    return out


def magnetic_t3(x0, d, one, centres=None):
    """Coefficients of t and t^3 of curl K_ij = eps_i^kl D_k K_lj (symmetrised) along the comoving
    geodesic observers (the magnetic Weyl tensor up to sign), at x0."""
    J, g, kap = bm.series_jets(x0, d, 3, one, centres=centres, W=7)
    Gi = bm.inverse3(J, g, one)
    Gam = christoffel(J, g, Gi, one)
    DK = cov_sym(J, kap, Gam)                      # D_k kappa_lj  as (l, j, k)
    # eps^{mkl} = [mkl] / sqrt(det gamma);  B_ij = gamma_im eps^{mkl} D_k kappa_lj
    dd = J.mul(J.mul(g[(0, 0)], g[(1, 1)]), g[(2, 2)])
    for a, b, c, sgn in (((0, 1), (1, 2), (0, 2), 2), ((0, 1), (0, 1), (2, 2), -1), ((0, 2), (0, 2), (1, 1), -1),
                         ((1, 2), (1, 2), (0, 0), -1)):
        dd = J.add_(dd, J.scale(J.mul(J.mul(g[a], g[b]), g[c]), sgn))
    det = dd
    isq = power(J, det, one, fmpq(-1, 2))
    perm = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}
    zero = 0 * one
    out = {}
    for (i, j) in PAIRS:
        sym = {}
        for (a, b) in ((i, j), (j, i)):
            acc = {}
            for (m, k, l), sg in perm.items():
                acc = J.add_(acc, J.scale(J.mul(bm.sym_get(g, a, m), J.mul(isq, DK[(l, b, k)])), sg))
            sym = J.add_(sym, acc)
        sym = J.scale(sym, one / 2)
        # K = t kappa(s): B = t [curl kappa]; the coefficient of t^(2n+1) is the s^n slice
        out[(i, j)] = [J.coef(sym, (0, 0, 0, n), zero) for n in range(2)]
    return out


def power(J, a, one, p):
    """a^p for a jet with positive constant term (binomial series), p rational."""
    a0 = a[0]
    x = J.scale({k: v for k, v in a.items() if k != 0}, one / a0)
    acc, term, coef = {0: one}, {0: one}, fmpq(1)
    for n in range(1, J.W + 1):
        coef = coef * (p - n + 1) / n
        term = J.mul(term, x)
        if not term:
            break
        acc = J.add_(acc, J.scale(term, one * coef))
    if isinstance(a0, arb):
        return J.scale(acc, a0 ** arb(p))
    raise TypeError("power() needs arb coefficients for a non-integer exponent")


# ---------------------------------------------------------------------------
# rational points
# ---------------------------------------------------------------------------
def plane_point(d, k):
    """A rational point (rho, 0, 0) on the plane z = 0 with rational distance r to both punctures:
    rho = p (k^2 - 1) / (2k), r = p (k^2 + 1) / (2k), p = d/2 (d, k rational)."""
    p = d / 2
    return (p * (k * k - 1) / (2 * k), fmpq(0), fmpq(0)), [p * (k * k + 1) / (2 * k)] * 2


def pade_poles(c, m, n):
    """Poles of the [m/n] Pade approximant of sum c_k s^k (numpy; exploratory)."""
    c = np.asarray(c, dtype=float)
    A = np.array([[c[m + i - j] if m + i - j >= 0 else 0.0 for j in range(1, n + 1)] for i in range(1, n + 1)])
    b = -c[m + 1: m + n + 1]
    q = np.linalg.solve(A, b)
    return np.roots(np.concatenate([q[::-1], [1.0]]))


def main(out=os.path.join(bm.DATA, "bl_structure.json")):
    res = {}
    t0 = time.time()
    # A: the t^4 identity, exactly at rational points (BL plane and axis, Schwarzschild), and in balls
    d = fmpq(3831, 5000)
    x0, rr = plane_point(d, fmpq(15, 8))
    g2, cand, _ = t4_identity(x0, d, fmpq(1), exact_r=rr)
    exact_plane = all(g2[pq] == cand[pq] for pq in PAIRS)
    za = fmpq(31, 50)
    xa = (fmpq(0), fmpq(0), za)
    g2a, canda, _ = t4_identity(xa, d, fmpq(1), exact_r=[za - d / 2, za + d / 2])
    exact_axis = all(g2a[pq] == canda[pq] for pq in PAIRS)
    xs = (fmpq(3, 5), fmpq(0), fmpq(4, 5))        # Schwarzschild, r = 1
    g2s, cands, _ = t4_identity(xs, None, fmpq(1), centres=[(fmpq(1, 2), (fmpq(0),) * 3, fmpq(1))])
    exact_schw = all(g2s[pq] == cands[pq] for pq in PAIRS)
    with bm.working_precision(256):
        xg = (arb("0.17"), arb("0.11"), arb("0.5"))
        g2g, candg, _ = t4_identity(xg, arb("0.7662"), arb(1))
        ball_generic = max(float(abs(g2g[pq] - candg[pq]).upper()) for pq in PAIRS)
    res["A_t4_identity"] = {"exact_plane_point": exact_plane, "exact_axis_point": exact_axis,
                            "exact_schwarzschild": exact_schw, "generic_point_ball_residual": ball_generic,
                            "plane_point": [str(v) for v in x0], "axis_point": str(za), "d": str(d)}
    print("A", res["A_t4_identity"], f"{time.time() - t0:.0f}s", flush=True)
    # C1, C2: obstructions in ball arithmetic at a generic point
    with bm.working_precision(256):
        xg = (arb("0.17"), arb("0.11"), arb("0.5"))
        C = cotton_s1(xg, arb("0.7662"), arb(1))
        c0 = max(float(abs(v[0]).upper()) for v in C.values())
        c1 = [(k, v[1]) for k, v in C.items()]
        nz = [(k, str(v)) for k, v in c1 if not v.contains(0)]
        B = magnetic_t3(xg, arb("0.7662"), arb(1))
        b0 = max(float(abs(v[0]).upper()) for v in B.values())
        bnz = [(k, str(v[1])) for k, v in B.items() if not v[1].contains(0)]
        # controls: Schwarzschild (both must vanish identically)
        cs = [(arb(1) / 2, (arb(0),) * 3, (arb("0.36") + arb("0.04") + arb("0.81")).sqrt())]
        xs2 = (arb("0.6"), arb("0.2"), arb("0.9"))
        Cs = cotton_s1(xs2, None, arb(1), centres=cs)
        Bs = magnetic_t3(xs2, None, arb(1), centres=cs)
        cs_max = max(float(abs(x).upper()) for v in Cs.values() for x in v)
        bs_max = max(float(abs(x).upper()) for v in Bs.values() for x in v)
    res["C1_conformally_flat_slices"] = {"order_t0_max": c0, "order_t2_nonzero_components": nz[:6],
                                         "excluded": bool(nz), "schwarzschild_control_max": cs_max}
    res["C2_silent_flow"] = {"order_t_max": b0, "order_t3_nonzero_components": bnz[:4],
                             "excluded": bool(bnz), "schwarzschild_control_max": bs_max}
    print("C1", res["C1_conformally_flat_slices"], flush=True)
    print("C2", res["C2_silent_flow"], f"{time.time() - t0:.0f}s", flush=True)
    # D: Pade (exploratory): coefficients from the Cartesian engine (ball arithmetic, midpoints used)
    dd = arb("0.7662")
    pts = {"plane, near the waist of S_c": (arb("0.26"), arb(0), arb(0)),
           "near the top of S_c": (arb("0.05"), arb(0), arb("0.6"))}
    D = {}
    NP = 8
    with bm.working_precision(128):
        for name, x0 in pts.items():
            J, g, _ = bm.series_jets(x0, dd, NP, arb(1))
            row = {}
            for comp, pq in (("xx", (0, 0)), ("yy", (1, 1)), ("zz", (2, 2))):
                c = [float(J.coef(g[pq], (0, 0, 0, k), arb(0)).mid()) for k in range(NP + 1)]
                poles = pade_poles(c, 4, 4)
                real_pos = sorted(float(p.real) for p in poles if abs(p.imag) < 1e-6 * abs(p) and p.real > 0)
                row[comp] = {"coefficients_s": c,
                             "pade_4_4_poles_s": [[float(p.real), float(p.imag)] for p in poles],
                             "smallest_positive_real_pole_as_t": (real_pos[0] ** 0.5) if real_pos else None,
                             "ratio_estimates_T": [abs(c[k] / c[k + 1]) ** 0.5 for k in range(2, NP) if c[k + 1]]}
            D[name] = row
            print("D", name, f"{time.time() - t0:.0f}s", flush=True)
        # Schwarzschild control (M = 1, isotropic r = 1, R0 = 9/4): crunch at t = pi sqrt(R0^3 / 8)
        cs = [(arb(1) / 2, (arb(0),) * 3, arb(1))]
        J, g, _ = bm.series_jets((arb(1), arb(0), arb(0)), None, NP, arb(1), centres=cs)
        c = [float(J.coef(g[(1, 1)], (0, 0, 0, k), arb(0)).mid()) for k in range(NP + 1)]
        poles = pade_poles(c, 4, 4)
        D["schwarzschild_control"] = {"exact_crunch_t": float(np.pi * np.sqrt(2.25 ** 3 / 8)),
                                      "coefficients_s": c,
                                      "pade_4_4_poles_s": [[float(p.real), float(p.imag)] for p in poles]}
    res["D_pade_exploratory"] = D
    print("D done", f"{time.time() - t0:.0f}s", flush=True)
    res["seconds"] = time.time() - t0
    res["replay"] = "python -m research.calc.bl_structure"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    main()
