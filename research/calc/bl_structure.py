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
from fractions import Fraction

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


# ---------------------------------------------------------------------------
# Session 13: the t^6 coefficient in a covariant basis, and geometric (Weyl) variables
# ---------------------------------------------------------------------------
T6_NAMES = ["Lap2Ric", "sym(Ric.LapRic)", "dRic.dRic(l,l)", "d_iRic.d_jRic", "sym(d_kR_il d^lR^k_j)",
            "Ric^kl d_k d_l Ric_ij", "Ric^kl d_(i d_j) Ric_kl", "sym(Ric^kl d_k d_i Ric_jl)",
            "g Ric.LapRic", "g |dRic|^2", "g d_mR_kl d^lR^km", "sym(d_iR_kl d^kR^l_j)",
            "|Ric|^2 Ric", "tr(Ric^3) g"]


def curvature_values(x0, d, one, exact_r=None, centres=None):
    """Values at x0 of G = gamma_0^-1, Ric, nabla Ric, nabla nabla Ric, Delta Ric, nabla nabla Delta Ric of
    the BL (or given) data, plus gamma_3 (the t^6 coefficient) from the recursion."""
    J, g, _ = bm.series_jets(x0, d, 3, one, exact_r=exact_r, centres=centres, W=8)
    zero = 0 * one
    at = lambda a: J.coef(a, (0, 0, 0, 0), zero)  # noqa: E731
    g3 = {pq: J.coef(g[pq], (0, 0, 0, 3), zero) for pq in PAIRS}
    g2 = {pq: J.coef(g[pq], (0, 0, 0, 2), zero) for pq in PAIRS}
    g0 = {pq: J.s_slice(g[pq], 0) for pq in PAIRS}
    Gi = bm.inverse3(J, g0, one)
    Ric = bm.ricci3(J, g0, one)
    Gam = christoffel(J, g0, Gi, one)
    dR = cov_sym(J, Ric, Gam)
    ddR = cov_3(J, dR, Gam)
    lap = {pq: {} for pq in PAIRS}
    for (i, j) in PAIRS:
        acc = {}
        for k in R3:
            for l in R3:
                acc = J.add_(acc, J.mul(bm.sym_get(Gi, l, k), ddR[(i, j, k, l)]))
        lap[(i, j)] = acc
    dL = cov_sym(J, lap, Gam)
    ddL = cov_3(J, dL, Gam)
    G = [[at(bm.sym_get(Gi, i, j)) for j in R3] for i in R3]
    Rv = [[at(bm.sym_get(Ric, i, j)) for j in R3] for i in R3]
    dRv = {k: at(v) for k, v in dR.items()}
    ddRv = {k: at(v) for k, v in ddR.items()}
    Lv = [[at(bm.sym_get(lap, i, j)) for j in R3] for i in R3]
    L2 = [[sum((G[l][k] * at(ddL[(i, j, k, l)]) for k in R3 for l in R3), zero) for j in R3] for i in R3]
    g0v = [[at(bm.sym_get(g0, i, j)) for j in R3] for i in R3]
    return {"G": G, "R": Rv, "dR": dRv, "ddR": ddRv, "L": Lv, "L2": L2, "g0": g0v, "g3": g3, "g2": g2}


def t6_basis(cv, zero):
    """The 14 basis tensors T6_NAMES at the point (dicts over PAIRS), from curvature_values."""
    G, R, dR, ddR, L, L2, g0 = (cv[k] for k in ("G", "R", "dR", "ddR", "L", "L2", "g0"))
    S = lambda f: sum(f, zero)  # noqa: E731
    up = lambda M: [[S(G[i][a] * M[a][b] * G[b][j] for a in R3 for b in R3) for j in R3] for i in R3]  # noqa: E731
    Ru = up(R)                                            # R^{ij}
    D = lambda i, j, k: dR[(i, j, k)]                     # nabla_k R_ij  # noqa: E731
    DD = lambda i, j, k, l: ddR[(i, j, k, l)]             # nabla_l nabla_k R_ij  # noqa: E731
    RR = S(R[i][j] * Ru[i][j] for i in R3 for j in R3)    # |Ric|^2
    RG = [[S(R[i][k] * G[k][j] for k in R3) for j in R3] for i in R3]   # R_i^j
    trR3 = S(RG[i][j] * RG[j][k] * RG[k][i] for i in R3 for j in R3 for k in R3)
    RLs = S(Ru[k][l] * L[k][l] for k in R3 for l in R3)
    dd2 = S(G[m][n] * G[k][a] * G[l][b] * D(k, l, m) * D(a, b, n) for m in R3 for n in R3 for k in R3
            for a in R3 for l in R3 for b in R3)
    dd3 = S(G[m][n] * G[k][a] * G[l][b] * D(k, l, m) * D(n, a, b) for m in R3 for n in R3 for k in R3
            for a in R3 for l in R3 for b in R3)
    out = []
    for name in T6_NAMES:
        T = {}
        for (i, j) in PAIRS:
            if name == "Lap2Ric":
                v = L2[i][j]
            elif name == "sym(Ric.LapRic)":
                v = (S(R[i][k] * G[k][l] * L[l][j] for k in R3 for l in R3)
                     + S(R[j][k] * G[k][l] * L[l][i] for k in R3 for l in R3)) / 2
            elif name == "dRic.dRic(l,l)":
                v = S(G[l][m] * G[k][n] * D(i, k, l) * D(n, j, m) for l in R3 for m in R3 for k in R3 for n in R3)
            elif name == "d_iRic.d_jRic":
                v = S(G[k][a] * G[l][b] * D(k, l, i) * D(a, b, j) for k in R3 for a in R3 for l in R3 for b in R3)
            elif name == "sym(d_kR_il d^lR^k_j)":
                v = (S(G[k][a] * G[l][b] * D(i, l, k) * D(a, j, b) for k in R3 for a in R3 for l in R3 for b in R3)
                     + S(G[k][a] * G[l][b] * D(j, l, k) * D(a, i, b) for k in R3 for a in R3 for l in R3 for b in R3)) / 2
            elif name == "Ric^kl d_k d_l Ric_ij":
                v = S(Ru[k][l] * DD(i, j, l, k) for k in R3 for l in R3)
            elif name == "Ric^kl d_(i d_j) Ric_kl":
                v = S(Ru[k][l] * (DD(k, l, j, i) + DD(k, l, i, j)) for k in R3 for l in R3) / 2
            elif name == "sym(Ric^kl d_k d_i Ric_jl)":
                v = S(Ru[k][l] * (DD(j, l, i, k) + DD(i, l, j, k)) for k in R3 for l in R3) / 2
            elif name == "g Ric.LapRic":
                v = g0[i][j] * RLs
            elif name == "g |dRic|^2":
                v = g0[i][j] * dd2
            elif name == "g d_mR_kl d^lR^km":
                v = g0[i][j] * dd3
            elif name == "sym(d_iR_kl d^kR^l_j)":
                v = (S(G[k][a] * G[l][b] * D(k, l, i) * D(b, j, a) for k in R3 for a in R3 for l in R3 for b in R3)
                     + S(G[k][a] * G[l][b] * D(k, l, j) * D(b, i, a) for k in R3 for a in R3 for l in R3 for b in R3)) / 2
            elif name == "|Ric|^2 Ric":
                v = RR * R[i][j]
            elif name == "tr(Ric^3) g":
                v = trR3 * g0[i][j]
            T[(i, j)] = v
        out.append(T)
    return out


def weyl_speciality(x0, d, one, centres=None):
    """Complex Weyl invariants along the comoving geodesic observer at x0, as series in s = t^2 to
    O(s): Q = E + iB, E_ij = R_ij(gamma) + K K_ij - K_ik K^k_j (vacuum Gauss equation), B = curl K
    (O(t^3) for conformally flat data), I = tr(Q^2)/2, J = tr(Q^3)/6, S = 27 J^2 / I^3 (S = 1 for
    Petrov D/II). Returns (I0, I1, J0, J1, S0, S1, commutator) where *1 are coefficients of t^2 and
    commutator = |[E0, E2]| (mixed indices with gamma_0): zero iff the E-eigenframe is frozen at order
    t^2. Uses E only (B contributes from t^3 on)."""
    J, g, kap = bm.series_jets(x0, d, 2, one, centres=centres, W=4)
    zero = 0 * one
    Gi = bm.inverse3(J, g, one)
    Ric = bm.ricci3(J, g, one)
    trK = bm.trace_prod(J, Gi, kap)
    KGK = bm.matmul_sym(J, kap, Gi, kap)
    E = {pq: J.add_(Ric[pq], J.times_s(J.sub(J.mul(trK, kap[pq]), KGK[pq]))) for pq in PAIRS}
    # mixed E^i_j = G^ik E_kj as jets; keep the point values of the s^0 and s^1 coefficients
    Em = [[{} for _ in R3] for _ in R3]
    for i in R3:
        for j in R3:
            acc = {}
            for k in R3:
                acc = J.add_(acc, J.mul(bm.sym_get(Gi, i, k), bm.sym_get(E, k, j)))
            Em[i][j] = acc
    c = lambda a, n: J.coef(a, (0, 0, 0, n), zero)  # noqa: E731
    E0 = [[c(Em[i][j], 0) for j in R3] for i in R3]
    E1 = [[c(Em[i][j], 1) for j in R3] for i in R3]

    def mm(A, Bm):
        return [[sum((A[i][k] * Bm[k][j] for k in R3), zero) for j in R3] for i in R3]

    def tr(A):
        return A[0][0] + A[1][1] + A[2][2]
    E0E0, E0E1 = mm(E0, E0), mm(E0, E1)
    I0 = tr(E0E0) / 2
    I1 = tr(E0E1)                                 # d/ds of tr(E^2)/2
    J0 = tr(mm(E0E0, E0)) / 6
    J1 = tr(mm(E0E0, E1)) / 2                     # d/ds of tr(E^3)/6
    S0 = 27 * J0 * J0 / I0 ** 3
    S1 = 27 * (2 * J0 * J1 * I0 - 3 * J0 * J0 * I1) / I0 ** 4
    com = mm(E0, E1)
    com2 = mm(E1, E0)
    commutator = sum(((com[i][j] - com2[i][j]) * (com[i][j] - com2[i][j]) for i in R3 for j in R3), zero)
    return {"I0": I0, "I1": I1, "J0": J0, "J1": J1, "S0": S0, "S1": S1, "commutator_sq": commutator}


# the t^6 coefficient for conformally flat time-symmetric vacuum data (fitted, held-out checked, and
# checked exactly at rational points; notes/17 §5): coefficients of the reduced basis
T6_REDUCED = {"|Ric|^2 Ric": Fraction(-1, 15), "tr(Ric^3) g": Fraction(-41, 540),
              "d_iRic.d_jRic": Fraction(1, 60), "g |dRic|^2": Fraction(-1, 180),
              "Ric^kl d_k d_l Ric_ij": Fraction(-1, 180)}


def cf_identities(cv, one):
    """Residuals (dicts over PAIRS) of, for conformally flat data with R = 0:
       lap:  Delta Ric - (3 Ric.Ric - |Ric|^2 gamma_0)                 (derived, notes/17 §5)
       t4:   gamma_2 - (7 Ric.Ric - |Ric|^2 gamma_0) / 12               (derived from the t^4 formula)
       t6:   gamma_3 - sum_c c T_c over T6_REDUCED                       (fitted; this checks it)"""
    zero = 0 * one
    G, R, L, g0 = cv["G"], cv["R"], cv["L"], cv["g0"]
    RR = [[sum((R[i][k] * G[k][l] * R[l][j] for k in R3 for l in R3), zero) for j in R3] for i in R3]
    R2 = sum((R[i][j] * G[j][k] * R[k][l] * G[l][i] for i in R3 for j in R3 for k in R3 for l in R3), zero)
    T = t6_basis(cv, zero)
    idx = {n: k for k, n in enumerate(T6_NAMES)}
    out = {"lap": {}, "t4": {}, "t6": {}}
    for (i, j) in PAIRS:
        out["lap"][(i, j)] = L[i][j] - (3 * RR[i][j] - R2 * g0[i][j])
        out["t4"][(i, j)] = cv["g2"][(i, j)] - (7 * RR[i][j] - R2 * g0[i][j]) / 12
        acc = cv["g3"][(i, j)]
        for name, c in T6_REDUCED.items():
            acc = acc - T[idx[name]][(i, j)] * (one * c.numerator) / c.denominator
        out["t6"][(i, j)] = acc
    return out


def algebraic_closure_minor(cv, one):
    """Largest-|.|-lower-bound 3x3 minor of the 6x3 matrix [|Ric|^2 Ric, tr(Ric^3) g_0, gamma_3] over the
    six components. If its ball excludes 0, gamma_3 is not a combination a |Ric|^2 Ric + b tr(Ric^3) g_0
    at this point for any a, b: the coefficients cannot be pointwise algebraic in Ric at order t^6
    (Ric^3 reduces to these two by Cayley-Hamilton for traceless Ric)."""
    from itertools import combinations
    zero = 0 * one
    T = t6_basis(cv, zero)
    idx = {n: k for k, n in enumerate(T6_NAMES)}
    cols = [T[idx["|Ric|^2 Ric"]], T[idx["tr(Ric^3) g"]], cv["g3"]]
    M = [[c[pq] for c in cols] for pq in PAIRS]
    best = None
    for rows in combinations(range(6), 3):
        a = [M[r] for r in rows]
        det = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
               + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
        if best is None or abs(det).lower() > abs(best).lower():
            best = det
    return best


def electric_series(x0, d, one, exact_r=None, centres=None, order=2):
    """Mixed electric Weyl tensor E^i_j(t) = G^ik E_kj along the comoving observer at x0, as values of
    its s = t^2 coefficients E_0..E_order (3x3 lists), with E_ij = R_ij(gamma) + K K_ij - K_ik K^k_j."""
    W = 2 * order + 2            # Ric of a weight-W jet is exact to weight W - 2 (s^order at the point)
    J, g, kap = bm.series_jets(x0, d, order + 1, one, exact_r=exact_r, centres=centres, W=W)
    zero = 0 * one
    Gi = bm.inverse3(J, g, one)
    Ric = bm.ricci3(J, g, one)
    trK = bm.trace_prod(J, Gi, kap)
    KGK = bm.matmul_sym(J, kap, Gi, kap)
    E = {pq: J.add_(Ric[pq], J.times_s(J.sub(J.mul(trK, kap[pq]), KGK[pq]))) for pq in PAIRS}
    out = []
    for n in range(order + 1):
        M = [[zero] * 3 for _ in R3]
        for i in R3:
            for j in R3:
                acc = {}
                for k in R3:
                    acc = J.add_(acc, J.mul(bm.sym_get(Gi, i, k), bm.sym_get(E, k, j)))
                M[i][j] = J.coef(acc, (0, 0, 0, n), zero)
        out.append(M)
    return out


def _mm(A, B, zero):
    return [[sum((A[i][k] * B[k][j] for k in R3), zero) for j in R3] for i in R3]


def eigenframe_tests(x0, d, one, exact_r=None, centres=None):
    """For conformally flat data: (i) E_1 (coefficient of t^2) - (-(5/2) E0^2 + (1/2) tr(E0^2) I) (derived,
    notes/17 §6); (ii) |[E0, E1]|^2 and |[E0, E2]|^2 (frozen eigenframe at orders t^2 and t^4)."""
    zero = 0 * one
    Es = electric_series(x0, d, one, exact_r=exact_r, centres=centres, order=2)
    E0, E1, E2 = Es
    E00 = _mm(E0, E0, zero)
    tr2 = E00[0][0] + E00[1][1] + E00[2][2]
    resid = [[E1[i][j] - (-(5 * one / 2) * E00[i][j] + (tr2 / 2 if i == j else zero)) for j in R3] for i in R3]

    def comm_sq(A, B):
        P, Q = _mm(A, B, zero), _mm(B, A, zero)
        return sum(((P[i][j] - Q[i][j]) * (P[i][j] - Q[i][j]) for i in R3 for j in R3), zero)
    return {"E1_formula_residual": resid, "comm_t2": comm_sq(E0, E1), "comm_t4": comm_sq(E0, E2)}


def schwarzschild_crunch(R0=2.25, M=1.0, n=2000):
    """Single-hole control for the meaning of the series' convergence time (notes/17 §7). Geodesic slicing
    from rest at areal radius R0: R = R0 F(tau), tau = t sqrt(M / R0^3), F = cos^2(eta/2),
    tau = (eta + sin eta)/sqrt(8). Returns the crunch time t_c = pi sqrt(R0^3 / 8M), the minimum over
    eta in (0, pi) of H = F - (3/2) tau dF/dtau (> 0: gamma_rr = psi^4 H^2 never degenerates, no caustic
    of the congruence before the crunch), and the Kretschmann scalar 48 M^2 / R^6 near the crunch
    (it diverges: the slicing ends on the r = 0 curvature singularity, not on a coordinate caustic).
    dF/dtau = -sqrt(2) tan(eta/2) in closed form; dtau/deta = 0 only at eta = pi (mod 2 pi), so the
    Taylor series of F in tau has radius exactly pi/sqrt(8), i.e. t_c."""
    import math
    t_c = math.pi * math.sqrt(R0 ** 3 / (8 * M))
    Hmin = math.inf
    for k in range(1, n):
        eta = math.pi * k / n
        F = math.cos(eta / 2) ** 2
        tau = (eta + math.sin(eta)) / math.sqrt(8)
        Fp = -math.sqrt(2) * math.tan(eta / 2)
        Hmin = min(Hmin, F - 1.5 * tau * Fp)
    eta = math.pi * (1 - 1e-4)
    R = R0 * math.cos(eta / 2) ** 2
    t = math.sqrt(R0 ** 3 / M) * (eta + math.sin(eta)) / math.sqrt(8)
    return {"t_crunch": t_c, "H_min_on_(0,pi)": Hmin, "t_near_crunch": t, "kretschmann_near_crunch": 48 * M * M / R ** 6,
            "kretschmann_initial": 48 * M * M / R0 ** 6}
