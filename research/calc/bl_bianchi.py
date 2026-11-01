"""The Einstein-Bianchi (E, H) form of the geodesic-slicing development, as jets (notes/19).

Gauss gauge: comoving coordinates on a geodesic, hypersurface-orthogonal congruence u = d_t, metric
g = -dt^2 + gamma_ij dx^i dx^j, K_ij = (1/2) d_t gamma_ij = nabla_i u_j (repo sign; K = sigma + theta/3
gamma). Electric and magnetic Weyl tensors E_ij = C_iajb u^a u^b, H = curl K (constraint C3 of
Ellis-van Elst 1998, eq. 34 with omega = 0), curl A_ij = sym(eps_i^kl D_k A_lj), eps^mkl = [mkl]/sqrt(det
gamma). In lower-index comoving components (L_u = d_t in this chart) the vacuum equations are

    d_t gamma = 2 K
    d_t K     = -E + K.K                                         (K.K)_ij = K_ik K^k_j
    d_t E     =  eps curl H - 2 theta E + (5/2)(K.E + E.K) - gamma tr(K.E)
    d_t H     = -eps curl E - 2 theta H + (5/2)(K.H + H.K) - gamma tr(K.H)

(eps = 1). These are Ellis-van Elst (40)-(41) with Edot_<ab> = nabla_u E (spatial and trace-free for a
geodesic u) and nabla_u E_ij = d_t E_ij - K_i^k E_kj - K_j^k E_ik; they coincide with Friedrich-Rendall's
P_ij = 0 (2000, §4.4) term by term. The parameter eps multiplies every curl: eps = 0 is the silent
system (notes/18 §6), and the coefficients are polynomials in eps (the curl degree, §4 of notes/19).

The variant 'FR-printed' replaces the H equation by the algebraic terms as extracted from the printed
Q_kl of Friedrich-Rendall (2000) eq. before (4.47): d_t H = -curl E + (1/2)(K.H + H.K) + gamma tr(K.H).
It is kept only as a control that the jet comparison can tell the two apart.

Initial data (time symmetric, conformally flat): gamma = psi^4 delta, K = 0, E = Ric(gamma_0), H = 0.
"""

import time

from flint import arb, fmpq

from research.calc import bl_metric as bm
from research.calc import bl_structure as st

PAIRS = bm.PAIRS
R3 = range(3)
PERM = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def _full(J, A, G, B):
    """(A G B)_ij for symmetric A, B (dicts over PAIRS), full 3x3 dict over (i, j)."""
    AG = {}
    for i in R3:
        for k in R3:
            acc = {}
            for m in R3:
                acc = J.add_(acc, J.mul(bm.sym_get(A, i, m), bm.sym_get(G, m, k)))
            AG[(i, k)] = acc
    out = {}
    for i in R3:
        for j in R3:
            acc = {}
            for k in R3:
                acc = J.add_(acc, J.mul(AG[(i, k)], bm.sym_get(B, k, j)))
            out[(i, j)] = acc
    return out


def _trace_full(J, G, X):
    """G^ij X_ji for a full 3x3 dict X."""
    acc = {}
    for i in R3:
        for j in R3:
            acc = J.add_(acc, J.mul(bm.sym_get(G, i, j), X[(j, i)]))
    return acc


def _det3(J, g):
    G = lambda i, j: bm.sym_get(g, i, j)  # noqa: E731
    m = J.mul
    d = J.sub(m(G(1, 1), G(2, 2)), m(G(1, 2), G(1, 2)))
    acc = m(G(0, 0), d)
    acc = J.sub(acc, m(G(0, 1), J.sub(m(G(0, 1), G(2, 2)), m(G(1, 2), G(0, 2)))))
    acc = J.add_(acc, m(G(0, 2), J.sub(m(G(0, 1), G(1, 2)), m(G(1, 1), G(0, 2)))))
    return acc


def _inv_sqrt_det(J, g, psi, one):
    """det(gamma)^(-1/2) as a jet, rational when psi is: psi^-6 (det gamma / psi^12)^(-1/2), the second
    factor by the binomial series (its constant term is 1 at t = 0 for gamma(0) = psi^4 delta)."""
    p2 = J.mul(psi, psi)
    p6 = J.mul(J.mul(p2, p2), p2)
    ip6 = J.inv(p6, one)
    q = J.mul(_det3(J, g), J.mul(ip6, ip6))
    q0 = q.get(0, 0 * one)
    x = J.scale({k: v for k, v in q.items() if k != 0}, one / q0)
    acc, term, c = {0: one}, {0: one}, fmpq(1)
    for n in range(1, J.W + 1):
        c = c * (fmpq(-1, 2) - n + 1) / n
        term = J.mul(term, x)
        if not term:
            break
        acc = J.add_(acc, J.scale(term, one * c))
    if q0 != 1:
        raise ValueError("normalised determinant must be 1 at the expansion point")
    return J.mul(ip6, acc)


def curl(J, A, g, Gam, isq, one):
    """curl A_ij = (1/2)(eps_i^kl D_k A_lj + eps_j^kl D_k A_li), eps^mkl = [mkl] isq (same orientation as
    bl_structure.magnetic_t3)."""
    DA = st.cov_sym(J, A, Gam)          # (l, j, k) -> D_k A_lj
    out = {}
    for (i, j) in PAIRS:
        acc = {}
        for (a, b) in ((i, j), (j, i)):
            for (m, k, l), sg in PERM.items():
                acc = J.add_(acc, J.scale(J.mul(bm.sym_get(g, a, m), DA[(l, b, k)]), sg))
        out[(i, j)] = J.scale(J.mul(acc, isq), one / 2)
    return out


def rhs(J, st_, psi, one, eps, variant="EvE"):
    g, K, E, H = st_["g"], st_["K"], st_["E"], st_["H"]
    Gi = bm.inverse3(J, g, one)
    Gam = st.christoffel(J, g, Gi, one)
    isq = _inv_sqrt_det(J, g, psi, one)
    th = bm.trace_prod(J, Gi, K)
    KK = bm.matmul_sym(J, K, Gi, K)
    KE = _full(J, K, Gi, E)
    KH = _full(J, K, Gi, H)
    trKE = _trace_full(J, Gi, KE)
    trKH = _trace_full(J, Gi, KH)
    cH = curl(J, H, g, Gam, isq, one)
    cE = curl(J, E, g, Gam, isq, one)
    dg = {pq: J.scale(K[pq], 2 * one) for pq in PAIRS}
    dK = {pq: J.sub(KK[pq], E[pq]) for pq in PAIRS}
    dE, dH = {}, {}
    for (i, j) in PAIRS:
        sKE = J.add_(KE[(i, j)], KE[(j, i)])        # (K.E + E.K)_ij
        sKH = J.add_(KH[(i, j)], KH[(j, i)])
        dE[(i, j)] = J.add_(J.add_(J.scale(cH[(i, j)], eps), J.scale(J.mul(th, E[(i, j)]), -2 * one)),
                            J.sub(J.scale(sKE, one * fmpq(5, 2)), J.mul(g[(i, j)], trKE)))
        if variant == "EvE":
            dH[(i, j)] = J.add_(J.add_(J.scale(cE[(i, j)], -eps), J.scale(J.mul(th, H[(i, j)]), -2 * one)),
                                J.sub(J.scale(sKH, one * fmpq(5, 2)), J.mul(g[(i, j)], trKH)))
        elif variant == "FR-printed":
            dH[(i, j)] = J.add_(J.scale(cE[(i, j)], -eps),
                                J.add_(J.scale(sKH, one / 2), J.mul(g[(i, j)], trKH)))
        else:
            raise ValueError(variant)
    return {"g": dg, "K": dK, "E": dE, "H": dH}


def bianchi_jets(psi_fn, W, one, eps=1, variant="EvE"):
    """Solve the E-H system as (X, Y, Z, t)-jets of total weight W (t of weight 1) from time-symmetric
    conformally flat data gamma_0 = psi^4 delta; psi_fn(J) returns the spatial jet of psi. The t^n slice
    of every field is exact to spatial weight W - n."""
    J = bm.Jets(W, sw=1)
    psi = psi_fn(J)
    p2 = J.mul(psi, psi)
    p4 = J.mul(p2, p2)
    g0 = {pq: (dict(p4) if pq[0] == pq[1] else {}) for pq in PAIRS}
    # E_0 = Ric(gamma_0) must be exact to spatial weight W: compute it on jets of weight W + 2 and copy
    J2 = bm.Jets(W + 2, sw=1)
    psi2 = psi_fn(J2)
    q2 = J2.mul(psi2, psi2)
    q4 = J2.mul(q2, q2)
    Ric2 = bm.ricci3(J2, {pq: (dict(q4) if pq[0] == pq[1] else {}) for pq in PAIRS}, one)
    E0 = {}
    for pq in PAIRS:
        E0[pq] = {J.index[J2.mons[i]]: v for i, v in Ric2[pq].items()
                  if J2.mons[i][3] == 0 and sum(J2.mons[i][:3]) <= W}
    state = {"g": g0, "K": {pq: {} for pq in PAIRS}, "E": E0, "H": {pq: {} for pq in PAIRS}}
    for n in range(W):
        F = rhs(J, state, psi, one, eps, variant)
        for key in state:
            for pq in PAIRS:
                sl = J.drop_above(J.s_slice(F[key][pq], n), W - 1 - n)
                J.put_s(state[key][pq], n + 1, J.scale(sl, one / (n + 1)), W)
    return J, state


def bl_psi(x0, d, one, exact_r=None, centres=None):
    cs = centres if centres is not None else bm.bl_centres(d, one, x0, exact_r)
    return lambda J: bm.psi_jet(J, x0, cs, one)


def point_coefficients(J, state, key, n_max, zero):
    """Values at the expansion point of the t^n coefficients (n <= n_max) of a field, as dicts over PAIRS."""
    return [{pq: J.coef(state[key][pq], (0, 0, 0, n), zero) for pq in PAIRS} for n in range(n_max + 1)]




def sectors(x0s, ds, N, eps_values=None, prec=256):
    """Curl-degree sectors of gamma at a point: gamma_{t^(2k)} = sum_j eps^(2j) G_{k,j}. Runs the E-H jets at
    eps = 0, 1, ..., N - 1 (exact interpolation in u = eps^2, degree <= N - 2 by the lag lemma) plus one
    extra eps as a consistency check. Returns mixed matrices G0^-1 G_{k,j} as nested lists of arb."""
    from flint import arb_mat
    eps_values = eps_values or list(range(N)) + [N]
    with bm.working_precision(prec):
        x0 = tuple(arb(v) for v in x0s)
        d = arb(ds)
        one = arb(1)
        vals = {}
        for e in eps_values:
            J, s = bianchi_jets(bl_psi(x0, d, one), 2 * N, one, eps=arb(e))
            vals[e] = point_coefficients(J, s, "g", 2 * N, arb(0))
        g0 = vals[eps_values[0]][0]
        G0 = arb_mat([[bm.sym_get(g0, i, j) for j in R3] for i in R3]).inv()
        fit = eps_values[:-1]
        u = [arb(e) ** 2 if e else arb(0) for e in fit]
        V = arb_mat([[u[r] ** c if c else arb(1) for c in range(len(fit))] for r in range(len(fit))]).inv()
        out, check = [], 0.0
        for k in range(N + 1):
            per = {}
            for pq in PAIRS:
                f = [vals[e][2 * k][pq] for e in fit]
                coef = [sum((V[c, r] * f[r] for r in range(len(fit))), arb(0)) for c in range(len(fit))]
                per[pq] = coef
                e = eps_values[-1]
                pred = sum((coef[c] * (arb(e) ** (2 * c) if c else arb(1)) for c in range(len(fit))), arb(0))
                check = max(check, float(abs(pred - vals[e][2 * k][pq]).upper()))
            secs = []
            for c in range(len(fit)):
                S = arb_mat([[bm.sym_get({pq: per[pq][c] for pq in PAIRS}, i, j) for j in R3] for i in R3])
                M = G0 * S
                secs.append([[M[i, j] for j in R3] for i in R3])
            out.append(secs)
        return out, check


SECTOR_POINTS = {"P1": (("0.17", "0.11", "0.5"), "0.7662"), "P2": (("0.26", "-0.09", "0.12"), "0.767"),
                 "P3": (("0.6", "0.3", "0.9"), "0.7662")}


def sector_record(name, N=6):
    """Frobenius norms (midpoints) of the curl-degree sectors G_{k,j} (mixed) at a named point."""
    import math
    x0s, ds = SECTOR_POINTS[name]
    t0 = time.time()
    secs, check = sectors(x0s, ds, N)
    norms = [[math.sqrt(sum(float(M[i][j].mid()) ** 2 for i in R3 for j in R3)) for M in secs[k]] for k in range(N + 1)]
    rad = max(float(M[i][j].rad()) for k in range(N + 1) for M in secs[k] for i in R3 for j in R3)
    return {"x0": x0s, "d": ds, "N": N, "sector_norms": norms, "extrapolation_check": check,
            "max_ball_radius": rad, "seconds": time.time() - t0}


def main_sectors(out=None, names=("P1", "P2", "P3"), N=6):
    import json
    import os
    from multiprocessing import Pool
    out = out or os.path.join(bm.DATA, "bl_bianchi_sectors.json")
    with Pool(len(names)) as pool:
        rows = pool.starmap(sector_record, [(n, N) for n in names])
    res = dict(zip(names, rows))
    res["replay"] = "python -m research.calc.bl_bianchi sectors"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "sectors":
        main_sectors()


def energy_identity_residual(psi_fn, W, one, n_max):
    """Check, at the expansion point and for t^0..t^n_max, the exact energy identity (notes/19 §5)
        d_t(sqrt(g) w) - 2 d_k(sqrt(g) F^k) = sqrt(g) (-theta w + 6 tr(sigma (E^2 + H^2))),
        w = |E|^2 + |H|^2,  F^k = eps_i^{kl} E^{ij} H_lj,  sigma = K - (theta/3) gamma.
    Returns the list of residual values (one per order)."""
    J, s = bianchi_jets(psi_fn, W, one)
    g, K, E, H = s["g"], s["K"], s["E"], s["H"]
    psi = psi_fn(J)
    Gi = bm.inverse3(J, g, one)
    isq = _inv_sqrt_det(J, g, psi, one)
    sq = J.mul(isq, _det3(J, g))                      # sqrt(det g) = det * det^(-1/2)
    th = bm.trace_prod(J, Gi, K)
    EE = bm.matmul_sym(J, E, Gi, E)
    HH = bm.matmul_sym(J, H, Gi, H)
    w = J.add_(bm.trace_prod(J, Gi, EE), bm.trace_prod(J, Gi, HH))
    S = {pq: J.sub(K[pq], J.scale(J.mul(th, g[pq]), one / 3)) for pq in PAIRS}
    Q = {pq: J.add_(EE[pq], HH[pq]) for pq in PAIRS}
    trSQ = _trace_full(J, Gi, _full(J, S, Gi, Q))
    rhs_ = J.mul(sq, J.add_(J.scale(J.mul(th, w), -one), J.scale(trSQ, 6 * one)))
    # F^k = eps^{m k l} g_mi ... : eps_i^{kl} E^{ij} H_lj = [mkl] isq g_mi E^{ij} H_lj = [mkl] isq E_m^j H_lj
    EH = _full(J, E, Gi, H)                           # E_m^j H_jl  -> (m, l)
    flux = []
    for k in R3:
        acc = {}
        for (m, kk, l), sg in PERM.items():
            if kk == k:
                acc = J.add_(acc, J.scale(EH[(m, l)], sg))
        flux.append(J.mul(J.mul(acc, isq), sq))      # sqrt(g) F^k
    # d_t of a (X, Y, Z, t)-jet: shift t^(n+1) -> (n+1) t^n
    sw = J.mul(sq, w)
    res = []
    zero = 0 * one
    for n in range(n_max + 1):
        dt = J.coef(sw, (0, 0, 0, n + 1), zero) * (n + 1)
        div = sum((J.coef(J.d(flux[k], k), (0, 0, 0, n), zero) for k in R3), zero)
        res.append(dt - 2 * div - J.coef(rhs_, (0, 0, 0, n), zero))
    return res
