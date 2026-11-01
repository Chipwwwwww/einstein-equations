"""Experiment H2 (notes/18 §4): can any slicing of the BL development with a bounded maximal-type initial
lapse have conformally flat leaves? (Session 12 excluded only the geodesic slicing, N = 1.)

Facts used. Zero-shift slicing with lapse N through the time-symmetric slice (K = 0 at t = 0): with the
repo sign d_t gamma = 2 N K, d_t K_ij = D_iD_j N - N R_ij + O(K), so the leaves are
    gamma(t) = gamma_0 + t^2 h + O(t^3),   h = N D_iD_j N - N^2 Ric(gamma_0),   N = N(0, x).
A shift only adds L_beta gamma_0 at this order. The 3D Cotton tensor C_ijk = nabla_k S_ij - nabla_j S_ik
(S the Schouten tensor) vanishes iff a metric is conformally flat; C[gamma_0] = 0 for BL data, so
    C[gamma(t)] = t^2 DC[h] + O(t^3),
and conformally flat leaves to O(t^2) require DC[h] = 0 (DC the linearised Cotton operator at gamma_0).
DC is invariant under h -> h + L_xi gamma_0 + phi gamma_0 (diffeomorphism and conformal invariance), so
with N DDN = (1/2) L_{N grad N} gamma_0 - dN (x) dN:   DC[h] = -DC[dN (x) dN + N^2 Ric]  (independent check).

The family. For gamma_0 = psi^4 delta, Delta_gamma(phi/psi) = psi^-5 Delta_flat phi, so the lapses with
D^2 N = 0 at t = 0 that are bounded on all three ends are exactly
    N_c = (c_0 + c_+ / r_+ + c_- / r_-) / psi,      c in R^3
(Bocher/Liouville: a flat-harmonic phi on R^3 minus the punctures with bounded phi/psi is of this
form). It contains the geodesic slicing (c = (1, 1/4, 1/4): N = 1) and the odd lapse
(1, -1/4, -1/4). DC[h_c] is a quadratic form in c with tensor values; its components Q^m(c) are conics.
Excluding the whole family at one point x amounts to certifying that the conics have no common zero
on RP^2: every face {c_i = 1, |c_j| <= 1} (i = 0, 1, 2; c and -c give the same h) is covered by boxes,
on each of which some component's ball enclosure excludes 0.

Run: python -m research.calc.bl_foliation   (writes research/calc/data/bl_foliation.json)
"""

import json
import os
import time

from flint import arb, fmpq

from research.calc import bl_metric as bm
from research.calc import bl_structure as st

PAIRS = bm.PAIRS
R3 = range(3)
W = 5          # spatial weight of the jets: Cotton (3 derivatives) of h (2 derivatives of N) needs N to weight 5


def _jets():
    """Jets in (X, Y, Z, eps) with eps of weight 1: eps marks the linear perturbation gamma_0 + eps h."""
    return bm.Jets(W, sw=1)


def inv_r(J, x0, p, r0, one):
    """1 / |x - p| as a jet at x0."""
    j = bm.psi_jet(J, x0, [(one, p, r0)], one)
    return J.sub(j, {0: one})


def hessian(J, N, Gam):
    """D_i D_j N (dict over PAIRS) with Christoffels Gam of gamma_0."""
    dN = [J.d(N, a) for a in R3]
    out = {}
    for (i, j) in PAIRS:
        v = J.d(dN[j], i)
        for k in R3:
            v = J.sub(v, J.mul(Gam[(k, i, j)], dN[k]))
        out[(i, j)] = v
    return out


def lin_cotton(J, g0, h, one):
    """DC[h] at the point: the eps^1 point coefficient of the Cotton tensor of gamma_0 + eps h, as a dict
    over (i, j, k) with j < k (C_ijk is antisymmetric in j, k)."""
    zero = 0 * one
    g = {pq: J.add_(g0[pq], {J.smap[i]: v for i, v in h[pq].items() if i in J.smap}) for pq in PAIRS}
    Gi = bm.inverse3(J, g, one)
    Ric = bm.ricci3(J, g, one)
    R = bm.trace_prod(J, Gi, Ric)
    S = {pq: J.sub(Ric[pq], J.scale(J.mul(R, g[pq]), one / 4)) for pq in PAIRS}
    Gam = st.christoffel(J, g, Gi, one)
    dS = st.cov_sym(J, S, Gam)
    out = {}
    for i in R3:
        for j in R3:
            for k in R3:
                if j < k:
                    C = J.sub(dS[(i, j, k)], dS[(i, k, j)])
                    out[(i, j, k)] = J.coef(C, (0, 0, 0, 1), zero)
    return out, J.coef(J.sub(dS[(0, 0, 1)], dS[(0, 1, 0)]), (0, 0, 0, 0), zero)


def setup(x0, d, one, exact_r=None, centres=None, lapse_centres=None):
    """gamma_0 = psi^4 delta, its Ricci tensor and Christoffels, and the lapse basis N_a = phi_a / psi with
    phi = 1 and 1/|x - q| for q in lapse_centres (default: the punctures)."""
    J = _jets()
    cs = centres if centres is not None else bm.bl_centres(d, one, x0, exact_r)
    psi = bm.psi_jet(J, x0, cs, one)
    p2 = J.mul(psi, psi)
    p4 = J.mul(p2, p2)
    g0 = {pq: (dict(p4) if pq[0] == pq[1] else {}) for pq in PAIRS}
    Gi = bm.inverse3(J, g0, one)
    Ric = bm.ricci3(J, g0, one)
    Gam = st.christoffel(J, g0, Gi, one)
    ipsi = J.inv(psi, one)
    lc = lapse_centres if lapse_centres is not None else [(c[1], c[2]) for c in cs]
    basis = [ipsi] + [J.mul(inv_r(J, x0, q, r0, one), ipsi) for q, r0 in lc]
    return J, g0, Ric, Gam, basis


def quadratic_forms(x0, d, one, exact_r=None, centres=None, lapse_centres=None, check=True):
    """Q[(a, b)] = DC[h_ab] (dict over Cotton components), a <= b, with
    h_ab = (1/2)(N_a DDN_b + N_b DDN_a) - N_a N_b Ric, so DC[h_c] = sum_{a,b} c_a c_b Q_ab (Q_ba = Q_ab).
    If check, also returns the maximal deviation from the independent form
    -DC[(1/2)(dN_a dN_b + dN_b dN_a) + N_a N_b Ric] and the zeroth-order Cotton of gamma_0 (must be 0)."""
    J, g0, Ric, Gam, basis = setup(x0, d, one, exact_r, centres, lapse_centres)
    H = [hessian(J, N, Gam) for N in basis]
    dN = [[J.d(N, a) for a in R3] for N in basis]
    n = len(basis)
    Q, dev, c0 = {}, 0.0, None
    for a in range(n):
        for b in range(a, n):
            NN = J.mul(basis[a], basis[b])
            h = {pq: J.sub(J.scale(J.add_(J.mul(basis[a], H[b][pq]), J.mul(basis[b], H[a][pq])), one / 2),
                           J.mul(NN, Ric[pq])) for pq in PAIRS}
            Q[(a, b)], c0 = lin_cotton(J, g0, h, one)
            if check:
                h2 = {(i, j): J.scale(J.add_(J.scale(J.add_(J.mul(dN[a][i], dN[b][j]), J.mul(dN[b][i], dN[a][j])),
                                                     one / 2), J.mul(NN, Ric[(i, j)])), -one)
                      for (i, j) in PAIRS}
                Q2, _ = lin_cotton(J, g0, h2, one)
                for key in Q2:
                    e = Q[(a, b)][key] - Q2[key]
                    dev = max(dev, float(abs(e).upper()) if hasattr(e, "upper") else (0.0 if e == 0 else float("inf")))
    return Q, dev, c0


def gauge_controls(x0, d, one):
    """DC of pure-gauge perturbations L_xi gamma_0 (polynomial xi) and phi gamma_0 (polynomial phi):
    both must vanish (ball radius). Returns the max |.| over components."""
    J, g0, Ric, Gam, basis = setup(x0, d, one)
    X = [J.var(i, one) for i in R3]
    xi = [J.add_(J.mul(X[1], X[2]), J.scale(X[0], one * 3)), J.add_(J.mul(X[0], X[0]), {0: one}),
          J.sub(J.mul(X[0], X[1]), J.scale(X[2], one * 2))]
    # (L_xi g)_ij = xi^k d_k g_ij + g_kj d_i xi^k + g_ik d_j xi^k
    L = {}
    for (i, j) in PAIRS:
        v = {}
        for k in R3:
            v = J.add_(v, J.mul(xi[k], J.d(bm.sym_get(g0, i, j), k)))
            v = J.add_(v, J.mul(bm.sym_get(g0, k, j), J.d(xi[k], i)))
            v = J.add_(v, J.mul(bm.sym_get(g0, i, k), J.d(xi[k], j)))
        L[(i, j)] = v
    phi = J.add_(J.mul(X[0], X[1]), J.add_(J.scale(X[2], one * 5), {0: one * 2}))
    C1, _ = lin_cotton(J, g0, L, one)
    C2, _ = lin_cotton(J, g0, {pq: J.mul(phi, g0[pq]) for pq in PAIRS}, one)
    # a non-gauge perturbation must give a nonzero DC (the tool is not identically zero)
    C3, _ = lin_cotton(J, g0, {(0, 0): J.mul(X[1], J.mul(X[1], X[2])), (0, 1): {}, (0, 2): {},
                               (1, 1): {}, (1, 2): {}, (2, 2): {}}, one)
    mx = lambda C: max(float(abs(v).upper()) for v in C.values())  # noqa: E731
    return {"lie": mx(C1), "conformal": mx(C2), "non_gauge_example": max(float(abs(v).mid()) for v in C3.values())}


def form_value(Qm, c):
    """sum_{a <= b} (2 - delta_ab) c_a c_b Q_ab for one component (Q dict over (a, b))."""
    n = len(c)
    acc = 0 * c[0]
    for a in range(n):
        for b in range(a, n):
            acc = acc + (Qm[(a, b)] * c[a] * c[b]) * (1 if a == b else 2)
    return acc


def certify_no_common_zero(Qs, max_depth=18, prec=128):
    """Branch and bound on the faces {c_i = 1, |c_j| <= 1} of the cube (i over the basis indices): each
    box must have a component, of the forms Q at any of the points in Qs (a list of Q, or one Q), whose
    ball enclosure excludes 0. Returns (ok, boxes, deepest, failures)."""
    if isinstance(Qs, dict):
        Qs = [Qs]
    n = 1 + max(b for (_, b) in Qs[0])
    comps = [{ab: Q[ab][k] for ab in Q} for Q in Qs for k in next(iter(Q.values()))]
    boxes, deepest, fails = 0, 0, []
    with bm.working_precision(prec):
        for i in range(n):
            free = [j for j in range(n) if j != i]
            stack = [([(-1.0, 1.0)] * len(free), 0)]
            while stack:
                box, depth = stack.pop()
                c = [None] * n
                c[i] = arb(1)
                for j, (lo, hi) in zip(free, box):
                    c[j] = arb((lo + hi) / 2, (hi - lo) / 2)
                boxes += 1
                if any(not form_value(Qm, c).contains(0) for Qm in comps):
                    deepest = max(deepest, depth)
                    continue
                if depth >= max_depth:
                    fails.append((i, box))
                    continue
                # bisect the widest side
                w = max(range(len(box)), key=lambda k: box[k][1] - box[k][0])
                lo, hi = box[w]
                m = (lo + hi) / 2
                for half in ((lo, m), (m, hi)):
                    nb = list(box)
                    nb[w] = half
                    stack.append((nb, depth + 1))
    return not fails, boxes, deepest, fails[:5]


def component_rank(Q, tol=1e-12):
    """Numerical rank of the (components x pairs) coefficient matrix: the number of independent conics.
    For the axisymmetric family it is 2 (the Cotton-York tensor of a twist-free axisymmetric metric
    has only the two components odd under phi -> -phi), so two conics in RP^2: a single point cannot
    exclude the family, and real common zeros exist on open sets of x."""
    import numpy as np
    keys = list(next(iter(Q.values())).keys())
    A = np.array([[float(Q[ab][k].mid()) for ab in sorted(Q)] for k in keys])
    sv = np.linalg.svd(A, compute_uv=False)
    return int((sv > tol * sv[0]).sum()), [float(v) for v in sv]


def pointwise_zeros(Q, tries=200, seed=1):
    """Numerical common zeros c* (normalised c_0 = 1) of the conics at one point (exploratory: BFGS on
    sum_m Q^m(c)^2 / |c|^4 from random starts). Returns sorted (residual, c*) for distinct minima."""
    import numpy as np
    from scipy.optimize import minimize
    n = 1 + max(b for (_, b) in Q)
    Ms = []
    for k in next(iter(Q.values())):
        M = np.zeros((n, n))
        for (a, b), F in Q.items():
            M[a, b] = M[b, a] = float(F[k].mid())
        Ms.append(M)
    f = lambda c: sum((c @ M @ c) ** 2 for M in Ms) / (c @ c) ** 2  # noqa: E731
    rng = np.random.default_rng(seed)
    found = []
    for _ in range(tries):
        r = minimize(f, rng.normal(size=n), method="BFGS", options={"gtol": 1e-14})
        if r.fun > 1e-16 or abs(r.x[0]) < 1e-9:
            continue
        c = r.x / r.x[0]
        if all(np.abs(c - g).max() > 1e-5 for _, g in found):
            found.append((float(r.fun), c))
    return sorted((f0, [float(v) for v in c]) for f0, c in found)


def geodesic_member(Q, one):
    """DC[h] for the geodesic lapse c = (1, 1/4, 1/4) (N = 1, h = -Ric): must equal -DC[Ric] (nonzero)."""
    c = [one, one / 4, one / 4]
    return {k: form_value({ab: Q[ab][k] for ab in Q}, c) for k in next(iter(Q.values()))}


def main(out=os.path.join(bm.DATA, "bl_foliation.json")):
    res = {}
    t0 = time.time()
    # exact rational point (BL plane): forms in Q, gauge-identity check exact
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    Q, dev, c0 = quadratic_forms(x0, d, fmpq(1), exact_r=rr)
    res["plane_point_exact"] = {"x0": [str(v) for v in x0], "d": str(d), "identity_deviation": dev,
                                "zeroth_order_cotton": str(c0),
                                "nonzero_components": sum(1 for k in next(iter(Q.values()))
                                                          if any(Q[ab][k] != 0 for ab in Q))}
    print("plane", res["plane_point_exact"], f"{time.time() - t0:.0f}s", flush=True)
    with bm.working_precision(256):
        # pairs of generic points at the same separation (a single point cannot decide: rank 2)
        groups = {"d = 0.7662": ("0.7662", [("0.17", "0.11", "0.5"), ("0.6", "0.3", "0.9")]),
                  "d = 0.767": ("0.767", [("0.26", "-0.09", "0.12"), ("0.05", "0.3", "-0.2")]),
                  "d = 1.5": ("1.5", [("0.4", "0.2", "1.1"), ("0.3", "0.5", "-0.6")])}
        for name, (ds, xs_list) in groups.items():
            Qs, info = [], []
            for xs in xs_list:
                x = tuple(arb(v) for v in xs)
                Q, dev, c0 = quadratic_forms(x, arb(ds), arb(1))
                ok1, nb1, _, _ = certify_no_common_zero(Q, max_depth=14)
                rank, sv = component_rank(Q)
                geo = geodesic_member(Q, arb(1))
                info.append({"x0": xs, "identity_deviation": dev, "zeroth_order_cotton_abs": float(abs(c0).upper()),
                             "independent_conics": rank, "singular_values": sv,
                             "single_point_certificate": ok1,
                             "pointwise_common_zeros_numerical": pointwise_zeros(Q)[:4],
                             "geodesic_member_max_abs": max(float(abs(v).mid()) for v in geo.values())})
                Qs.append(Q)
            ok, nb, depth, fails = certify_no_common_zero(Qs)
            res[name] = {"points": info, "no_common_zero_certified": ok, "boxes": nb, "max_depth": depth,
                         "failures": [str(f) for f in fails]}
            print(name, ok, nb, depth, [(i["independent_conics"], i["single_point_certificate"],
                                         i["pointwise_common_zeros_numerical"][:2]) for i in info],
                  f"{time.time() - t0:.0f}s", flush=True)
        pts = {"G1": (("0.17", "0.11", "0.5"), "0.7662")}
        res["gauge_controls_G1"] = gauge_controls(tuple(arb(v) for v in pts["G1"][0]), arb("0.7662"), arb(1))
        # Schwarzschild controls: lapses (c_0 + c_1/r)/psi are spherical, so DC = 0 identically;
        # an off-centre lapse 1/|x - q| breaks the symmetry and must give DC != 0
        xs2 = (arb("0.6"), arb("0.2"), arb("0.9"))
        cs = [(arb(1) / 2, (arb(0),) * 3, (arb("0.36") + arb("0.04") + arb("0.81")).sqrt())]
        Qs, devs, _ = quadratic_forms(xs2, None, arb(1), centres=cs)
        q = (arb(0), arb(0), arb("0.3"))
        rq = (arb("0.36") + arb("0.04") + arb("0.36")).sqrt()
        Qo, _, _ = quadratic_forms(xs2, None, arb(1), centres=cs, lapse_centres=[((arb(0),) * 3, cs[0][2]), (q, rq)])
        res["schwarzschild_control"] = {
            "spherical_family_max_abs": max(float(abs(v).upper()) for F in Qs.values() for v in F.values()),
            "identity_deviation": devs,
            "off_centre_lapse_max_abs": max(float(abs(v).mid()) for F in Qo.values() for v in F.values())}
        print("controls", res["gauge_controls_G1"], res["schwarzschild_control"], flush=True)
    res["seconds"] = time.time() - t0
    res["replay"] = "python -m research.calc.bl_foliation"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    main()
