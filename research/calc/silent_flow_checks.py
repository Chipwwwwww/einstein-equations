"""Check for research/nonlinear_silent_flow_sprint_02.md: where the metric variation of curl enters.

The E-H system (bl_bianchi) uses curl_gamma with the FULL metric. Replace, inside both curls only,
gamma(t) by the silent metric gamma_sil(t) (the eps = 0 solution): the 'frozen-metric' system. The
difference is exactly the metric-variation part eps * D_gamma(curl)[gamma - gamma_sil] (all orders:
volume form, index lowering, Christoffels, inverse metric). Predicted from the lag lemma
(gamma - gamma_sil = O(eps^2 t^6)) and Ebar != 0:
    H differs first at t^7 (pure eps^3),  E first at t^8 (pure eps^4),  gamma identical through t^8
(first possible difference t^10). Checked in 256-bit balls at a generic BL point (W = 8).

Run: python -m research.calc.silent_flow_checks
"""

from flint import arb, fmpq

from research.calc import bl_bianchi as bb
from research.calc import bl_metric as bm
from research.calc import bl_structure as st

PAIRS = bm.PAIRS


def frozen_jets(psi_fn, W, one, eps, g_frozen):
    """bianchi_jets with curl taken in the metric g_frozen (jets) instead of the evolving metric."""
    J = bm.Jets(W, sw=1)
    _, ref = bb.bianchi_jets(psi_fn, W, one, eps=0)       # only for the initial state layout
    state = {k: {pq: {i: v for i, v in ref[k][pq].items() if J.mons[i][3] == 0} for pq in PAIRS} for k in ref}
    psi = psi_fn(J)
    Gi_f = bm.inverse3(J, g_frozen, one)
    Gam_f = st.christoffel(J, g_frozen, Gi_f, one)
    isq_f = bb._inv_sqrt_det(J, g_frozen, psi, one)
    for n in range(W):
        g, K, E, H = state["g"], state["K"], state["E"], state["H"]
        Gi = bm.inverse3(J, g, one)
        th = bm.trace_prod(J, Gi, K)
        KK = bm.matmul_sym(J, K, Gi, K)
        KE, KH = bb._full(J, K, Gi, E), bb._full(J, K, Gi, H)
        trKE, trKH = bb._trace_full(J, Gi, KE), bb._trace_full(J, Gi, KH)
        cH = bb.curl(J, H, g_frozen, Gam_f, isq_f, one)
        cE = bb.curl(J, E, g_frozen, Gam_f, isq_f, one)
        F = {"g": {pq: J.scale(K[pq], 2 * one) for pq in PAIRS}, "K": {pq: J.sub(KK[pq], E[pq]) for pq in PAIRS},
             "E": {}, "H": {}}
        for (i, j) in PAIRS:
            sKE, sKH = J.add_(KE[(i, j)], KE[(j, i)]), J.add_(KH[(i, j)], KH[(j, i)])
            F["E"][(i, j)] = J.add_(J.add_(J.scale(cH[(i, j)], eps), J.scale(J.mul(th, E[(i, j)]), -2 * one)),
                                    J.sub(J.scale(sKE, one * fmpq(5, 2)), J.mul(g[(i, j)], trKE)))
            F["H"][(i, j)] = J.add_(J.add_(J.scale(cE[(i, j)], -eps), J.scale(J.mul(th, H[(i, j)]), -2 * one)),
                                    J.sub(J.scale(sKH, one * fmpq(5, 2)), J.mul(g[(i, j)], trKH)))
        for key in state:
            for pq in PAIRS:
                sl = J.drop_above(J.s_slice(F[key][pq], n), W - 1 - n)
                J.put_s(state[key][pq], n + 1, J.scale(sl, one / (n + 1)), W)
    return J, state


def main(xs=("0.17", "0.11", "0.5"), ds="0.7662", W=8):
    with bm.working_precision(256):
        x = tuple(arb(v) for v in xs)
        one = arb(1)
        psi_fn = bb.bl_psi(x, arb(ds), one)
        J0, s0 = bb.bianchi_jets(psi_fn, W, one, eps=0)
        diffs = {}
        for e in (1, 2):
            J1, full = bb.bianchi_jets(psi_fn, W, one, eps=arb(e))
            J2, froz = frozen_jets(psi_fn, W, one, arb(e), s0["g"])
            row = {}
            for key in ("H", "E", "g"):
                row[key] = [max(float(abs(J1.coef(full[key][pq], (0, 0, 0, n), arb(0))
                                          - J2.coef(froz[key][pq], (0, 0, 0, n), arb(0))).mid()) for pq in PAIRS)
                            for n in range(W + 1)]
            diffs[e] = row
            print("eps", e, {k: [f"{v:.2e}" for v in row[k]] for k in row}, flush=True)
        r7 = diffs[2]["H"][7] / diffs[1]["H"][7]
        r8 = diffs[2]["E"][8] / diffs[1]["E"][8]
        print(f"H t^7 diff ratio eps=2/eps=1: {r7:.6f} (eps^3 -> 8);  E t^8 ratio: {r8:.6f} (eps^4 -> 16)")
        return diffs


if __name__ == "__main__":
    main()
