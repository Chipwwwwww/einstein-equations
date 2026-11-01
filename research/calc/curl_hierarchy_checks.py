"""Exact checks for research/curl_hierarchy_convergence_sprint_01.md.

(1) CGR (59) against the repository's (curl H)_3 law. Clifton-Gregoris-Rosquist 2017, eq. (59):
    d^3/dt^3 (curl H)|_0 = 4 curl curl (E^2) for time-symmetric vacuum data with curl E = 0 (their (46)-(47)).
    With the Taylor coefficient (curl H)_3 = (1/3!) d^3/dt^3 (curl H)|_0 this reads (curl H)_3 = (2/3) curl curl(Ric^2).
    Checked: (a) exactly at BL rational points against the exact chain (curl H)_3 = 4 gamma_0 (E_4 - E_4^sil);
    (b) for ALL conformally flat time-symmetric vacuum data by the harmonic-jet reduction of notes/19 §3.3
    (same 72 evaluation points; weight 6), against the four-term law.
(2) Principal symbol of the coordinate E-H form (bianchi_symbol.coordinate_symbol) in exact rational
    arithmetic at gamma = I, E = diag(-1, -1, 2), H = 0, for xi = e_3 and xi = (3/5, 0, 4/5): algebraic and
    geometric multiplicity of the eigenvalue 0.

Run: python -m research.calc.curl_hierarchy_checks
"""

import itertools

import sympy as sp
from flint import fmpq

from research.calc import bl_bianchi as bb
from research.calc import bl_metric as bm
from research.calc import bl_silent as sl
from research.calc import bl_structure as st
from research.calc import cf_universal as cu

PAIRS = bm.PAIRS


def curlcurl_ric2(psi_fn, one):
    """(2/3) curl curl (Ric.Ric) at the expansion point (lower indices), gamma_0 = psi^4 delta."""
    J = bm.Jets(6, sw=7)                       # spatial jets of weight 6
    psi = psi_fn(J)
    p2 = J.mul(psi, psi)
    p4 = J.mul(p2, p2)
    g = {pq: (dict(p4) if pq[0] == pq[1] else {}) for pq in PAIRS}
    Gi = bm.inverse3(J, g, one)
    Gam = st.christoffel(J, g, Gi, one)
    isq = J.inv(J.mul(p4, p2), one)            # det(psi^4 delta)^(-1/2) = psi^-6
    Ric = bm.ricci3(J, g, one)
    R2 = bm.matmul_sym(J, Ric, Gi, Ric)
    C2 = bb.curl(J, bb.curl(J, R2, g, Gam, isq, one), g, Gam, isq, one)
    zero = 0 * one
    return {pq: J.coef(C2[pq], (0, 0, 0, 0), zero) * one * 2 / 3 for pq in PAIRS}


def cgr_bl_points():
    out = {}
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    za = fmpq(31, 50)
    for name, x, kw in (("plane", x0, dict(exact_r=rr)), ("axis", (fmpq(0), fmpq(0), za), dict(exact_r=[za - d / 2, za + d / 2]))):
        lhs, _ = sl.curl_h3_law(x, d, fmpq(1), **kw)
        cc = curlcurl_ric2(bb.bl_psi(x, d, fmpq(1), **kw), fmpq(1))
        out[name] = {"equal": all(lhs[p] == cc[p] for p in PAIRS), "nonzero": any(lhs[p] != 0 for p in PAIRS)}
    return out


def cgr_universal():
    """(2/3) curl curl Ric^2 - law(curl H)_3 at the 72 harmonic-jet points (weight 6 subset)."""
    ok, n = True, 0
    for label, weights, terms in cu.evaluation_set():
        if 6 not in weights:
            continue
        with cu.custom_data(cu.phi_fn(terms)):
            _, law = sl.curl_h3_law((fmpq(0),) * 3, None, fmpq(1), centres=[])
        cc = curlcurl_ric2(cu.phi_fn(terms), fmpq(1))
        ok &= all(law[p] == cc[p] for p in PAIRS)
        n += 1
    return {"all_equal": ok, "evaluations": n}


def symbol_exact(xi, E):
    """Exact 24x24 symbol of the coordinate form at gamma = I, H = 0 (as bianchi_symbol.coordinate_symbol)."""
    eps = sp.zeros(3, 3)
    epsT = {}
    for (i, j, k), s in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
        epsT[(i, j, k)] = s
    e = lambda i, j, k: epsT.get((i, j, k), 0)  # noqa: E731
    SYM = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    def to_mat(v):
        A = sp.zeros(3, 3)
        for c, (i, j) in enumerate(SYM):
            A[i, j] = A[j, i] = v[c]
        return A

    def vec(A):
        return [A[i, j] for (i, j) in SYM]

    def curl_sym(A):
        T = sp.Matrix(3, 3, lambda a, b: sum(e(a, c, d) * xi[c] * A[d, b] for c in range(3) for d in range(3)))
        return (T + T.T) / 2

    def curl_gamma(dg, A):
        G = {(m, k, l): sp.Rational(1, 2) * (xi[k] * dg[m, l] + xi[l] * dg[m, k] - xi[m] * dg[k, l])
             for m, k, l in itertools.product(range(3), repeat=3)}
        T = sp.Matrix(3, 3, lambda i, j: -sum(e(i, k, l) * G[(m, k, j)] * A[l, m]
                                               for k in range(3) for l in range(3) for m in range(3)))
        return (T + T.T) / 2

    P = sp.zeros(24, 24)
    for c in range(6):
        v = [0] * 6
        v[c] = 1
        A = to_mat(v)
        cs = vec(curl_sym(A))
        cg = vec(curl_gamma(A, E))
        for r in range(6):
            P[12 + r, 18 + c] = cs[r]
            P[18 + r, 12 + c] = -cs[r]
            P[18 + r, c] = -cg[r]
    lam = sp.Symbol("lam")
    cp = P.charpoly(lam).as_expr()
    alg = sp.Poly(cp, lam).as_expr().as_poly(lam)
    mult = 0
    while alg.eval(0) == 0:
        alg = sp.Poly(sp.quo(alg.as_expr(), lam), lam)
        mult += 1
    geo = 24 - P.rank()
    return {"algebraic": mult, "geometric": geo, "rank_P": P.rank(), "rank_P2": (P * P).rank()}


def main():
    res = {}
    res["cgr59_bl_exact"] = cgr_bl_points()
    print("CGR (59) at BL rational points:", res["cgr59_bl_exact"], flush=True)
    res["cgr59_universal"] = cgr_universal()
    print("CGR (59) vs four-term law, all conformally flat vacuum data:", res["cgr59_universal"], flush=True)
    E = sp.diag(-1, -1, 2)
    for name, xi in (("xi=e3", [0, 0, 1]), ("xi=(3/5,0,4/5)", [sp.Rational(3, 5), 0, sp.Rational(4, 5)])):
        res[name] = symbol_exact(xi, E)
        print("symbol", name, res[name], flush=True)
    return res


if __name__ == "__main__":
    main()
