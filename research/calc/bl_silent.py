"""Experiment H1 (notes/18): the silent (local tidal) part of the BL development and its remainder.

Exact facts used (geodesic slicing, comoving coordinates, vacuum, K_ij = (1/2) d_t gamma_ij, mixed
indices K = K^i_j, E = E^i_j with E_ij = R_ij(gamma) + K K_ij - K_ik K^k_j the electric Weyl tensor):

    d_t gamma = 2 gamma K                       (exact)
    d_t K     = -E - K^2                        (exact: the evolution equation rewritten with E)
    d_t E     = -theta E + 3 (sigma E)^TF + [E, K] + curl H     (Bianchi identities; curl H is the
                                                only term with spatial derivatives)

with theta = tr K, sigma = K - theta/3. Dropping curl H gives the *silent* system: a closed ODE along
each comoving worldline, started from K = 0, E = E_0 = Ric(gamma_0) (mixed), gamma = gamma_0. Its
solution keeps every matrix a polynomial in E_0 (they commute), so gamma_silent = gamma_0 M(t), and
the whole silent metric is determined pointwise by (gamma_0, E_0). For one hole curl H = 0 and the
silent solution is the exact (cycloid) free fall.

For conformally flat time-symmetric data H = O(t^3), so curl H perturbs E from t^4, K from t^5 and
gamma from t^6: gamma - gamma_silent = O(t^6). The experiments below test this, identify the t^6
remainder, and compare the silent and true series to high order (analytic structure).

Run: python -m research.calc.bl_silent   (writes research/calc/data/bl_silent.json)
"""

import json
import math
import os
import time

import numpy as np
from flint import arb, fmpq

from research.calc import bl_metric as bm

R3 = range(3)


def _mm(A, B, zero):
    return [[sum((A[i][k] * B[k][j] for k in R3), zero) for j in R3] for i in R3]


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in R3] for i in R3]


def _scale(A, c):
    return [[A[i][j] * c for j in R3] for i in R3]


def _tr(A):
    return A[0][0] + A[1][1] + A[2][2]


def silent_series(E0, n_max, one):
    """Taylor coefficients (in t) of the silent system with K(0) = 0, E(0) = E0, M(0) = I:
    K' = -E - K^2, E' = -theta E + 3 (sigma E - tr(sigma E) I / 3), M' = 2 M K (all matrices commute).
    Exact in the arithmetic of `one` (fmpq or arb). Returns lists K[n], E[n], M[n] (3x3), n <= n_max."""
    zero = 0 * one
    Id = [[one if i == j else zero for j in R3] for i in R3]
    Z = [[zero] * 3 for _ in R3]
    K, E, M = [Z], [E0], [Id]
    for n in range(n_max):
        # coefficients of t^n of the right-hand sides
        KK = Z
        for i in range(n + 1):
            KK = _add(KK, _mm(K[i], K[n - i], zero))
        rhsK = _scale(_add(E[n], KK), -one)
        rhsE = Z
        for i in range(n + 1):
            Ki, Ej = K[i], E[n - i]
            th = _tr(Ki)
            sig = _add(Ki, _scale(Id, -th / 3))
            sE = _mm(sig, Ej, zero)
            term = _add(_add(_scale(Ej, -th), _scale(sE, 3 * one)), _scale(Id, -_tr(sE)))
            rhsE = _add(rhsE, term)
        rhsM = Z
        for i in range(n + 1):
            rhsM = _add(rhsM, _scale(_mm(M[i], K[n - i], zero), 2 * one))
        inv = one / (n + 1)
        K.append(_scale(rhsK, inv))
        E.append(_scale(rhsE, inv))
        M.append(_scale(rhsM, inv))
    return K, E, M


def true_and_silent(x0, d, N, one, exact_r=None, centres=None):
    """At x0: the true coefficients gamma_k (lower, t^(2k), k <= N) from the jet recursion, and the
    silent ones gamma_0 M_(2k), with E_0 = -G_0 gamma_1 (mixed Ric). Returns dicts of 3x3 lists."""
    J, g, _ = bm.series_jets(x0, d, N, one, exact_r=exact_r, centres=centres)
    zero = 0 * one
    vals, _, _ = bm.coefficients(J, g, zero, N)
    true = [[[bm.sym_get(vals[k], i, j) for j in R3] for i in R3] for k in range(N + 1)]
    g0 = true[0]
    G0 = _inv3(g0, one)
    E0 = _scale(_mm(G0, true[1], zero), -one)            # mixed Ric = -G0 gamma_1
    _, _, M = silent_series(E0, 2 * N, one)
    silent = [_mm(g0, M[2 * k], zero) for k in range(N + 1)]
    odd = [M[2 * k + 1] for k in range(N)]               # must vanish (time reflection)
    return true, silent, E0, odd


def _inv3(A, one):
    a, b, c = A[0]
    d_, e, f = A[1]
    g, h, i = A[2]
    det = a * (e * i - f * h) - b * (d_ * i - f * g) + c * (d_ * h - e * g)
    adj = [[e * i - f * h, c * h - b * i, b * f - c * e],
           [f * g - d_ * i, a * i - c * g, c * d_ - a * f],
           [d_ * h - e * g, b * g - a * h, a * e - b * d_]]
    return [[adj[r][s] / det for s in R3] for r in R3]


def silent_t6_coefficients():
    """gamma_3^silent = gamma_0 M_6 in the basis {tr(E^2) E, tr(E^3) I} (E^3 removed by Cayley-Hamilton
    for traceless E): exact rational coefficients (alpha, beta), from two generic diagonal E_0."""
    out = []
    for lam in ((fmpq(1), fmpq(2)), (fmpq(3), fmpq(-5))):
        l1, l2 = lam
        E0 = [[l1, 0, 0], [0, l2, 0], [0, 0, -l1 - l2]]
        E0 = [[fmpq(v) for v in row] for row in E0]
        _, _, M = silent_series(E0, 6, fmpq(1))
        z = fmpq(0)
        E2 = _mm(E0, E0, z)
        E3 = _mm(E2, E0, z)
        t2, t3 = _tr(E2), _tr(E3)
        # M6 = alpha t2 E + beta t3 I  (diagonal): two equations from entries 0, 1 (and check entry 2)
        a11, a12, b1 = t2 * E0[0][0], t3, M[6][0][0]
        a21, a22, b2 = t2 * E0[1][1], t3, M[6][1][1]
        det = a11 * a22 - a12 * a21
        alpha = (b1 * a22 - a12 * b2) / det
        beta = (a11 * b2 - b1 * a21) / det
        assert M[6][2][2] == alpha * t2 * E0[2][2] + beta * t3
        out.append((alpha, beta))
    assert out[0] == out[1]
    return out[0]


def silent_crunch(E0_eigs, t_max=20.0):
    """Integrate the silent ODE in the eigenframe of E_0 (float, scipy, rtol 1e-12) until a scale
    factor reaches 1e-6: the silent collapse time. Exploratory (numerical)."""
    from scipy.integrate import solve_ivp
    lam = np.array(E0_eigs, dtype=float)

    def rhs(t, y):
        a, K, E = y[0:3], y[3:6], y[6:9]
        th = K.sum()
        sig = K - th / 3
        sE = sig * E
        dE = -th * E + 3 * sE - sE.sum()
        return np.concatenate([a * K, -E - K * K, dE])

    def hit(t, y):
        return y[0:3].min() - 1e-6
    hit.terminal = True
    sol = solve_ivp(rhs, (0, t_max), np.concatenate([np.ones(3), np.zeros(3), lam]), events=hit,
                    rtol=1e-12, atol=1e-14, method="DOP853")
    return float(sol.t[-1]) if sol.status == 1 else None


def ratio_radius(c):
    """Ratio estimates |c_k / c_(k+1)|^(1/2) of a series in s = t^2 (exploratory)."""
    return [abs(c[k] / c[k + 1]) ** 0.5 for k in range(2, len(c) - 1) if c[k + 1] != 0]


def frob(A):
    """Frobenius norm upper bound (float) of a 3x3 list of arb (for conformally flat gamma_0 the mixed
    matrices are symmetric in Cartesian coordinates, so this is the gamma_0-orthonormal norm)."""
    return math.sqrt(sum(float(A[i][j].mid()) ** 2 for i in R3 for j in R3))


def mixed(G0, S, zero):
    return _mm(G0, S, zero)


def remainder_chain(x0, d, one, exact_r=None, centres=None):
    """The first deviation from the silent system is a single traceless source c3 = (curl H)_3 (mixed,
    coefficient of t^3 in dE/dt). It enters E at t^4 (c3/4), K at t^5 (-c3/20), gamma at t^6 (-c3/60),
    so the predicted relation between the two independently computed deviations is
        G0 gamma_3 - M_6^silent = -(1/15) (E_true,4 - E_silent,4).
    Returns (dM6, dE4) as 3x3 lists (mixed), and E0."""
    from research.calc import bl_structure as st
    zero = 0 * one
    true, silent, E0, odd = true_and_silent(x0, d, 3, one, exact_r=exact_r, centres=centres)
    G0 = _inv3(true[0], one)
    dM6 = _add(mixed(G0, true[3], zero), _scale(mixed(G0, silent[3], zero), -one))
    Et = st.electric_series(x0, d, one, exact_r=exact_r, centres=centres, order=2)
    _, Es, _ = silent_series(E0, 4, one)
    dE4 = _add(Et[2], _scale(Es[4], -one))
    return dM6, dE4, E0, Et, Es


def eig_sym3(E0):
    """Eigenvalues (float) of the mixed tidal matrix (symmetric for conformally flat gamma_0)."""
    A = np.array([[float(E0[i][j].mid()) if hasattr(E0[i][j], "mid") else float(E0[i][j]) for j in R3]
                  for i in R3])
    return sorted(np.linalg.eigvalsh((A + A.T) / 2))


def high_order(x0, d, N=8, n_sil=40, prec=256):
    """True coefficients gamma_k (k <= N, s = t^2) at x0 (decimal strings, converted at the working
    precision) from the jets, the silent ones (to order n_sil in s, cheap), and the remainder
    Delta_k = G0 gamma_k - M_(2k)^silent. Ball arithmetic; the reported norms are midpoints of
    Frobenius norms (exploratory diagnostics), with the largest ball radius reported alongside."""
    with bm.working_precision(prec):
        x0 = tuple(arb(v) for v in x0)
        d = arb(d)
        one = arb(1)
        zero = arb(0)
        t0 = time.time()
        J, g, _ = bm.series_jets(x0, d, N, one)
        vals = [{pq: J.coef(g[pq], (0, 0, 0, k), zero) for pq in bm.PAIRS} for k in range(N + 1)]
        true = [[[bm.sym_get(v, i, j) for j in R3] for i in R3] for v in vals]
        G0 = _inv3(true[0], one)
        A = [mixed(G0, T, zero) for T in true]
        E0 = _scale(A[1], -one)
        _, _, M = silent_series(E0, 2 * n_sil, one)
        Ms = [M[2 * k] for k in range(n_sil + 1)]
        D = [_add(A[k], _scale(Ms[k], -one)) for k in range(N + 1)]
        # largest ball radius among the true coefficients (to check precision is adequate)
        rad = max(float(A[k][i][j].rad()) for k in range(N + 1) for i in R3 for j in R3)
        out = {"x0": [str(v) for v in x0], "d": str(d), "N": N, "jet_seconds": time.time() - t0,
               "max_ball_radius_true": rad,
               "E0_eigenvalues": eig_sym3(E0),
               "norm_true": [frob(a) for a in A], "norm_silent": [frob(m) for m in Ms],
               "norm_remainder": [frob(x) for x in D],
               "components_true_xx": [float(A[k][0][0].mid()) for k in range(N + 1)],
               "components_silent_xx": [float(Ms[k][0][0].mid()) for k in range(n_sil + 1)],
               "components_remainder_xx": [float(D[k][0][0].mid()) for k in range(N + 1)]}
        out["relative_remainder"] = [out["norm_remainder"][k] / out["norm_true"][k] if out["norm_true"][k] else None
                                     for k in range(N + 1)]
        out["ratio_T_true"] = ratio_radius(out["norm_true"])
        out["ratio_T_silent"] = ratio_radius(out["norm_silent"])
        out["ratio_T_remainder"] = [abs(out["norm_remainder"][k] / out["norm_remainder"][k + 1]) ** 0.5
                                    for k in range(3, N) if out["norm_remainder"][k + 1]]
        out["silent_crunch_t"] = silent_crunch(out["E0_eigenvalues"])
        return out


def curl_h3_law(x0, d, one, exact_r=None, centres=None):
    """The source of the first post-silent deviation, (curl H)_3 = 4 (E_true,4 - E_silent,4) (lowered with
    gamma_0), against the covariant law implied by the silent t^6 coefficients (-13/180, -41/540) and
    the session-13 t^6 law (notes/17 §5):
        (curl H)_3 = -(1/3)|Ric|^2 Ric - nabla_i Ric . nabla_j Ric + (1/3)|nabla Ric|^2 gamma_0
                     + (1/3) Ric^kl nabla_k nabla_l Ric_ij.
    Returns (lhs, rhs) as dicts over PAIRS."""
    from research.calc import bl_structure as st
    zero = 0 * one
    dM6, dE4, E0, _, _ = remainder_chain(x0, d, one, exact_r=exact_r, centres=centres)
    cv = st.curvature_values(x0, d, one, exact_r=exact_r, centres=centres)
    g0 = cv["g0"]
    low = _mm(g0, dE4, zero)
    lhs = {pq: 4 * low[pq[0]][pq[1]] for pq in bm.PAIRS}
    B = st.t6_basis(cv, zero)
    nm = st.T6_NAMES
    coef = {"|Ric|^2 Ric": -one / 3, "d_iRic.d_jRic": -one, "g |dRic|^2": one / 3, "Ric^kl d_k d_l Ric_ij": one / 3}
    rhs = {pq: sum((c * B[nm.index(k)][pq] for k, c in coef.items()), zero) for pq in bm.PAIRS}
    return lhs, rhs


POINTS = {"P1": (("0.17", "0.11", "0.5"), "0.7662"), "P2": (("0.26", "-0.09", "0.12"), "0.767"),
          "P3": (("0.6", "0.3", "0.9"), "0.7662"), "P4": (("0.05", "0.02", "0.6"), "0.7662"),
          "P5": (("0.33", "0.21", "0.05"), "0.7662"), "P6": (("1.1", "-0.6", "1.4"), "0.7662")}


def main(out=os.path.join(bm.DATA, "bl_silent.json"), N=8, names=("P1", "P2", "P3", "P4", "P5", "P6")):
    """Exact checks (rational points), the covariant (curl H)_3 law (exact and held-out balls), and the
    high-order comparison at generic points (exploratory)."""
    from research.calc import bl_structure as st
    res = {}
    t0 = time.time()
    a, b = silent_t6_coefficients()
    res["silent_t6_coefficients"] = {"tr(E^2)E": str(a), "tr(E^3)I": str(b),
                                     "fitted_true_t6 (notes/17)": {"|Ric|^2 Ric": "-1/15", "tr(Ric^3) g": "-41/540"}}
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    za = fmpq(31, 50)
    exact_pts = {"plane": (x0, d, dict(exact_r=rr)),
                 "axis": ((fmpq(0), fmpq(0), za), d, dict(exact_r=[za - d / 2, za + d / 2])),
                 "schwarzschild": ((fmpq(3, 5), fmpq(0), fmpq(4, 5)), None,
                                   dict(centres=[(fmpq(1, 2), (fmpq(0),) * 3, fmpq(1))]))}
    ex = {}
    for name, (x, dd, kw) in exact_pts.items():
        true, sil, E0, odd = true_and_silent(x, dd, 3, fmpq(1), **kw)
        agree = [all(true[k][i][j] == sil[k][i][j] for i in R3 for j in R3) for k in range(4)]
        dM6, dE4, _, _, _ = remainder_chain(x, dd, fmpq(1), **kw)
        chain = all(dM6[i][j] == -fmpq(1, 15) * dE4[i][j] for i in R3 for j in R3)
        lhs, rhs = curl_h3_law(x, dd, fmpq(1), **kw)
        ex[name] = {"true_equals_silent_t^0..t^6": agree, "odd_silent_coefficients_vanish":
                    all(v == 0 for M in odd for row in M for v in row),
                    "chain_dM6_eq_-dE4/15": chain, "curlH3_law_exact": all(lhs[p] == rhs[p] for p in lhs),
                    "curlH3_nonzero": any(lhs[p] != 0 for p in lhs)}
    res["exact_points"] = ex
    print("exact", ex, f"{time.time() - t0:.0f}s", flush=True)
    held = {}
    with bm.working_precision(256):
        for xs, ds in ((("0.13", "0.07", "0.46"), "0.7662"), (("0.31", "0.2", "-0.7"), "0.7"),
                       (("1.3", "-0.4", "0.25"), "1.5")):
            lhs, rhs = curl_h3_law(tuple(arb(v) for v in xs), arb(ds), arb(1))
            held[f"{xs} d={ds}"] = {"max_abs_residual": max(float(abs(lhs[p] - rhs[p]).upper()) for p in lhs),
                                    "max_abs_curlH3": max(float(abs(lhs[p]).mid()) for p in lhs)}
    res["curlH3_law_held_out"] = held
    print("held-out", held, f"{time.time() - t0:.0f}s", flush=True)
    hi = {}
    for name in names:
        xs, ds = POINTS[name]
        hi[name] = high_order(xs, ds, N=N)
        print(name, f"{time.time() - t0:.0f}s", flush=True)
    res["high_order"] = hi
    res["seconds"] = time.time() - t0
    res["replay"] = "python -m research.calc.bl_silent"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    main()
