"""Axisymmetric structure sprint 01 (research/axisym_structure_sprint_01.md).

Reduction (standard: Geroch 1971; Brink-Zimmerman-Hinderer 2013, eqs. 2.3-2.10, twist-free case):
g = g3 + Lambda^2 dphi^2 with Lambda = |d_phi| and g3 the orbit-space metric. Vacuum is equivalent to
    Box_{g3} Lambda = 0,      Ric(g3) = Lambda^-1 Hess_{g3} Lambda.

CONTROL 1  Schwarzschild (isotropic, static): both reduced equations hold identically (SymPy).
CONTROL 2  Teukolsky's even-parity m = 0 linearised quadrupole wave on Minkowski (PRD 26, 745 (1982)):
           linearised Ricci = 0 for the given F, linearised Riemann != 0 and O(1/r) (radiative), no dphi
           cross terms (twist-free), delta Lambda != 0 (the wave lives in the reduced scalar).
A3 TEST    Hypothesis KFA (frozen, see the note): the first nonlocal correction S = (curl H)_3 of the
           BL geodesic development has the unit gradient n of the initial axial norm Lambda_0 = psi^2 rho
           as an eigenvector, i.e. S(n, m) = 0 with m the third leg of {n, m, e_phi}.

Run: python -m research.calc.axisym_sprint01
"""

import sympy as sp
from flint import arb

from research.calc import bl_metric as bm
from research.calc import bl_silent as sl

R3 = range(3)


# ---------------------------------------------------------------------------------------------
# generic symbolic helpers
# ---------------------------------------------------------------------------------------------
def christoffel(g, X):
    gi = g.inv()
    n = len(X)
    return [[[sp.simplify(sum(gi[a, e] * (sp.diff(g[e, b], X[c]) + sp.diff(g[e, c], X[b]) - sp.diff(g[b, c], X[e]))
                              for e in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]


def ricci(g, X, G=None):
    n = len(X)
    G = G or christoffel(g, X)
    R = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            R[b, d] = sp.simplify(sum(sp.diff(G[a][b][d], X[a]) - sp.diff(G[a][a][b], X[d]) for a in range(n))
                                  + sum(G[a][a][e] * G[e][b][d] - G[a][d][e] * G[e][a][b] for a in range(n) for e in range(n)))
    return R


def control_schwarzschild():
    """Reduced equations for Schwarzschild in isotropic coordinates, meridian polar (r, th)."""
    t, r, th, M = sp.symbols("t r theta M", positive=True)
    X = [t, r, th]
    psi = 1 + M / (2 * r)
    N = (1 - M / (2 * r)) / (1 + M / (2 * r))
    g3 = sp.diag(-N ** 2, psi ** 4, psi ** 4 * r ** 2)
    Lam = psi ** 2 * r * sp.sin(th)
    G = christoffel(g3, X)
    hess = sp.Matrix(3, 3, lambda a, b: sp.diff(Lam, X[a], X[b]) - sum(G[c][a][b] * sp.diff(Lam, X[c]) for c in range(3)))
    gi = g3.inv()
    box = sp.simplify(sum(gi[a, b] * hess[a, b] for a in range(3) for b in range(3)))
    Ric = ricci(g3, X, G)
    res = (Ric - hess / Lam).applyfunc(sp.simplify)
    return {"box_Lambda": box, "Ric_minus_Hess_over_Lambda_zero": res == sp.zeros(3, 3)}


def control_teukolsky(profile_power=7):
    """Even-parity m = 0 Teukolsky wave with F(x) = x^p (x = t - r), linearised about Minkowski."""
    t, r, th, ph, eps = sp.symbols("t r theta phi epsilon")
    X = [t, r, th, ph]
    x = sp.Symbol("x")
    Fx = x ** profile_power
    Fd = [sp.diff(Fx, x, k).subs(x, t - r) for k in range(5)]
    A = 3 * (Fd[2] / r ** 3 + 3 * Fd[1] / r ** 4 + 3 * Fd[0] / r ** 5)
    B = -(Fd[3] / r ** 2 + 3 * Fd[2] / r ** 3 + 6 * Fd[1] / r ** 4 + 6 * Fd[0] / r ** 5)
    C = sp.Rational(1, 4) * (Fd[4] / r + 2 * Fd[3] / r ** 2 + 9 * Fd[2] / r ** 3 + 21 * Fd[1] / r ** 4 + 21 * Fd[0] / r ** 5)
    s, c = sp.sin(th), sp.cos(th)
    frr, frth = 2 - 3 * s ** 2, -3 * s * c
    f1thth, f2thth = 3 * s ** 2, -1
    f1phph, f2phph = -3 * s ** 2, 3 * s ** 2 - 1
    h = sp.zeros(4, 4)
    h[1, 1] = A * frr
    h[1, 2] = h[2, 1] = B * frth * r
    h[2, 2] = (C * f1thth + A * f2thth) * r ** 2
    h[3, 3] = (C * f1phph + A * f2phph) * r ** 2 * s ** 2
    eta = sp.diag(-1, 1, r ** 2, r ** 2 * s ** 2)
    g = eta + eps * h
    # linearised Christoffels and Ricci (first order in eps)
    G0 = christoffel(eta, X)
    gi0 = eta.inv()
    G1 = [[[sp.simplify(sum(gi0[a, e] * (sp.diff(h[e, b], X[cc]) + sp.diff(h[e, cc], X[b]) - sp.diff(h[b, cc], X[e]))
                            for e in range(4)) / 2 - sum(gi0[a, e] * h[e, f] * G0[f][b][cc] for e in range(4) for f in range(4)))
            for cc in range(4)] for b in range(4)] for a in range(4)]
    R1 = sp.zeros(4, 4)
    for b in range(4):
        for d in range(4):
            R1[b, d] = sp.simplify(sum(sp.diff(G1[a][b][d], X[a]) - sp.diff(G1[a][a][b], X[d]) for a in range(4))
                                   + sum(G0[a][a][e] * G1[e][b][d] + G1[a][a][e] * G0[e][b][d]
                                         - G0[a][d][e] * G1[e][a][b] - G1[a][d][e] * G0[e][a][b]
                                         for a in range(4) for e in range(4)))
    # a transverse linearised Riemann component R1^t_{theta t theta} (orthonormal: divide by r^2); the
    # longitudinal R^t_{rtr} falls off faster than 1/r and does not show the radiation
    a_, b_, c_, d_ = 0, 2, 0, 2
    Riem = sp.simplify(sp.diff(G1[a_][d_][b_], X[c_]) - sp.diff(G1[a_][c_][b_], X[d_])
                       + sum(G0[a_][c_][e] * G1[e][d_][b_] + G1[a_][c_][e] * G0[e][d_][b_]
                             - G0[a_][d_][e] * G1[e][c_][b_] - G1[a_][d_][e] * G0[e][c_][b_] for e in range(4)))
    twist_free = all(g[3, k] == 0 for k in range(3))
    dLam_over_Lam = sp.simplify(h[3, 3] / (2 * r ** 2 * s ** 2))
    # large-r behaviour of the Riemann component at fixed retarded time u = t - r
    u = sp.Symbol("u")
    Riem_u = sp.simplify(Riem.subs(t, u + r))
    lead = sp.simplify(sp.limit(Riem_u / r ** 2 * r, r, sp.oo))
    return {"linearised_Ricci_zero": R1 == sp.zeros(4, 4), "Riemann_t_rtr": Riem, "r*Riemann_at_fixed_u_limit": lead,
            "twist_free": twist_free, "deltaLambda_over_Lambda": dLam_over_Lam}


# ---------------------------------------------------------------------------------------------
# A3: Killing-frame alignment of the first nonlocal source
# ---------------------------------------------------------------------------------------------
def kfa_test(xs, ds, prec=256):
    """S = (curl H)_3 (lower, Cartesian) from the exact remainder chain S = 4 gamma_0 (E_4 - E_4^sil); the
    frame {n, m, e_phi} from Lambda_0 = psi^2 rho; returns S(n, m) and its parts, gamma_0-orthonormal."""
    with bm.working_precision(prec):
        x = tuple(arb(v) for v in xs)
        d = arb(ds)
        one = arb(1)
        dM6, dE4, E0, Et, Es = sl.remainder_chain(x, d, one)
        lhs, rhs = sl.curl_h3_law(x, d, one)                       # S two ways (chain vs proved law)
        S = [[bm.sym_get(lhs, i, j) for j in R3] for i in R3]
        chain_vs_law = max(float(abs(lhs[p] - rhs[p]).upper()) for p in lhs)
        # psi and grad psi (Euclidean), punctures at z = +-d/2, coefficient 1/4
        psi, dpsi = one, [arb(0)] * 3
        for zc in (d / 2, -d / 2):
            v = [x[0], x[1], x[2] - zc]
            r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
            psi = psi + one / (4 * r)
            dpsi = [dpsi[k] - v[k] / (4 * r * r * r) for k in R3]
        rho = (x[0] * x[0] + x[1] * x[1]).sqrt()
        drho = [x[0] / rho, x[1] / rho, arb(0)]
        dLam = [2 * psi * rho * dpsi[k] + psi * psi * drho[k] for k in R3]   # grad Lambda_0 (Cartesian)
        nrm = lambda v: (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()  # noqa: E731
        n = [c / nrm(dLam) for c in dLam]
        ephi = [-x[1] / rho, x[0] / rho, arb(0)]
        m = [ephi[1] * n[2] - ephi[2] * n[1], ephi[2] * n[0] - ephi[0] * n[2], ephi[0] * n[1] - ephi[1] * n[0]]
        p4 = psi * psi * psi * psi
        comp = lambda a, b, T: sum((a[i] * b[j] * T[i][j] for i in R3 for j in R3), arb(0)) / p4  # noqa: E731
        Snm, Snphi, Smphi = comp(n, m, S), comp(n, ephi, S), comp(m, ephi, S)
        normS = (sum((S[i][j] * S[i][j] for i in R3 for j in R3), arb(0))).sqrt() / p4
        # which law terms break the alignment: Ric(n, m) and the four terms of the law separately
        from research.calc import bl_structure as st
        cv = st.curvature_values(x, d, one)
        B = st.t6_basis(cv, arb(0))
        nm = st.T6_NAMES
        terms = {"-(1/3)|Ric|^2 Ric": (-one / 3, "|Ric|^2 Ric"), "-d_iRic.d_jRic": (-one, "d_iRic.d_jRic"),
                 "+(1/3)|dRic|^2 g": (one / 3, "g |dRic|^2"), "+(1/3)Ric^kl d_k d_l Ric_ij": (one / 3, "Ric^kl d_k d_l Ric_ij")}
        parts = {}
        for label, (cf, key) in terms.items():
            T = [[cf * bm.sym_get(B[nm.index(key)], i, j) for j in R3] for i in R3]
            parts[label] = comp(n, m, T)
        Ric = [[cv["R"][i][j] for j in R3] for i in R3]
        return {"point": xs, "d": ds, "S(n,m)": Snm, "S(n,e_phi)": Snphi, "S(m,e_phi)": Smphi, "|S|": normS,
                "chain_vs_law_max_residual": chain_vs_law, "Ric(n,m)": comp(n, m, Ric), "parts_(n,m)": parts,
                "excludes_zero": not Snm.contains(0)}


def main():
    out = {}
    out["control_schwarzschild"] = control_schwarzschild()
    print("CONTROL 1", out["control_schwarzschild"], flush=True)
    out["control_teukolsky"] = control_teukolsky()
    print("CONTROL 2", {k: (v if k != "Riemann_t_rtr" else "computed") for k, v in out["control_teukolsky"].items()}, flush=True)
    for xs, ds in ((("0.17", "0.11", "0.5"), "0.7662"), (("0.31", "-0.2", "0.27"), "0.7662"),
                   (("0.4", "0.2", "1.1"), "1.5"), (("0.26", "-0.09", "0"), "0.767")):
        r = kfa_test(xs, ds)
        out[f"{xs} d={ds}"] = r
        mr = lambda b: f"{float(b.mid()):.10g} (rad {float(b.rad()):.1e})"  # noqa: E731
        print("A3", xs, ds, "S(n,m) =", mr(r["S(n,m)"]), "|S| =", mr(r["|S|"]), "excl0:", r["excludes_zero"],
              "S(n,ephi) =", mr(r["S(n,e_phi)"]), "chain-law", f"{r['chain_vs_law_max_residual']:.1e}",
              "Ric(n,m) =", mr(r["Ric(n,m)"]), {k: mr(v) for k, v in r["parts_(n,m)"].items()}, flush=True)
    return out


if __name__ == "__main__":
    main()
