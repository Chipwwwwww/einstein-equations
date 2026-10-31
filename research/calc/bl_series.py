"""Route E (research/route_map.md): the vacuum development of Brill-Lindquist data as a time series.

Gauge: geodesic slicing, lapse N = 1, shift 0, so the spatial coordinates (rho, z, phi) are comoving
with the observers at rest on the initial slice. With the repository's sign K_ij = (1/2N) d_t gamma_ij
(notes/10, tests/test_validation.py::test_adm_evolution_equation_sign), the vacuum ADM equations are
    d_t gamma_ij = 2 K_ij,
    d_t K_ij     = -R_ij(gamma) + 2 K_ik K^k_j - K K_ij,
    R(gamma) + K^2 - K_ij K^ij = 0,       D_j (K^j_i - delta^j_i K) = 0.
Time-symmetric data (K = 0 at t = 0) give gamma even and K odd in t. With s = t^2,
gamma = sum_k g_k s^k and K = t kappa(s), kappa = sum_k k_k s^k, the evolution equations become the
explicit recursion
    g_{k+1} = k_k / (k + 1),
    (2k + 1) k_k = [-Ric(gamma)]_k + [2 kappa gamma^-1 kappa - tr(gamma^-1 kappa) kappa]_{k-1},
starting from g_0 = psi^4 (d rho^2 + dz^2 + rho^2 d phi^2). So g_1 = k_0 = -Ric(g_0), and every
coefficient is an explicit differential expression in psi.

Why the series is an exact construction (category B, local). Written second order in time,
d_t^2 gamma = F(gamma, d_x gamma, d_x^2 gamma, d_t gamma) is in Kovalevskaya form, and the data are
real-analytic away from the punctures. By the Cauchy-Kovalevskaya theorem the series converges for
|t| < T(x) near each point and gives an analytic solution of the evolution equations. The constraint
quantities satisfy a linear homogeneous system in Kovalevskaya form (contracted Bianchi identities) and
vanish at t = 0, so they vanish identically. Hence the sum solves the full vacuum Einstein equations
on its domain. T(x) is not computed here: the ratio estimates below are numerical evidence only.

Implementation: everything is a truncated Taylor polynomial (a "jet") in (u, w, s) around a point
(rho0, z0) off the axis, u = rho - rho0, w = z - z0, with weight (deg_u + deg_w) + 2 deg_s <= W. Each
spatial derivative lowers the valid weight by one. The recursion then determines g_k at the point
for 2k <= W (the s^k coefficient needs 2k spatial derivatives). Floating point (numpy float64,
direct convolution): an exploration tool, validated against Schwarzschild free fall
(tests/test_bl_series.py), not a certificate.

Run: python -m research.calc.bl_series
"""

import math

import numpy as np
from scipy.signal import convolve


class JetSpace:
    """Truncated polynomials in (u, w, s) with weight i + j + 2k <= W."""

    def __init__(self, W):
        self.W = W
        self.N = W // 2
        i, j, k = np.indices((W + 1, W + 1, self.N + 1))
        self.mask = (i + j + 2 * k) <= W

    def zero(self):
        return np.zeros((self.W + 1, self.W + 1, self.N + 1))

    def const(self, c):
        a = self.zero()
        a[0, 0, 0] = c
        return a

    def var(self, which):
        a = self.zero()
        a[(1, 0, 0) if which == "u" else (0, 1, 0)] = 1.0
        return a

    def mul(self, a, b):
        c = convolve(a, b, method="direct")[: self.W + 1, : self.W + 1, : self.N + 1]
        return np.where(self.mask, c, 0.0)

    def d(self, a, axis):
        """d/du (axis 0) or d/dw (axis 1)."""
        out = self.zero()
        n = self.W
        idx = np.arange(1, n + 1).reshape((-1, 1, 1) if axis == 0 else (1, -1, 1))
        if axis == 0:
            out[:-1] = a[1:] * idx
        else:
            out[:, :-1] = a[:, 1:] * idx
        return np.where(self.mask, out, 0.0)

    def inv(self, a):
        a0 = a[0, 0, 0]
        dl = a.copy()
        dl[0, 0, 0] = 0.0
        x = -dl / a0
        term, acc = self.const(1.0), self.const(1.0)
        for _ in range(self.W):
            term = self.mul(term, x)
            acc = acc + term
        return acc / a0

    def powr(self, a, p):
        """a^p for a with a[0,0,0] > 0 (binomial series)."""
        a0 = a[0, 0, 0]
        dl = a.copy()
        dl[0, 0, 0] = 0.0
        x = dl / a0
        term, acc, coef = self.const(1.0), self.const(1.0), 1.0
        for n in range(1, self.W + 1):
            coef *= (p - n + 1) / n
            term = self.mul(term, x)
            acc = acc + coef * term
        return acc * a0**p


def psi_jet(J, rho0, z0, centres):
    """psi = 1 + sum m / (2 |x - p|) as a jet at (rho0, z0); centres = [(m, z_c), ...] on the axis."""
    u, w = J.var("u"), J.var("w")
    psi = J.const(1.0)
    for m, zc in centres:
        X, Y = u + J.const(rho0), w + J.const(z0 - zc)  # constants enter the s^0 u^0 w^0 slot only
        q = J.mul(X, X) + J.mul(Y, Y)
        psi = psi + (m / 2) * J.powr(q, -0.5)
    return psi


def initial_metric(J, rho0, z0, centres):
    """g_0 in (rho, z, phi): psi^4 (1, 1, rho^2), as components (rr, rz, zz, pp)."""
    psi = psi_jet(J, rho0, z0, centres)
    p4 = J.mul(J.mul(psi, psi), J.mul(psi, psi))
    rho = J.var("u") + J.const(rho0)
    # separate arrays: the recursion fills the s-slices in place (rr and zz must not alias)
    return {"rr": p4.copy(), "rz": J.zero(), "zz": p4.copy(), "pp": J.mul(p4, J.mul(rho, rho))}


def _full(g):
    """3x3 nested list of jets from components (rr, rz, zz, pp); phi-mixed components vanish."""
    Z = np.zeros_like(g["rr"])
    return [[g["rr"], g["rz"], Z], [g["rz"], g["zz"], Z], [Z, Z, g["pp"]]]


def inverse(J, g):
    det = J.mul(g["rr"], g["zz"]) - J.mul(g["rz"], g["rz"])
    idet = J.inv(det)
    Z = J.zero()
    gi = {"rr": J.mul(g["zz"], idet), "rz": -J.mul(g["rz"], idet), "zz": J.mul(g["rr"], idet),
          "pp": J.inv(g["pp"])}
    return [[gi["rr"], gi["rz"], Z], [gi["rz"], gi["zz"], Z], [Z, Z, gi["pp"]]]


def ricci(J, g):
    """R_bd = d_a G^a_bd - d_d G^a_ab + G^a_ae G^e_bd - G^a_de G^e_ab (repository convention), for a
    metric independent of phi (index 2). Returns components (rr, rz, zz, pp)."""
    G = _full(g)
    Gi = inverse(J, g)

    def dd(a, k):
        return J.d(a, k) if k < 2 else np.zeros_like(a)
    dG = [[[dd(G[i][j], k) for k in range(3)] for j in range(3)] for i in range(3)]
    Gam = [[[None] * 3 for _ in range(3)] for _ in range(3)]
    for a in range(3):
        for b in range(3):
            for c in range(b, 3):
                acc = J.zero()
                for e in range(3):
                    if not Gi[a][e].any():
                        continue
                    t = dG[e][c][b] + dG[e][b][c] - dG[b][c][e]
                    if t.any():
                        acc = acc + J.mul(Gi[a][e], t)
                Gam[a][b][c] = Gam[a][c][b] = 0.5 * acc
    out = {}
    for name, (b, dI) in {"rr": (0, 0), "rz": (0, 1), "zz": (1, 1), "pp": (2, 2)}.items():
        acc = J.zero()
        for a in range(3):
            acc = acc + dd(Gam[a][b][dI], a) - dd(Gam[a][a][b], dI)
            for e in range(3):
                acc = acc + J.mul(Gam[a][a][e], Gam[e][b][dI]) - J.mul(Gam[a][dI][e], Gam[e][a][b])
        out[name] = acc
    return out


def _mat(J, c):
    return _full(c)


def _matmul(J, A, B):
    return [[sum((J.mul(A[i][k], B[k][j]) for k in range(3)), J.zero()) for j in range(3)]
            for i in range(3)]


def series(rho0, z0, centres, W):
    """The coefficients g_k(rho0, z0), k_k(rho0, z0) of the geodesic-slicing development, k <= W/2,
    with the constraint residuals at each order. Returns (g, kappa, ham, mom): dicts of component ->
    list over k of values at the point, and lists of constraint residual values."""
    J = JetSpace(W)
    g = initial_metric(J, rho0, z0, centres)
    kap = {c: J.zero() for c in g}
    shift = lambda a: np.concatenate([np.zeros_like(a[:, :, :1]), a[:, :, :-1]], axis=2)  # noqa: E731
    for k in range(J.N):
        Ric = ricci(J, g)
        Gi = inverse(J, g)
        K = _mat(J, kap)
        KK = _matmul(J, _matmul(J, K, Gi), K)          # kappa gamma^-1 kappa
        trK = sum((J.mul(Gi[i][j], K[j][i]) for i in range(3) for j in range(3)), J.zero())
        quad = {c: 2 * KK[i][j] - J.mul(trK, K[i][j])
                for c, (i, j) in {"rr": (0, 0), "rz": (0, 1), "zz": (1, 1), "pp": (2, 2)}.items()}
        for c in g:
            rhs = -Ric[c] + shift(quad[c])               # s * quad shifts the s-index by one
            kap[c][:, :, k] = rhs[:, :, k] / (2 * k + 1)
            if k + 1 <= J.N:
                g[c][:, :, k + 1] = kap[c][:, :, k] / (k + 1)
        for c in g:
            g[c] = np.where(J.mask, g[c], 0.0)
            kap[c] = np.where(J.mask, kap[c], 0.0)
    # constraints, at the point, order by order in s: valid for orders k with 2k + 2 <= W
    Ric = ricci(J, g)
    Gi = inverse(J, g)
    K = _mat(J, kap)
    R = sum((J.mul(Gi[i][j], Ric[c]) * (2 if i != j else 1)
             for c, (i, j) in {"rr": (0, 0), "rz": (0, 1), "zz": (1, 1), "pp": (2, 2)}.items()),
            J.zero())
    trK = sum((J.mul(Gi[i][j], K[j][i]) for i in range(3) for j in range(3)), J.zero())
    Kup = _matmul(J, _matmul(J, Gi, K), Gi)
    KK = sum((J.mul(Kup[i][j], K[i][j]) for i in range(3) for j in range(3)), J.zero())
    ham = R + shift(J.mul(trK, trK) - KK)
    vals = {c: [g[c][0, 0, k] for k in range(J.N + 1)] for c in g}
    kv = {c: [kap[c][0, 0, k] for k in range(J.N)] for c in g}
    hres = [ham[0, 0, k] for k in range(J.N)]
    return vals, kv, hres


def schwarzschild_ratio(R0, M, N):
    """Coefficients of (r(t)/R0)^2 in s = t^2 for radial free fall from rest at areal radius R0:
    r'' = -M / r^2, r(0) = R0, r'(0) = 0 (exact recursion for the Taylor coefficients)."""
    # r = sum a_n t^n; r^2 r'' = -M
    n_max = 2 * N + 2
    a = [0.0] * (n_max + 1)
    a[0] = R0

    def coef_r2(n):  # coefficient of t^n in r^2
        return sum(a[i] * a[n - i] for i in range(n + 1))
    for n in range(0, n_max - 1):
        # coefficient of t^n in r^2 r'' : sum_{i} [r^2]_i (n - i + 2)(n - i + 1) a_{n - i + 2}
        known = sum(coef_r2(i) * (n - i + 2) * (n - i + 1) * a[n - i + 2] for i in range(1, n + 1))
        rhs = (-M if n == 0 else 0.0) - known
        a[n + 2] = rhs / (coef_r2(0) * (n + 2) * (n + 1))
    r2 = [coef_r2(2 * k) / R0**2 for k in range(N + 1)]
    return r2


def main():
    W = 12
    # validation: Schwarzschild (one centre, M = 1) at isotropic radius 1 on the plane z = 0
    vals, _, ham = series(1.0, 0.0, [(1.0, 0.0)], W)
    R0 = 1.0 * (1 + 1 / 2) ** 2
    exact = schwarzschild_ratio(R0, 1.0, W // 2)
    got = [v / vals["pp"][0] for v in vals["pp"]]
    print("Schwarzschild, isotropic r = 1 (areal R0 = 2.25): (r(t)/R0)^2 in s = t^2")
    for k, (a, b) in enumerate(zip(got, exact)):
        print(f"  s^{k}: series {a:+.12e}  exact {b:+.12e}")
    print("  Hamiltonian constraint residuals:", [f"{h:.1e}" for h in ham])
    # BL near-critical: a point near the critical surface's equator (d = d_c, rho ~ 0.26)
    d = 0.766197428240
    vals, kv, ham = series(0.26, 0.0, [(0.5, d / 2), (0.5, -d / 2)], W)
    print("BL d = d_c, point (rho, z) = (0.26, 0): g_pp coefficients", [f"{v:+.4e}" for v in vals["pp"]])
    print("  Hamiltonian constraint residuals:", [f"{h:.1e}" for h in ham])
    ratios = [abs(vals["pp"][k] / vals["pp"][k + 1]) for k in range(1, len(vals["pp"]) - 1)
              if vals["pp"][k + 1] != 0]
    print("  |g_k / g_{k+1}| (estimates of T^2):", [f"{r:.3f}" for r in ratios])
    print("  Schwarzschild check: T = pi sqrt(R0^3 / 8M), at the throat R0 = 2M: T = pi M")


if __name__ == "__main__":
    main()


def magnetic_weyl(rho0, z0, centres, W):
    """B_{rho phi} and B_{z phi} of the development at the point, by order in t (valid for s-orders
    k <= (W - 3) / 2: one more spatial derivative than the recursion uses). B_ij =
    eps_i^{kl} D_k K_lj (symmetrised), with K = t kappa(s). Returns lists over s-orders of the
    coefficients of t * s^k. In the axisymmetric twist-free case these are the only components
    (B is odd under phi -> -phi). At order t (s^0) B = -curl Ric(gamma_0) = -Cotton(gamma_0), which
    vanishes for conformally flat data: B starts at t^3."""
    J = JetSpace(W)
    g = initial_metric(J, rho0, z0, centres)
    kap = {c: J.zero() for c in g}
    shift = lambda a: np.concatenate([np.zeros_like(a[:, :, :1]), a[:, :, :-1]], axis=2)  # noqa: E731
    for k in range(J.N):
        Ric = ricci(J, g)
        Gi = inverse(J, g)
        K = _mat(J, kap)
        KK = _matmul(J, _matmul(J, K, Gi), K)
        trK = sum((J.mul(Gi[i][j], K[j][i]) for i in range(3) for j in range(3)), J.zero())
        quad = {c: 2 * KK[i][j] - J.mul(trK, K[i][j])
                for c, (i, j) in {"rr": (0, 0), "rz": (0, 1), "zz": (1, 1), "pp": (2, 2)}.items()}
        for c in g:
            rhs = -Ric[c] + shift(quad[c])
            kap[c][:, :, k] = rhs[:, :, k] / (2 * k + 1)
            if k + 1 <= J.N:
                g[c][:, :, k + 1] = kap[c][:, :, k] / (k + 1)
        for c in g:
            g[c] = np.where(J.mask, g[c], 0.0)
            kap[c] = np.where(J.mask, kap[c], 0.0)
    G, Gi, K = _full(g), inverse(J, g), _mat(J, kap)

    def dd(a, k):
        return J.d(a, k) if k < 2 else np.zeros_like(a)
    dG = [[[dd(G[i][j], k) for k in range(3)] for j in range(3)] for i in range(3)]
    Gam = [[[0.5 * sum((J.mul(Gi[a][e], dG[e][c][b] + dG[e][b][c] - dG[b][c][e]) for e in range(3)),
                       J.zero()) for c in range(3)] for b in range(3)] for a in range(3)]
    # D_k K_lj
    DK = [[[dd(K[l][j], k) - sum((J.mul(Gam[m][k][l], K[m][j]) + J.mul(Gam[m][k][j], K[l][m])
                                  for m in range(3)), J.zero())
            for j in range(3)] for l in range(3)] for k in range(3)]
    det = J.mul(J.mul(G[0][0], G[1][1]) - J.mul(G[0][1], G[0][1]), G[2][2])
    isq = J.powr(det, -0.5)
    perm = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}

    def B(i, j):  # eps_i^{kl} D_k K_lj = g_im eps^{mkl} D_k K_lj
        acc = J.zero()
        for (m, k, l), sgn in perm.items():
            if G[i][m].any():
                acc = acc + sgn * J.mul(G[i][m], J.mul(isq, DK[k][l][j]))
        return acc
    out = {}
    for name, (i, j) in {"rp": (0, 2), "zp": (1, 2)}.items():
        sym = 0.5 * (B(i, j) + B(j, i))
        out[name] = [sym[0, 0, k] for k in range((W - 3) // 2 + 1)]
    return out
