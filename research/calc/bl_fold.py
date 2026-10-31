"""H4 (notes/09): how fast a common apparent horizon forms in the vacuum evolution of
near-critical Brill-Lindquist data.

Data. Time-symmetric Brill-Lindquist data, two punctures of bare mass m = 1/2 (ADM mass
M = 1) at z = +-d/2:  gamma = psi^4 delta, psi = 1 + m/(2 r1) + m/(2 r2), K = 0. Vacuum,
Lambda = 0, asymptotically flat: the primary target of notes/07.

Evolution to first order. With lapse N and zero shift, the vacuum development has
gamma(t) = gamma + O(t^2) and K(t) = t Kdot + O(t^2) with (tests/test_validation.py)

    Kdot_ij = D_i D_j N - N R_ij.

So on the slice Sigma_t a closed surface S has theta_+ = H_S + t p + O(t^2), with
p = tr_S Kdot. A common MOTS on Sigma_t solves H_S + t p = 0 to this order.

The fold. As d grows past d_c the two common minimal surfaces (outer, inner) of the t = 0
slice merge and disappear: a degenerate minimal surface S_c with a positive Jacobi field
phi_0 (L phi_0 = 0, L = -Lap_S - Ric(nu, nu) - |A|^2). Lyapunov-Schmidt reduction of the
MOTS equation near S_c in the parameters (d, t) gives, to leading order,

    t_*(d) = (d - d_c) <phi_0, dH/dd> / Int N phi_0 |A|^2 dA,          (*)

where <f, g> = Int_S f g dA. The denominator is minus the AMMS coefficient
<phi_0, p> = -Int N phi_0 |A|^2 dA (vacuum, time symmetric; Andersson-Mars-Metzger-Simon
2009, Prop. 5.1, with W = |chi+|^2 = |A|^2): positive for every lapse N > 0, so the
common MOTS always appears to the FUTURE when d slightly exceeds d_c.

This module computes everything in (*) at the fold, checks the identity
<phi_0, p> = -Int N phi_0 |A|^2 numerically, and compares (*) with t_*(d) found directly
by solving H_S + t p = 0 near the fold. Two lapses: N = 1 and N = psi^-2.

Run: python -m research.calc.bl_fold   (a few minutes)
"""

import math
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, simpson
from scipy.optimize import brentq, minimize_scalar

MASS = 0.5  # bare mass of each puncture; ADM mass M = 1


@lru_cache(maxsize=None)
def kdot_functions(lapse):
    """Numerical functions (rho, z, d) -> Kdot_ij at the point (x, y, z) = (rho, 0, z),
    for lapse 'one' (N = 1) or 'psi-2' (N = psi^-2). Built from the conformally flat
    Ricci formula of tests/test_validation.py with w = 2 log psi:
        R_ij = -(w_ij - w_i w_j) - (Lap w + |dw|^2) delta_ij,
        D_i D_j N = N_ij - w_i N_j - w_j N_i + delta_ij w_k N_k."""
    X, Y, Z, D = sp.symbols("x y z d", real=True)
    m = sp.Rational(1, 2)
    r1 = sp.sqrt(X**2 + Y**2 + (Z - D / 2) ** 2)
    r2 = sp.sqrt(X**2 + Y**2 + (Z + D / 2) ** 2)
    psi = 1 + m / (2 * r1) + m / (2 * r2)
    w = 2 * sp.log(psi)
    co = (X, Y, Z)
    dw = [sp.diff(w, c) for c in co]
    lap = sum(sp.diff(w, c, 2) for c in co)
    grad2 = sum(v**2 for v in dw)
    N = sp.Integer(1) if lapse == "one" else psi**-2
    dN = [sp.diff(N, c) for c in co]
    wN = sum(a * b for a, b in zip(dw, dN))
    comps = {}
    for (i, j) in ((0, 0), (1, 1), (2, 2), (0, 2)):
        Ric = -(sp.diff(w, co[i], co[j]) - dw[i] * dw[j]) - (lap + grad2) * (1 if i == j else 0)
        DDN = sp.diff(N, co[i], co[j]) - dw[i] * dN[j] - dw[j] * dN[i] + (wN if i == j else 0)
        comps[(i, j)] = (DDN - N * Ric).subs(Y, 0)
    f = sp.lambdify((X, Z, D), [comps[(0, 0)], comps[(1, 1)], comps[(2, 2)], comps[(0, 2)]],
                    "math", cse=True)
    return f


class BL:
    """Brill-Lindquist data with separation d, evolved to first order in t with a lapse."""

    def __init__(self, d, t=0.0, lapse="one", rtol=1e-11):
        self.d, self.t, self.lapse, self.rtol = d, t, lapse, rtol
        self._kdot = kdot_functions(lapse)

    def psi(self, rho, z, d=None):
        d = self.d if d is None else d
        return 1 + MASS / (2 * math.hypot(rho, z - d / 2)) + MASS / (2 * math.hypot(rho, z + d / 2))

    def grad_psi(self, rho, z, d=None):
        d = self.d if d is None else d
        r1, r2 = math.hypot(rho, z - d / 2), math.hypot(rho, z + d / 2)
        c1, c2 = MASS / (2 * r1**3), MASS / (2 * r2**3)
        return -(c1 + c2) * rho, -c1 * (z - d / 2) - c2 * (z + d / 2)

    def N(self, rho, z):
        return 1.0 if self.lapse == "one" else self.psi(rho, z) ** -2

    def p(self, rho, z, alpha):
        """tr_S Kdot = psi^-4 (Kdot(t,t) + Kdot(phi,phi)) in the flat orthonormal frame."""
        kxx, kyy, kzz, kxz = self._kdot(rho, z, self.d)
        c, s = math.cos(alpha), math.sin(alpha)
        ktt = c * c * kxx + 2 * c * s * kxz + s * s * kzz
        return (ktt + kyy) / self.psi(rho, z) ** 4

    def n_dlogpsi(self, rho, z, alpha, d=None):
        pr, pz = self.grad_psi(rho, z, d)
        return (-math.sin(alpha) * pr + math.cos(alpha) * pz) / self.psi(rho, z, d)

    # theta_+ = psi^-2 (H_flat + 4 n.grad psi / psi) + t p = 0, H_flat = -alpha' - sin(alpha)/rho
    def _rhs(self, s, y):
        rho, z, alpha = y
        extra = self.t * self.psi(rho, z) ** 2 * self.p(rho, z, alpha) if self.t else 0.0
        return [math.cos(alpha), math.sin(alpha),
                -math.sin(alpha) / rho + 4 * self.n_dlogpsi(rho, z, alpha) + extra]

    def shoot(self, z0):
        k0 = (4 * self.n_dlogpsi(0.0, z0, 0.0)
              + (self.t * self.psi(0.0, z0) ** 2 * self.p(0.0, z0, 0.0) if self.t else 0.0)) / 2
        ds = 1e-7
        y0 = [ds, z0 + k0 * ds**2 / 2, k0 * ds]

        def hit_axis(s, y):
            return y[0] - 1e-10
        hit_axis.terminal, hit_axis.direction = True, -1

        def hit_plane(s, y):
            return y[1]
        hit_plane.terminal, hit_plane.direction = True, -1
        return solve_ivp(self._rhs, (ds, 20.0), y0, method="DOP853", rtol=self.rtol,
                         atol=self.rtol * 1e-2, events=[hit_axis, hit_plane], dense_output=True)

    def mismatch(self, z0):
        sol = self.shoot(z0)
        return sol.y_events[1][0][2] + math.pi / 2 if sol.t_events[1].size else math.nan

    def min_mismatch(self, lo, hi):
        zs = np.linspace(lo, hi, 40)
        vals = np.array([self.mismatch(z) for z in zs])
        k = int(np.nanargmin(vals))
        a, b = zs[max(k - 1, 0)], zs[min(k + 1, len(zs) - 1)]
        res = minimize_scalar(self.mismatch, bounds=(a, b), method="bounded",
                              options={"xatol": 1e-12})
        return res.fun, res.x


ZLO, ZHI = 0.45, 1.2  # search window for the axis height of common surfaces (d near 0.77)


def critical_separation(lapse="one"):
    """d_c: the minimum of the t = 0 mismatch touches zero (degenerate minimal surface)."""
    def g(d):
        return BL(d, 0.0, lapse).min_mismatch(ZLO, ZHI)[0]
    return brentq(g, 0.70, 0.85, xtol=1e-12)


def onset_time(d, lapse, t_hi=0.05):
    """Smallest t >= 0 at which H_S + t p = 0 has a symmetric common solution (d > d_c)."""
    def g(t):
        return BL(d, t, lapse).min_mismatch(ZLO, ZHI)[0]
    return brentq(g, 0.0, t_hi, xtol=1e-13)


def fold_coefficients(dc, lapse, h=1e-5, npts=4001):
    """Terms of (*) on the degenerate surface S_c, by quadrature in flat arclength s.

    phi_0: proper normal displacement psi^2 n.(dX/dz0) at fixed s (finite differences at
    fixed d = d_c; at the fold this variation is a Jacobi field). dA = 2 pi rho psi^4 ds,
    doubled for the lower half (everything is symmetric under z -> -z).
    """
    bl = BL(dc, 0.0, lapse)
    z0c = bl.min_mismatch(ZLO, ZHI)[1]
    sol0, solp, solm = bl.shoot(z0c), bl.shoot(z0c + h), bl.shoot(z0c - h)
    s_end = sol0.t_events[1][0]
    s = np.linspace(1e-6, s_end, npts)
    rho, z, alpha = sol0.sol(s)
    dX = (solp.sol(s)[:2] - solm.sol(s)[:2]) / (2 * h)
    n = np.array([-np.sin(alpha), np.cos(alpha)])
    psi = np.array([bl.psi(a, b) for a, b in zip(rho, z)])
    phi = psi**2 * np.sum(n * dX, axis=0)
    dA = 2 * 2 * math.pi * rho * psi**4  # per unit s, both halves
    # |A|^2 from gamma-principal curvatures psi^-2 (kappa_flat + 2 dpsi/dn / psi)
    alpha_s = np.gradient(alpha, s)
    ndl = np.array([bl.n_dlogpsi(a, b, c) for a, b, c in zip(rho, z, alpha)])
    k_m = (-alpha_s + 2 * ndl) / psi**2
    k_p = (-np.sin(alpha) / rho + 2 * ndl) / psi**2
    A2 = k_m**2 + k_p**2
    Nvals = np.array([bl.N(a, b) for a, b in zip(rho, z)])
    pvals = np.array([bl.p(a, b, c) for a, b, c in zip(rho, z, alpha)])
    # dH/dd on the fixed surface: H = psi^-2 (H_flat + 4 n.grad psi/psi), H = 0 on S_c
    eps = 1e-6
    dndl = np.array([(bl.n_dlogpsi(a, b, c, dc + eps) - bl.n_dlogpsi(a, b, c, dc - eps)) / (2 * eps)
                     for a, b, c in zip(rho, z, alpha)])
    dHdd = 4 * dndl / psi**2
    H = (-alpha_s - np.sin(alpha) / rho + 4 * ndl) / psi**2
    integ = lambda f: simpson(f * dA, x=s)  # noqa: E731
    return {
        "z0c": z0c,
        "area": integ(np.ones_like(s)),
        "max|H| on S_c": float(np.max(np.abs(H[5:-5]))),
        "phi0 min/max": (float(phi.min()), float(phi.max())),
        "<phi0,dH/dd>": integ(phi * dHdd),
        "<phi0,p>": integ(phi * pvals),
        "-Int N phi0 |A|^2": -integ(Nvals * phi * A2),
        "Int N phi0 |A|^2": integ(Nvals * phi * A2),
    }


def certify_individual(d, R_in=0.01, R_out=0.3, pieces=4000):
    """Interval-arithmetic barrier certificate for the minimal surface around the upper
    puncture at t = 0 (theta_+ = H_S since K = 0): H_S < 0 on the sphere of radius R_in
    about z = d/2 and H_S > 0 on the sphere of radius R_out, with H_S = psi^-2 (2/R +
    4 dpsi/dn / psi). By Andersson-Metzger Thm 3.1 a MOTS (here a minimal surface)
    enclosing only the upper puncture lies between them; the lower one is its mirror image.
    d may be a pair (lo, hi) of decimal strings: the check then holds for every separation
    in [lo, hi] at once. Returns the two certified bounds, or raises if a check fails."""
    import mpmath as mp
    mp.iv.dps = 30
    iv = mp.iv
    dd = iv.mpf([str(d[0]), str(d[1])]) if isinstance(d, tuple) else iv.mpf(repr(d))
    m = iv.mpf("0.5")

    def H(R, th):
        R = iv.mpf(repr(R))
        c, s = iv.cos(th), iv.sin(th)
        psi, dpsi = iv.mpf(1), iv.mpf(0)
        # displacement from each puncture, written without cancellation so that an
        # interval d does not inflate it: own puncture (R s, R c), other (R s, R c + d)
        for dz in (R * c, R * c + dd):
            r = iv.sqrt((R * s) ** 2 + dz**2)
            psi = psi + m / (2 * r)
            dpsi = dpsi - m * (R * s * s + dz * c) / (2 * r**3)
        return (2 / R + 4 * dpsi / psi) / psi**2

    bounds = []
    for R, sign in ((R_in, -1), (R_out, +1)):
        worst = math.inf
        for k in range(pieces):
            th = iv.mpf([mp.mpf(k) / pieces, mp.mpf(k + 1) / pieces]) * iv.pi
            lower = float((H(R, th) * sign).a)
            if lower <= 0:
                raise AssertionError(f"barrier at radius {R} not certified (bound {lower})")
            worst = min(worst, lower)
        bounds.append(sign * worst)
    return bounds


def main():
    print("Brill-Lindquist, bare masses 1/2 each (ADM mass 1), punctures at z = +-d/2")
    dc = critical_separation()
    print(f"critical separation (common minimal surface born): d_c = {dc:.10f} M"
          f"  [= {2 * dc:.8f} in units of the bare mass; Hui & Lin 2025: 1.5323948, 1.5323949]")
    for lapse in ("one", "psi-2"):
        c = fold_coefficients(dc, lapse)
        slope = c["<phi0,dH/dd>"] / c["Int N phi0 |A|^2"]
        print(f"\nlapse N = {'1' if lapse == 'one' else 'psi^-2'}")
        print(f"  area of S_c: {c['area']:.6f} M^2 [x4 = {4 * c['area']:.5f}; Hui & Lin: 196.40795]")
        print(f"  axis height of S_c: {c['z0c']:.8f}; max |H| on S_c: {c['max|H| on S_c']:.1e}; "
              f"phi_0 in [{c['phi0 min/max'][0]:.4f}, {c['phi0 min/max'][1]:.4f}] (one sign)")
        print(f"  AMMS identity: <phi0, p> = {c['<phi0,p>']:.6e}  vs  -Int N phi0 |A|^2 = "
              f"{c['-Int N phi0 |A|^2']:.6e}")
        print(f"  <phi0, dH/dd> = {c['<phi0,dH/dd>']:.6e}  =>  predicted slope dt*/dd = {slope:.6f}")
        for delta in (2e-4, 1e-4, 5e-5):
            ts = onset_time(dc + delta, lapse)
            print(f"  direct: d = d_c + {delta:.0e}: t_* = {ts:.6e}, t_*/delta = {ts/delta:.6f}")


if __name__ == "__main__":
    main()
