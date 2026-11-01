"""Apparent horizons of the two-centre Kastor-Traschen spacetime (notes/08, section 6).

Setting: H = -1 (contracting branch, tau < 0), equal masses m at z = +d and z = -d,
U = |tau| + m/r1 + m/r2. On a slice tau = const, K_ij = H gamma_ij, so a closed surface
is a MOTS (theta_+ = 0) iff its mean curvature in gamma = U^2 delta is H_S = 2|H| = 2.

Two calculations:

1. Numerical (evidence): axisymmetric common MOTSs by shooting. The generating curve
   (rho(s), z(s)), s = flat arclength, tangent angle alpha, outward normal
   n = (-sin alpha, cos alpha), satisfies

       alpha' = -sin(alpha)/rho - 2|H| U + 2 (n . grad U) / U,

   which is H_S = U^-1 (H_flat + 2 n . grad U / U) = 2|H| with
   H_flat = -alpha' - sin(alpha)/rho. A common MOTS is symmetric under z -> -z, so we
   shoot from the top of the axis, z = z0 > d, and require a perpendicular crossing of
   z = 0. The first time tau_c at which one exists is where the minimum over z0 of the
   mismatch touches zero (a fold: an inner and an outer common MOTS are born together).

   Individual MOTSs sit inside the near-cylindrical extremal throats, where shooting is
   ill-conditioned; their existence is settled by barriers instead.

2. Computer-assisted (rigorous, given Andersson-Metzger 2009, Thm 3.1): coordinate
   spheres as barriers. theta_+ of a sphere is an explicit function of the polar angle;
   mpmath interval arithmetic encloses it on a subdivision of [0, pi]. A compact region
   whose inner boundary has theta_+ < 0 (normal pointing into the region) and whose
   outer boundary has theta_+ > 0 (normal pointing out) contains a smooth embedded MOTS
   homologous to the outer boundary.

Run: python -m research.calc.kt_horizons   (about two minutes)
"""

import math

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

H_ABS = 1.0  # |H|; by the homothety of notes/08 this fixes the unit of length


def _exact(iv):
    """Interval conversion of a parameter given as a decimal: repr() recovers the shortest
    decimal string of a float (e.g. '0.05'), and iv.mpf of that string is an interval
    containing the exact decimal value, so certificates hold for the stated decimals."""
    return lambda v: v if isinstance(v, type(iv.mpf(0))) else iv.mpf(repr(v) if isinstance(v, float) else str(v))


def s_bh(M):
    """|H| r |tau| on the black-hole horizon of a single centre of mass M (notes/08)."""
    return ((1 - 2 * M) - math.sqrt(1 - 4 * M)) / 2


class KT:
    """Equal-mass, two-centre Kastor-Traschen slice at time tau < 0 (H = -1)."""

    def __init__(self, m, d, tau, rtol=1e-11):
        self.m, self.d, self.tau, self.rtol = m, d, tau, rtol

    def U(self, rho, z):
        r1, r2 = math.hypot(rho, z - self.d), math.hypot(rho, z + self.d)
        return H_ABS * abs(self.tau) + self.m / r1 + self.m / r2

    def grad_U(self, rho, z):
        r1, r2 = math.hypot(rho, z - self.d), math.hypot(rho, z + self.d)
        c1, c2 = self.m / r1**3, self.m / r2**3
        return -(c1 + c2) * rho, -c1 * (z - self.d) - c2 * (z + self.d)

    # ---- shooting for symmetric common MOTSs ----
    def _rhs(self, s, y):
        rho, z, alpha = y
        U = self.U(rho, z)
        Ur, Uz = self.grad_U(rho, z)
        n_grad = -math.sin(alpha) * Ur + math.cos(alpha) * Uz
        return [math.cos(alpha), math.sin(alpha),
                -math.sin(alpha) / rho - 2 * H_ABS * U + 2 * n_grad / U]

    def shoot(self, z0):
        """Integrate from the axis point (0, z0), heading right, to the plane z = 0."""
        U = self.U(0.0, z0)
        _, Uz = self.grad_U(0.0, z0)
        k0 = (-2 * H_ABS * U + 2 * Uz / U) / 2  # alpha'(0), from regularity on the axis
        ds = 1e-6 * max(1.0, z0)
        y0 = [ds, z0 + k0 * ds**2 / 2, k0 * ds]

        def hit_axis(s, y):
            return y[0] - 1e-9
        hit_axis.terminal, hit_axis.direction = True, -1

        def hit_plane(s, y):
            return y[1]
        hit_plane.terminal, hit_plane.direction = True, -1

        return solve_ivp(self._rhs, (ds, 50.0 * max(1.0, z0)), y0, method="DOP853",
                         rtol=self.rtol, atol=self.rtol * 1e-2, events=[hit_axis, hit_plane],
                         dense_output=True)

    def mismatch(self, z0):
        """alpha + pi/2 at the first crossing of z = 0 (zero for a symmetric MOTS);
        nan if the curve returns to the axis first or never crosses."""
        sol = self.shoot(z0)
        if sol.t_events[1].size:
            return sol.y_events[1][0][2] + math.pi / 2
        return math.nan

    def common_mots(self, lo, hi, n=200):
        """Top heights z0 of symmetric common MOTSs with z0 in (lo, hi). A bracket counts
        only if f is finite at both ends and brentq converges to |f| < 1e-8 (so jumps of
        alpha by 2 pi are not mistaken for roots)."""
        zs = np.geomspace(lo, hi, n)
        vals = [self.mismatch(z) for z in zs]
        out = []
        for a, b, fa, fb in zip(zs, zs[1:], vals, vals[1:]):
            if np.isfinite(fa) and np.isfinite(fb) and fa * fb < 0:
                try:
                    root = brentq(self.mismatch, a, b, xtol=1e-13)
                except ValueError:
                    continue
                if abs(self.mismatch(root)) < 1e-8:
                    out.append(root)
        return out

    def min_mismatch(self, lo, hi):
        """Minimum of the mismatch over z0 in [lo, hi] (the fold diagnostic)."""
        zs = np.linspace(lo, hi, 60)
        vals = np.array([self.mismatch(z) for z in zs])
        k = int(np.nanargmin(vals))
        a, b = zs[max(k - 1, 0)], zs[min(k + 1, len(zs) - 1)]
        res = minimize_scalar(self.mismatch, bounds=(a, b), method="bounded",
                              options={"xatol": 1e-10})
        return res.fun, res.x

    # ---- theta_+ of coordinate spheres ----
    def theta_plus_sphere(self, centre_z, R, theta, iv=None):
        """theta_+ = 2H + U^-1 (2/R + 2 dU/dn / U) on the sphere |x - (0,0,centre_z)| = R
        at polar angle theta. With iv=mpmath.iv the same formula runs in interval
        arithmetic (inputs converted exactly from their binary floating-point values)."""
        lib, conv = (math, float) if iv is None else (iv, _exact(iv))
        m, d, tau_abs, R, cz = (conv(v) for v in (self.m, self.d, abs(self.tau), R, centre_z))
        c, s = lib.cos(theta), lib.sin(theta)
        px, pz = R * s, cz + R * c  # point on the sphere: (rho, z)
        U = tau_abs * conv(H_ABS)
        dU = conv(0)
        for zc in (d, -d):
            r = lib.sqrt(px**2 + (pz - zc) ** 2)
            U = U + m / r
            dU = dU - m * (px * s + (pz - zc) * c) / r**3  # n . grad(m/r), n = (s, c)
        return -2 * conv(H_ABS) + (2 / R + 2 * dU / U) / U


def theta_plus_polar(kt, coeffs, theta, iv=None):
    """theta_+ on the axisymmetric surface r = f(theta) = sum_k a_k cos(2k theta) about the
    origin (symmetric under z -> -z), in the Kastor-Traschen slice kt.

    Flat geometry of the meridian curve (rho, z) = f (sin, cos):
      kappa_meridian = (f^2 + 2 f'^2 - f f'') / v^3,   v = sqrt(f^2 + f'^2),
      kappa_parallel = -(z'/sin theta) / (v f),        z'/sin = (f'/sin) cos - f,
    with f'/sin theta = -sum 2k a_k U_{2k-1}(cos theta) (Chebyshev U), so nothing is
    divided by sin theta and the formula is valid on the axis. Outward normal
    n = (-z', rho') / v. Then theta_+ = 2H + U^-1 (kappa_m + kappa_p + 2 n.grad U / U).
    """
    lib, conv = (math, float) if iv is None else (iv, _exact(iv))
    m, d, tau_abs = (conv(v) for v in (kt.m, kt.d, abs(kt.tau)))
    a = [conv(c) for c in coeffs]
    c, s = lib.cos(theta), lib.sin(theta)
    # cos(2k theta), sin(2k theta)/sin(theta) = U_{2k-1}(c), via recurrences in c
    cheb_T = [conv(1), c]  # T_n(c) = cos(n theta)
    cheb_U = [conv(1), 2 * c]  # U_n(c) = sin((n+1) theta)/sin(theta)
    for n in range(2, 2 * len(a) + 1):
        cheb_T.append(2 * c * cheb_T[-1] - cheb_T[-2])
        cheb_U.append(2 * c * cheb_U[-1] - cheb_U[-2])
    f = sum(ak * cheb_T[2 * k] for k, ak in enumerate(a))
    fp_over_s = -sum(2 * k * ak * cheb_U[2 * k - 1] for k, ak in enumerate(a) if k > 0)
    fp = fp_over_s * s
    fpp = -sum(4 * k * k * ak * cheb_T[2 * k] for k, ak in enumerate(a))
    v = lib.sqrt(f**2 + fp**2)
    kappa_m = (f**2 + 2 * fp**2 - f * fpp) / v**3
    zp_over_s = fp_over_s * c - f
    kappa_p = -zp_over_s / (v * f)
    rho_p, z_p = fp * s + f * c, zp_over_s * s
    n_rho, n_z = -z_p / v, rho_p / v
    px, pz = f * s, f * c
    U = tau_abs * conv(H_ABS)
    n_grad = conv(0)
    for zc in (d, -d):
        r = lib.sqrt(px**2 + (pz - zc) ** 2)
        U = U + m / r
        n_grad = n_grad - m * (n_rho * px + n_z * (pz - zc)) / r**3
    return -2 * conv(H_ABS) + (kappa_m + kappa_p + 2 * n_grad / U) / U


def certify_polar(kt, coeffs, sign, pieces=4000):
    """Interval check of sign * theta_+ > 0 on r = f(theta); theta in [0, pi/2] suffices
    by the z -> -z symmetry of both the surface and the slice."""
    mp.iv.dps = 30
    worst = math.inf
    for k in range(pieces):
        th = mp.iv.mpf([mp.mpf(k) / pieces, mp.mpf(k + 1) / pieces]) * (mp.iv.pi / 2)
        lower = float((theta_plus_polar(kt, coeffs, th, iv=mp.iv) * sign).a)
        worst = min(worst, lower)
        if lower <= 0:
            return False, worst
    return True, worst


def fit_polar(kt, z0s, K=4):
    """Least-squares fit of r(theta) = sum_{k<K} a_k cos(2k theta) to the average of the
    common MOTSs with top heights z0s (inner and outer)."""
    th = np.linspace(0, math.pi / 2, 400)
    rs = []
    for z0 in z0s:
        sol = kt.shoot(z0)
        y = sol.sol(np.linspace(sol.t[0], sol.t[-1], 20000))  # dense output, not raw steps
        r, ang = np.hypot(y[0], y[1]), np.arctan2(y[0], y[1])
        rs.append(np.interp(th, np.concatenate([[0.0], ang]), np.concatenate([[z0], r])))
    target = np.mean(rs, axis=0)
    A = np.stack([np.cos(2 * k * th) for k in range(K)], axis=1)
    return np.linalg.lstsq(A, target, rcond=None)[0]


def certify_sign(kt, centre_z, R, sign, pieces=4000):
    """Rigorously check sign * theta_+ > 0 on the whole sphere by interval arithmetic on
    `pieces` subintervals of theta in [0, pi]. Returns (certified, lowest lower bound)."""
    mp.iv.dps = 30
    worst = math.inf
    for k in range(pieces):
        th = mp.iv.mpf([mp.mpf(k) / pieces, mp.mpf(k + 1) / pieces]) * mp.iv.pi
        lower = float((kt.theta_plus_sphere(centre_z, R, th, iv=mp.iv) * sign).a)
        worst = min(worst, lower)
        if lower <= 0:
            return False, worst
    return True, worst


def fold_time(m, d, tau_no, tau_yes, rtol=1e-11, zlo=None, zhi=None):
    """tau at which the minimum mismatch crosses zero, i.e. the first common MOTS."""
    zlo = zlo or 1.2 * d
    zhi = zhi or 6 * d

    def g(tau):
        return KT(m, d, tau, rtol=rtol).min_mismatch(zlo, zhi)[0]
    return brentq(g, tau_no, tau_yes, xtol=1e-10)


# Inner barrier for the common MOTS at tau = -0.009: the outer common MOTS of the tau = -0.01
# slice, fitted with eight cos(2k theta) modes and rounded to six significant figures. The
# certificate is about this explicit surface, however it was obtained.
COMMON_BARRIER = ["1.24825", "0.0959577", "0.0024235", "-5.71413e-05", "-6.30472e-06",
                  "6.79142e-08", "2.53883e-08", "2.40738e-10"]


def main():
    m, d = 0.05, 0.5
    print(f"Equal masses m = {m} at z = +-{d}; H = -1; total M = {2*m} < M_ext = 0.25")
    print(f"s_bh(m) = {s_bh(m):.6f}, s_bh(2m) = {s_bh(2*m):.6f}; BHKT order-of-magnitude "
          f"event-horizon topology change |tau| ~ s_bh(M)/(2d) = {s_bh(2*m)/(2*d):.5f}")

    print("\n1. Numerical: symmetric common MOTSs (top height z0 in (1.2 d, 6 d))")
    for tau in (-0.05, -0.02, -0.015, -0.0125, -0.01, -0.008):
        kt = KT(m, d, tau)
        print(f"   tau = {tau:8.4f}: z0 = {[round(z, 6) for z in kt.common_mots(1.2*d, 6*d)]}, "
              f"min mismatch = {kt.min_mismatch(1.2*d, 6*d)[0]:+.3e}")
    tc = fold_time(m, d, -0.02, -0.0125)
    tc8 = fold_time(m, d, -0.02, -0.0125, rtol=1e-8)
    print(f"   fold (first common MOTS): tau_c = {tc:.9f}  (rtol 1e-8: {tc8:.9f}; "
          f"difference {abs(tc - tc8):.1e})")

    print("\n2. Cross-check: theta_+ from the polar formula on fits of the shooting solutions")
    kt = KT(m, d, -0.01)
    zin, zout = kt.common_mots(1.2*d, 6*d)
    th = np.linspace(0, math.pi / 2, 1001)
    for name, z0, K in (("inner", zin, 14), ("outer", zout, 10)):
        a = fit_polar(kt, [z0], K)
        worst = max(abs(theta_plus_polar(kt, a, t)) for t in th)
        print(f"   {name} common MOTS at tau = -0.01, {K} modes: max |theta_+| = {worst:.2e}")

    print("\n3. Computer-assisted barriers (interval arithmetic; Andersson-Metzger Thm 3.1)")
    kt = KT(m, d, -0.05)
    for name, R, sign in (("inner sphere about upper centre", 0.001, -1),
                          ("outer sphere about upper centre", 0.2, +1)):
        ok, w = certify_sign(kt, d, R, sign)
        print(f"   tau = -0.05, {name}, radius {R}: {'theta_+ < 0' if sign < 0 else 'theta_+ > 0'} "
              f"certified: {ok} (sign*theta_+ >= {w:.4f})")
    kt = KT(m, d, -0.009)
    ok, w = certify_polar(kt, COMMON_BARRIER, -1)
    print(f"   tau = -0.009, inner barrier r = sum a_k cos(2k theta) (COMMON_BARRIER): "
          f"theta_+ < 0 certified: {ok} (-theta_+ >= {w:.4f})")
    ok, w = certify_sign(kt, 0.0, 3.0, +1)
    print(f"   tau = -0.009, outer sphere about origin, radius 3: theta_+ > 0 certified: {ok} "
          f"(theta_+ >= {w:.4f})")


if __name__ == "__main__":
    main()
