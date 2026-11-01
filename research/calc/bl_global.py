"""Gap 3 (research/status.md): toward a proof that time-symmetric Brill-Lindquist data with
d > d_c have NO common minimal surface (no common apparent horizon at t = 0).

Data and conventions as in research/calc/bl_fold.py: two punctures of bare mass m = 1/2 at
z = +-d/2, gamma = psi^4 delta, psi = 1 + m/(2 r1) + m/(2 r2), ADM mass M = 1.

What this module establishes (rigorously, by interval arithmetic in mpmath.iv):

  Lemma G1 (throat barrier). For s in (0, s_*], every coordinate sphere |x - p+| = s has
    H < 0 for the normal pointing away from p+. Hence (maximum principle, notes on gap 3 in
    research/idea_log.md) a closed minimal surface enclosing p+ stays outside B(p+, s_*); by
    symmetry the same holds at p-.
  Lemma G2 (far barrier). For R >= R_*, every sphere |x| = R has H > 0. Hence every closed
    minimal surface lies in B(0, R_*). For R >= 2 this is the analytic bound
    psi^2 H >= 2/R - 2/(R - d/2)^2 > 0, valid when (R - d/2)^2 > R.

Consequence. A common minimal surface S* that is G-invariant (as the outermost one is: it is
unique, hence invariant under rotations about the axis and z -> -z) has an embedded profile
curve from the axis point (0, z0) to the plane z = 0, meeting the plane once and at a right
angle, with d/2 + s_* <= z0 < R_* and the curve inside B(0, R_*) minus the two throat balls.

What it does NOT establish: the sweep over z0. landscape() classifies the shooting curves
(floating point) and shows the obstruction: just above the top of the individual horizon S_+,
the curves hug S_+, pass through a neck of radius -> 0 near the axis below it, and flare out.
In the arclength parameter the axis is singular (alpha' = -sin(alpha)/rho + G). With
ds = rho dsigma and (u, v) = rho (cos alpha, sin alpha) the system becomes

    rho' = u,  z' = v,  u' = rho - v P,  v' = u P,   P := rho G = 4 (u psi_z - v psi_rho) / psi,

which is regular (polynomial in u, v) up to the axis; u^2 + v^2 = rho^2 is invariant. The
axis {rho = u = v = 0} is a line of fixed points with eigenvalues +-1 in the (rho, u) plane:
the necked curves pass near a hyperbolic fixed point. A proof needs an isolating-block
(Shilnikov-type) argument there; see research/status.md. (Both identities: tests/test_bl_global.py.)

Run: python -m research.calc.bl_global
"""

import math

import mpmath as mp
import numpy as np

from research.calc.bl_fold import BL, MASS

D_DEFAULT = "0.78"


def _iv(x):
    return mp.iv.mpf(x) if isinstance(x, (list, tuple)) else mp.iv.mpf(str(x))


def throat_quantity(s, th, d):
    """s psi^2 H for the sphere |x - p+| = s, normal pointing away from p+, written so that
    it stays regular at s = 0 (where it tends to 2 - 4 = -2):

        2 + 4 (-m/2 + s^2 dn(m / (2 r2))) / (s + m/2 + s m / (2 r2)).
    Arguments may be mpmath intervals."""
    m = MASS if not isinstance(s, mp.ctx_iv.ivmpf) else _iv("0.5")
    sn, cs = (mp.iv.sin(th), mp.iv.cos(th)) if isinstance(th, mp.ctx_iv.ivmpf) else (math.sin(th), math.cos(th))
    x = s * sn
    zz = s * cs + d  # z + d/2 measured from p- ; p+ at the origin of this local frame
    r2 = (mp.iv.sqrt(x * x + zz * zz) if isinstance(x, mp.ctx_iv.ivmpf) else math.hypot(x, zz))
    dn2 = -m / (2 * r2**3) * (x * sn + zz * cs)
    return 2 + 4 * (-m / 2 + s * s * dn2) / (s + m / 2 + s * m / (2 * r2))


def far_quantity(R, th, d):
    """psi^2 H for the sphere |x| = R (outward normal). Arguments may be intervals."""
    m = _iv("0.5") if isinstance(R, mp.ctx_iv.ivmpf) else MASS
    iv = isinstance(R, mp.ctx_iv.ivmpf) or isinstance(th, mp.ctx_iv.ivmpf)
    sqrt = mp.iv.sqrt if iv else math.sqrt
    sn, cs = (mp.iv.sin(th), mp.iv.cos(th)) if iv else (math.sin(th), math.cos(th))
    psi, dpsi = 1, 0
    for c in (d / 2, -d / 2):
        dz = R * cs - c
        r = sqrt((R * sn) ** 2 + dz * dz)
        psi = psi + m / (2 * r)
        dpsi = dpsi - m * (R * sn * sn + dz * cs) / (2 * r**3)
    return 2 / R + 4 * dpsi / psi


def certify_throat(d=D_DEFAULT, s_star="0.155", n_s=64, n_th=256):
    """Lemma G1: throat_quantity < 0 on (0, s_*] x [0, pi]. Returns the largest upper bound."""
    mp.iv.dps = 30
    dd, S = _iv(d), mp.mpf(s_star)
    worst = -mp.inf
    for i in range(n_s):
        s = mp.iv.mpf([S * i / n_s, S * (i + 1) / n_s])
        for k in range(n_th):
            th = mp.iv.mpf([mp.mpf(k) / n_th, mp.mpf(k + 1) / n_th]) * mp.iv.pi
            up = throat_quantity(s, th, dd).b
            if up >= 0:
                raise AssertionError(f"throat barrier fails near s={float(s.a):.4f}")
            worst = max(worst, up)
    return float(worst)


def certify_far(d=D_DEFAULT, R_star="1.05", R_switch=2, n_R=256, n_th=256):
    """Lemma G2: far_quantity > 0 for R in [R_*, R_switch] (interval boxes) and for R >= R_switch
    (analytic bound). Returns the smallest certified lower bound on the boxes."""
    mp.iv.dps = 30
    dd = _iv(d)
    lo, hi = mp.mpf(R_star), mp.mpf(R_switch)
    a = (hi - dd.b / 2) ** 2 - hi  # analytic tail: (R - d/2)^2 - R increasing for R > d/2 + 1/2
    if not (a > 0 and hi > dd.b / 2 + mp.mpf(1) / 2):
        raise AssertionError("analytic tail bound does not apply")
    worst = mp.inf
    for i in range(n_R):
        R = mp.iv.mpf([lo + (hi - lo) * i / n_R, lo + (hi - lo) * (i + 1) / n_R])
        for k in range(n_th):
            th = mp.iv.mpf([mp.mpf(k) / n_th, mp.mpf(k + 1) / n_th]) * mp.iv.pi
            low = far_quantity(R, th, dd).a
            if low <= 0:
                raise AssertionError(f"far barrier fails near R={float(R.a):.4f}")
            worst = min(worst, low)
    return float(worst)


def landscape(d=0.78, s_star=0.155, R_star=1.05, n=60):
    """Floating-point classification of the shooting curves from (0, z0), z0 in the window
    [d/2 + s_*, R_*]: the first of 'plane' (with F = alpha + pi/2 there), 'throat' (enters a
    throat ball), 'leaves' (leaves B(0, R_*)). Also the minimal distance to the axis after the
    start (the neck)."""
    bl = BL(d)
    rows = []
    for z0 in np.linspace(d / 2 + s_star, R_star, n):
        sol = bl.shoot(z0)
        s = np.linspace(sol.t[0], sol.t[-1], 20000)
        rho, z, a = sol.sol(s)
        r = np.hypot(rho, z)
        inner = (np.hypot(rho, z - d / 2) < s_star) | (np.hypot(rho, z + d / 2) < s_star)
        events = []
        if (r > R_star).any():
            events.append((s[np.argmax(r > R_star)], "leaves"))
        if inner.any():
            events.append((s[np.argmax(inner)], "throat"))
        if sol.t_events[1].size:
            events.append((sol.t_events[1][0], "plane"))
        first = min(events) if events else (None, "none")
        F = sol.y_events[1][0][2] + math.pi / 2 if first[1] == "plane" else math.nan
        back = a < -math.pi / 2  # after the curve has turned back toward the axis
        neck = float(rho[np.argmax(back):].min()) if back.any() else math.nan
        rows.append((float(z0), first[1], F, neck))
    return rows


def main():
    print("Lemma G1, d = 0.78, s_* = 0.155: max of s psi^2 H on the throat spheres <=",
          f"{certify_throat():.4f}")
    print("Lemma G2, d = 0.78, R_* = 1.05: min of psi^2 H on the far spheres >=",
          f"{certify_far():.5f}")
    print("\nShooting landscape (floating point):  z0, first event, F, neck radius")
    for z0, ev, F, neck in landscape():
        print(f"  {z0:.4f}  {ev:7s}  {F:8.4f}  {neck:.4f}")


if __name__ == "__main__":
    main()
