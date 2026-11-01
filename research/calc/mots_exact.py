"""Exact evolved-data shooting identities for note 22.

These algebraic functions accept floats or Arb balls. They are NOT a global
MOTS solver or a new spacetime certificate. The caller must supply valid
enclosures of the TRUE metric, its derivatives, and K on its whole domain.
See bl_barrier.theta_general for the component layout and sign conventions.
"""

from research.calc.bl_barrier import theta_general


def alpha_prime(sa, ca, k2, gam, dg, K, sqrt_=None):
    """Solve theta_+ = 0 for flat-arclength curvature alpha'.

    sa=sin(alpha), ca=cos(alpha), k2=sin(alpha)/rho. Assumes a positive
    definite axisymmetric metric with meridional reflection symmetry.
    Interval arguments enclose every corresponding exact RHS value; they
    do not assert that arbitrary interval selections form a smooth metric.
    """
    sq = sqrt_ or (lambda x: x.sqrt())
    xx, yy, zz, xz = (gam[k] for k in ("xx", "yy", "zz", "xz"))
    det = xx * zz - xz * xz
    n2 = (zz * sa * sa + 2 * xz * sa * ca + xx * ca * ca) / det
    gtt = xx * ca * ca + 2 * xz * ca * sa + zz * sa * sa
    # theta(alpha') = theta(0) - alpha'/(sqrt(n2)*gtt).
    return theta_general(sa, ca, 0, k2, gam, dg, K, sqrt_=sq) * sq(n2) * gtt


def upper_axis_curvature(gam, dg, K, sqrt_=None):
    """Regular upper-pole limit alpha'(0), with tangent toward +rho.

    Rotational smoothness at rho=0 gives gxx=gyy, gxz=0,
    Gamma_{z,xx}=Gamma_{z,yy}, and sin(alpha)/rho -> alpha'(0).
    This is a necessary axis-start coefficient, not a validated axis collar.
    """
    sq = sqrt_ or (lambda x: x.sqrt())
    return K["xx"] / sq(gam["zz"]) - (
        dg[("x", "xz")] - dg[("z", "xx")] / 2
    ) / gam["zz"]


def fold_time_derivatives(Fd, Ft, Fzz, Fzd, Fzt, Fdd, Fdt, Ftt):
    """(t'(dc), z'(dc), t''(dc)) for F=F_z=0 at a simple fold.

    All arguments are derivatives of the TRUE shooting/reduced function at
    the critical point, not just its first-order-in-time truncation.
    No numerical BL value for the second derivative is asserted here.
    """
    c = -Fd / Ft
    z1 = -(Fzd + c * Fzt) / Fzz
    t2 = -(Fdd + 2 * c * Fdt + c * c * Ftt
           - (Fzd + c * Fzt) ** 2 / Fzz) / Ft
    return c, z1, t2
