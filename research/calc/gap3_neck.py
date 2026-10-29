"""Gap 3, thin necks: a meridian curve with a thin waist cannot belong to a stable surface.

Rescale around a waist point (rho, z) = (lam, z_w) where the curve is vertical (u = 0,
v = s lam, s = -1 going down, +1 going up): rho = lam R, z = z_w + lam Z, u = lam U,
v = lam V. In sigma-time the meridian system (research/calc/gap3_sweep.py) becomes

    R' = U,  Z' = V,  U' = R - V P,  V' = U P,   P = 4 lam (U psi_z - V psi_rho_hat) / psi,

with psi_rho_hat = psi_rho / lam regular, all evaluated at (lam R, z_w + lam Z). This is regular
in lam down to lam = 0, where it is the catenoid R = cosh(sigma), Z = s sigma.

Jacobi field. Varying lam at fixed z_w moves the curve through minimal surfaces, so the normal
component of d/dlam (lam R, z_w + lam Z) is a rotationally symmetric Jacobi field. Up to the
positive factor 1/R (and psi^2) it is

    phi = -V (R + lam R_lam) + U (Z + lam Z_lam),

equal to -s (cosh(sigma) - sigma sinh(sigma)) on the catenoid, with zeros at sigma = +-1.1997.
If, for every lam in [0, lam0] and z_w in a range, phi(+-SIGMA_1) has the sign opposite to
phi(0), then on every such neck phi has a zero on each side of the waist. The annulus between
the two zeros nearest the waist then has a Dirichlet Jacobi eigenvalue 0, so the closed
surface containing it has a negative eigenvalue: it is unstable. The outermost minimal
surface (the boundary of the trapped region) is stable (Andersson-Metzger 2009), so it has no
such neck.
"""

from flint import arb

from research.calc import fold_proof as fp
from research.calc.fold_proof import J

QUARTER = arb(1) / 4
SIGMA_1 = 1.5


def make_rhs(lam, zw, d, sign=1):
    def rhs(y):
        _, R, Z, U, V = y
        rho, z = lam * R, zw + lam * Z
        psi, prh, pz = 1, 0, 0
        for c in (d * 0.5, -(d * 0.5)):
            Zc = z - c
            ir = fp.inv(fp.sqrt(rho * rho + Zc * Zc))
            ir3 = ir * ir * ir
            psi = psi + QUARTER * ir
            prh = prh - QUARTER * R * ir3  # psi_rho / lam
            pz = pz - QUARTER * Zc * ir3
        P = 4 * lam * (U * pz - V * lam * prh) / psi
        one = R * 0 + 1
        f = [one, U, V, R - V * P, U * P]
        return f if sign > 0 else [one] + [-x for x in f[1:]]
    return rhs


def phi(y, lam):
    _, R, Z, U, V = y
    return (-V.c[0] * (R.c[0] + lam * R.c[1]) + U.c[0] * (Z.c[0] + lam * Z.c[1]))


# The certified region: waists of radius lam in [0, LAM0] at heights z_w in [ZW_LO, ZW_HI]
# (dyadic numbers, so that every piece endpoint below is exact in binary).
LAM0 = 5 * 2.0**-10  # 0.0048828125
ZW_LO, ZW_HI = 23 * 2.0**-7, 31 * 2.0**-7  # [0.1796875, 0.2421875]


def _ball(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)


def certify(lam_lo, lam_hi, zw_lo, zw_hi, s=-1, d="0.78", h=2.0**-5, L=16):
    """Returns (phi(0), phi(-SIGMA_1), phi(+SIGMA_1)) enclosures over lam in [lam_lo, lam_hi],
    z_w in [zw_lo, zw_hi], for a waist moving down (s = -1) or up (s = +1)."""
    lam_b, zw_b = _ball(lam_lo, lam_hi), _ball(zw_lo, zw_hi)
    lam = J([lam_b, arb(1), arb(0), arb(0), arb(0), arb(0)])
    zw = J([zw_b, arb(0), arb(0), arb(0), arb(0), arb(0)])
    dj = J([arb(d), arb(0), arb(0), arb(0), arb(0), arb(0)])
    one, zero = J([arb(1)] + [arb(0)] * 5), J([arb(0)] * 6)
    y0 = [zero, one, zero, zero, one * s]
    out = [phi(y0, lam_b)]
    for direction in (-1, 1):
        y, _ = fp.integrate(make_rhs(lam, zw, dj, direction), y0, 0.0, SIGMA_1, h, L)
        out.append(phi(y, lam_b))
    return tuple(out)


def check_piece(piece):
    """piece = (lam_lo, lam_hi, zw_lo, zw_hi, s). Returns the margin (min of -s phi(0) and
    s phi(+-SIGMA_1), lower bounds) if the sign pattern is certified, else None."""
    *box, s = piece
    try:
        p0, pm, pp = certify(*box, s=s)
    except (ValueError, ZeroDivisionError, RuntimeError):
        return None
    # phi(0) = -s; certified: -s phi(0) > 0 and s phi(+-SIGMA_1) > 0
    vals = [-s * p0, s * pm, s * pp]
    if all(v > 0 for v in vals):
        return min(float(v.lower()) for v in vals)
    return None


def _split(piece):
    l0, l1, z0, z1, s = piece
    lm, zm = (l0 + l1) / 2, (z0 + z1) / 2
    return [(a, b, c, e, s) for a, b in ((l0, lm), (lm, l1)) for c, e in ((z0, zm), (zm, z1))]


def certify_region(lam0=LAM0, zw_lo=ZW_LO, zw_hi=ZW_HI, dlam=2.0**-11, dzw=2.0**-9,
                   min_width=2.0**-15, processes=4, log=print):
    """The thin-neck lemma on [0, lam0] x [zw_lo, zw_hi], both directions: start from pieces of
    size dlam x dzw and quarter every piece whose sign pattern is not certified. Returns
    (number of certified pieces, worst margin). Raises if a piece gets narrower than min_width."""
    from multiprocessing import Pool
    nl, nz = round(lam0 / dlam), round((zw_hi - zw_lo) / dzw)
    todo = [(i * dlam, (i + 1) * dlam, zw_lo + k * dzw, zw_lo + (k + 1) * dzw, s)
            for i in range(nl) for k in range(nz) for s in (-1, 1)]
    done, worst = 0, float("inf")
    with Pool(processes) as pool:
        while todo:
            results = pool.map(check_piece, todo, chunksize=1)
            failed = [p for p, m in zip(todo, results) if m is None]
            ok = [m for m in results if m is not None]
            done += len(ok)
            worst = min([worst] + ok)
            if log:
                log(f"  round: {len(todo)} pieces, {len(failed)} split; certified so far {done}, "
                    f"worst margin {worst:.4f}")
            todo = []
            for p in failed:
                if p[1] - p[0] < min_width or p[3] - p[2] < min_width:
                    raise AssertionError(f"neck lemma fails on {p}")
                todo += _split(p)
    return done, worst


def main():
    n, w = certify_region()
    print(f"Thin-neck lemma: lam in [0, {LAM0}], z_w in [{ZW_LO}, {ZW_HI}], both directions, "
          f"{n} pieces: phi changes sign on both sides within |sigma| <= {SIGMA_1}, "
          f"margin >= {w:.4f}")


if __name__ == "__main__":
    main()
