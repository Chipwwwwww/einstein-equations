"""Gap 3: the thin-neck Jacobi field computed a second way, with Lohner enclosures
(research/calc/lohner.py). Used as an independent cross-check of research/calc/gap3_neck.py at
points (tests/test_gap3.py); NOT part of the proof.

Status (sessions 7-8). Meant to extend the thin-neck lemma to necks of radius up to ~0.03,
where gap3_neck's plain interval steps need impractically small parameter pieces.
  * Session 7: with (lam, z_w) as state variables, pieces of width ~1e-3 stay tight for about
    0.7 in sigma and then blow up (the constant parameters were mixed into the QR frame).
  * Session 8: the parameters are kept out of the frame (lohner.PSet, x = m + A r + P q). At
    lam ~ 0.03, pieces of +-1e-4 now certify phi(+-1.5) with radius < 1e-2 against a margin
    >= 0.2; pieces of +-5e-4 still fail (the excess is second order in the piece width and
    grows fast in sigma). Covering the waists that the approach lemma allows for the medium
    necks would take ~2e4 pieces, so the region certificate has not been run.

Same lemma as research/calc/gap3_neck.py (read its docstring first): in the rescaled variables
around a waist of radius lam at height z_w, the lam-variation Jacobi field
    phi = -V (R + lam R_lam) + U (Z + lam Z_lam)
must have the sign opposite to phi(0) = -s at sigma = +-sigma1. There, lam and z_w are balls
carried through plain interval Taylor steps, which wrap badly once lam ~ 0.01 (pieces of width
2^-10 give phi(+-1.5) only to +-1). Here (lam, z_w) are state variables with zero derivative,
and the log-derivatives Rc = lam R_lam, Zc = lam Z_lam, Uc = lam U_lam, Vc = lam V_lam satisfy
the variational equations, so the whole 10-dimensional system goes through Lohner's method,
which follows the dependence on the parameters instead of re-boxing it at every step. Then
phi = -V (R + Rc) + U (Z + Zc). (With plain lam-derivatives the variational coefficients are
of order psi_rho / (lam psi) ~ 100, and the QR frames mix parameter and state directions
badly; the log-derivative lam d/dlam keeps them of order one and needs no division by lam.)

The variational right-hand side is not written out by hand: f_neck evaluates the rescaled
system on dual numbers (class Dl) seeded with lam d/dlam, whose derivative parts are exactly
the variational equations. Checked against gap3_neck.certify at points (tests/test_gap3.py).
"""

from flint import arb

from research.calc import lohner
from research.calc.lohner import Dl

QUARTER = arb(1) / 4
D78 = arb("0.78")
SIGMA_2 = 2.0


def f_neck(y, direction=1, d=D78):
    """Rescaled meridian system plus its variation lam d/dlam, state
    (R, Z, U, V, Rc, Zc, Uc, Vc, lam, z_w); direction -1 reverses sigma."""
    R, Z, U, V, Rl, Zl, Ul, Vl, lam, zw = y
    zero = lam * 0
    Rd, Zd, Ud, Vd = Dl(R, Rl), Dl(Z, Zl), Dl(U, Ul), Dl(V, Vl)
    lamd, zwd = Dl(lam, lam), Dl(zw, zero)  # lam d/dlam (lam) = lam
    rho, z = lamd * Rd, zwd + lamd * Zd
    psi, prh, pz = 1, 0, 0
    for c in (d / 2, -d / 2):
        Zc = z - c
        ir = 1 / (rho * rho + Zc * Zc).sqrt()
        ir3 = ir * ir * ir
        psi = psi + QUARTER * ir
        prh = prh - QUARTER * Rd * ir3  # psi_rho / lam
        pz = pz - QUARTER * Zc * ir3
    P = 4 * lamd * (Ud * pz - Vd * lamd * prh) / psi
    F = [Ud, Vd, Rd - Vd * P, Ud * P]
    out = [x.v for x in F] + [x.g for x in F]
    if direction < 0:
        out = [-x for x in out]
    return out + [zero, zero]


def _phi_parts(y):
    """phi and its gradient with respect to the state."""
    R, Z, U, V, Rl, Zl, Ul, Vl, lam, zw = y
    phi = -V * (R + Rl) + U * (Z + Zl)
    grad = [-V, U, Z + Zl, -(R + Rl), -V, U, 0, 0, 0, 0]
    return phi, grad


def phi_enclosure(S):
    """phi over the set S = m + A r + P q (lohner.PSet; phi depends on the state only), in
    mean-value form around m."""
    n, k = len(S.m), len(S.p0)
    phi_m, _ = _phi_parts(S.m + S.p0)
    _, grad = _phi_parts(S.hull() + S.params())
    acc = phi_m
    for j in range(n):
        acc += sum((grad[i] * S.A[i, j] for i in range(n)), arb(0)) * S.r[j]
    for l in range(k):
        acc += sum((grad[i] * S.P[i, l] for i in range(n)), arb(0)) * S.q[l]
    return acc


def certify(lam_lo, lam_hi, zw_lo, zw_hi, s=-1, sigma1=SIGMA_2, h=2.0**-4, L=14, d=D78):
    """(phi(0), phi(-sigma1), phi(+sigma1)) enclosures over lam in [lam_lo, lam_hi], z_w in
    [zw_lo, zw_hi], for a waist moving down (s = -1) or up (s = +1). The parameters (lam, z_w)
    are kept out of the moving frame (lohner.step_param)."""
    def ball(lo, hi):
        lo, hi = arb(lo), arb(hi)
        return (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
    lam, zw = ball(lam_lo, lam_hi), ball(zw_lo, zw_hi)
    x0 = [arb(1), arb(0), arb(0), arb(s), arb(0), arb(0), arb(0), arb(0)]
    out = [arb(-s)]  # phi(0) = -V R = -s exactly
    nsteps = round(sigma1 / h)
    for direction in (-1, 1):
        S = lohner.PSet.from_box(x0, [lam, zw])
        for _ in range(nsteps):
            S, _ = lohner.step_param(lambda y: f_neck(y, direction, d), S, h, L)
        out.append(phi_enclosure(S))
    return tuple(out)
