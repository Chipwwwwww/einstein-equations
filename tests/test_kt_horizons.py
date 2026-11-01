"""Apparent horizons of two-centre Kastor-Traschen (notes/08 section 6, research/calc/kt_horizons.py).

The certificates are rigorous up to mpmath's interval arithmetic and Andersson-Metzger
(2009) Theorem 3.1; the fold location is numerical evidence.
"""

import math

import pytest

from research.calc import kt_horizons as kh

M, D = 0.05, 0.5


@pytest.mark.parametrize("R", [0.3, 1.0, 2.5])
def test_polar_formula_reduces_to_sphere_formula(R):
    kt = kh.KT(M, D, -0.01)
    for th in (0.0, 0.2, 0.9, math.pi / 2):
        assert math.isclose(kh.theta_plus_polar(kt, [R], th), kt.theta_plus_sphere(0.0, R, th),
                            rel_tol=1e-12, abs_tol=1e-12)


def test_two_disjoint_mots_at_early_time():
    """tau = -0.05: in the annulus 0.001 <= |x - x_1| <= 0.2 the inner sphere is outer trapped
    and the outer one outer untrapped, so it contains a MOTS enclosing x_1 only. The slice is
    symmetric under z -> -z, so the same holds about x_2; the two balls of radius 0.2 about
    z = +-0.5 are disjoint."""
    kt = kh.KT(M, D, -0.05)
    assert kh.certify_sign(kt, D, 0.001, -1)[0]
    assert kh.certify_sign(kt, D, 0.2, +1)[0]


def test_common_mots_at_late_time():
    """tau = -0.009: the explicit surface COMMON_BARRIER (enclosing both centres: its axis
    radius 1.35 exceeds d = 0.5) is outer trapped, the sphere |x| = 3 is outer untrapped, so a
    MOTS enclosing both centres lies between them, far inside the cosmological horizon."""
    kt = kh.KT(M, D, -0.009)
    assert sum(float(c) for c in kh.COMMON_BARRIER) > D
    assert kh.certify_polar(kt, kh.COMMON_BARRIER, -1)[0]
    assert kh.certify_sign(kt, 0.0, 3.0, +1)[0]
    # the individual MOTSs persist inside it
    assert kh.certify_sign(kt, D, 0.001, -1)[0] and kh.certify_sign(kt, D, 0.1, +1)[0]


def test_fold_is_bracketed_numerically():
    """Shooting finds no symmetric common MOTS at tau = -0.015 and two (inner and outer) at
    tau = -0.0125: the fold tau_c = -0.012884 lies between."""
    assert kh.KT(M, D, -0.015).common_mots(1.2 * D, 6 * D, n=120) == []
    assert len(kh.KT(M, D, -0.0125).common_mots(1.2 * D, 6 * D, n=120)) == 2
