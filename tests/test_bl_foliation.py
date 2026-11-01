"""Experiment H2 (notes/18 §4): no slicing of the BL development whose initial lapse is bounded and
D^2-harmonic (N = (c_0 + c_+/r_+ + c_-/r_-)/psi, which contains the geodesic N = 1) has leaves that are
conformally flat at O(t^2) (research/calc/bl_foliation.py)."""

from flint import arb, fmpq

from research.calc import bl_foliation as fo
from research.calc import bl_metric as bm
from research.calc import bl_structure as st


def test_linearised_cotton_controls():
    """Pure-gauge perturbations L_xi gamma_0 and phi gamma_0 give DC = 0; a non-gauge one does not; the two
    forms of h agree; the geodesic member reproduces the session-12 Cotton coefficient (independent
    engine, notes/16 §5)."""
    with bm.working_precision(256):
        x, d = (arb("0.17"), arb("0.11"), arb("0.5")), arb("0.7662")
        g = fo.gauge_controls(x, d, arb(1))
        assert g["lie"] < 1e-60 and g["conformal"] < 1e-60 and g["non_gauge_example"] > 1e-3
        Q, dev, c0 = fo.quadratic_forms(x, d, arb(1))
        assert dev < 1e-60 and abs(c0) < 1e-60
        C = st.cotton_s1(x, d, arb(1))
        geo = fo.geodesic_member(Q, arb(1))
        assert abs(geo[(0, 0, 1)] - C[(0, 0, 1)][1]) < 1e-60 and abs(geo[(0, 1, 2)] - C[(0, 1, 2)][1]) < 1e-60
    x0, rr = st.plane_point(fmpq(3831, 5000), fmpq(15, 8))
    Q, dev, c0 = fo.quadratic_forms(x0, fmpq(3831, 5000), fmpq(1), exact_r=rr)
    assert dev == 0.0 and c0 == 0


def test_two_points_exclude_the_family_one_point_cannot():
    """Axisymmetry leaves two independent conics per point, which have real common zeros (a pointwise
    'conformally flat' lapse exists); the conics at two points have none (box covering of RP^2)."""
    with bm.working_precision(256):
        Qs = []
        for xs in (("0.17", "0.11", "0.5"), ("0.6", "0.3", "0.9")):
            Q, _, _ = fo.quadratic_forms(tuple(arb(v) for v in xs), arb("0.7662"), arb(1), check=False)
            assert fo.component_rank(Q)[0] == 2
            assert not fo.certify_no_common_zero(Q, max_depth=12)[0]
            assert fo.pointwise_zeros(Q, tries=60)[0][0] < 1e-14
            Qs.append(Q)
        ok, boxes, depth, fails = fo.certify_no_common_zero(Qs)
        assert ok and not fails


def test_schwarzschild_control():
    """One hole: every lapse (c_0 + c_1/r)/psi is spherical, so DC = 0; an off-centre 1/|x - q| is not."""
    with bm.working_precision(256):
        xs = (arb("0.6"), arb("0.2"), arb("0.9"))
        r0 = (arb("0.36") + arb("0.04") + arb("0.81")).sqrt()
        cs = [(arb(1) / 2, (arb(0),) * 3, r0)]
        Q, dev, _ = fo.quadratic_forms(xs, None, arb(1), centres=cs)
        assert max(abs(v) for F in Q.values() for v in F.values()) < 1e-60
        q = (arb(0), arb(0), arb("0.3"))
        rq = (arb("0.36") + arb("0.04") + arb("0.36")).sqrt()
        Qo, _, _ = fo.quadratic_forms(xs, None, arb(1), centres=cs, lapse_centres=[((arb(0),) * 3, r0), (q, rq)])
        assert max(abs(v) for F in Qo.values() for v in F.values()) > 1e-3
