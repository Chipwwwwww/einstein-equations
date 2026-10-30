"""Gap 3 over an interval of separations (research/calc/gap3_interval.py): consistency checks of
the machinery. The certificate itself is `python -m research.calc.gap3_interval 4` (hours)."""

import math

from flint import arb

from research.calc import gap3_interval as gi
from research.calc import gap3_sweep as gs
from research.calc.gap3_neck import check_piece_rel, D_HI, D_LO, LAM0


def test_relative_height_system_matches_the_original():
    """f_rel at (rho, z - d/2, u, v, d) equals f_bl at (rho, z, u, v) with separation d."""
    for rho, z, u, v, d in [(0.1, 0.3, 0.05, -0.08, 0.77), (0.4, -0.1, 0.3, 0.2, 0.7801)]:
        a = gi.f_rel([arb(rho), arb(z) - arb(d) / 2, arb(u), arb(v), arb(d)])
        b = gs.f_bl([arb(rho), arb(z), arb(u), arb(v)], arb(d))
        for x, y in zip(a[:4], b):
            assert x.overlaps(y) and abs(float(x.mid() - y.mid())) < 1e-25


def test_start_set_contains_floating_point_starts():
    """The band-aligned start set contains the states computed at sample (delta, d)."""
    from research.calc.bl_fold import BL
    box = (0.001, 0.0011, 0.7795, 0.78)
    S = gi.start_pset(*box)
    X = S.hull()
    for delta, d in [(0.001, 0.7795), (0.00105, 0.7797), (0.0011, 0.78)]:
        z0 = gi.Z0_REF + delta + gi.SLOPE * (d - gi.D_REF)
        k0 = 2 * BL(d).n_dlogpsi(0.0, z0, 0.0)
        s1 = 2.0**-26
        zt = z0 - d / 2 + k0 * s1**2 / 2
        assert X[1].lower() - 1e-14 <= zt <= X[1].upper() + 1e-14
        assert abs(float(X[0].mid()) - s1) < 1e-15


def test_small_neck_certificate_over_the_interval():
    """One corner piece of certify_region_rel: lam near LAM0, zeta near the neck, all d."""
    piece = (LAM0 - 2.0**-11, LAM0, -0.1796875, -0.1796875 + 2.0**-9, D_LO, D_HI, -1)
    assert check_piece_rel(piece) > 0.5


def test_medium_neck_certificate_over_the_interval():
    """certify_rel on a grid cell at lam ~ 0.03 for all d in [D_LO, D_HI]."""
    w = 2.0**-12
    cell = (0.03, 0.03 + w, -0.2535, -0.2535 + w, -1)
    assert gi._certify_cell(cell) > 0


def test_interval_boxes_end_in_N_and_W():
    """A box with a tiny neck ends in N; one with a medium neck reaches a waist (W) whose box
    contains the floating-point waist."""
    ev = gi.run_box((0.0001, 0.00016, 0.7802734375, 0.78125))[1]
    assert ev == "N"
    r = gi.run_box((0.0029296875, 0.002990722656, 0.7802734375, 0.78125))
    assert r[1] == "W"
    lam_lo, lam_hi, ze_lo, ze_hi, s = r[5]
    # floating point: waists between (0.01355, -0.22238) and (0.01388, -0.22322) over the box
    assert s == -1 and lam_lo < 0.01355 and 0.01388 < lam_hi
    assert ze_lo < -0.22322 and -0.22238 < ze_hi


def test_band_cover_check_is_dyadic_and_recursive():
    """covered_by: a waist box is covered by a certified GRID cell, or by certified quarters
    of it that meet the box; a missing quarter or the other direction s is not covered."""
    G = gi.GRID
    whole = [(0.0, G, -G, 0.0, -1)]
    assert gi.covered_by((1e-4, 5e-4, -5e-4, -1e-4, -1), whole)
    assert not gi.covered_by((1e-4, 5e-4, -5e-4, -1e-4, 1), whole)
    assert not gi.covered_by((1e-4, 1.2e-3, -5e-4, -1e-4, -1), whole)
    h = G / 2
    quarters = [(0.0, h, -G, -h, -1), (h, G, -G, -h, -1), (0.0, h, -h, 0.0, -1)]
    assert gi.covered_by((1e-4, 2e-4, -2e-4, -1e-4, -1), quarters)
    assert not gi.covered_by((6e-4, 8e-4, -2e-4, -1e-4, -1), quarters)


def test_slab_leaves_out_exactly_the_fold_window():
    """The slab boxes cover [DELTA_MIN, DELTA_MAX] minus HOLE, the hole's z0-image lies inside
    the fold lemma's window, and the z0-window [d/2 + s*, R*] is covered for d in [D0, D1]."""
    from research.calc import gap3_fold as gf
    boxes = gi.initial_boxes(gi.D0, gi.D1, gi.HOLE)
    edges = sorted({(b[0], b[1]) for b in boxes})
    x = gi.DELTA_MIN
    for lo, hi in edges:
        if lo != x:
            assert (x, lo) == gi.HOLE
        x = hi
    assert x == gi.DELTA_MAX
    assert all(b[2] == gi.D0 and b[3] == gi.D1 for b in boxes)
    lo, hi = gi.hole_z0()
    a, b = gf.fold_window()
    assert a <= lo and hi <= b and a <= gf.Z_C_LO and gf.Z_C_HI <= b
    assert gi.z0_window_covered(gi.D0, gi.D1)
    assert gi.D0 <= gf.D_C_LO and gf.D_C_HI < gi.D1


def test_falsification_below_d_c():
    """At d = 0.76595 < d_c two common minimal surfaces exist, at z0 = 0.630956 and 0.634135
    (floating point; delta = 0.048929 and 0.052108 in the sweep's coordinates). Boxes
    containing them must not be excluded; a box between them, where F ~ -3e-4, is (E1)."""
    d = (0.76594, 0.76596)
    for delta in (0.048929, 0.052108):
        r = gi.run_box((delta - 2.0**-16, delta + 2.0**-16) + d)
        assert r[1] not in gi.EXCLUDED, r
    r = gi.run_box((0.0504, 0.0504 + 2.0**-16) + d)
    assert r[1] == "E1", r
