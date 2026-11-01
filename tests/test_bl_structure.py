"""B-to-A structure of the BL development's time series (research/calc/bl_structure.py, notes/16 §5,
notes/17 §§5-7). Identities are checked in exact rational arithmetic at rational-distance points
(BL plane and axis, Schwarzschild); exclusions are certified by ball enclosures that exclude 0, with
Schwarzschild controls where the excluded property holds."""

from flint import arb, fmpq

from research.calc import bl_metric as bm
from research.calc import bl_structure as st

P = bm.PAIRS


def _exact_points():
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    za = fmpq(31, 50)
    return [(x0, d, dict(exact_r=rr)),
            ((fmpq(0), fmpq(0), za), d, dict(exact_r=[za - d / 2, za + d / 2])),
            ((fmpq(3, 5), fmpq(0), fmpq(4, 5)), None, dict(centres=[(fmpq(1, 2), (fmpq(0),) * 3, fmpq(1))]))]


def test_conformally_flat_identities_exactly():
    """For conformally flat time-symmetric vacuum data (R = 0, Cotton = 0):
    Delta Ric = 3 Ric.Ric - |Ric|^2 g, gamma_2 = (7 Ric.Ric - |Ric|^2 g)/12, and the t^6 coefficient
    gamma_3 = -(1/15)|Ric|^2 Ric - (41/540) tr(Ric^3) g + (1/60) d_iRic.d_jRic - (1/180)|dRic|^2 g
              - (1/180) Ric^kl d_k d_l Ric_ij,
    exactly at BL plane and axis points and at a Schwarzschild point."""
    for x0, d, kw in _exact_points():
        cv = st.curvature_values(x0, d, fmpq(1), **kw)
        r = st.cf_identities(cv, fmpq(1))
        for key in ("lap", "t4", "t6"):
            assert all(v == 0 for v in r[key].values()), key
        assert any(cv["g3"][pq] != 0 for pq in P)


def test_t6_formula_on_held_out_generic_points_and_algebraic_closure_excluded():
    """Held-out generic points (not used in the fit): the t^6 formula holds to the ball radius; and the
    6x3 matrix [|Ric|^2 Ric, tr(Ric^3) g, gamma_3] has a 3x3 minor excluding 0, so gamma_3 is not
    pointwise algebraic in Ric (the class gamma(t) = F(t, Ric) is excluded at t^6)."""
    with bm.working_precision(256):
        for x, d in (((arb("0.13"), arb("0.07"), arb("0.46")), arb("0.7662")),
                     ((arb("0.26"), arb("-0.09"), arb("0.12")), arb("0.767"))):
            cv = st.curvature_values(x, d, arb(1))
            r = st.cf_identities(cv, arb(1))
            assert all(abs(v) < 1e-60 for k in r for v in r[k].values())
            assert not st.algebraic_closure_minor(cv, arb(1)).contains(0)


def test_electric_weyl_t2_law_and_eigenframe_rotation_at_t4():
    """E^i_j(t) = E_0 + t^2 (-(3/2) E_0^2 + (1/2) tr(E_0^2) I) + O(t^4) exactly (BL plane point), so the
    tidal eigenframe is frozen at t^2; at t^4 it rotates: |[E_0, E_2]|^2 > 0 at generic BL points
    (ball excluding 0), = 0 for Schwarzschild. This excludes every closed pointwise isotropic
    evolution law E' = f(E) along the comoving observers."""
    z = fmpq(0)
    x0, d, kw = _exact_points()[0]
    E0, E1 = st.electric_series(x0, d, fmpq(1), order=1, **kw)
    E00 = st._mm(E0, E0, z)
    tr2 = E00[0][0] + E00[1][1] + E00[2][2]
    assert all(E1[i][j] == -fmpq(3, 2) * E00[i][j] + (tr2 / 2 if i == j else 0) for i in range(3) for j in range(3))
    with bm.working_precision(256):
        r = st.eigenframe_tests((arb("0.17"), arb("0.11"), arb("0.5")), arb("0.7662"), arb(1))
        assert r["comm_t2"].contains(0) and not r["comm_t4"].contains(0)
        cs = [(arb(1) / 2, (arb(0),) * 3, arb("1.21").sqrt())]
        rs = st.eigenframe_tests((arb("0.6"), arb("0.2"), arb("0.9")), None, arb(1), centres=cs)
        assert rs["comm_t4"].contains(0)


def test_weyl_speciality_is_not_frozen():
    """S = 27 J^2 / I^3 along the comoving observer: S = 1 and dS/d(t^2) = 0 for Schwarzschild (type D);
    for BL at a generic point dS/d(t^2) != 0 (ball excluding 0): the eigenvalue ratios of the tidal
    tensor drift, so no 'isospectral up to scale' tidal evolution."""
    with bm.working_precision(256):
        w = st.weyl_speciality((arb("0.26"), arb("-0.09"), arb("0.12")), arb("0.767"), arb(1))
        assert not w["S1"].contains(0)
        cs = [(arb(1) / 2, (arb(0),) * 3, arb("1.21").sqrt())]
        ws = st.weyl_speciality((arb("0.6"), arb("0.2"), arb("0.9")), None, arb(1), centres=cs)
        assert ws["S1"].contains(0) and abs(ws["S0"] - 1) < 1e-10


def test_single_hole_series_ends_at_the_curvature_singularity():
    """Single-hole control (notes/17 §7): the comoving congruence has no caustic before the crunch
    (H = F - (3/2) tau F' >= 1 > 0), and the Kretschmann scalar diverges at the crunch time
    pi sqrt(R0^3 / 8M), which is also the Taylor radius in t at that comoving point."""
    c = st.schwarzschild_crunch(R0=2.25, M=1.0)
    assert c["H_min_on_(0,pi)"] > 0.99
    assert c["kretschmann_near_crunch"] > 1e40 * c["kretschmann_initial"]
    assert abs(c["t_crunch"] - 3.7486824790711215) < 1e-12
