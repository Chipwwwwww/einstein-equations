"""The evaluable metric of the near-critical BL development (research/calc/bl_metric.py, notes/16).
Independent checks: exact Schwarzschild free fall in every Cartesian component, the independent
cylindrical floating-point code, the first-order system of the certificate (with its own Ricci
reduction), the constraints, symmetries, overlapping enclosures (precision, order, patches) and
refusals outside the certified domain."""

import math

import pytest
from flint import arb, ctx

from research.calc import bl_metric as bm
from research.calc import bl_series as bs

P = bm.PAIRS


def _free_fall_series(n):
    """Taylor coefficients a_m of F(tau): F'' = -1/F^2, F(0) = 1, F'(0) = 0 (radial free fall from rest,
    R = R0 F(t sqrt(M / R0^3)))."""
    a = [arb(0)] * (2 * n + 3)
    a[0] = arb(1)

    def c2(m):
        return sum((a[i] * a[m - i] for i in range(m + 1)), arb(0))
    for m in range(2 * n + 1):
        known = sum((c2(i) * (m - i + 2) * (m - i + 1) * a[m - i + 2] for i in range(1, m + 1)), arb(0))
        a[m + 2] = ((arb(-1) if m == 0 else arb(0)) - known) / (c2(0) * (m + 2) * (m + 1))
    return a


def test_schwarzschild_free_fall_in_every_component():
    """One centre of mass M = 1. Geodesic slicing from rest is exactly
    gamma_ij = psi^4 [F^2 (delta_ij - n_i n_j) + (F - 3/2 tau F')^2 n_i n_j], tau = t sqrt(M / R0^3),
    R0 = r psi^2 (Novikov-type comoving coordinates). The jet engine must reproduce its coefficients
    through t^6 at a generic point, in all six Cartesian components."""
    with bm.working_precision(128):
        N = 3
        a = _free_fall_series(N + 1)
        L = 2 * N + 3
        F = a[:L]
        H = [F[m] - arb(3) / 2 * (m * F[m]) for m in range(L)]
        sq = lambda c: [sum((c[i] * c[m - i] for i in range(m + 1)), arb(0)) for m in range(L)]  # noqa: E731
        F2, H2 = sq(F), sq(H)
        x0 = (arb("0.6"), arb("0.3"), arb("0.9"))
        r = (x0[0] ** 2 + x0[1] ** 2 + x0[2] ** 2).sqrt()
        psi = 1 + 1 / (2 * r)
        R0 = r * psi ** 2
        n = [v / r for v in x0]
        J, g, _ = bm.series_jets(x0, None, N, arb(1), centres=[(arb(1) / 2, (arb(0),) * 3, r)])
        vals, _, _ = bm.coefficients(J, g, arb(0), N)
        for k in range(N + 1):
            fac = psi ** 4 / R0 ** (3 * k)
            for (i, j) in P:
                ex = fac * (F2[2 * k] * ((1 if i == j else 0) - n[i] * n[j]) + H2[2 * k] * n[i] * n[j])
                assert abs(vals[k][(i, j)] - ex) < 1e-25


def test_cartesian_engine_matches_the_cylindrical_float_code():
    """bl_series (floating point, cylindrical (rho, z, phi) components, independent code) and the
    Cartesian ball-arithmetic engine agree at an off-axis BL point through t^8."""
    d, rho0, z0, phi = 0.7662, 0.26, 0.1, 0.7
    with bm.working_precision(128):
        x0 = (arb(rho0) * arb(phi).cos(), arb(rho0) * arb(phi).sin(), arb(z0))
        J, g, _ = bm.series_jets(x0, arb(d), 4, arb(1))
        vals, _, _ = bm.coefficients(J, g, arb(0), 4)
    vf, _, _ = bs.series(rho0, z0, [(0.5, d / 2), (0.5, -d / 2)], 8)
    c, s = math.cos(phi), math.sin(phi)
    er, ep, ez = (c, s, 0), (-s, c, 0), (0, 0, 1)

    def con(k, u, v):
        return sum(float(bm.sym_get(vals[k], i, j).mid()) * u[i] * v[j] for i in range(3) for j in range(3))
    for k in range(5):
        for name, val in (("rr", con(k, er, er)), ("rz", con(k, er, ez)), ("zz", con(k, ez, ez)),
                          ("pp", rho0 ** 2 * con(k, ep, ep))):
            assert abs(val - vf[name][k]) <= 1e-9 * max(1.0, abs(vf[name][k]))
        assert abs(con(k, er, ep)) < 1e-15 and abs(con(k, ez, ep)) < 1e-15   # axisymmetry


def test_first_order_system_of_the_certificate_has_the_same_coefficients():
    """The Taylor coefficients of the first-order system the CK certificate is proved for (built from
    ck_bound.ricci_parts and quad_K, no parity assumed) equal those of the second-order engine; the
    odd ones vanish (time reflection)."""
    with bm.working_precision(128):
        x0, d = (arb("0.17"), arb("0.11"), arb("0.5")), arb("0.7662")
        gv, kv = bm.first_order_coefficients(x0, d, 5, arb(1))
        J, g, kap = bm.series_jets(x0, d, 2, arb(1))
        vals, _, _ = bm.coefficients(J, g, arb(0), 2)
        for n in range(5):
            for pq in P:
                ref = vals[n // 2][pq] if n % 2 == 0 else arb(0)
                assert abs(gv[n][pq] - ref) < 1e-28
                # K_n = (n + 1)/2 gamma_(n+1)
                if n + 1 < 5:
                    assert abs(kv[n][pq] - (n + 1) * gv[n + 1][pq] / 2) < 1e-28


def test_constraints_vanish_order_by_order():
    with bm.working_precision(128):
        J, g, kap = bm.series_jets((arb("0.2"), arb("-0.05"), arb("0.45")), arb("0.767"), 3, arb(1))
        ham, mom = bm.constraint_residuals(J, g, kap, arb(1))
    assert max(float(abs(h).upper()) for h in ham) < 1e-25
    assert max(float(abs(v).upper()) for m in mom for v in m) < 1e-25


def _tube():
    """The whole-tube certificate with the session-12 constants (times are fractions of its T)."""
    return bm.load_certificate("ck_bound_sharp.json")


def test_axis_is_regular_and_the_metric_is_axisymmetric():
    """On the axis the Cartesian chart is regular: gamma is finite, diagonal and gamma_xx = gamma_yy.
    Off the axis, two points related by a rotation have rotated metrics (equivariance)."""
    cert = _tube()
    t = arb("0.2") * arb(cert["T_min"])
    r = bm.evaluate_metric(t, "0", "0", "0.62", "0.7662", order=2)
    g = r["g"]
    assert r["status"] == "certified" and r["lorentzian"]
    assert g[1][2].contains(0) and g[1][3].contains(0) and g[2][3].contains(0)
    assert (g[1][1] - g[2][2]).contains(0)
    phi = arb("0.9")
    rho = arb("0.2")
    a = bm.evaluate_metric(t, rho, "0", "0.5", "0.7662", order=2)["g"]
    b = bm.evaluate_metric(t, rho * phi.cos(), rho * phi.sin(), "0.5", "0.7662", order=2)["g"]
    c, s = phi.cos(), phi.sin()
    Rm = [[c, -s, arb(0)], [s, c, arb(0)], [arb(0), arb(0), arb(1)]]
    for i in range(3):
        for j in range(3):
            rot = sum((Rm[i][k] * Rm[j][l] * a[k + 1][l + 1] for k in range(3) for l in range(3)), arb(0))
            assert rot.overlaps(b[i + 1][j + 1])


def test_refuses_outside_the_certified_domain():
    T = arb(max(c["T_min"] for c in bm.load_certificates()))
    with pytest.raises(bm.DomainError):
        bm.evaluate_metric(T * arb("1.01"), "0.26", "0", "0", "0.7662", order=2)    # |t| >= T_min
    with pytest.raises(bm.DomainError):
        bm.evaluate_metric(T / 4, "0.1", "0", "0", "0.7662", order=2)              # far from S_c
    with pytest.raises(bm.DomainError):
        bm.evaluate_metric(T / 4, "0.26", "0", "0", "0.78", order=2)               # d not in I
    u = bm.evaluate_metric(T * 3, "0.26", "0", "0", "0.7662", order=2, allow_uncertified=True)
    assert u["status"] == "uncertified" and u["error_budget"]["time_truncation"] is None


def test_unreachable_tolerance_fails_informatively():
    cert = _tube()
    with pytest.raises(bm.DomainError, match="not reached"):
        bm.evaluate_metric(arb(cert["T_min"]) / 2, "0.26", "0", "0", "0.7662", tolerance=1e-40,
                           max_order=2)


def test_overlapping_enclosures_precision_order_patch_and_time_reflection():
    """Certified enclosures computed differently must overlap: 128 vs 192 bits, order 2 vs 3, a
    spatial patch vs points inside it; and gamma is even, d_t gamma odd in t."""
    cert = _tube()
    t = arb("0.3") * arb(cert["T_min"])
    x, y, z, d = "0.17", "0.11", "0.5", "0.7662"
    r1 = bm.evaluate_metric(t, x, y, z, d, order=2, precision=128)
    r2 = bm.evaluate_metric(t, x, y, z, d, order=3, precision=192)
    rm = bm.evaluate_metric(-t, x, y, z, d, order=2)
    patch = bm.evaluate_metric(arb(0, float(t.upper())), arb(x, "0.0005"), arb(y, "0.0005"),
                               arb(z, "0.0005"), d, order=2)
    q = bm.evaluate_metric(t, arb(x) + arb("0.0003"), arb(y) - arb("0.0004"), z, d, order=2)
    for i in range(4):
        for j in range(4):
            assert r1["g"][i][j].overlaps(r2["g"][i][j])
            assert r1["g"][i][j].overlaps(rm["g"][i][j])
            assert (r1["dg"][0][i][j] + rm["dg"][0][i][j]).contains(0)
            assert patch["g"][i][j].contains(r1["g"][i][j]) and patch["g"][i][j].contains(q["g"][i][j])
            for l in range(4):
                assert r1["dg"][l][i][j].overlaps(r2["dg"][l][i][j])
    assert patch["lorentzian"] and r1["error_budget"]["spatial_approximation"] == 0.0


def test_lorentzian_on_the_whole_certified_domain():
    """|gamma(t) - gamma_0|_F <= R_gamma B on U_c for every certified pair (a, B) of every certificate,
    and gamma_0 = psi^4 delta >= psi_min^4: so gamma(t) is positive definite (g Lorentzian) wherever
    the evaluator certifies."""
    for cert in bm.load_certificates():
        p = arb(cert["psi_abs_lower"])
        for a, B in bm.certified_pairs(cert):
            assert p * p * p * p - arb(cert["radii"][0]) * arb(B) > 15
    assert ctx.prec >= 53


def test_universal_t4_coefficient_exactly_at_rational_points():
    """gamma_2 = (5/6) Ric.Ric - (1/6) |Ric|^2 gamma_0 - (1/12) Delta Ric (notes/16 §4), as an exact
    rational identity at a rational-distance point of BL data and at a Schwarzschild point."""
    from flint import fmpq
    from research.calc import bl_structure as st
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    g2, cand, _ = st.t4_identity(x0, d, fmpq(1), exact_r=rr)
    assert all(g2[pq] == cand[pq] for pq in P)
    g2s, cs, _ = st.t4_identity((fmpq(3, 5), fmpq(0), fmpq(4, 5)), None, fmpq(1),
                                centres=[(fmpq(1, 2), (fmpq(0),) * 3, fmpq(1))])
    assert all(g2s[pq] == cs[pq] for pq in P)
    assert any(g2[pq] != 0 for pq in P)


def test_conformally_flat_and_silent_ansatze_are_excluded():
    """C1: the t^2 coefficient of the Cotton tensor of gamma(t) is nonzero for BL (zero at t^0 and
    for Schwarzschild). C2: the t^3 coefficient of the magnetic Weyl tensor is nonzero for BL (zero
    for Schwarzschild). Nonzero = a ball excluding 0."""
    from research.calc import bl_structure as st
    with bm.working_precision(128):
        xg, d = (arb("0.17"), arb("0.11"), arb("0.5")), arb("0.7662")
        C = st.cotton_s1(xg, d, arb(1))
        assert all(v[0].contains(0) for v in C.values())
        assert any(not v[1].contains(0) for v in C.values())
        B = st.magnetic_t3(xg, d, arb(1))
        assert all(v[0].contains(0) for v in B.values())
        assert any(not v[1].contains(0) for v in B.values())
        cs = [(arb(1) / 2, (arb(0),) * 3, arb("1.21").sqrt())]
        xs = (arb("0.6"), arb("0.2"), arb("0.9"))
        assert all(x.contains(0) for v in st.cotton_s1(xs, None, arb(1), centres=cs).values() for x in v)
        assert all(x.contains(0) for v in st.magnetic_t3(xs, None, arb(1), centres=cs).values() for x in v)


def test_certificates_combine_and_respect_their_regions():
    """A waist point uses the waist certificate and the tube certificates; a point near the top of S_c
    only the tube ones; the waist certificate alone refuses the top point; tails from the combination
    are no larger than from any single certificate; the corrected Baouendi-Goulaouic time is below
    the CK' time for the same constants."""
    certs = {c["file"]: c for c in bm.load_certificates()}
    waist, tube = certs["ck_bound_waist.json"], certs["ck_bound_sharp.json"]
    t = arb("0.3") * arb(tube["T_min"])
    r = bm.evaluate_metric(t, "0.26", "0", "0.05", "0.7662", order=3)
    assert "ck_bound_waist.json" in r["error_budget"]["certificate"] and len(r["error_budget"]["certificate"]) >= 2
    r_w = bm.evaluate_metric(t, "0.26", "0", "0.05", "0.7662", order=3, certificate=waist)
    r_t = bm.evaluate_metric(t, "0.26", "0", "0.05", "0.7662", order=3, certificate=tube)
    for k, v in r["error_budget"]["time_truncation"].items():
        assert v <= r_w["error_budget"]["time_truncation"][k] and v <= r_t["error_budget"]["time_truncation"][k]
    top = bm.evaluate_metric(t, "0", "0", "0.62", "0.7662", order=2)
    assert "ck_bound_waist.json" not in top["error_budget"]["certificate"]
    with pytest.raises(bm.DomainError):
        bm.evaluate_metric(t, "0", "0", "0.62", "0.7662", order=2, certificate=waist)
    # beyond the tube's time but inside the waist's: certified only on the waist region
    tw = arb("0.9") * arb(waist["T_min"])
    assert bm.evaluate_metric(tw, "0.26", "0", "0.05", "0.7662", order=3)["status"] == "certified"
    with pytest.raises(bm.DomainError):
        bm.evaluate_metric(tw, "0.17", "0.11", "0.5", "0.7662", order=3)
    assert tube["T_bg"] < tube["T_min"]
