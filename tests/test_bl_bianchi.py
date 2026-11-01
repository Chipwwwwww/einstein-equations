"""Session 15 (notes/19): the Einstein-Bianchi (E, H) form of the geodesic-slicing development, the
curl-degree grading, the universal t^4/t^6 laws, the principal symbols, and the energy identity."""

from flint import fmpq

from research.calc import bianchi_symbol as bs
from research.calc import bl_bianchi as bb
from research.calc import bl_metric as bm
from research.calc import bl_silent as sl
from research.calc import bl_structure as st
from research.calc import cf_universal as cu


def _plane():
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    return x0, d, rr


def _gamma(variant="EvE", eps=1, N=4):
    x0, d, rr = _plane()
    J, s = bb.bianchi_jets(bb.bl_psi(x0, d, fmpq(1), exact_r=rr), 2 * N, fmpq(1), eps=fmpq(eps), variant=variant)
    return bb.point_coefficients(J, s, "g", 2 * N, fmpq(0)), (J, s)


def test_bianchi_system_reproduces_adm_coefficients_exactly():
    """The E-H system (Ellis-van Elst 40-41 / Friedrich-Rendall P, Q converted to (-,+,+,+)) gives the same
    gamma as the ADM recursion through t^8 at a rational BL point; the unconverted transcription of
    Friedrich-Rendall's Q_kl (signature +---) fails first at t^8, so the comparison is sensitive."""
    x0, d, rr = _plane()
    JA, gA, _ = bm.series_jets(x0, d, 4, fmpq(1), exact_r=rr)
    ref = [{pq: JA.coef(gA[pq], (0, 0, 0, k), fmpq(0)) for pq in bm.PAIRS} for k in range(5)]
    g, _ = _gamma()
    assert all(g[2 * k][pq] == ref[k][pq] for k in range(5) for pq in bm.PAIRS)
    assert all(g[2 * k + 1][pq] == 0 for k in range(4) for pq in bm.PAIRS)
    gw, _ = _gamma("FR-printed")
    assert [all(gw[2 * k][pq] == ref[k][pq] for pq in bm.PAIRS) for k in range(5)] == [True] * 4 + [False]


def test_eps_zero_is_silent_and_lag_lemma():
    """eps = 0 is the silent ODE exactly; the eps^2 sector of gamma starts at t^6 and the eps^4 sector at t^8."""
    x0, d, rr = _plane()
    true, sil, _, _ = sl.true_and_silent(x0, d, 4, fmpq(1), exact_r=rr)
    g = {e: _gamma(eps=e)[0] for e in (0, 1, 2)}
    assert all(g[0][2 * k][pq] == sil[k][pq[0]][pq[1]] for k in range(5) for pq in bm.PAIRS)
    pq = (0, 0)
    for k in range(5):
        f0, f1, f2 = (g[e][2 * k][pq] for e in (0, 1, 2))
        c2 = (f2 - 4 * f1 + 3 * f0) / 12          # eps^4 coefficient (degree <= 2 in eps^2 through t^8)
        c1 = f1 - f0 - c2
        assert (c1 != 0) == (k >= 3) and (c2 != 0) == (k >= 4)


def test_t4_t6_laws_hold_for_all_conformally_flat_vacuum_data():
    """The finite exact verification of notes/19 §3 (72 evaluations on harmonic jets), plus controls: a
    non-harmonic phi breaks the identities."""
    for label, weights, terms in cu.evaluation_set():
        r, _ = cu.residuals(terms)
        keys = (["lap", "t4"] if 4 in weights else []) + (["t6"] if 6 in weights else [])
        for key in keys:
            assert all(v == 0 for v in r[key].values()), (label, key)
    r, _ = cu.residuals([((4, 0, 0), fmpq(1))])
    assert not all(v == 0 for v in r["lap"].values())
    r, _ = cu.residuals([((3, 0, 0), fmpq(1))])
    assert not all(v == 0 for v in r["t6"].values())


def test_principal_symbols():
    """Frame form: C_xi antisymmetric (symmetric hyperbolic), speeds 0, +-1/2, +-1. Coordinate form
    (Christoffels from gamma): at the sampled state and direction there is a Jordan block at eigenvalue 0,
    so not strongly hyperbolic there (exact example in research/calc/curl_hierarchy_checks.py; no claim
    for every state or direction)."""
    out = bs.main()
    assert out["C_antisymmetric"] < 1e-12
    assert sorted(abs(v) for v in out["frame_speeds"]) == [0, 0, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1]
    alg, geo, _, _ = out["coordinate_H0"]
    assert alg == 16 and geo == 15
    alg0, geo0, _, _ = out["coordinate_flat"]
    assert alg0 == geo0 == 16


def test_energy_identity_exact():
    """d_t(sqrt g w) - 2 d_k(sqrt g F^k) = sqrt g (-theta w + 6 tr(sigma(E^2 + H^2))), exactly, t^0..t^3."""
    x0, d, rr = _plane()
    res = bb.energy_identity_residual(bb.bl_psi(x0, d, fmpq(1), exact_r=rr), 6, fmpq(1), 3)
    assert all(v == 0 for v in res)
