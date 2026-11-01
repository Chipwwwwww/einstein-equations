"""Experiment H1 (notes/18 §3): the silent (local tidal) part of the BL geodesic development and its first
nonlocal correction (research/calc/bl_silent.py). Identities are exact in rational arithmetic at
rational-distance points; held-out generic points use 256-bit balls."""

import math

from flint import arb, fmpq

from research.calc import bl_metric as bm
from research.calc import bl_silent as sl
from research.calc import bl_structure as st

R3 = range(3)


def _plane():
    d = fmpq(3831, 5000)
    x0, rr = st.plane_point(d, fmpq(15, 8))
    return x0, d, dict(exact_r=rr)


def test_silent_system_is_exact_through_t4_and_fails_at_t6():
    """gamma = gamma_0 M_silent(t) + O(t^6) exactly (BL plane and axis points), with the t^6 coefficient
    different; the time-odd silent coefficients vanish."""
    za = fmpq(31, 50)
    d = fmpq(3831, 5000)
    for x0, dd, kw in (_plane(), ((fmpq(0), fmpq(0), za), d, dict(exact_r=[za - d / 2, za + d / 2]))):
        true, sil, _, odd = sl.true_and_silent(x0, dd, 3, fmpq(1), **kw)
        agree = [all(true[k][i][j] == sil[k][i][j] for i in R3 for j in R3) for k in range(4)]
        assert agree == [True, True, True, False]
        assert all(v == 0 for M in odd for row in M for v in row)


def test_schwarzschild_is_silent_to_t8():
    """Control: for one hole the magnetic part vanishes and the silent ODE is the exact free fall."""
    cs = [(fmpq(1, 2), (fmpq(0),) * 3, fmpq(1))]
    true, sil, _, _ = sl.true_and_silent((fmpq(3, 5), fmpq(0), fmpq(4, 5)), None, 4, fmpq(1), centres=cs)
    assert all(true[k][i][j] == sil[k][i][j] for k in range(5) for i in R3 for j in R3)


def test_silent_t6_coefficients():
    """gamma_0 M_6 = -(13/180) tr(E^2) E - (41/540) tr(E^3) I for any traceless E_0 (commuting system);
    the true t^6 coefficient has -(1/15) = -(12/180) and the same -(41/540) (notes/17 §5)."""
    assert sl.silent_t6_coefficients() == (fmpq(-13, 180), fmpq(-41, 540))


def test_single_source_chain_and_covariant_law():
    """The first deviation is one traceless source (curl H)_3: G0 gamma_3 - M_6 = -(E_4 - E_4^silent)/15,
    and (curl H)_3 = -(1/3)|Ric|^2 Ric - nabla_i Ric.nabla_j Ric + (1/3)|nabla Ric|^2 g
    + (1/3) Ric^kl nabla_k nabla_l Ric_ij: exact at the plane point, to the ball radius at a held-out
    point with another separation."""
    x0, d, kw = _plane()
    dM6, dE4, _, _, _ = sl.remainder_chain(x0, d, fmpq(1), **kw)
    assert all(dM6[i][j] == -fmpq(1, 15) * dE4[i][j] for i in R3 for j in R3)
    assert dE4[0][0] + dE4[1][1] + dE4[2][2] == 0 and any(v != 0 for r in dE4 for v in r)
    lhs, rhs = sl.curl_h3_law(x0, d, fmpq(1), **kw)
    assert all(lhs[p] == rhs[p] for p in lhs)
    with bm.working_precision(256):
        lhs, rhs = sl.curl_h3_law((arb("0.31"), arb("0.2"), arb("-0.7")), arb("0.7"), arb(1))
        assert all(abs(lhs[p] - rhs[p]) < 1e-60 for p in lhs) and max(abs(lhs[p]) for p in lhs) > 0.01


def test_silent_crunch_type_d_control():
    """For E_0 = mu diag(-2, 1, 1) the silent ODE is Schwarzschild free fall: collapse at pi / sqrt(8 mu)."""
    for mu in (0.29, 0.05):
        tc = sl.silent_crunch([-2 * mu, mu, mu])
        assert abs(tc - math.pi / math.sqrt(8 * mu)) < 1e-6
