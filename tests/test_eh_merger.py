"""Background facts used in notes/11 (event horizons of extreme-mass-ratio mergers).

Kerr-de Sitter in Boyer-Lindquist form, coordinates (t, r, chi = cos theta, phi):

- it solves Ric = Lambda g, and reduces to Kerr (Lambda = 0) and Schwarzschild-de Sitter (a = 0);
- g^{rr} = Delta_r / rho^2, so r is a time function wherever Delta_r < 0 (black-hole interior and
  beyond the cosmological horizon);
- the spheres {t, r = const} have area density E / (Xi^2 Delta_chi) with dE/dr > 0 for every r > 0,
  so inside the black hole (where future null normals have dr < 0) they are trapped.
"""

import pytest
import sympy as sp

from einstein import metrics
from einstein.metrics import Lambda, M, a_kerr, chi, r
from einstein.tensors import Spacetime

ZERO4 = sp.zeros(4, 4)
RHO2 = r**2 + a_kerr**2 * chi**2
DELTA_R = (r**2 + a_kerr**2) * (1 - Lambda * r**2 / 3) - 2 * M * r
DELTA_CHI = 1 + Lambda * a_kerr**2 * chi**2 / 3
XI = 1 + Lambda * a_kerr**2 / 3
# Area density squared of {t = const, r = const} in the coordinates (chi, phi), times Xi^2 Delta_chi
E = DELTA_CHI * (r**2 + a_kerr**2) ** 2 - DELTA_R * a_kerr**2 * (1 - chi**2)


def test_kerr_de_sitter_reduces_to_kerr():
    kds = metrics.kerr_de_sitter().metric.subs(Lambda, 0)
    assert (kds - metrics.kerr().metric).applyfunc(sp.cancel) == ZERO4


def test_schwarzschild_de_sitter_solves_einstein_lambda():
    """a = 0, M and Lambda symbolic."""
    st = metrics.kerr_de_sitter()
    g = st.metric.subs(a_kerr, 0)
    ric = Spacetime(g, st.coords, simplify=sp.cancel).ricci
    assert (ric - Lambda * g).applyfunc(sp.cancel) == ZERO4


@pytest.mark.slow
def test_kerr_de_sitter_solves_einstein_lambda():
    """M = 1, a = 1/10, Lambda = 3/100: every component of R_ab - Lambda g_ab is a rational
    function of (r, chi) that cancel() reduces to 0 (~2.5 min)."""
    st = metrics.kerr_de_sitter()
    vals = {M: 1, a_kerr: sp.Rational(1, 10), Lambda: sp.Rational(3, 100)}
    g = st.metric.subs(vals)
    ric = Spacetime(g, st.coords, simplify=sp.cancel).ricci
    assert (ric - sp.Rational(3, 100) * g).applyfunc(sp.cancel) == ZERO4


def test_r_is_a_time_function_where_delta_r_is_negative():
    """The metric is block diagonal in (t, phi), r, chi, so g^{rr} = 1 / g_rr = Delta_r / rho^2."""
    g = metrics.kerr_de_sitter().metric
    assert all(g[1, j] == 0 for j in (0, 2, 3))
    assert sp.cancel(1 / g[1, 1] - DELTA_R / RHO2) == 0


def test_area_density_of_the_r_spheres():
    """{t, r = const} has induced metric diag(g_chichi, g_phiphi); its determinant is
    E / (Xi^2 Delta_chi). Delta_chi does not depend on r, so the first variation of area
    along a normal V is V^r d_r log sqrt(E)."""
    g = metrics.kerr_de_sitter().metric
    assert sp.cancel(g[2, 2] * g[3, 3] - E / (XI**2 * DELTA_CHI)) == 0
    assert sp.diff(DELTA_CHI, r) == 0


def test_spheres_inside_the_black_hole_are_trapped():
    """dE/dr = 4r^3 + 2a^2 r (1 + chi^2) + 2 M a^2 (1 - chi^2) + Lambda * (nonnegative terms):
    positive for every r > 0, every spin and every Lambda >= 0 (no slow-rotation assumption).
    With dr(l_pm) < 0 for both future null normals inside the black hole, theta_pm < 0."""
    dE = sp.expand(sp.diff(E, r))
    kerr_part = 4 * r**3 + 2 * a_kerr**2 * r * (1 + chi**2) + 2 * M * a_kerr**2 * (1 - chi**2)
    lambda_part = sp.expand((dE - kerr_part) / Lambda)
    expected = sp.Rational(4, 3) * r * a_kerr**2 * chi**2 * (r**2 + a_kerr**2) + a_kerr**2 * (
        1 - chi**2
    ) * (sp.Rational(4, 3) * r**3 + sp.Rational(2, 3) * a_kerr**2 * r)
    assert sp.expand(lambda_part - expected) == 0
    # every term is a product of nonnegative factors for r > 0, chi^2 <= 1
    assert sp.expand(dE.subs(Lambda, 0) - kerr_part) == 0


def test_spheres_beyond_the_cosmological_horizon_are_untrapped():
    """Same E: beyond r_C both future null normals have dr > 0 (expanding region), so both
    expansions are positive. Check the sign of Delta_r on either side of the horizons for one
    slowly rotating example (M = 1, a = 1/10, Lambda = 3/100)."""
    vals = {M: 1, a_kerr: sp.Rational(1, 10), Lambda: sp.Rational(3, 100)}
    roots = sorted(sp.Poly(DELTA_R.subs(vals) * 100, r).nroots(n=30), key=sp.re)
    real = [float(sp.re(x)) for x in roots if abs(sp.im(x)) < 1e-20]
    assert len(real) == 4  # r_neg < 0 < r_CH < r_H < r_C
    r_neg, r_ch, r_h, r_c = real
    assert r_neg < 0 < r_ch < r_h < r_c
    d = sp.lambdify(r, DELTA_R.subs(vals))
    assert d((r_ch + r_h) / 2) < 0 < d((r_h + r_c) / 2) and d(r_c + 1) < 0


def _kerr_inverse_tt_rr():
    """g^{tt} and g^{rr} of Kerr in Boyer-Lindquist form; g^{tr} = 0 and the (t, phi) block
    is inverted directly."""
    g = metrics.kerr().metric
    det_t_phi = g[0, 0] * g[3, 3] - g[0, 3] ** 2
    return sp.cancel(g[3, 3] / det_t_phi), sp.cancel(1 / g[1, 1])


def test_corrected_retarded_time_is_a_time_function_in_the_kerr_exterior():
    """notes/11, Lemma 5: u~ = t - r_* + a arctan(r / a), with dr_*/dr = (r^2 + a^2) / Delta,
    has rho^2 g^{-1}(du~, du~) = -a^2 r^2 / (r^2 + a^2) - (manifestly nonnegative) <= 0 where
    Delta > 0. Along outgoing principal null geodesics u = t - r_* is constant, so u~ changes
    by at most pi |a| / 2 there."""
    gtt, grr = _kerr_inverse_tt_rr()
    rho2 = r**2 + a_kerr**2 * chi**2
    delta = r**2 - 2 * M * r + a_kerr**2
    du_dr = -(r**2 + a_kerr**2) / delta + a_kerr**2 / (r**2 + a_kerr**2)
    norm = sp.cancel(rho2 * (gtt + du_dr**2 * grr))
    excess = a_kerr**2 * (
        2 * M * a_kerr**2 * r + a_kerr**4 * chi**2 + 2 * a_kerr**2 * chi**2 * r**2 + chi**2 * r**4
    ) / (r**2 + a_kerr**2) ** 2
    assert sp.cancel(norm - (-(a_kerr**2) * r**2 / (r**2 + a_kerr**2) - excess)) == 0


def test_outgoing_field_is_uniformly_timelike_far_out():
    """W = d_t + d_r / 2 has g(W, W) <= -1/4 for r >= 10 M and |a| <= M (Kerr), with bounded
    components: the steering field used at large r in notes/11, Lemma 2(b).
    4 rho^2 Delta (g(W, W) + 1/4) = -3 rho^2 Delta + 8 M r Delta + rho^4; with M = 1,
    r = 10 + x, a^2 = q, chi^2 = c, every coefficient of this polynomial in x is negative
    for q, c in [0, 1]."""
    x, q, c = sp.symbols("x q c", nonnegative=True)
    g = metrics.kerr().metric
    gWW = g[0, 0] + g[1, 1] / 4
    rho2 = r**2 + a_kerr**2 * chi**2
    delta = r**2 - 2 * M * r + a_kerr**2
    poly = sp.expand(-3 * rho2 * delta + 8 * M * r * delta + rho2**2)
    assert sp.cancel(4 * rho2 * delta * (gWW + sp.Rational(1, 4)) - poly) == 0
    poly = sp.expand(poly.subs({M: 1, a_kerr: sp.sqrt(q), chi: sp.sqrt(c)}).subs(r, 10 + x))
    coeffs = sp.Poly(poly, x).all_coeffs()
    assert coeffs[:2] == [-2, -66]
    # bound each coefficient on the box q, c in [0, 1]: positive monomials are at most 1 there
    for coeff in coeffs[2:]:
        terms = sp.Add.make_args(sp.expand(coeff))
        positive = sum(t for t in terms if t.as_coeff_Mul()[0] > 0 and t.free_symbols)
        constant = sum(t for t in terms if not t.free_symbols)
        assert all(t.as_coeff_Mul()[0] <= 1 for t in sp.Add.make_args(positive) if t != 0)
        assert constant + len(sp.Add.make_args(positive)) < 0


def test_outgoing_null_component_of_kerr_decays_like_r_minus_2():
    """notes/11, Lemma 6. In the null frame of the mass-M Schwarzschild model used at null
    infinity (Hintz 2026, (3.30)-(3.32)), x^1 = t - r with t = t_BL - r_*^{Schw} + r and
    d_0 = (1 - 2M/r)/2 (d_t + d_r): g(d_0, d_0) vanishes for Schwarzschild and is
    -a^2 (1 - chi^2) / (4 r^2) + O(r^-3) for Kerr. Hence the background's outgoing null
    component is integrable along outgoing curves."""
    g = metrics.kerr().metric
    V = sp.Matrix([r / (r - 2 * M), 1, 0, 0])  # d_t + d_r in the new coordinates, in BL components
    g00 = sp.cancel((1 - 2 * M / r) ** 2 / 4 * (V.T * g * V)[0, 0])
    assert sp.cancel(g00.subs(a_kerr, 0)) == 0
    w = sp.Symbol("w", positive=True)  # w = 1/r
    lead = sp.limit(sp.cancel(g00.subs(r, 1 / w) / w**2), w, 0)
    assert sp.simplify(lead + a_kerr**2 * (1 - chi**2) / 4) == 0
