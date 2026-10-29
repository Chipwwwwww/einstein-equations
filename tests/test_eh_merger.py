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
