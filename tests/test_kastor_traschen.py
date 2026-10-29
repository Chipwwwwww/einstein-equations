"""Kastor-Traschen benchmark (notes/08): local verification of every claim the note makes.

Metric and potential (Brill-Horowitz-Kastor-Traschen 1994 form, G = c = 1, Gaussian units):

    ds^2 = -U^-2 dtau^2 + U^2 (dx^2 + dy^2 + dz^2),   A = U^-1 dtau,
    U = H tau + sum_i m_i / |x - x_i|,                H^2 = Lambda / 3.
"""

import sympy as sp

from einstein import Spacetime, metrics
from einstein.fields import field_strength, maxwell_divergence, maxwell_stress_energy
from einstein.metrics import H, tau, x, y, z

Lam = sp.Symbol("Lambda")
X = (x, y, z)


def residuals(U, simplify=sp.simplify):
    """E_ab = G_ab + Lambda g_ab - 8 pi T_ab and nabla_a F^ab for the ansatz with this U."""
    st = metrics.conformastatic(U, simplify=simplify)
    F = field_strength(metrics.conformastatic_potential(U), st.coords)
    E = (st.einstein + Lam * st.metric - 8 * sp.pi * maxwell_stress_energy(st, F)).applyfunc(simplify)
    return st, E, [simplify(e) for e in maxwell_divergence(st, F)]


def test_reduced_system_for_arbitrary_U():
    """For ANY U(tau, x, y, z) the field equations reduce to three scalar expressions.
    This is the identity the classification in notes/08 rests on."""
    U = sp.Function("U")(tau, *X)
    _, E, div = residuals(U)
    Ut = U.diff(tau)
    lap = sum(U.diff(v, 2) for v in X)
    expected = sp.zeros(4, 4)
    expected[0, 0] = -(Lam * U**3 - 3 * U**3 * Ut**2 + 2 * lap) / U**5
    for i, v in enumerate(X, start=1):
        expected[0, i] = expected[i, 0] = -2 * U.diff(tau, v) / U
        expected[i, i] = -U**2 * (2 * U * U.diff(tau, 2) + 3 * Ut**2 - Lam)
    assert (E - expected).applyfunc(sp.simplify) == sp.zeros(4, 4)
    assert [sp.simplify(d - e) for d, e in zip(div, [lap / U**2] + [-U.diff(tau, v) / U**2 for v in X])] == [0] * 4


def test_kastor_traschen_two_centres_reduced():
    """Proof for any number of centres, via the identity above: U = H tau + V gives
    U_tau = H, U_tau tau = 0, U_tau x = 0, so E_ab = 0 <=> Lambda = 3 H^2 and lap V = 0,
    and the Maxwell equations <=> lap V = 0. Here lap V = 0 for two centres off the centres."""
    m1, m2, d = sp.symbols("m1 m2 d", positive=True)
    U = metrics.kastor_traschen_U([m1, m2], [(0, 0, d), (0, 0, -d)])
    assert sp.simplify(sum(U.diff(v, 2) for v in X)) == 0
    assert U.diff(tau) == H and U.diff(tau, 2) == 0 and all(U.diff(tau, v) == 0 for v in X)


def test_kastor_traschen_two_centres_direct_exact_points():
    """Independent of the reduction: the full unsimplified residual, evaluated in exact
    rational arithmetic at points where both distances are rational, is exactly zero.
    (A pointwise check, so a diagnostic; the proof is the reduced-system test.)"""
    m1, m2, d = sp.symbols("m1 m2 d", positive=True)
    U = metrics.kastor_traschen_U([m1, m2], [(0, 0, d), (0, 0, -d)])
    _, E, div = residuals(U, simplify=lambda e: e)
    # (x, y, z, d) with r1 = |x - (0,0,d)| and r2 = |x + (0,0,d)| integers: 13 & 15, 25 & 29, 15 & 37
    points = [(12, 0, 7, 2), (0, 20, 18, 3), (0, 12, 22, 13)]
    for px, py, pz, pd in points:
        values = {x: px, y: py, z: pz, d: pd, tau: sp.Rational(-3, 7), H: sp.Rational(-2, 5),
                  m1: sp.Rational(1, 3), m2: sp.Rational(1, 11), Lam: 3 * sp.Rational(-2, 5) ** 2}
        assert E.subs(values) == sp.zeros(4, 4)
        assert [e.subs(values) for e in div] == [0] * 4


def test_lambda_zero_forces_static():
    """With Lambda = 0, U = h(tau) + V(x) solves the equations only if h is constant
    (Majumdar-Papapetrou): the E_ij equation becomes 2 U h'' + 3 h'^2 = 0 for every x."""
    h = sp.Function("h")(tau)
    V = sp.Function("V")(*X)
    _, E, _ = residuals(h + V)
    space = sp.simplify(E[1, 1].subs(Lam, 0) / (-(h + V) ** 2))
    assert sp.simplify(space - (2 * (h + V) * h.diff(tau, 2) + 3 * h.diff(tau) ** 2)) == 0
    # d/dx of the bracket is 2 V_x h'': nonconstant V forces h'' = 0, then 3 h'^2 = 0.
    assert sp.simplify(space.diff(x) - 2 * V.diff(x) * h.diff(tau, 2)) == 0


def test_homothety_between_kastor_traschen_spacetimes():
    """(tau, x, m, H) -> (lam tau, lam x, lam m, H / lam) multiplies g by lam^2, so only
    H m_i and H |x_i - x_j| are physical. Lambda -> 0 at fixed masses is the same as
    m_i -> 0 at fixed Lambda: tiny holes, not an asymptotically flat merger."""
    lam, m, d = sp.symbols("lambda m d", positive=True)
    U = metrics.kastor_traschen_U([m, m], [(0, 0, d), (0, 0, -d)])
    scaled = U.subs(
        {tau: lam * tau, x: lam * x, y: lam * y, z: lam * z, d: lam * d, m: lam * m, H: H / lam},
        simultaneous=True,
    )
    assert sp.simplify(scaled - U) == 0  # U is invariant, so g = diag(-U^-2 dtau^2, U^2 dx^2) scales by lam^2


def test_single_centre_is_charge_equals_mass_reissner_nordstrom_de_sitter():
    """One centre: the areal radius R = r U satisfies |dR|^2 = (1 - m/R)^2 - H^2 R^2,
    the RNdS function 1 - 2m/R + Q^2/R^2 - Lambda R^2 / 3 with Q = m, Lambda = 3 H^2.
    |dR|^2 is a scalar, so this identification is coordinate independent."""
    m, rr = sp.Symbol("m", positive=True), sp.Symbol("r", positive=True)
    th, ph = sp.symbols("theta phi")
    U = H * tau + m / rr
    st = Spacetime(sp.diag(-U**-2, U**2, U**2 * rr**2, U**2 * rr**2 * sp.sin(th) ** 2), (tau, rr, th, ph))
    R = rr * U
    grad2 = sum(st.inverse[a, b] * R.diff(st.coords[a]) * R.diff(st.coords[b]) for a in range(4) for b in range(4))
    assert sp.simplify(grad2 - ((1 - m / R) ** 2 - H**2 * R**2)) == 0


def test_extremal_mass_bound():
    """Black-hole and cosmological horizons, roots of |H| R^2 - R + M = 0, coincide at
    M = 1 / (4 |H|): the bound of Brill et al. (1994) for coalescence."""
    R, Mt, h = sp.symbols("R M h", positive=True)  # h = |H|
    assert sp.solve(sp.discriminant(h * R**2 - R + Mt, R), Mt) == [1 / (4 * h)]


def test_slices_are_umbilic_with_constant_mean_curvature():
    """On tau = const, K_ij = (1/2N) d_tau gamma_ij = U_tau gamma_ij (sign: K > 0 expands).
    For U = H tau + V(x): K_ij = H gamma_ij, K = 3H, constant on the slice. Hence a 2-surface
    is marginally outer trapped (theta_+ = tr_S K + H_S = 0) iff its mean curvature is
    H_S = -2H: MOTSs in Kastor-Traschen slices are constant-mean-curvature surfaces."""
    V = sp.Function("V")(*X)
    U = H * tau + V
    lapse, gamma = 1 / U, U**2 * sp.eye(3)
    K = (gamma.diff(tau) / (2 * lapse)).applyfunc(sp.simplify)
    assert (K - H * gamma).applyfunc(sp.simplify) == sp.zeros(3, 3)


def test_single_centre_mots_are_the_rnds_horizons():
    """Coordinate spheres about one centre: theta_+ = 2H + H_S with
    H_S = U^-1 (2/r + 2 U_r / U). Setting theta_+ = 0 gives exactly 1 - m/R = -H R
    (R = r U), the black-hole and cosmological horizons of RNdS for H < 0."""
    m, rr = sp.Symbol("m", positive=True), sp.Symbol("r", positive=True)
    U = H * tau + m / rr
    theta_plus = 2 * H + (2 / rr + 2 * U.diff(rr) / U) / U
    R = rr * U
    # theta_+ * R / 2 = (1 - m/R) + H R identically, so theta_+ = 0 <=> 1 - m/R = -H R.
    assert sp.simplify(theta_plus * R / 2 - ((1 - m / R) + H * R)) == 0


def test_original_cosmological_form_matches():
    """Kastor-Traschen (1993) write ds^2 = -Omega^-2 dt^2 + a^2 Omega^2 dx^2 with a = e^{Ht},
    Omega = 1 + sum m_i / (a r_i). With tau = a / H (so dtau = a dt) and U = a Omega this is
    the BHKT form: check U = H tau + m / r and g_tt dt^2 = -U^-2 dtau^2."""
    t, m, rr = sp.symbols("t m r", positive=True)
    a = sp.exp(H * t)
    Omega = 1 + m / (a * rr)
    tau_of_t = a / H
    U = (a * Omega).subs(t, sp.log(H * tau) / H)
    assert sp.simplify(U - (H * tau + m / rr)) == 0
    dtau_dt = tau_of_t.diff(t)
    assert sp.simplify(-(Omega**-2) - (-(a * Omega) ** -2) * dtau_dt**2) == 0
