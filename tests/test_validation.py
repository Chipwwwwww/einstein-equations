"""Validation of the curvature code against independent calculations (research/README.md).

Each check here would catch a different class of bug in einstein/tensors.py:

- an independent all-lower-index Riemann formula (index-order and sign errors);
- the algebraic symmetries of R_abcd (a wrong contraction or transposition);
- the closed-form Ricci tensor of a conformally flat metric (overall sign);
- Reissner-Nordstrom: a non-vacuum solution with known curvature (matter coupling);
- Kerr: a vacuum solution with off-diagonal terms (exact rational identity in chi = cos theta).
"""

from itertools import product

import pytest
import sympy as sp

from einstein import Spacetime, metrics
from einstein.fields import field_strength, maxwell_divergence, maxwell_stress_energy
from einstein.metrics import M, Q, a_kerr, r

ZERO4 = sp.zeros(4, 4)


def riemann_all_lower(g, coords):
    """R_abcd from second derivatives of g_ab, written independently of tensors.py:

    R_abcd = 1/2 (g_ad,bc + g_bc,ad - g_ac,bd - g_bd,ac)
             + g_ef (Gamma^e_bc Gamma^f_ad - Gamma^e_bd Gamma^f_ac),

    with the Christoffel symbols of the first kind built here from scratch.
    """
    n = len(coords)
    gi = g.inv()
    d = lambda expr, i: sp.diff(expr, coords[i])  # noqa: E731
    first_kind = {
        (e, b, c): (d(g[e, b], c) + d(g[e, c], b) - d(g[b, c], e)) / 2
        for e, b, c in product(range(n), repeat=3)
    }
    gamma = {
        (a, b, c): sum(gi[a, e] * first_kind[(e, b, c)] for e in range(n))
        for a, b, c in product(range(n), repeat=3)
    }
    out = {}
    for a, b, c, dd in product(range(n), repeat=4):
        second = (
            d(d(g[a, dd], b), c) + d(d(g[b, c], a), dd) - d(d(g[a, c], b), dd) - d(d(g[b, dd], a), c)
        ) / 2
        quad = sum(
            g[e, f] * (gamma[(e, b, c)] * gamma[(f, a, dd)] - gamma[(e, b, dd)] * gamma[(f, a, c)])
            for e, f in product(range(n), repeat=2)
        )
        out[(a, b, c, dd)] = second + quad
    return out


def lowered_riemann(st):
    """R_abcd = g_ae R^e_bcd from the library."""
    n = st.dim
    return {
        idx: sum(st.metric[idx[0], e] * st.riemann[e, idx[1], idx[2], idx[3]] for e in range(n))
        for idx in product(range(n), repeat=4)
    }


def generic_metric():
    """A 3-metric with three free functions of two variables and an off-diagonal term."""
    u, v, w = sp.symbols("u v w")
    A, B, C = (sp.Function(name)(u, v) for name in "ABC")
    g = sp.Matrix([[A, B, 0], [B, C, 0], [0, 0, A * C]])
    return g, (u, v, w)


def test_riemann_matches_independent_formula():
    g, coords = generic_metric()
    st = Spacetime(g, coords, simplify=False)
    ours, theirs = lowered_riemann(st), riemann_all_lower(g, coords)
    assert all(sp.cancel(ours[i] - theirs[i]) == 0 for i in ours)


def test_riemann_algebraic_symmetries():
    g, coords = generic_metric()
    R = lowered_riemann(Spacetime(g, coords, simplify=False))
    n = len(coords)
    for a, b, c, d in product(range(n), repeat=4):
        assert sp.cancel(R[(a, b, c, d)] + R[(b, a, c, d)]) == 0
        assert sp.cancel(R[(a, b, c, d)] - R[(c, d, a, b)]) == 0
        assert sp.cancel(R[(a, b, c, d)] + R[(a, c, d, b)] + R[(a, d, b, c)]) == 0


def test_conformally_flat_ricci_formula():
    """For g = e^{2w} delta in 3 dimensions,
    R_ab = -(d_a d_b w - d_a w d_b w) - (lap w + |dw|^2) delta_ab   (Wald, appendix D)."""
    X = sp.symbols("x y z")
    w = sp.Function("w")(*X)
    st = Spacetime(sp.exp(2 * w) * sp.eye(3), X, simplify=False)
    grad2 = sum(sp.diff(w, v) ** 2 for v in X)
    lap = sum(sp.diff(w, v, 2) for v in X)
    expected = sp.Matrix(
        3,
        3,
        lambda a, b: -(sp.diff(w, X[a], X[b]) - sp.diff(w, X[a]) * sp.diff(w, X[b]))
        - (lap + grad2) * (1 if a == b else 0),
    )
    assert (st.ricci - expected).applyfunc(sp.simplify) == sp.zeros(3, 3)


def test_minkowski_cartesian_is_flat():
    st = Spacetime(sp.diag(-1, 1, 1, 1), sp.symbols("t x y z"))
    assert st.ricci == ZERO4 and st.kretschmann == 0


def test_reissner_nordstrom_solves_einstein_maxwell():
    st = metrics.reissner_nordstrom()
    F = field_strength(metrics.reissner_nordstrom_potential(), st.coords)
    T = maxwell_stress_energy(st, F)
    assert (st.einstein - 8 * sp.pi * T).applyfunc(sp.simplify) == ZERO4
    assert [sp.simplify(e) for e in maxwell_divergence(st, F)] == [0] * 4


def test_reissner_nordstrom_kretschmann():
    K = metrics.reissner_nordstrom().kretschmann
    assert sp.simplify(K - 8 * (6 * M**2 * r**2 - 12 * M * Q**2 * r + 7 * Q**4) / r**8) == 0


def test_kerr_is_vacuum_at_fixed_parameters():
    """M = 1, a = 3/5. Every component of R_ab is then a rational function of
    (r, chi) and cancel() reduces it to 0: an exact identity, not a spot check (~35 s)."""
    st = metrics.kerr()
    g = st.metric.subs({M: 1, a_kerr: sp.Rational(3, 5)})
    assert Spacetime(g, st.coords, simplify=sp.cancel).ricci == ZERO4


@pytest.mark.slow
def test_kerr_is_vacuum():
    """The same identity with M and a symbolic (~3 min)."""
    assert metrics.kerr().ricci == ZERO4


def test_adm_evolution_equation_sign():
    """The 3+1 identity used in notes/07 and notes/09, with K_ij = (1/2N) d_t gamma_ij:

        d_t K_ij = D_i D_j N - N (R3_ij + K K_ij - 2 K_ik K^k_j) + N R4_ij

    checked for a zero-shift metric with four free functions of (t, x)."""
    t, X1, X2, X3 = sp.symbols("t x1 x2 x3")
    A, B, C, N = (sp.Function(name)(t, X1) for name in ("A", "B", "C", "N"))
    gamma = sp.diag(A, B, C)
    four = Spacetime(sp.diag(-(N**2), A, B, C), (t, X1, X2, X3), simplify=False)
    three = Spacetime(gamma, (X1, X2, X3), simplify=False)
    K = gamma.diff(t) / (2 * N)
    Kmix = three.inverse * K  # K^i_j
    trK = Kmix.trace()
    gam = three.christoffel
    X = (X1, X2, X3)
    for i, j in product(range(3), repeat=2):
        DDN = sp.diff(N, X[i], X[j]) - sum(gam[k, i, j] * sp.diff(N, X[k]) for k in range(3))
        KK = sum(K[i, k] * Kmix[k, j] for k in range(3))
        rhs = DDN - N * (three.ricci[i, j] + trK * K[i, j] - 2 * KK) + N * four.ricci[i + 1, j + 1]
        assert sp.simplify(K[i, j].diff(t) - rhs) == 0


def test_conformal_mean_curvature_formula():
    """The formula behind every horizon calculation in research/calc: for gamma = W^2 delta
    and a level set of F with flat unit normal n = grad F / |grad F|, the gamma-mean
    curvature H = div_gamma(nu) (nu the gamma-unit normal) equals
    W^-1 (H_flat + 2 n . grad W / W), H_flat = div_flat(n). Checked from first principles:
    H = (1/sqrt det gamma) d_i (sqrt det gamma nu^i), with arbitrary F and W."""
    X = sp.symbols("x y z")
    F = sp.Function("F")(*X)
    W = sp.Function("W")(*X)
    grad = [sp.diff(F, v) for v in X]
    norm = sp.sqrt(sum(c**2 for c in grad))
    n = [c / norm for c in grad]
    nu_up = [W**-1 * c for c in n]  # gamma^{ij} d_j F / |dF|_gamma = W^-1 n^i
    H = sum(sp.diff(W**3 * nu_up[i], X[i]) for i in range(3)) / W**3
    H_flat = sum(sp.diff(n[i], X[i]) for i in range(3))
    expected = (H_flat + 2 * sum(n[i] * sp.diff(W, X[i]) for i in range(3)) / W) / W
    assert sp.simplify(H - expected) == 0
