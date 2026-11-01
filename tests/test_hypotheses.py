"""Symbolic steps behind the hypothesis tests in notes/09 (H1: Lemma A; H2: hyperboloidal
Kastor-Traschen). Each test is an exact identity; the logical conclusions drawn from them
are stated, with their remaining hypotheses, in the note."""

import sympy as sp

from einstein import Spacetime
from einstein.fields import field_strength, maxwell_stress_energy

t, x, y, z = sp.symbols("t x y z", real=True)
X = (x, y, z)
Lam = sp.Symbol("Lambda")


def test_lemma_A_momentum_constraint_for_conformally_flat_zero_shift():
    """H1 / Lemma A, step 1. For g = -N^2 dt^2 + Psi^4 delta (any N, Psi of all four
    coordinates) the slices are umbilic, K_ij = (K/3) gamma_ij with K = 6 Psi_t / (N Psi),
    and the momentum constraint is G_{ti} = -(2/3) N d_i K. So wherever the momentum density
    vanishes (vacuum, or a purely electric field), K is constant on each slice."""
    N = sp.Function("N")(t, *X)
    Psi = sp.Function("Psi")(t, *X)
    st = Spacetime(sp.diag(-(N**2), Psi**4, Psi**4, Psi**4), (t, *X), simplify=False)
    K = 6 * Psi.diff(t) / (N * Psi)
    for i, v in enumerate(X, start=1):
        assert sp.simplify(st.einstein[0, i] + sp.Rational(2, 3) * N * K.diff(v)) == 0


def test_lemma_A_electric_field_carries_no_momentum():
    """A purely electric potential A = phi dt on a zero-shift metric has T_{ti} = 0, so
    Lemma A applies to Kastor-Traschen-type matter as well as to vacuum."""
    N, Psi, phi = (sp.Function(name)(t, *X) for name in ("N", "Psi", "phi"))
    st = Spacetime(sp.diag(-(N**2), Psi**4, Psi**4, Psi**4), (t, *X), simplify=False)
    T = maxwell_stress_energy(st, field_strength([phi, 0, 0, 0], st.coords))
    assert all(sp.simplify(T[0, i]) == 0 for i in (1, 2, 3))


def hyperboloidal_residual(W):
    """E_ab = G_ab + Lambda g_ab - 8 pi T_ab for g = -W^-2 dt^2 + W^2 h_{-1}, A = W^-1 dt,
    with h_{-1} = 4 delta / (1 - |x|^2)^2 the hyperbolic metric on the Poincare ball."""
    Om2 = 4 / (1 - x**2 - y**2 - z**2) ** 2
    st = Spacetime(sp.diag(-(W**-2), *(3 * [W**2 * Om2])), (t, *X), simplify=sp.simplify)
    T = maxwell_stress_energy(st, field_strength([1 / W, 0, 0, 0], st.coords))
    return st, (st.einstein + Lam * st.metric - 8 * sp.pi * T).applyfunc(sp.simplify)


def test_hyperboloidal_kastor_traschen_reduced_system():
    """H2. For arbitrary W(t, x): E_ti = -2 W_ti / W and
    E_ij = -(2 W W_tt + 3 W_t^2 - Lambda - W^-2) g_ij. Compared with Kastor-Traschen
    (notes/08) the only new term is -W^-2, from the curvature of the hyperbolic slices."""
    W = sp.Function("W")(t, *X)
    st, E = hyperboloidal_residual(W)
    bracket = 2 * W * W.diff(t, 2) + 3 * W.diff(t) ** 2 - Lam - W**-2
    for i, v in enumerate(X, start=1):
        assert sp.simplify(E[0, i] + 2 * W.diff(t, v) / W) == 0
        assert sp.simplify(E[i, i] + bracket * st.metric[i, i]) == 0
    assert all(E[i, j] == 0 for i in (1, 2, 3) for j in (1, 2, 3) if i != j)


def test_hyperboloidal_kastor_traschen_has_no_centres():
    """With W = h(t) + V(x) (forced by E_ti = 0), the x-derivative of the E_ij bracket is
    2 V_x (h'' + W^-3). Where V_x != 0 this forces W^-3 = -h''(t), i.e. W independent of x:
    a contradiction. So V is constant: no black holes, for any Lambda."""
    h = sp.Function("h")(t)
    V = sp.Function("V")(*X)
    W = h + V
    bracket = 2 * W * W.diff(t, 2) + 3 * W.diff(t) ** 2 - Lam - W**-2
    assert sp.simplify(bracket.diff(x) - 2 * V.diff(x) * (h.diff(t, 2) + W**-3)) == 0


def test_hyperboloidal_background_is_milne():
    """Sanity check of the ansatz: with no centres, W = sqrt(2t) solves the equations with
    Lambda = 0; proper time T = sqrt(2t) gives ds^2 = -dT^2 + T^2 h_{-1}, the Milne wedge of
    Minkowski space, whose slices are umbilic with K = 3/T != 0."""
    st, E = hyperboloidal_residual(sp.sqrt(2 * t))
    assert E.subs(Lam, 0).applyfunc(sp.simplify) == sp.zeros(4, 4)
