"""The Lohner integrator (research/calc/lohner.py) against flows with exact solutions, including
the two cases that defeat plain interval boxes: rotation and a hyperbolic saddle."""

import math

from flint import arb

from research.calc import lohner


def _run(f, X, h, nsteps, L=14):
    S = lohner.Set.from_box(X)
    for _ in range(nsteps):
        S, _ = lohner.step(f, S, h, L)
    return S.hull()


def test_rotation_keeps_small_enclosures():
    """x' = -y, y' = x for time 20 (about three turns): exact solution (cos 20, sin 20) from
    (1, 0). Plain boxes would grow by a factor of order e^20; Lohner's stay small."""
    f = lambda y: [-y[1], y[0]]  # noqa: E731
    out = _run(f, [arb(1, 1e-12), arb(0, 1e-12)], 0.25, 80)
    assert out[0].contains(arb(20).cos()) and out[1].contains(arb(20).sin())
    assert max(float(o.rad()) for o in out) < 1e-10


def test_hyperbolic_saddle():
    """x' = y, y' = x from (1, -1 + 1e-6) for time 8: x = cosh t - (1 - 1e-6) sinh t."""
    f = lambda y: [y[1], y[0]]  # noqa: E731
    out = _run(f, [arb(1), arb("-0.999999")], 0.25, 32)
    t = arb(8)
    x = t.cosh() + arb("-0.999999") * t.sinh()
    assert out[0].contains(x)
    assert float(out[0].rad()) < 1e-10  # errors grow like the true expansion e^8 only


def test_nonlinear_with_known_solution():
    """Riccati y' = -y^2 coupled to t' = 1: y(t) = 1/(1 + t) from y(0) = 1."""
    f = lambda y: [y[0] * 0 + 1, -(y[1] * y[1])]  # noqa: E731
    out = _run(f, [arb(0), arb(1)], 0.125, 16)
    assert out[1].contains(arb(1) / 3)
    assert float(out[1].rad()) < 1e-8


def test_parameters_kept_out_of_the_frame():
    """lohner.step_param: x' = -p y, y' = p x with p in [0.999, 1.001], from (1, 0) for time 20.
    The enclosure must contain (cos 20p, sin 20p) for p at both ends and the middle, and stay
    close to the true spread 20 * 0.001 (the parameter error is carried by P q, never
    re-boxed; the excess is second order)."""
    f = lambda y: [-y[2] * y[1], y[2] * y[0], y[0] * 0]  # noqa: E731
    S = lohner.PSet.from_box([arb(1), arb(0)], [arb(1, 0.001)])
    for _ in range(80):
        S, _ = lohner.step_param(f, S, 0.25, 14)
    X = S.hull()
    for p in ("0.999", "1", "1.001"):
        t = arb(p) * 20
        assert X[0].contains(t.cos()) and X[1].contains(t.sin())
    assert max(float(x.rad()) for x in X) < 0.025
    assert max(float(r.rad()) for r in S.r) < 1e-3
