"""Gap 3 near the fold (research/calc/gap3_fold.py): checks of the jet/Lohner machinery behind
the local fold lemma. The lemma itself is `python -m research.calc.gap3_fold 4` (~2 CPU-hours)."""

import pytest
from flint import arb

from research.calc import gap3_fold as gf
from research.calc import lohner


def test_tau_system_is_the_reparametrised_sigma_system():
    """make_f_tau = -z_sw f / v componentwise, with zt' = -z_sw, on the base slot."""
    y = [arb(0.25), arb(-0.3), arb(0.01), arb(-0.24), arb(0.05)]
    d = arb(0.7662)
    F = gf.f_state(*y[:4], d)
    G = gf.make_f_tau((0,))(y + [d, arb(0.63)])
    for i in (0, 2, 3):
        assert G[i].overlaps(-y[4] * F[i] / y[3])
    assert G[1].overlaps(-y[4]) and G[4] == 0 and G[5] == 0 and G[6] == 0


def test_switch_appends_the_height_with_its_jets():
    """to_tau: the new z_sw rows of a sigma-phase set equal z = zt + d/2 (slot 0, through the
    set's parts), zt_z (slot 1) and zt_d + 1/2 (slot 2), and the old rows are unchanged."""
    slots = (0, 1, 2)
    S = gf.start(0.6326, 0.6327, gf.D_C_LO, gf.D_C_LO + 1e-5, slots)
    f = gf.make_f_sigma(slots)
    blocks = [list(range(4 * a, 4 * a + 4)) for a in range(3)]
    for _ in range(8):
        S, _ = lohner.step_param(f, S, 0.25, 14, blocks)
    T = gf.to_tau(S, slots)
    X, Y = S.hull(), T.hull()
    for a in range(3):
        for i in range(4):
            assert Y[5 * a + i].overlaps(X[4 * a + i])
            assert abs(float(Y[5 * a + i].rad() - X[4 * a + i].rad())) < 1e-12
    z = gf.z_abs(S)
    assert Y[4].overlaps(z) and float(Y[4].rad()) < float(z.rad()) + 1e-12
    assert Y[9].overlaps(X[5])
    assert Y[14].overlaps(X[9] + arb(1) / 2)


@pytest.mark.slow
def test_jets_reproduce_the_fold_certificate():
    """A thin run at the fold point (z0 = z* +- 1e-9, d = d*) reproduces Theorem C's
    independent enclosures (research/calc/fold_proof.py, a different integrator and
    parametrisation): F_zz in [240.87, 242.70], F_d in [1.2315, 1.2364], F_z ~ 0."""
    Fz, Fzz, Fd, _ = gf.run(gf.Z_C - 1e-9, gf.Z_C + 1e-9, gf.D_C_LO, gf.D_C_HI, (0, 1, 2, 4))
    assert Fzz.overlaps(arb("241.785", "0.915")) and float(Fzz.rad()) < 1
    assert Fd.overlaps(arb("1.23395", "0.00245")) and float(Fd.rad()) < 0.01
    assert abs(float(Fz.mid())) < 1e-3 and float(Fz.rad()) < 1e-3


@pytest.mark.slow
def test_fold_lemma_piece():
    """One piece of the fold window: F_zz > 0 and F_d > 0 on z0 in [0.6366, 0.6368],
    d in [D_C_LO, D1] (measured: F_zz >= 150, F_d >= 0.657)."""
    _, fzz, fd, n, _ = gf.check_piece((0.6366, 0.6368))
    assert fzz > 100 and fd > 0.5
