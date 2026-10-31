"""H4 (notes/09): onset of the common apparent horizon for near-critical Brill-Lindquist data.

The barrier certificates (rigorous, interval arithmetic), the location of the
critical separation (numerical), and the onset-time formula (*) of
research/calc/bl_fold.py against a direct solution, and the AMMS identity.
"""

from research.calc import bl_fold as bf


def test_individual_horizons_certified_for_all_near_critical_separations():
    """For every d in [0.76, 0.78] (d_c = 0.7662), each puncture is enclosed by a minimal
    surface lying in the shell 0.01 <= |x - x_i| <= 0.3; the two shells are disjoint."""
    inner, outer = bf.certify_individual(("0.76", "0.78"))
    assert inner < 0 < outer


def test_critical_separation_is_bracketed():
    """The t = 0 slice has a symmetric common minimal surface at d = 0.765 (mismatch
    changes sign) and none in the search window at d = 0.77."""
    assert bf.BL(0.765).min_mismatch(bf.ZLO, bf.ZHI)[0] < 0
    assert bf.BL(0.770).min_mismatch(bf.ZLO, bf.ZHI)[0] > 0


def test_onset_formula_matches_direct_solution():
    dc = bf.critical_separation()
    assert abs(dc - 0.7661974282) < 1e-8
    c = bf.fold_coefficients(dc, "one")
    # literature (Hui & Lin 2025, bare masses 1, so lengths x2 and areas x4):
    # d_crit = 1.5323948 (1.5323949 in their ref. [42]), A_crit = 196.40795
    assert abs(2 * dc - 1.53239485) < 1e-7
    assert abs(4 * c["area"] - 196.40795) < 1e-4
    # AMMS (2009) Prop. 5.1 in the time-symmetric vacuum case: <phi0, p> = -Int N phi0 |A|^2
    assert abs(c["<phi0,p>"] / c["-Int N phi0 |A|^2"] - 1) < 1e-5
    slope = c["<phi0,dH/dd>"] / c["Int N phi0 |A|^2"]
    delta = 5e-5
    assert abs(bf.onset_time(dc + delta, "one") / delta / slope - 1) < 1e-4
