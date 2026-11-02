from fractions import Fraction as Q
import pytest
from research.calc.gh_ks_geometry_tube import certificate


def test_nonflat_schwarzschild_geometric_tube():
    c=certificate()
    assert c["all_metric_tube_members_Lorentzian"]
    assert c["reference_g_inverse_operator_bound"]==Q(11,3)
    assert c["neumann_multiplier"]==Q(11,50)
    assert c["spatial_metric_eig_min"]==Q(47,50)
    assert c["spatial_metric_eig_max"]==Q(359,150)
    assert c["GH_symmetrizer_50_min_eigen"]==Q(150,359)
    assert c["minus_g_inverse_00_min"]==Q(448,585)
    assert c["full_H3_GH_source_residual_and_vacuum"] is False


def test_fail_closed_on_wide_tube():
    with pytest.raises(ValueError,match="Neumann"):
        certificate(R=Q(1,2))
    with pytest.raises(ValueError,match="positive"):
        certificate(rmin=Q(0))
