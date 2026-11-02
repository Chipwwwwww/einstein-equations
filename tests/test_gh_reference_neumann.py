from fractions import Fraction as Q
import pytest
from research.calc.gh_reference_neumann import validate,periodic_variable_reference_example

def test_nonconstant_reference_full_H3_inverse_tube():
    a=periodic_variable_reference_example()
    assert a['Neumann_ratio']==Q(36,3125)
    assert a['pointwise_metric_operator_error']==Q(3,500000)
    assert a['inverse_H3_bound']==Q(5)/(1-Q(36,3125))
    assert a['inverse_H3_Lipschitz']==2304*a['inverse_H3_bound']**2
    assert a['Lorentzian_inertia_and_spatial_positivity']
    assert not a['binary_black_hole_exterior_H3_certificate']

def test_inverse_and_lorentzian_failure_modes():
    with pytest.raises(ValueError,match='Neumann'):
        validate(5,2,1,Q(1,100))
    with pytest.raises(ValueError,match='Lorentzian'):
        validate(Q(1,1000),100,1,Q(1,100))
    with pytest.raises(ValueError,match='spatial'):
        validate(Q(1,1000),Q(1,100),Q(1,1000),Q(1,100))
    with pytest.raises(ValueError,match='positive'):
        validate(5,2,1,Q(0))
