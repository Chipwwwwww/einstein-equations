"""Validated H3 variable-reference 4x4 inverse matrix and Lorentzian tube.

A standard background-centered Neumann series; not a new PDE theorem.
The inputs B,J,Lambda are certified reference inverse/spatial bounds.
All estimates concern the ENTIRE H3 ball, not points sampled from it.
"""
from fractions import Fraction as Q
import json
P=Q(48)
C0=Q(6)

def validate(B,J,spatial_lower,R,product=P,embedding=C0):
    B,J,spatial_lower,R,product,embedding=map(Q,(B,J,spatial_lower,R,product,embedding))
    if min(B,J,spatial_lower,R,product,embedding)<=0:
        raise ValueError('all input certificate bounds must be positive')
    z=product**2*B*R
    if z>=1:
        raise ValueError('H3 matrix-inverse Neumann series does not converge')
    pointwise=embedding*R
    if J*pointwise>=1:
        raise ValueError('pointwise Lorentzian homotopy cannot be certified')
    if pointwise>=spatial_lower:
        raise ValueError('spatial positive-definiteness cannot be certified')
    M=B/(1-z)
    return {'H3_product_constant':product,'H3_embedding_constant':embedding,
      'reference_inverse_H3_bound':B,'reference_inverse_Linf_operator_bound':J,
      'reference_spatial_metric_eigen_lower':spatial_lower,'H3_tube_radius':R,
      'Neumann_ratio':z,'inverse_H3_bound':M,
      'inverse_H3_Lipschitz':product**2*M*M,
      'pointwise_metric_operator_error':pointwise,
      'Lorentzian_inertia_and_spatial_positivity':True,
      'reference_may_vary_in_space':True,
      'binary_black_hole_exterior_H3_certificate':False}

def periodic_variable_reference_example():
    # G(x)=diag(-(1+a cos x1),(1+a cos x1),1,1), a=1/1000.
    # The exact Fourier H3 norm ||a cos x1||3 = 1/500.
    a=Q(1,1000);d=2*a
    inverse_scalar_deviation=d/(1-P*d)
    B=4+2*inverse_scalar_deviation
    if not B<5:
        raise ValueError("reference inverse norm enclosure failed")
    return validate(B=Q(5),J=Q(2),spatial_lower=Q(999,1000),R=Q(1,10**6))

def serial(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    return x

if __name__=='__main__':
    print(json.dumps(serial(periodic_variable_reference_example()),indent=2))
