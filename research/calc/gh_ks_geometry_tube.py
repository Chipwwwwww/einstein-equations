"""Exact rational *geometric* GH hyperbolicity tube on a Schwarzschild KS shell.

Not a full GH evolution certificate. Uniform pointwise metric error is an
INPUT condition, not inferred from a shell boundary Sobolev theorem.
"""
from fractions import Fraction as Q
import json


def certificate(R=Q(1,100),M=Q(1),rmin=Q(3,2),rmax=Q(5,2),C0=Q(6)):
    if min(R,M,rmin,rmax,C0)<=0 or rmin>=rmax:
        raise ValueError("positive radius, mass, shell, and Sobolev bound required")
    hmax=2*M/rmin
    hmin=2*M/rmax
    B=1+2*hmax  # gKS inverse op norm <= 1+2h
    z=C0*R*B
    if z>=1:
        raise ValueError("Neumann nondegenerate Lorentzian tube inequality fails")
    max_inverse_error=B*B*C0*R/(1-z)
    if max_inverse_error>=1+hmin:
        raise ValueError("lapse positivity not certified")
    spatial_min=1-C0*R
    spatial_max=1+hmax+C0*R
    if spatial_min<=0:
        raise ValueError("spatial metric positivity not certified")
    sym_min=min(Q(1),1/spatial_max)
    lapse_inv00_lower=1+hmin-max_inverse_error
    lapse_inv00_upper=1+hmax+max_inverse_error
    return {
      "classification":"KS_GEOMETRIC_HYPERBOLICITY_ONLY",
      "shell":(rmin,rmax),"mass":M,"H3_metric_error_radius":R,
      "C0_embedding":C0,"hmin":hmin,"hmax":hmax,
      "reference_g_inverse_operator_bound":B,
      "neumann_multiplier":z,"inverse_error_operator_bound":max_inverse_error,
      "spatial_metric_eig_min":spatial_min,
      "spatial_metric_eig_max":spatial_max,
      "GH_symmetrizer_50_min_eigen":sym_min,
      "minus_g_inverse_00_min":lapse_inv00_lower,
      "minus_g_inverse_00_max":lapse_inv00_upper,
      "all_metric_tube_members_Lorentzian":True,
      "full_H3_GH_source_residual_and_vacuum":False,
    }


def exact_json(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:exact_json(v) for k,v in x.items()}
    if isinstance(x,tuple):return [exact_json(v) for v in x]
    return x


if __name__=="__main__":
    print(json.dumps(exact_json(certificate()),indent=2))
