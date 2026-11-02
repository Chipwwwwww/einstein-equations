"""Exact-rational closed H^3 bootstrap: nonlinear 2-field symmetric hyperbolic system.

The GENERAL mathematical theorem is in notes/25-full-hs-a-posteriori.md.
This implementation certifies the stated coupled example, NOT Einstein vacuum.
"""
from fractions import Fraction as Q
import json

C0, C1, CP, CC = 6, 12, 48, 144


def verify(a=Q(1,10000), T=Q(1,100), R=Q(1,1000),
           residual_override=None, require_exact_initial=True, nslabs=1):
    """Certify u_t+A1(u)u_x+A2(u)u_y=0 near v=a*t*(cos(x),sin(x)).

    A1(u)=[[u0,1/20],[1/20,u1]], A2(u)=[[u1,0],[0,-u0]], A3=0;
    S=I; u0=0 is exact initial data. The true comparison u=0 is known
    independently but not used in the bootstrap.
    """
    if not (a>0 and T>0 and R>0 and isinstance(nslabs,int) and nslabs>=1):
        raise ValueError("positive amplitude, time, tube and nslabs required")
    if not require_exact_initial:
        raise ValueError("initial matching is mandatory, not inferred from residual")
    V=3*a*T
    if V+R>Q(1,100):
        raise ValueError("cannot certify the fixed A-matrix H3 bound 1/10")
    MA=Q(1,10)
    rho=3*a+CP*MA*V
    if residual_override is not None:
        if residual_override<rho:
            raise ValueError("claimed residual bound is below proved analytic majorant")
        rho=residual_override
    D=2*C1*(V+R)
    C=CC*2*(V+R)+CP*V
    K=D/2+C
    h=T/nslabs
    if K*h>=1:
        raise ValueError("rational exponential majorant cannot certify slab")
    bound=Q(0)
    slab_bounds=[]
    for _ in range(nslabs):
        bound=(bound+h*rho)/(1-K*h)
        if bound>=R:
            raise ValueError("nonlinear H3 error bootstrap does not close")
        slab_bounds.append(bound)
    return {
      "classification":"VALIDATED_COUPLED_SYMMETRIC_HYPERBOLIC_EXAMPLE_NOT_EINSTEIN",
      "Sobolev":"H^3(T^3;R^2)","a":a,"T":T,"tube":R,
      "C0_H3_to_C0":C0,"C1_H3_to_C1":C1,
      "product_H3":CP,"commutator_H3":CC,
      "v_H3_and_vx_H3_upper":V,"A1_A2_H3_upper":MA,
      "residual_H3_upper":rho,"energy_divergence_bound":D,
      "nonlinear_coefficient_bound":C,"exponent_K":K,
      "K_times_T":K*T,"slab_count":nslabs,"slab_errors":slab_bounds,
      "closed_H3_error":bound,"independent_actual_H3_error_upper":V,
      "bootstrap_closed":True,
      "full_50_field_GH_vacuum_certified":False,
    }


def rational_json(v):
    if isinstance(v,Q): return str(v)
    if isinstance(v,dict): return {k:rational_json(x) for k,x in v.items()}
    if isinstance(v,list): return [rational_json(x) for x in v]
    return v


if __name__=="__main__":
    print(json.dumps(rational_json(verify()),indent=2))
