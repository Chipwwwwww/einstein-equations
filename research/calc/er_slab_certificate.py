"""Verified rational posteriori bounds for polarized Einstein--Rosen vacuum waves.

This is a symmetry-reduced, curvature-carrying cylindrical Einstein spacetime,
NOT an asymptotically flat binary black-hole merger or a full GH theorem.

    psi_tt - psi_rr - psi_r/r = 0
    gamma_r = r (psi_t**2 + psi_r**2)
    gamma_t = 2r psi_t psi_r

Analytic rational estimates are used on a shrinking annular causal domain.
"""

from fractions import Fraction as F
import json


def er_certificate(*, a=F(1, 10), eps=F(1, 100), T=F(1, 4),
                   rmin=F(1), rmax=F(2), slabs=2):
    """Manufactured approximate Einstein--Rosen metric with rigorous bounds.

    psi_app = (a+eps)t^2 + (a/2)r^2
    gamma_app = 2(a+eps)^2*t^2*r^2 + a^2*r^4/4
    The independent exact vacuum comparator has eps=0.
    """
    if not (a > 0 and eps > 0 and T > 0 and rmin > 0 and rmax-rmin >= 2*T
            and slabs >= 1 and isinstance(slabs, int)):
        raise ValueError('Invalid positive annular causal diamond')

    # Exact pointwise domain-wide residual and approximate derivative bounds.
    # rho=2eps; q_r=0; q_t=4eps(a+eps)t*r^2.
    delta = 2*eps
    q = 4*eps*(a+eps)*T*rmax**2
    B = max(2*(a+eps)*T, a*rmax)
    psi_min = a*rmin**2/2
    psi_max = (a+eps)*T**2 + a*rmax**2/2
    gamma_max = 2*(a+eps)**2*T**2*rmax**2 + a**2*rmax**4/4

    P = E = G = F(0)
    h = T / slabs
    if h >= rmin:
        raise ValueError('Require each characteristic step h < rmin')
    rows=[]
    for j in range(slabs):
        # Characteristic radial error: Q <= P + h*(Q/rmin+delta).
        Q = (P + h*delta)/(1-h/rmin)
        E += h*Q
        # gamma(t,r) is integrated at fixed r from gamma_t=2r psi_t psi_r.
        G += h*(q + 2*rmax*(2*B+Q)*Q)
        P = Q
        rows.append(dict(slab=j+1, time=(j+1)*h,
                         radial_null_error=P, psi_error=E, gamma_error=G))

    # Majorant exp(z) <= 1/(1-z) valid for rational z in [0,1).
    max_chi_app = 2*(gamma_max - psi_min)
    if max_chi_app >= 1 or 2*psi_max >= 1:
        raise ValueError('Elementary rational exponential majorants unavailable')
    Amax = F(1) if max_chi_app <= 0 else 1/(1-max_chi_app)
    Zmax = 1/(1-2*psi_max)
    err_chi = 2*(G+E)
    err_psi_exp = 2*E
    if err_chi >= 1 or err_psi_exp >= 1:
        raise ValueError('Error tube exceeds chosen exponential majorant')
    err_tt_rr = Amax*err_chi/(1-err_chi)
    err_zz = Zmax*err_psi_exp/(1-err_psi_exp)
    err_phiphi = rmax**2*err_psi_exp/(1-err_psi_exp)
    actual_psi = eps*T**2
    actual_gamma = 2*((a+eps)**2-a**2)*T**2*rmax**2
    assert E >= actual_psi and G >= actual_gamma
    return {
        'classification': 'VALIDATED_NONFLAT_CYLINDRICAL_VACUUM_SUBSECTOR',
        'domain': '0<=t<=T, rmin+t<=r<=rmax-t, z and phi arbitrary',
        'a':a, 'eps':eps, 'T':T, 'rmin':rmin, 'rmax':rmax,
        'radial_wave_residual_rho':delta, 'gamma_r_constraint_residual':F(0),
        'gamma_t_residual_bound':q, 'psi_gradient_app_bound':B,
        'psi_error_bound':E, 'gamma_error_bound':G,
        'null_derivative_error_bound':P,
        'metric_tt_rr_error_bound':err_tt_rr,
        'metric_zz_error_bound':err_zz,
        'metric_phiphi_error_bound':err_phiphi,
        'metric_component_error_bound':max(err_tt_rr, err_zz, err_phiphi),
        'independent_exact_psi_error':actual_psi,
        'independent_exact_gamma_error':actual_gamma,
        'slabs': rows,
        'asymptotically_flat_binary_merger_certified':False,
    }


def fmt(v):
    if isinstance(v, F): return f'{v.numerator}/{v.denominator}'
    if isinstance(v, dict): return {k:fmt(x) for k,x in v.items()}
    if isinstance(v, list): return [fmt(x) for x in v]
    return v


if __name__ == '__main__':
    print(json.dumps(fmt(er_certificate()),indent=2))
