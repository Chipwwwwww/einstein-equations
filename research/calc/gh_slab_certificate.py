"""Exact-rational certificate for a RESTRICTED harmonic gauge-wave Einstein test.

This is not a generic Einstein solver, a validated discretisation, or a merger
certificate.  It checks explicitly proved bounds on a 2+2 product metric.

Run: python -m research.calc.gh_slab_certificate
"""

from dataclasses import dataclass
from fractions import Fraction as F
import json


@dataclass(frozen=True)
class NullSlab:
    width: F
    residual: F


@dataclass(frozen=True)
class NullError:
    value: F  # || log(H_true) - log(H_app) ||_infty
    null_derivative: F  # max ||(dt +/- dx)(log H_true-log H_app)||_infty


def propagate(slab: NullSlab, incoming: NullError) -> NullError:
    """Rigorous null-characteristic recursion on a shrinking causal diamond."""
    h, d = slab.width, slab.residual
    if h < 0 or d < 0 or incoming.value < 0 or incoming.null_derivative < 0:
        raise ValueError("All certificate quantities must be nonnegative")
    return NullError(
        incoming.value + h * incoming.null_derivative + h * h * d / 2,
        incoming.null_derivative + h * d,
    )


def gauge_wave_certificate(
    *, a: F = F(1, 2), epsilon: F = F(1, 100),
    T: F = F(1, 4), X: F = F(1), nslabs: int = 2,
) -> dict:
    """Certify H_app=1+a(x-t)+epsilon*t^2 on 0<=t<=T, |x|+t<=X.

    Relies on EXACT inequalities x-t>=-X and x+t<=X on the diamond.
    R=(dt^2-dx^2) log(H_app)
     =2*epsilon*(1+a*(x+t)-epsilon*t^2)/H_app**2.
    The wave initial data coincide with the exact vacuum solution H=1+a(x-t).
    """
    if not (a >= 0 and epsilon >= 0 and X >= T > 0 and nslabs >= 1):
        raise ValueError("Need a>=0, epsilon>=0, X>=T>0, nslabs>=1")
    lower_H = 1 - a * X
    upper_H = 1 + a * X + epsilon * T * T
    lower_numerator = 1 - a * X - epsilon * T * T
    upper_numerator = 1 + a * X
    if lower_H <= 0 or lower_numerator < 0:
        raise ValueError("Cannot use the nonnegative-residual numerator bound")
    residual = 2 * epsilon * upper_numerator / (lower_H * lower_H)
    # rho = |(dt^2-dx^2)log H_app| bound; |Ric_tt|=|Ric_xx|<=rho/2.
    error = NullError(F(0), F(0))
    h = T / nslabs
    rows = []
    for j in range(nslabs):
        error = propagate(NullSlab(h, residual), error)
        rows.append({"slab": j + 1, "through_t": (j + 1) * h,
                     "log_error": error.value, "null_error": error.null_derivative})
    if error.value >= 1:
        raise ValueError("Basic rational exponential majorant requires log error < 1")
    # For 0<=z<1: exp(z)<=1/(1-z), exp(z)-1<=z/(1-z).
    c0_error = upper_H * error.value / (1 - error.value)
    inverse_c0_error = error.value / (lower_H * (1 - error.value))
    # |d log Happ| <= (a+2 epsilon T)/lower_H for t and x.
    log_gradient_app = (a + 2 * epsilon * T) / lower_H
    c1_error = upper_H * (error.value * log_gradient_app + error.null_derivative) / (1 - error.value)
    # Exact independent comparator: H_true=1+a(x-t), |H_true-H_app|=epsilon*t^2.
    exact_c0_error = epsilon * T * T
    assert c0_error >= exact_c0_error
    assert c0_error < lower_H, "Error ball must preserve the elementary positivity bound"
    assert error.value == residual * T * T / 2
    assert error.null_derivative == residual * T
    return {
        "classification": "VALIDATED_RESTRICTED_GAUGE_WAVE_ONLY",
        "initial_data": "exactly matched at t=0 in logH and dt(logH)",
        "domain": "0<=t<=T, |x|+t<=X, y,z arbitrary",
        "a": a, "epsilon": epsilon, "T": T, "X": X,
        "Happ_lower": lower_H, "Happ_upper": upper_H,
        "R_wave_bound": residual,
        "Ric_tt_Ric_xx_bound": residual / 2,
        "log_error_bound": error.value,
        "null_derivative_bound": error.null_derivative,
        "metric_C0_error_bound": c0_error,
        "metric_inverse_C0_error_bound": inverse_c0_error,
        "metric_C1_derivative_error_bound": c1_error,
        "independently_known_actual_metric_C0_error": exact_c0_error,
        "slabs": rows,
        "full_einstein_merger_certified": False,
    }


def _stringify(obj):
    if isinstance(obj, F):
        return f"{obj.numerator}/{obj.denominator}"
    if isinstance(obj, dict):
        return {k: _stringify(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_stringify(v) for v in obj]
    return obj


if __name__ == "__main__":
    print(json.dumps(_stringify(gauge_wave_certificate()), indent=2))
