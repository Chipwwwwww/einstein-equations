"""Independent identities/controls for the all-mode horizon-gap argument."""

import pytest
import sympy as sp
from flint import arb

from research.calc.horizon_gap import geometric_weights, krawczyk_from_jets
from research.calc.fold_proof import J


def test_radial_graph_metric_and_energy_weights():
    """Induced metric computed from the embedding, including both angles."""
    th, ph = sp.symbols("theta phi", real=True)
    h = sp.Function("h", positive=True)(th)
    psi, beta = sp.symbols("psi beta", positive=True)
    X = sp.Matrix([h*sp.sin(th)*sp.cos(ph), h*sp.sin(th)*sp.sin(ph), h*sp.cos(th)])
    derivatives = [X.diff(th), X.diff(ph)]
    q = sp.simplify(psi**4 * sp.Matrix(2, 2, lambda i,j: derivatives[i].dot(derivatives[j])))
    assert sp.simplify(q[0,0] - psi**4*(h*h+sp.diff(h,th)**2)) == 0
    assert sp.simplify(q[1,1] - psi**4*h*h*sp.sin(th)**2) == 0
    assert q[0,1] == 0
    # On the open polar chart sin(theta)>0, extended smoothly at the poles.
    c = h/sp.sqrt(h*h+sp.diff(h,th)**2)
    density = psi**4*h*h/c
    assert sp.simplify(beta**2*density/q[0,0] - beta**2*c) == 0
    assert sp.simplify(beta**2*density*sp.sin(th)**2/q[1,1] - beta**2/c) == 0


def test_ground_state_operator_identity_in_two_angles():
    th, ph = sp.symbols("theta phi", real=True)
    P, c, beta = [sp.Function(n, positive=True)(th) for n in ("P", "c", "beta")]
    v = sp.Function("v")(th, ph)
    def lap(f):
        return c/(P*sp.sin(th))*sp.diff(c*sp.sin(th)*sp.diff(f,th),th) \
            + sp.diff(f,ph,2)/(P*sp.sin(th)**2)
    # L=-Delta+V; V=Delta(beta)/beta, so L beta=0 exactly.
    transformed = (-lap(beta*v) + lap(beta)*v)/beta
    expected = -lap(v) - 2*c*c/P * sp.diff(beta,th)/beta * sp.diff(v,th)
    assert sp.simplify(transformed-expected) == 0


def test_round_sphere_calibration_and_physical_normal_factor():
    # d=0 gives Schwarzschild conformal factor. A unit coordinate sphere
    # psi=3/2 has physical radius 9/4. Unit Euclidean normal variation gives
    # physical Jacobi amplitude psi^2=9/4 (the synthetic operator has V=0).
    w = geometric_weights(arb(1), arb(0), arb(-1), arb(0), arb(1), arb(0), arb(0))
    assert w["phi"].contains(arb(9)/4)
    assert w["c"].contains(1)
    assert (2*w["energy_weight"]/w["mass_weight"]).contains(arb(32)/81)
    # This calibration is NOT a claim that the Schwarzschild minimal-sphere
    # Jacobi operator has a zero eigenvalue.


def test_approximate_ground_state_residual_identity():
    x = sp.symbols("x", real=True)
    beta = 2 + x  # positive on the round-sphere coordinate x=cos(theta)
    v = 1 + x**2
    V = 3*x
    lap = lambda f: sp.diff((1-x*x)*sp.diff(f,x),x)
    L = lambda f: -lap(f)+V*f
    r = L(beta)/beta
    # Q[beta v] equals weighted gradient energy plus r|beta v|^2.
    left = sp.integrate(beta*v*L(beta*v),(x,-1,1))
    right = sp.integrate(beta**2*(1-x*x)*sp.diff(v,x)**2 + r*beta**2*v*v,(x,-1,1))
    assert sp.simplify(left-right)==0


def test_reparametrization_does_not_change_weights():
    rho,z,sa,ca,d = map(arb, (3,4,0,1,1))
    base = geometric_weights(rho,z,sa,ca,arb(0),arb(1),d)
    changed = geometric_weights(rho,z,sa,ca,arb(13),arb(1),d)
    assert all(base[k].overlaps(changed[k]) for k in base)


def test_invalid_star_shape_and_nonfinite_bounds_are_rejected():
    with pytest.raises(ArithmeticError):
        geometric_weights(arb(1),arb(0),arb(1),arb(0),arb(-1),arb(0),arb(0))
    with pytest.raises(ArithmeticError):
        geometric_weights(arb(0),arb(1),arb(0),arb(1),arb(0),arb(1),arb(2))


def test_krawczyk_on_exact_fold_and_nonroot_box():
    # F=z^2+d, Psi=(F,F_z), centre (0,0), D Psi exact on the box.
    r = 0.125
    Fp = J([arb(0),arb(0),arb(1),arb(0),arb(1),arb(0)])
    Fb = J([arb(0,r+r*r),arb(0,2*r),arb(1),arb(0),arb(1),arb(0)])
    ok,K,_ = krawczyk_from_jets(Fp,Fb,(0.,0.),r)
    assert ok and all(x.contains(0) for x in K)
    wrong = J([arb(2),arb(0),arb(1),arb(0),arb(1),arb(0)])
    ok,_,_ = krawczyk_from_jets(wrong,Fb,(0.,0.),r)
    assert not ok


@pytest.mark.slow
def test_full_horizon_gap_certificate():
    from research.calc.horizon_gap import certify
    result = certify()
    assert result["krawczyk_verified"]
    assert all(result["assertions"].values())
