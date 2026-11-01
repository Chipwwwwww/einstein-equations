"""Algebraic and falsification checks for the repaired global-exclusion route."""

import sympy as sp
from flint import arb

from research.calc import bl_barrier as bb
from research.calc import bl_metric as bm
from research.calc.mots_exact import alpha_prime, upper_axis_curvature, fold_time_derivatives


def test_decimal_time_bridge_covers_the_moore_parameter(monkeypatch):
    from research.calc import bl_foliate as bf
    from research.calc.gap3_interval import _ball
    with bm.working_precision(128):
        for text in ("6.5e-6", "1.2e-5"):
            # This falsifies the old conversion, rather than merely checking
            # that a new helper agrees with its own implementation.
            assert not arb(float(text)).overlaps(arb(text))
            lo, hi = bb.parameter_endpoints(text)
            assert _ball(lo, hi).contains(arb(text))
        seen = {}
        def fake_run(zlo, zhi, d, thi, h, **kw):
            seen.update(z=_ball(zlo, zhi), t=_ball(kw["t_lo"], thi))
            return {}
        monkeypatch.setattr(bf, "run", fake_run)
        z, radius = 0.632632142085469, 1e-10
        bf.barrier_geometry("0.7662", "6.5e-6", -4e-8, z, radius)
        assert seen["t"].contains(arb("6.5e-6"))
        assert seen["z"].contains(arb(z, radius))


def test_radial_bins_cover_a_profile_at_bin_boundaries():
    # Boxes overlap every bin; a boundary must not be dropped by a rounded
    # quotient. This checks enclosure/coverage, not a particular bin count.
    with bm.working_precision(128):
        geo = []
        for k in range(33):
            theta = arb.pi()*k/64
            rho, z = theta.sin(), theta.cos()
            geo.append((bb._down(rho.lower()),bb._up(rho.upper()),
                        bb._down(z.lower()),bb._up(z.upper())))
        # A continuous enclosing box for each intervening arc.
        arcs = [(min(a[0],b[0]),max(a[1],b[1]),min(a[2],b[2]),max(a[3],b[3]))
                for a,b in zip(geo,geo[1:])]
        lo,hi,width = bb.radial_profile({"_geo":arcs,"_z0ball":(1.0,1.0)},nbins=32)
        assert arb(width)*32 > arb.pi()/2
        assert all(a <= 1 <= b for a,b in zip(lo,hi))


def test_exact_rhs_conformal_limit_with_general_K():
    p, px, pz = sp.symbols("p px pz", positive=True)
    k2 = sp.Symbol("k2", real=True)
    # Exact Pythagorean pair avoids any trigonometric simplifier assumptions.
    sa, ca = sp.Rational(-3, 5), sp.Rational(4, 5)
    gam = dict(xx=p**4, yy=p**4, zz=p**4, xz=0)
    dg = {(a, k): (4*p**3*v if k != "xz" else 0)
          for a, v in (("x", px), ("z", pz)) for k in gam}
    dg.update({("y", "xy"): 0, ("y", "yz"): 0})
    K = dict(zip(gam, sp.symbols("Kxx Kyy Kzz Kxz")))
    expected = -k2 + 4*(-sa*px+ca*pz)/p + (
        ca*ca*K["xx"] + 2*ca*sa*K["xz"] + sa*sa*K["zz"] + K["yy"])/p**2
    assert sp.simplify(alpha_prime(sa, ca, k2, gam, dg, K, sp.sqrt)-expected) == 0


def test_exact_rhs_anisotropic_cross_metric_and_nonzero_K():
    gam = dict(xx=3, yy=2, zz=4, xz=sp.Rational(1, 2))
    dg = {(a, k): sp.Rational(i+1, 11)
          for i, (a, k) in enumerate((a, k) for a in ("x", "z") for k in gam)}
    dg.update({("y", "xy"): sp.Rational(2, 7), ("y", "yz"): sp.Rational(-1, 3)})
    K = dict(zip(gam, map(sp.Rational, (1, 2, -1, 3))))
    sa, ca, k2 = -sp.Rational(3, 5), sp.Rational(4, 5), -sp.Rational(6, 7)
    k1 = alpha_prime(sa, ca, k2, gam, dg, K, sp.sqrt)
    assert sp.simplify(bb.theta_general(sa, ca, k1, k2, gam, dg, K, sp.sqrt)) == 0
    # A changed sign in the K term must be detected, not absorbed in a tolerance.
    wrong = alpha_prime(sa, ca, k2, gam, dg, {k: -v for k, v in K.items()}, sp.sqrt)
    assert sp.simplify(bb.theta_general(sa, ca, wrong, k2, gam, dg, K, sp.sqrt)) != 0


def test_regular_pole_limit():
    a, b = sp.symbols("a b", positive=True)
    q, dz, k = sp.symbols("q dz k", real=True)
    gam = dict(xx=a, yy=a, zz=b, xz=0)
    dg = {(x, y): 0 for x in ("x", "z") for y in gam}
    dg.update({("x", "xz"): q, ("z", "xx"): dz, ("z", "yy"): dz,
               ("y", "xy"): 0, ("y", "yz"): q})
    K = dict(xx=k, yy=k, zz=0, xz=0)
    c = upper_axis_curvature(gam, dg, K, sp.sqrt)
    assert sp.simplify(alpha_prime(0, 1, c, gam, dg, K, sp.sqrt)-c) == 0


def test_exact_rhs_encloses_evolved_metric_point():
    # Tail bounds are inherited; this is an enclosure-composition check,
    # not independent verification of the CK certificate or a MOTS location.
    with bm.working_precision(128):
        result = bm.evaluate_metric("0.0000065", "0.25", "0", "0", "0.7662", order=2)
        assert result["status"] == "certified"
        keys = {k: tuple(i+1 for i in ij) for k, ij in bb.KEYS.items()}
        gam = {k: result["g"][i][j] for k, (i,j) in keys.items() if k in ("xx","yy","zz","xz")}
        K = {k: result["dg"][0][i][j]/2 for k, (i,j) in keys.items() if k in gam}
        dg = {(a,k): result["dg"][ax][keys[k][0]][keys[k][1]]
              for a,ax in (("x",1),("z",3)) for k in gam}
        dg.update({("y",k): result["dg"][2][keys[k][0]][keys[k][1]] for k in ("xy","yz")})
        sa, ca, k2 = -arb(3)/5, arb(4)/5, -arb(12)/5
        curvature = alpha_prime(sa, ca, k2, gam, dg, K)
        assert curvature.is_finite()
        assert bb.theta_general(sa, ca, curvature, k2, gam, dg, K).contains(0)


def test_constant_parameter_sweep_does_not_enclose_variable_residual():
    s, e, h = sp.symbols("s e h", positive=True)
    # y''+y=h(s), y(0)=0, y'(0)=pi*e. Every constant h gives
    # y'(2*pi)=pi*e>0, even for arbitrary constant h. The smooth admissible
    # forcing h(s)=-e*cos(s), |h(s)|<=e, instead gives y'(2*pi)=0.
    constant = sp.pi*e*sp.sin(s) + h*(1-sp.cos(s))
    variable = sp.pi*e*sp.sin(s) - e*s*sp.sin(s)/2
    assert sp.simplify(constant.diff(s, 2)+constant-h) == 0
    assert sp.simplify(variable.diff(s, 2)+variable+e*sp.cos(s)) == 0
    assert constant.diff(s).subs(s, 2*sp.pi) == sp.pi*e
    assert variable.diff(s).subs(s, 2*sp.pi) == 0
    assert variable.subs(s, 0) == 0 and variable.diff(s).subs(s, 0) == sp.pi*e


def test_fold_second_derivative_includes_moving_critical_point():
    # An exact fold with nonzero mixed and temporal quadratic coefficients.
    z, d, t = sp.symbols("z d t")
    F = 3*d-2*t+5*z*z/2+7*z*d+11*z*t+13*d*d/2+17*d*t+19*t*t/2
    origin = {z: 0, d: 0, t: 0}
    derivs = [F.diff(*args).subs(origin) for args in
              [(d,), (t,), (z,z), (z,d), (z,t), (d,d), (d,t), (t,t)]]
    c, z1, t2 = fold_time_derivatives(*derivs)
    z2 = -(11*t2)/5
    path = {t: c*d+t2*d*d/2, z: z1*d+z2*d*d/2}
    assert sp.expand(F.subs(path)).coeff(d, 1) == 0
    assert sp.expand(F.subs(path)).coeff(d, 2) == 0
    assert sp.expand(F.diff(z).subs(path)) == 0
    frozen_z_answer = -(derivs[5]+2*c*derivs[6]+c*c*derivs[7])/derivs[1]
    assert frozen_z_answer != t2


def test_exact_drift_conjugation_on_general_axisymmetric_metric():
    x, y = sp.symbols("x y", real=True)
    E, G = sp.Function("E")(x), sp.Function("G")(x)
    w = sp.sqrt(E*G)
    chi, V = sp.Function("chi")(x), sp.Function("V")(x)
    f = sp.Function("f")(x, y)
    def lap(v):
        return sp.diff(w/E*sp.diff(v, x), x)/w + sp.diff(v, y, 2)/G
    def L(v):
        return -lap(v)+2*sp.diff(chi,x)*sp.diff(v,x)/E + (
            V+lap(chi)-sp.diff(chi,x)**2/E)*v
    assert sp.simplify(sp.exp(-chi)*L(sp.exp(chi)*f)+lap(f)-V*f) == 0
    # Rotation invariance alone permits omega=sin(theta)^2 dphi; its
    # exterior derivative is nonzero. Meridional reflection is essential.
    th = sp.Symbol("theta", real=True)
    assert sp.diff(sp.sin(th)**2, th) != 0
