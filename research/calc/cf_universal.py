"""Universal identities for conformally flat, time-symmetric vacuum data (notes/19 §3).

Claim checked: for EVERY conformally flat scalar-flat 3-metric, the geodesic-slicing development obeys
    (lap)  Delta Ric = 3 Ric.Ric - |Ric|^2 gamma_0
    (t4)   gamma_2   = (7 Ric.Ric - |Ric|^2 gamma_0) / 12
    (t6)   gamma_3   = -(1/15)|Ric|^2 Ric - (41/540) tr(Ric^3) gamma_0 + (1/60) nabla_i Ric . nabla_j Ric
                       - (1/180)|nabla Ric|^2 gamma_0 - (1/180) Ric^kl nabla_k nabla_l Ric_ij
(bl_structure.cf_identities computes the residuals).

Reduction to a finite exact computation (proof in notes/19 §3):
 1. Locally gamma_0 = phi^4 delta with Delta phi = 0 (R = -8 phi^-5 Delta phi).
 2. Both sides are natural tensors (the development in Gauss gauge is diffeomorphism-natural), so it
    suffices to check them at x = 0 after an isometry. Translations, rotations and special conformal
    maps x -> (x + b|x|^2)/(1 + 2b.x + |b|^2|x|^2) (with the Kelvin weight |DF|^(1/2) for phi) preserve the
    form phi^4 delta and harmonicity; a homothety rescales both sides by the same power. Hence one may
    assume phi(0) = 1, grad phi(0) = 0 (b = grad phi(0)/phi(0)), and a diagonal Hessian.
 3. Then phi = 1 + q2 + q3 + ... with q_k harmonic homogeneous of degree k, q2 = a(x^2 - z^2) + b(y^2 - z^2).
    The value at 0 of a residual of derivative order w (w = 4 for lap and t4, w = 6 for t6) is a
    polynomial in the coefficients, weighted-homogeneous of weight w (coefficient of q_k has weight k),
    because the recursion only differentiates and multiplies. So:
      weight 4: P = P4(q4) + P22(a, b);   weight 6: P = P6(q6) + P24(q2, q4) + P33(q3) + P222(a, b).
    (q5 cannot occur: no weight-1 coefficient remains.)
 4. P vanishes identically iff it vanishes on: each basis vector of q4 (q6); q2 = e_a with each basis
    vector of q4; q3 = e_i and e_i + e_j (polarisation); (a, b) in {(1,0), (0,1), (1,1)} for the binary
    quadratic, and {(1,0), (0,1), (1,1), (1,-1)} for the binary cubic (a nonzero binary cubic has at
    most three projective zeros). All evaluations are exact rational arithmetic.

Run: python -m research.calc.cf_universal   (writes research/calc/data/cf_universal.json)
"""

import contextlib
import json
import os
import time

import sympy as sp
from flint import fmpq

from research.calc import bl_metric as bm
from research.calc import bl_structure as st

X, Y, Z = sp.symbols("x y z")


def harmonic_basis(k):
    """A rational basis of the homogeneous harmonic polynomials of degree k (dimension 2k + 1), as lists
    of ((a, b, c), coefficient)."""
    mons = [(a, b, k - a - b) for a in range(k, -1, -1) for b in range(k - a, -1, -1)]
    lap_rows = {}
    cols = []
    for (a, b, c) in mons:
        p = X ** a * Y ** b * Z ** c
        lp = sp.Poly(sp.diff(p, X, 2) + sp.diff(p, Y, 2) + sp.diff(p, Z, 2), X, Y, Z) if k >= 2 else None
        cols.append(dict(lp.terms()) if lp is not None else {})
        if lp is not None:
            for m in lp.monoms():
                lap_rows.setdefault(m, len(lap_rows))
    M = sp.zeros(max(len(lap_rows), 1), len(mons))
    for j, col in enumerate(cols):
        for m, v in col.items():
            M[lap_rows[m], j] = v
    ns = M.nullspace() if k >= 2 else [sp.eye(len(mons))[:, j] for j in range(len(mons))]
    out = []
    for v in ns:
        den = sp.ilcm(*[sp.fraction(e)[1] for e in v]) if len(v) else 1
        out.append([(mons[j], fmpq(int(v[j] * den), 1)) for j in range(len(mons)) if v[j] != 0])
    assert len(out) == 2 * k + 1, (k, len(out))
    return out


def phi_fn(terms):
    """psi_fn for phi = 1 + sum of the given ((a, b, c), coefficient) terms."""
    def f(J):
        jet = {0: fmpq(1)}
        for (a, b, c), v in terms:
            if a + b + c <= J.W:
                i = J.index[(a, b, c, 0)]
                jet[i] = jet.get(i, fmpq(0)) + v
        return jet
    return f


@contextlib.contextmanager
def custom_data(psi_fn):
    """Temporarily make bl_metric.series_jets use gamma_0 = phi^4 delta for the given phi (expansion at 0)."""
    orig = bm.psi_jet
    bm.psi_jet = lambda J, x0, cs, one: psi_fn(J)
    try:
        yield
    finally:
        bm.psi_jet = orig


def residuals(terms):
    """Residuals (lap, t4, t6) at 0 for phi = 1 + terms: dicts over PAIRS of exact rationals."""
    with custom_data(phi_fn(terms)):
        cv = st.curvature_values((fmpq(0),) * 3, None, fmpq(1), centres=[])
    return st.cf_identities(cv, fmpq(1)), cv


def q2(a, b):
    return [((2, 0, 0), fmpq(a)), ((0, 2, 0), fmpq(b)), ((0, 0, 2), fmpq(-a - b))]


def scale(basis_vec, c=1):
    return [(m, v * c) for m, v in basis_vec]


def evaluation_set():
    """The evaluation points of step 4: (label, weight, terms)."""
    H = {k: harmonic_basis(k) for k in (3, 4, 6)}
    pts = []
    for (a, b) in ((1, 0), (0, 1), (1, 1), (1, -1)):
        pts.append((f"q2=({a},{b})", (4, 6), q2(a, b)))
    for j, v in enumerate(H[4]):
        pts.append((f"q4=e{j}", (4,), v))
        for (a, b) in ((1, 0), (0, 1)):
            pts.append((f"q2=({a},{b}),q4=e{j}", (6,), q2(a, b) + v))
    for i, v in enumerate(H[3]):
        pts.append((f"q3=e{i}", (6,), v))
        for j in range(i + 1, len(H[3])):
            pts.append((f"q3=e{i}+e{j}", (6,), v + H[3][j]))
    for j, v in enumerate(H[6]):
        pts.append((f"q6=e{j}", (6,), v))
    return pts


def main(out=os.path.join(bm.DATA, "cf_universal.json")):
    t0 = time.time()
    res = {"points": {}, "all_zero": {"lap": True, "t4": True, "t6": True}}
    for label, weights, terms in evaluation_set():
        r, cv = residuals(terms)
        z = {k: all(v == 0 for v in r[k].values()) for k in ("lap", "t4", "t6")}
        nontriv = any(cv["g3"][pq] != 0 for pq in bm.PAIRS)
        res["points"][label] = {"zero": z, "gamma3_nonzero": nontriv}
        for k in z:
            if (k in ("lap", "t4") and 4 in weights) or (k == "t6" and 6 in weights):
                res["all_zero"][k] &= z[k]
        print(label, z, nontriv, f"{time.time() - t0:.0f}s", flush=True)
    res["evaluations"] = len(res["points"])
    res["seconds"] = time.time() - t0
    res["replay"] = "python -m research.calc.cf_universal"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    main()
