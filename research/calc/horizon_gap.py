"""All-mode Jacobi gap at the critical Brill--Lindquist common minimal sphere.

See notes/20-horizon-spectral-gap.md for the analytic proof. This replays the
fold's Krawczyk test and checks geometric inequalities on EVERY validated ODE
step box, including the pole collar. No eigenvalue discretisation is used.

Run: python -m research.calc.horizon_gap [--output PATH] [--fine]
The --fine run halves every maximal regular step and the near-axis ratio.
Existing fold code supplies the validated ODE integrator, not stored numbers.
"""

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import flint
from flint import arb, ctx

from research.calc import fold_proof as fp


def geometric_weights(rho, z, sa, ca, rho_z, z_z, d):
    """Ground-state energy/mass weights relative to the unit round sphere.

    rho_z,z_z are derivatives with respect to the initial pole height, at any
    of the integrator's parameters: tangential reparametrisations cancel in j.
    phi = psi**2*j is the PHYSICAL normal Jacobi field, not just j.
    """
    r2 = rho**2 + z**2
    psi = 1 + arb(1) / (4 * (rho**2 + (z - d/2)**2).sqrt()) \
        + arb(1) / (4 * (rho**2 + (z + d/2)**2).sqrt())
    c = (-rho*sa + z*ca) / r2.sqrt()
    phi = psi**2 * (-sa*rho_z + ca*z_z)
    if not all(v.is_finite() for v in (psi, c, phi)):
        raise ArithmeticError("nonfinite geometric quantity")
    if not (c > 0 and phi > 0):
        raise ArithmeticError("star-shapedness or positive Jacobi field not certified")
    return {"c": c, "phi": phi,
            "energy_weight": phi**2*c,
            "mass_weight": phi**2*psi**4*r2/c,
            "axial_radius_squared": psi**4*rho**2}


class SpectralChecks(fp.Checks):
    """Enclose global extrema without converting a proof bound to a float."""

    def __init__(self, d):
        super().__init__()
        self.d = d
        self.hulls = {}
        self.counts = {"axis": 0, "mu": 0, "alpha": 0, "tau": 0}

    def observe(self, phase, rho, z, sa, ca, rho_z, z_z):
        if phase != "axis" and not (rho > 0 and sa < 0):
            raise ArithmeticError("regular embedded meridian not certified")
        if phase in ("mu", "alpha") and not z > 0:
            raise ArithmeticError("profile may cross equator before final segment")
        values = geometric_weights(rho, z, sa, ca, rho_z, z_z, self.d)
        for key, val in values.items():
            self.hulls[key] = self.hulls[key].union(val) if key in self.hulls else val
        self.counts[phase] += 1

    def axis(self, z0):
        zero, one = arb(0), arb(1)
        zj = fp.J([z0, one, zero, zero, zero, zero])
        dj = fp.J([self.d, zero, one, zero, zero, zero])
        A, Q, Z = fp.axis_enclosure(zj, dj, fp.S1)
        s = arb(fp.S1/2, fp.S1/2)  # the entire closed collar, including s=0
        if not (A.c[0] < 0 and 1 + Q.c[0]*s**2 > 0):
            raise ArithmeticError("regular axis collar not certified")
        a = A.c[0]*s
        self.observe("axis", s + Q.c[0]*s**3, z0 + Z.c[0]*s**2,
                     a.sin(), a.cos(), Q.c[1]*s**3, 1 + Z.c[1]*s**2)

    def mu(self, B):
        super().mu(B)
        _, rho, z, mu = B
        sa = mu.c[0]/rho.c[0]
        self.observe("mu", rho.c[0], z.c[0], sa, (1-sa**2).sqrt(), rho.c[1], z.c[1])

    def alpha(self, B):
        super().alpha(B)
        _, rho, z, a = B
        self.observe("alpha", rho.c[0], z.c[0], a.c[0].sin(), a.c[0].cos(), rho.c[1], z.c[1])

    def tau(self, zsw):
        if not zsw.c[0] > 0:
            raise ArithmeticError("height switch must be above equator")
        parent = super().tau(zsw)
        def check(B):
            parent(B)
            tau, rho, a = B
            self.observe("tau", rho.c[0], zsw.c[0]*(1-tau.c[0]),
                         a.c[0].sin(), a.c[0].cos(), rho.c[1], zsw.c[1]*(1-tau.c[0]))
        return check


def krawczyk_from_jets(Fp, Fb, centre, radius):
    """Krawczyk with a fixed, exactly represented dyadic preconditioner."""
    D = [[Fb.c[1], Fb.c[2]], [2*Fb.c[4], Fb.c[5]]]
    m = [[float(v.mid()) for v in row] for row in D]
    det = m[0][0]*m[1][1] - m[0][1]*m[1][0]
    C = [[arb(m[1][1]/det), arb(-m[0][1]/det)],
         [arb(-m[1][0]/det), arb(m[0][0]/det)]]
    if (C[0][0]*C[1][1]-C[0][1]*C[1][0]).contains(0):
        raise ArithmeticError("singular preconditioner")
    X = [arb(v, radius) for v in centre]
    rad, K = arb(0, radius), []
    contraction_rows = []
    for i in range(2):
        val = arb(centre[i]) - C[i][0]*Fp.c[0] - C[i][1]*Fp.c[1]
        row_norm = arb(0)
        for j in range(2):
            e = (1 if i == j else 0) - C[i][0]*D[0][j] - C[i][1]*D[1][j]
            val += e*rad
            row_norm += e.abs_upper()
        K.append(val)
        contraction_rows.append(row_norm)
    verified = all(X[i].contains_interior(K[i]) for i in range(2)) \
        and all(row < 1 for row in contraction_rows)
    return verified, K, contraction_rows


def certify(fine=False):
    start = time.monotonic()
    source_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in (Path(__file__), Path(fp.__file__))}
    old_prec = ctx.prec
    ctx.prec = 256
    try:
        schedules = {}
        if fine:
            schedules = {"mu_schedule": [(end, L, ratio*2, h/2)
                                          for end, L, ratio, h in fp.MU_SCHEDULE],
                         "alpha_schedule": (fp.ALPHA_SCHEDULE[0], fp.ALPHA_SCHEDULE[1]/2),
                         "tau_schedule": (fp.TAU_SCHEDULE[0], fp.TAU_SCHEDULE[1]/2)}
        z, d = fp.FOLD_CENTRE
        print("Validated point shooting ...", flush=True)
        Fp, steps_p = fp.shoot(arb(z), arb(d), **schedules)
        chk = SpectralChecks(arb(d, fp.BOX_RADIUS))
        chk.axis(arb(z, fp.BOX_RADIUS))
        print("Validated box shooting and full-surface weights ...", flush=True)
        Fb, steps_b = fp.shoot(arb(z, fp.BOX_RADIUS), chk.d, checks=chk, **schedules)
        ok, K, contraction = krawczyk_from_jets(Fp, Fb, fp.FOLD_CENTRE, fp.BOX_RADIUS)
        if not ok:
            raise ArithmeticError("fold root not certified")
        if not (Fb.c[4] > 0 and Fb.c[2] > 0 and Fb.c[3] < 0):
            raise ArithmeticError("fold transversality not certified")
        if not all(chk.counts[k] > 0 for k in chk.counts):
            raise ArithmeticError("missing profile phase")
        h = chk.hulls
        # These are exact rational thresholds. Their truth is tested on Arb
        # enclosures, NOT on the rounded diagnostic decimals in the report.
        assertions = {
            "energy_weight_ge_24": bool(h["energy_weight"] > 24),
            "mass_weight_le_768": bool(h["mass_weight"] < 768),
            "axial_radius_squared_le_2": bool(h["axial_radius_squared"] < 2),
            "positive_phi": bool(h["phi"] > 0),
            "positive_c": bool(h["c"] > 0),
        }
        result = {
            "kind": "computer-assisted; analytic ground-state comparison plus validated fold ODE",
            "precision_bits": ctx.prec, "series_cap": ctx.cap, "fine_schedule": fine,
            "python": platform.python_version(), "python_flint": flint.__version__,
            "centre_exact_binary64_hex": [v.hex() for v in fp.FOLD_CENTRE],
            "box_radius_exact_binary64_hex": fp.BOX_RADIUS.hex(),
            "krawczyk_verified": ok,
            "krawczyk_image": [v.str(40, radius=True) for v in K],
            "contraction_rows": [v.str(30, radius=True) for v in contraction],
            "point_steps": steps_p, "box_steps": steps_b, "phase_counts": chk.counts,
            "weight_enclosures": {k: {"lower_decimal_diagnostic": v.lower().str(30, radius=False),
                                      "upper_decimal_diagnostic": v.upper().str(30, radius=False)}
                                  for k, v in h.items()},
            "assertions": assertions,
            "fold_jet": {k: v.str(30, radius=True) for k,v in zip(fp.NAMES, Fb.c)},
            "elapsed_seconds": time.monotonic()-start,
            "source_sha256": source_hashes,
        }
        if all(assertions.values()):
            result["proved_bounds_M_ADM_1"] = {
                "lambda_gap_ge": "1/16", "inverse_L2_norm_le": "16",
                "nonaxisymmetric_m_gap_ge": "m^2/2",
                "axisymmetric_equator_even_gap_ge": "3/16",
            }
        else:
            result["proved_bounds_M_ADM_1"] = None
        return result
    finally:
        ctx.prec = old_prec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fine", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = certify(args.fine)
    formatted = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(formatted)
    print(formatted)
    if result["proved_bounds_M_ADM_1"] is None:
        raise SystemExit("Exact target bounds not verified; inspect diagnostic enclosures.")


if __name__ == "__main__":
    main()
