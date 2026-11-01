"""Matter fields on a Spacetime: the electromagnetic field (Gaussian units, G = c = 1).

Conventions follow notes/00-conventions.md. With a potential A_a,

    F_ab = d_a A_b - d_b A_a
    T_ab = (1 / 4 pi) (F_ac F_b^c - 1/4 g_ab F_cd F^cd)

so that the Einstein-Maxwell-Lambda equations read

    G_ab + Lambda g_ab = 8 pi T_ab,      nabla_a F^ab = 0.

T_ab is the Maxwell stress tensor in Gaussian units: for a point charge Q in flat
space it gives energy density E^2 / 8 pi with E = Q / r^2, the convention in
which Reissner-Nordstrom has g_tt = -(1 - 2M/r + Q^2/r^2).
"""

from itertools import product

import sympy as sp


def field_strength(potential, coords):
    """F_ab = d_a A_b - d_b A_a for a covector potential A_a given as a list."""
    n = len(coords)
    return sp.Matrix(
        n, n, lambda a, b: sp.diff(potential[b], coords[a]) - sp.diff(potential[a], coords[b])
    )


def maxwell_stress_energy(st, F):
    """T_ab of the field strength F_ab (a SymPy matrix) on the Spacetime st."""
    n, g, gi = st.dim, st.metric, st.inverse
    F_mixed = F * gi  # F_a^c = F_ad g^dc
    F_up = gi * F * gi  # F^ab
    invariant = sum(F[a, b] * F_up[a, b] for a, b in product(range(n), repeat=2))
    return sp.Matrix(
        n,
        n,
        lambda a, b: (
            sum(F[a, c] * F_mixed[b, c] for c in range(n)) - g[a, b] * invariant / 4
        )
        / (4 * sp.pi),
    )


def maxwell_divergence(st, F):
    """nabla_a F^ab = (1 / sqrt|g|) d_a (sqrt|g| F^ab), one entry per b.

    Returned without the 1/sqrt|g| factor removed, i.e. as the full divergence; it
    vanishes exactly where the source-free Maxwell equations hold.
    """
    n, gi, x = st.dim, st.inverse, st.coords
    sqrt_g = sp.sqrt(-st.metric.det())
    F_up = gi * F * gi
    return [
        sum(sp.diff(sqrt_g * F_up[a, b], x[a]) for a in range(n)) / sqrt_g for b in range(n)
    ]
