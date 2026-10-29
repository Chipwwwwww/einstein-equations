"""Metrics discussed in the notes, in geometrized units (G = c = 1).

The symbols are module-level so that tests and notebooks can refer to them,
e.g. `metrics.M` or `metrics.a`.
"""

import sympy as sp

from .tensors import Spacetime

t, r, x, y, z = sp.symbols("t r x y z", real=True)
theta, phi = sp.symbols("theta phi", real=True)
M = sp.Symbol("M", positive=True)            # mass
Lambda = sp.Symbol("Lambda", positive=True)  # cosmological constant
k = sp.Symbol("k", real=True)                # FLRW spatial curvature: -1, 0 or +1
a = sp.Function("a")(t)                      # FLRW scale factor
alpha = sp.Function("alpha")(r)              # general static spherical metric
beta = sp.Function("beta")(r)
Phi = sp.Function("Phi")(x, y, z)            # Newtonian potential
f = sp.Function("f")                         # gravitational-wave profile
epsilon = sp.Symbol("epsilon")               # bookkeeping parameter: expand to first order
Q = sp.Symbol("Q", real=True)                # electric charge (Gaussian units)
a_kerr = sp.Symbol("a", real=True)           # Kerr spin parameter, J = a M
chi = sp.Symbol("chi", real=True)            # chi = cos(theta), keeps Kerr rational
tau = sp.Symbol("tau", real=True)            # Kastor-Traschen time coordinate
H = sp.Symbol("H", real=True)                # Hubble rate, H**2 = Lambda / 3; H < 0 contracts

SPHERICAL = (t, r, theta, phi)
CARTESIAN = (t, x, y, z)


def _sphere(radius_squared):
    return radius_squared, radius_squared * sp.sin(theta) ** 2


def minkowski_spherical():
    """Flat spacetime in spherical coordinates: nonzero Christoffels, zero curvature."""
    return Spacetime(sp.diag(-1, 1, *_sphere(r**2)), SPHERICAL)


def schwarzschild():
    """Vacuum outside a spherical mass M (notes/04)."""
    h = 1 - 2 * M / r
    return Spacetime(sp.diag(-h, 1 / h, *_sphere(r**2)), SPHERICAL)


def de_sitter_static():
    """Vacuum with cosmological constant Lambda, static patch (notes/04)."""
    h = 1 - Lambda * r**2 / 3
    return Spacetime(sp.diag(-h, 1 / h, *_sphere(r**2)), SPHERICAL)


def flrw():
    """Homogeneous, isotropic cosmology with scale factor a(t) and curvature k (notes/04)."""
    return Spacetime(sp.diag(-1, a**2 / (1 - k * r**2), *_sphere(a**2 * r**2)), SPHERICAL)


def reissner_nordstrom():
    """Charged black hole: Einstein-Maxwell with A = -Q/r dt (fields.py conventions)."""
    h = 1 - 2 * M / r + Q**2 / r**2
    return Spacetime(sp.diag(-h, 1 / h, *_sphere(r**2)), SPHERICAL)


def reissner_nordstrom_potential():
    """A_a for reissner_nordstrom(): electric potential Q/r, A_t = -Q/r."""
    return [-Q / r, 0, 0, 0]


def kerr(simplify=sp.cancel):
    """Kerr in Boyer-Lindquist coordinates (t, r, chi = cos theta, phi) (notes/04).

    Using chi instead of theta makes every component a rational function, so
    sp.cancel decides identities exactly and quickly.
    """
    sigma = r**2 + a_kerr**2 * chi**2
    delta = r**2 - 2 * M * r + a_kerr**2
    s2 = 1 - chi**2  # sin^2 theta
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M * r / sigma)
    g[0, 3] = g[3, 0] = -2 * M * a_kerr * r * s2 / sigma
    g[1, 1] = sigma / delta
    g[2, 2] = sigma / s2
    g[3, 3] = (r**2 + a_kerr**2 + 2 * M * a_kerr**2 * r * s2 / sigma) * s2
    return Spacetime(g, (t, r, chi, phi), simplify=simplify)


def conformastatic(U, simplify=False):
    """ds^2 = -U^-2 dtau^2 + U^2 (dx^2 + dy^2 + dz^2) for any function U(tau, x, y, z).

    U = sum m_i / r_i is Majumdar-Papapetrou; U = H tau + sum m_i / r_i is
    Kastor-Traschen (notes/08). Pair with conformastatic_potential(U).
    """
    return Spacetime(sp.diag(-U**-2, U**2, U**2, U**2), (tau, x, y, z), simplify=simplify)


def conformastatic_potential(U):
    """A = U^-1 dtau: the electric potential that goes with conformastatic(U)."""
    return [1 / U, 0, 0, 0]


def kastor_traschen_U(masses, centers):
    """U = H tau + sum_i m_i / |x - x_i| (Brill-Horowitz-Kastor-Traschen 1994 form)."""
    return H * tau + sum(
        m / sp.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        for m, (cx, cy, cz) in zip(masses, centers)
    )


def static_spherical():
    """Most general static, spherically symmetric metric: two free functions of r."""
    return Spacetime(sp.diag(-sp.exp(2 * alpha), sp.exp(2 * beta), *_sphere(r**2)), SPHERICAL)


def weak_field():
    """Static weak field with potential Phi(x, y, z) (notes/03). Expand in epsilon."""
    return Spacetime(
        sp.diag(-(1 + 2 * epsilon * Phi), *[1 - 2 * epsilon * Phi] * 3),
        CARTESIAN,
        simplify=False,
    )


def plane_wave(polarization):
    """Linearized gravitational wave along z with profile f(t - z) (notes/05). Expand in epsilon."""
    h = epsilon * f(t - z)
    g = sp.diag(-1, 1, 1, 1)
    if polarization == "plus":
        g[1, 1] += h
        g[2, 2] -= h
    elif polarization == "cross":
        g[1, 2] = g[2, 1] = h
    else:
        raise ValueError(f"polarization must be 'plus' or 'cross', not {polarization!r}")
    return Spacetime(g, CARTESIAN, simplify=False)


# Metrics that `python -m einstein` can print: exact, fully simplified ones.
REGISTRY = {
    "minkowski": minkowski_spherical,
    "schwarzschild": schwarzschild,
    "de-sitter": de_sitter_static,
    "flrw": flrw,
    "reissner-nordstrom": reissner_nordstrom,
    "static-spherical": static_spherical,
}
