"""Curvature of a metric, computed symbolically with SymPy.

Conventions are those of notes/00-conventions.md (MTW / Carroll), with
indices written as array positions:

    Gamma[a, b, c]   = Gamma^a_{bc}
                     = 1/2 g^{ad} (d_b g_{dc} + d_c g_{db} - d_d g_{bc})
    riemann[a, b, c, d] = R^a_{bcd}
                     = d_c Gamma^a_{db} - d_d Gamma^a_{cb}
                       + Gamma^a_{ce} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{cb}
    ricci[b, d]      = R_{bd} = R^a_{bad}
    einstein[a, b]   = G_{ab} = R_{ab} - 1/2 R g_{ab}
"""

from functools import cached_property
from itertools import product

import sympy as sp


class Spacetime:
    """A metric g_{ab} in coordinates x^a, with curvature computed on first use.

    Pass simplify=False when you will post-process the expressions yourself
    (for example, expanding to first order in a small parameter): simplifying
    the full nonlinear expressions first is slow and unnecessary.
    """

    def __init__(self, metric, coords, simplify=True):
        self.metric = sp.Matrix(metric)
        self.coords = tuple(coords)
        self.dim = len(self.coords)
        if self.metric.shape != (self.dim, self.dim):
            raise ValueError(f"metric is {self.metric.shape}, expected {self.dim}x{self.dim}")
        if self.metric != self.metric.T:
            raise ValueError("metric must be symmetric")
        self._simplify = sp.simplify if simplify else (lambda expr: expr)

    @cached_property
    def inverse(self):
        """g^{ab}."""
        return self.metric.inv().applyfunc(self._simplify)

    @cached_property
    def christoffel(self):
        """Gamma^a_{bc}, symmetric in b and c."""
        n, g, gi, x = self.dim, self.metric, self.inverse, self.coords
        dg = [[[sp.diff(g[i, j], x[k]) for k in range(n)] for j in range(n)] for i in range(n)]
        gamma = sp.MutableDenseNDimArray.zeros(n, n, n)
        for a, b in product(range(n), repeat=2):
            for c in range(b, n):
                value = sum(
                    gi[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) for d in range(n)
                ) / 2
                gamma[a, b, c] = gamma[a, c, b] = self._simplify(value)
        return gamma.as_immutable()

    @cached_property
    def riemann(self):
        """R^a_{bcd}, antisymmetric in c and d."""
        n, x, gamma = self.dim, self.coords, self.christoffel
        riem = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b in product(range(n), repeat=2):
            for c in range(n):
                for d in range(c + 1, n):
                    value = sp.diff(gamma[a, d, b], x[c]) - sp.diff(gamma[a, c, b], x[d]) + sum(
                        gamma[a, c, e] * gamma[e, d, b] - gamma[a, d, e] * gamma[e, c, b]
                        for e in range(n)
                    )
                    value = self._simplify(value)
                    riem[a, b, c, d] = value
                    riem[a, b, d, c] = -value
        return riem.as_immutable()

    @cached_property
    def ricci(self):
        """R_{bd} = R^a_{bad}."""
        n, riem = self.dim, self.riemann
        return sp.Matrix(n, n, lambda b, d: self._simplify(sum(riem[a, b, a, d] for a in range(n))))

    @cached_property
    def ricci_scalar(self):
        """R = g^{ab} R_{ab}."""
        n, gi, ric = self.dim, self.inverse, self.ricci
        return self._simplify(sum(gi[a, b] * ric[a, b] for a, b in product(range(n), repeat=2)))

    @cached_property
    def einstein(self):
        """G_{ab} = R_{ab} - 1/2 R g_{ab}."""
        return (self.ricci - self.ricci_scalar * self.metric / 2).applyfunc(self._simplify)

    @cached_property
    def einstein_mixed(self):
        """G^a_b = g^{ac} G_{cb}."""
        return (self.inverse * self.einstein).applyfunc(self._simplify)

    @cached_property
    def kretschmann(self):
        """K = R_{abcd} R^{abcd}, a curvature invariant independent of coordinates."""
        n = self.dim
        lowered = _contract(self.riemann, self.metric, 0)  # R_{abcd}
        raised = self.riemann
        for position in (1, 2, 3):
            raised = _contract(raised, self.inverse, position)  # R^{abcd}
        return self._simplify(
            sum(lowered[i] * raised[i] for i in product(range(n), repeat=4))
        )

    def divergence(self, tensor):
        """nabla_a T^a_b of a mixed tensor given as a matrix T[a, b] = T^a_b.

        Returns the n components (one per b) as a list.
        """
        n, x, gamma = self.dim, self.coords, self.christoffel
        return [
            self._simplify(
                sum(sp.diff(tensor[a, b], x[a]) for a in range(n))
                + sum(
                    gamma[a, a, e] * tensor[e, b] - gamma[e, a, b] * tensor[a, e]
                    for a, e in product(range(n), repeat=2)
                )
            )
            for b in range(n)
        ]


def _contract(array, matrix, position):
    """Contract `matrix` into one index of `array`: out[..i..] = matrix[i, j] array[..j..]."""
    n = array.shape[0]
    out = sp.MutableDenseNDimArray.zeros(*array.shape)
    for index in product(range(n), repeat=array.rank()):
        i = index[position]
        out[index] = sum(
            matrix[i, j] * array[index[:position] + (j,) + index[position + 1:]]
            for j in range(n)
            if matrix[i, j] != 0
        )
    return out


def nonzero(array):
    """(index, value) pairs of the nonzero components of a SymPy array or matrix."""
    shape = array.shape
    return [(i, array[i]) for i in product(*(range(s) for s in shape)) if array[i] != 0]
