# 00 · Conventions

Most disagreements about "the" Einstein equations are really disagreements
about signs. Fix them once here; every note and every line of code follows them.

## Signature and coordinates

Metric signature $(-,+,+,+)$. In flat spacetime

```math
ds^2 = \eta_{\mu\nu}\,dx^\mu dx^\nu = -c^2 dt^2 + dx^2 + dy^2 + dz^2 ,
\qquad x^0 = ct .
```

Timelike intervals have $ds^2 < 0$; proper time is $c^2 d\tau^2 = -ds^2$.

## Indices

- Greek $\mu, \nu, \ldots \in \{0,1,2,3\}$ run over spacetime; Latin $i, j, \ldots \in \{1,2,3\}$ over space.
- Repeated upper/lower indices are summed (Einstein summation).
- $\partial_\mu = \partial / \partial x^\mu$; $\nabla_\mu$ is the Levi-Civita covariant derivative.

## Curvature (MTW / Carroll)

```math
\Gamma^\lambda{}_{\mu\nu} = \tfrac12\, g^{\lambda\sigma}
  \left( \partial_\mu g_{\sigma\nu} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu} \right)
```

```math
R^\rho{}_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho{}_{\nu\sigma} - \partial_\nu \Gamma^\rho{}_{\mu\sigma}
  + \Gamma^\rho{}_{\mu\lambda} \Gamma^\lambda{}_{\nu\sigma} - \Gamma^\rho{}_{\nu\lambda} \Gamma^\lambda{}_{\mu\sigma}
```

```math
R_{\mu\nu} = R^\lambda{}_{\mu\lambda\nu}, \qquad
R = g^{\mu\nu} R_{\mu\nu}, \qquad
G_{\mu\nu} = R_{\mu\nu} - \tfrac12 R\, g_{\mu\nu}
```

With these choices a sphere has $R > 0$, and the field equations carry a
**plus** sign: $G_{\mu\nu} + \Lambda g_{\mu\nu} = +8\pi G\, T_{\mu\nu}/c^4$.

The code in [`einstein/tensors.py`](../einstein/tensors.py) implements exactly
these formulas; its docstring repeats them.

## Units

- Notes keep $G$ and $c$ where the physics needs them.
- The code uses **geometrized units** $G = c = 1$. To restore units, a mass $M$
  becomes the length $GM/c^2$, time $t$ becomes $ct$, and the cosmological
  constant $\Lambda$ has units of length$^{-2}$ either way.
- Useful numbers: $GM_\odot/c^2 \approx 1.48\ \text{km}$,
  $8\pi G/c^4 \approx 2.08 \times 10^{-43}\ \text{N}^{-1}$.

## Comparing with other books

- Landau & Lifshitz and much of the particle-physics literature use signature
  $(+,-,-,-)$. Metric components flip sign; $G_{\mu\nu} = 8\pi G T_{\mu\nu}/c^4$
  keeps its form, but formulas with explicit $g_{\mu\nu}$ (like the $\Lambda$
  term or $T_{\mu\nu}$ of a perfect fluid) change sign.
- Some authors (e.g. Weinberg 1972) define the Riemann tensor with the opposite
  overall sign, so the field equations appear with $-8\pi G$.
- MTW's endpapers tabulate the choices of many classic texts. Check before you
  argue that a sign is wrong.

## Discussion prompts

- Which signature do you find more natural, and does the choice ever matter physically?
- Find a formula in your favourite textbook that differs from ours only by convention. What is the conversion?
