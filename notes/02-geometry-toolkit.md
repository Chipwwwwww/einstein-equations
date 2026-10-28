# 02 · The geometry toolkit

Everything on the left side of the field equations is built from the metric in
a fixed chain. Each step below names the attribute in
[`einstein/tensors.py`](../einstein/tensors.py) that computes it.

```text
g_{μν}  →  Γ^λ_{μν}  →  R^ρ_{σμν}  →  R_{μν}, R  →  G_{μν}
metric     christoffel   riemann       ricci,       einstein
                                       ricci_scalar
```

## Metric `metric`, `inverse`

$g_{\mu\nu}$ turns coordinate differences into physical intervals,
$ds^2 = g_{\mu\nu} dx^\mu dx^\nu$. Its inverse $g^{\mu\nu}$ raises indices.

## Christoffel symbols `christoffel`

```math
\Gamma^\lambda{}_{\mu\nu} = \tfrac12\, g^{\lambda\sigma}
  \left( \partial_\mu g_{\sigma\nu} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu} \right)
```

They are **not** a tensor. They vanish at any chosen point in suitable
coordinates (the equivalence principle in math form) and are nonzero in flat
space in curvilinear coordinates. Free fall follows the geodesic equation

```math
\frac{d^2 x^\lambda}{d\tau^2} + \Gamma^\lambda{}_{\mu\nu} \frac{dx^\mu}{d\tau} \frac{dx^\nu}{d\tau} = 0 .
```

*Checked in code:* flat space in spherical coordinates has nonzero $\Gamma$ but
zero Riemann (`test_flat_space_in_spherical_coordinates`).

## Riemann tensor `riemann`

```math
R^\rho{}_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho{}_{\nu\sigma} - \partial_\nu \Gamma^\rho{}_{\mu\sigma}
  + \Gamma^\rho{}_{\mu\lambda} \Gamma^\lambda{}_{\nu\sigma} - \Gamma^\rho{}_{\nu\lambda} \Gamma^\lambda{}_{\mu\sigma}
```

Two ways to read it:

- **Parallel transport.** Carry a vector around a small loop; it comes back rotated by an amount proportional to $R^\rho{}_{\sigma\mu\nu}$.
- **Tidal acceleration.** Two nearby geodesics with separation $\xi^\mu$ and velocity $u^\mu$ obey the geodesic deviation equation

```math
\frac{D^2 \xi^\mu}{d\tau^2} = R^\mu{}_{\nu\rho\sigma}\, u^\nu u^\rho \xi^\sigma .
```

Symmetries $R_{\rho\sigma\mu\nu} = -R_{\sigma\rho\mu\nu} = -R_{\rho\sigma\nu\mu} = R_{\mu\nu\rho\sigma}$
and $R_{\rho[\sigma\mu\nu]} = 0$ cut $4^4 = 256$ components to
$n^2(n^2-1)/12 = 20$ in four dimensions: 10 in Ricci, 10 in the Weyl tensor.

## Kretschmann scalar `kretschmann`

$K = R_{\rho\sigma\mu\nu} R^{\rho\sigma\mu\nu}$ is coordinate-independent. If it
blows up, the curvature really is singular; if it stays finite, a divergence in
the metric may only be a bad coordinate choice. For Schwarzschild
$K = 48 G^2 M^2 / (c^4 r^6)$: finite at the horizon, infinite at $r = 0$.

## Ricci tensor and scalar `ricci`, `ricci_scalar`

```math
R_{\mu\nu} = R^\lambda{}_{\mu\lambda\nu}, \qquad R = g^{\mu\nu} R_{\mu\nu} .
```

Ricci measures how a small ball of freely falling test particles changes
**volume**; Weyl measures how it changes **shape**.

## Einstein tensor `einstein`

```math
G_{\mu\nu} = R_{\mu\nu} - \tfrac12 R\, g_{\mu\nu}
```

## Bianchi identities `divergence`

The second Bianchi identity

```math
\nabla_\lambda R_{\rho\sigma\mu\nu} + \nabla_\rho R_{\sigma\lambda\mu\nu} + \nabla_\sigma R_{\lambda\rho\mu\nu} = 0
```

contracted twice gives

```math
\nabla^\mu G_{\mu\nu} = 0
```

for **any** metric. [`tests/test_identities.py`](../tests/test_identities.py)
checks this on a general static spherically symmetric metric with two arbitrary
functions, and on FLRW with an arbitrary scale factor.

## Discussion prompts

- The Christoffel symbols are not a tensor, yet the geodesic equation is coordinate-independent. How?
- Why is $G_{\mu\nu}$, and not Riemann itself, the thing that equals matter?
- In 2+1 dimensions Weyl vanishes identically. What does that imply for gravity in a 3-dimensional spacetime?
