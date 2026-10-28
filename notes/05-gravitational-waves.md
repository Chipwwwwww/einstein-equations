# 05 · Gravitational waves

Far from sources the metric is nearly flat, and the field equations become a
wave equation. This note derives it, shows the two polarizations, and gives the
formula that sets how strong the waves are.

## Linearize

Write $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ with $|h_{\mu\nu}| \ll 1$ and keep
first order in $h$. Define the trace-reversed perturbation

```math
\bar h_{\mu\nu} = h_{\mu\nu} - \tfrac12 \eta_{\mu\nu} h, \qquad h = \eta^{\mu\nu} h_{\mu\nu} .
```

Use the coordinate freedom to impose the Lorenz (harmonic) gauge
$\partial^\mu \bar h_{\mu\nu} = 0$. Then the linearized Einstein tensor is
$G_{\mu\nu} = -\tfrac12 \Box \bar h_{\mu\nu}$, and the field equations read

```math
\Box \bar h_{\mu\nu} = -\frac{16\pi G}{c^4}\, T_{\mu\nu},
\qquad \Box = -\frac{1}{c^2}\partial_t^2 + \nabla^2 .
```

In vacuum, $\Box \bar h_{\mu\nu} = 0$: waves moving at exactly $c$.
GW170817 (a neutron-star merger seen in gravitational waves and gamma rays)
confirmed this to about one part in $10^{15}$.

## Two polarizations

Residual gauge freedom allows the transverse-traceless (TT) gauge:
$h_{0\mu} = 0$, $h^i{}_i = 0$, $\partial^i h_{ij} = 0$. A wave moving along $z$ then has only two components:

```math
h_{ij}^{TT} =
\begin{pmatrix}
h_+ & h_\times & 0 \\
h_\times & -h_+ & 0 \\
0 & 0 & 0
\end{pmatrix}
f(t - z/c) .
```

A ring of free particles in the $xy$ plane is stretched along $x$ and squeezed
along $y$ by $h_+$, and the same rotated by 45° by $h_\times$. These are the "2
degrees of freedom" counted in [01](01-the-field-equations.md).

*Checked in code:* `test_linearized_plane_wave_is_vacuum` verifies that both
polarizations, with an arbitrary profile $f$, have vanishing Ricci tensor to
first order.

## How strong? The quadrupole formula

Far from a slow source at distance $r$:

```math
h_{ij}^{TT}(t) = \frac{2G}{c^4 r}\, \ddot I_{ij}^{TT}(t - r/c),
\qquad I_{ij} = \int \rho\, x_i x_j\, d^3x ,
```

and the radiated power (Einstein 1918, fixing a factor of 2 in his 1916 paper) is

```math
P = \frac{G}{5c^5} \left\langle \dddot Q_{ij}\, \dddot Q_{ij} \right\rangle,
\qquad Q_{ij} = I_{ij} - \tfrac13 \delta_{ij} I_{kk} .
```

No dipole radiation: momentum conservation kills it. The factor $G/c^5$ is why
only astrophysical sources (merging black holes and neutron stars) are
detectable. GW150914 reached a peak strain of about $10^{-21}$, a length change
much smaller than a proton over LIGO's 4 km arms.

## Discussion prompts

- Gravitational waves carry energy, yet there is no local gravitational energy density. How do we make sense of the energy they carry?
- Why is TT gauge allowed, and which statements about $h_{ij}^{TT}$ are gauge-invariant?
- The Hulse–Taylor pulsar showed orbital decay matching the quadrupole formula decades before LIGO. Which was the stronger evidence for gravitational waves?
