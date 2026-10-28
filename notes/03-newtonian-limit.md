# 03 · The Newtonian limit

Any theory of gravity must reduce to Newton for weak fields and slow matter.
This note shows how, and in doing so fixes the constant $\kappa = 8\pi G/c^4$.

## Assumptions

1. **Weak field.** $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ with $|h_{\mu\nu}| \ll 1$; keep first order in $h$.
2. **Static.** No time dependence.
3. **Slow matter.** Velocities $\ll c$, pressure $\ll \rho c^2$: dust with $T_{00} = \rho c^2$.

Write the metric as

```math
ds^2 = -\left(1 + \frac{2\Phi}{c^2}\right) c^2 dt^2 + \left(1 - \frac{2\Phi}{c^2}\right) \delta_{ij}\, dx^i dx^j .
```

## Motion: Newton's second law

For a slow particle, $dx^0/d\tau \approx c$ dominates the geodesic equation:

```math
\frac{d^2 x^i}{dt^2} \approx -c^2\, \Gamma^i{}_{00}, \qquad
\Gamma^i{}_{00} = -\tfrac12 \partial_i g_{00} = \frac{\partial_i \Phi}{c^2}
\quad\Longrightarrow\quad
\frac{d^2 x^i}{dt^2} = -\partial_i \Phi .
```

So $\Phi$ is the Newtonian potential, and only $g_{00}$ matters for slow motion.

## Field: Poisson's equation

To first order and static, $R_{00} \approx \partial_i \Gamma^i{}_{00} = \nabla^2 \Phi / c^2$.
Use the trace-reversed field equations ([01](01-the-field-equations.md)) with
dust, $T_{00} = \rho c^2$ and $T = -\rho c^2$:

```math
R_{00} = \kappa \left( T_{00} - \tfrac12 T g_{00} \right) = \tfrac12 \kappa \rho c^2 .
```

Equate the two:

```math
\nabla^2 \Phi = \tfrac12 \kappa c^4 \rho .
```

This is Poisson's equation $\nabla^2 \Phi = 4\pi G \rho$ exactly when
$\kappa = 8\pi G / c^4$. That is where the $8\pi$ comes from.

*Checked in code:* `test_newtonian_limit` in
[`tests/test_solutions.py`](../tests/test_solutions.py) expands the curvature of
this metric to first order and finds $R_{00} = \nabla^2 \Phi$ and
$G_{00} = 2 \nabla^2 \Phi$ (with $c = 1$), so $G_{00} = 8\pi G \rho$ is Poisson's equation.

## What Newton misses

- **Clocks.** $d\tau/dt = \sqrt{-g_{00}} \approx 1 + \Phi/c^2$: gravitational time dilation (GPS corrects for it).
- **Light.** For light, the spatial part $(1 - 2\Phi/c^2)$ counts as much as $g_{00}$.
  Deflection past the Sun is $4GM/(c^2 b) \approx 1.75''$ at the limb, twice the Newtonian value (Eddington 1919).
- **Perihelion precession.** Mercury's extra 43″ per century needs second order in $\Phi/c^2$.

## How weak is weak?

| Object | $\lvert\Phi\rvert/c^2 \approx GM/(Rc^2)$ |
| --- | --- |
| Earth's surface | $7 \times 10^{-10}$ |
| Sun's surface | $2 \times 10^{-6}$ |
| White dwarf | $\sim 10^{-4}$ |
| Neutron star | $\sim 0.2$ |
| Black hole horizon | $0.5$ |

## Discussion prompts

- The Newtonian limit only fixes $g_{00}$. Why does the spatial part come out as $(1 - 2\Phi/c^2)$ and not $1$?
- Where in the derivation did we use $p \ll \rho c^2$? What changes for a radiation-dominated gas?
- Is there a sense in which Newtonian gravity is "curved time" rather than a force?
