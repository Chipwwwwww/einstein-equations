# 01 · The field equations

```math
G_{\mu\nu} + \Lambda g_{\mu\nu} = \kappa\, T_{\mu\nu},
\qquad \kappa = \frac{8\pi G}{c^4} \approx 2.08 \times 10^{-43}\ \text{N}^{-1}
```

Einstein presented this form (without $\Lambda$) to the Prussian Academy on
25 November 1915. He added $\Lambda$ in 1917 to allow a static universe.

## Reading it term by term

| Symbol | What it is | Built from |
| --- | --- | --- |
| $g_{\mu\nu}$ | the metric: lengths, times, causal structure | the unknown, 10 functions of 4 coordinates |
| $G_{\mu\nu}$ | Einstein tensor: a particular average of curvature | $g$ and its first and second derivatives |
| $\Lambda$ | cosmological constant | a number, units of length$^{-2}$ |
| $T_{\mu\nu}$ | stress-energy tensor of matter and radiation | $T_{00}$ energy density, $T_{0i}$ momentum density / energy flux, $T_{ij}$ stress |

As differential equations for $g_{\mu\nu}$ they are second order, quasi-linear
and highly nonlinear. That nonlinearity is gravity gravitating: the field's own
energy curves spacetime too.

Wheeler's summary: "Spacetime tells matter how to move; matter tells spacetime how to curve."

## Ten equations, two degrees of freedom

$G_{\mu\nu}$ and $T_{\mu\nu}$ are symmetric, so there are **10** component equations. They are not independent:

1. **Bianchi identity.** $\nabla^\mu G_{\mu\nu} = 0$ holds for *every* metric (see
   [02](02-geometry-toolkit.md)). That is 4 differential relations among the 10 equations.
2. **Coordinate freedom.** Any 4 functions $x^\mu \to x'^\mu(x)$ give a physically
   identical metric, so 4 of the 10 metric functions are gauge.
3. **Constraints.** The 4 components $G_{0\mu}$ contain no second time derivatives.
   They constrain the initial data (Hamiltonian and momentum constraints) rather than evolve it.

That leaves $10 - 4 - 4 = 2$ physical degrees of freedom per point: the two
polarizations of a gravitational wave ([05](05-gravitational-waves.md)).

## Conservation comes for free

Take the divergence of both sides. $\nabla^\mu G_{\mu\nu} = 0$ identically and
$\nabla^\mu g_{\mu\nu} = 0$, so

```math
\nabla^\mu T_{\mu\nu} = 0 .
```

Local energy-momentum conservation is not an extra assumption; the field
equations would be inconsistent without it. This is why the left side must be
$G_{\mu\nu}$ and not simply $R_{\mu\nu}$: $\nabla^\mu R_{\mu\nu} = \tfrac12 \partial_\nu R$,
which is not zero in general. Einstein's earlier November 1915 versions lacked
the trace term; the final version fixed this.

*Checked in code:* `test_contracted_bianchi_identity` in
[`tests/test_identities.py`](../tests/test_identities.py).

## Trace-reversed form

Contract with $g^{\mu\nu}$ (in 4 dimensions $g^{\mu\nu} g_{\mu\nu} = 4$):

```math
-R + 4\Lambda = \kappa T
\quad\Longrightarrow\quad
R = 4\Lambda - \kappa T .
```

Substitute back:

```math
R_{\mu\nu} = \kappa \left( T_{\mu\nu} - \tfrac12 T\, g_{\mu\nu} \right) + \Lambda g_{\mu\nu} .
```

In vacuum ($T_{\mu\nu} = 0$) this is $R_{\mu\nu} = \Lambda g_{\mu\nu}$, and with
$\Lambda = 0$ simply

```math
R_{\mu\nu} = 0 .
```

Ricci-flat does **not** mean flat. The Riemann tensor has 20 components; Ricci
fixes only 10. The other 10 (the Weyl tensor) carry tidal forces outside a star
and gravitational waves.

## From an action

In $c = 1$ units, the Einstein–Hilbert action plus matter is

```math
S = \frac{1}{16\pi G} \int (R - 2\Lambda) \sqrt{-g}\, d^4x \;+\; S_M ,
\qquad
T_{\mu\nu} \equiv -\frac{2}{\sqrt{-g}} \frac{\delta S_M}{\delta g^{\mu\nu}} .
```

Using $\delta\sqrt{-g} = -\tfrac12 \sqrt{-g}\, g_{\mu\nu}\, \delta g^{\mu\nu}$ and
$\delta(\sqrt{-g} R) = \sqrt{-g}\, G_{\mu\nu}\, \delta g^{\mu\nu} + (\text{total derivative})$,
setting $\delta S = 0$ gives the field equations. The total derivative needs a
boundary term (Gibbons–Hawking–York) for a well-posed variational problem.

## Why this equation? Lovelock's theorem

Lovelock (1971–72): in **four** dimensions, the only symmetric, divergence-free
tensor built locally from $g_{\mu\nu}$ and its first two derivatives is

```math
a\, G_{\mu\nu} + b\, g_{\mu\nu} .
```

Matching Newton fixes $a$ ([03](03-newtonian-limit.md)); $b$ is $\Lambda$. To
get a different theory you must break an assumption: add fields (scalar-tensor
theories), allow higher derivatives ($f(R)$ gravity), go to more dimensions
(Gauss–Bonnet), or give up locality.

## Discussion prompts

- Is $\Lambda$ part of the geometry (left side) or a form of energy (right side)? Does it matter?
- The equations determine the metric only up to coordinates. What, then, is a "solution"?
- Lovelock says GR is unique under its assumptions. Which assumption would you drop first, and why?
- Why is there no $T_{\mu\nu}$ for the gravitational field itself? (See also [06](06-open-questions.md).)
