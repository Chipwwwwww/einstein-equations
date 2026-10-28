# 04 · Exact solutions

The field equations are too nonlinear to solve in general, so nearly all exact
solutions assume a lot of symmetry. The handful below carry most of
astrophysics and cosmology. Units: $G = c = 1$ unless stated; $d\Omega^2 = d\theta^2 + \sin^2\theta\, d\phi^2$.

| Solution | Year | Symmetry | Matter | Test |
| --- | --- | --- | --- | --- |
| Minkowski | 1908 | maximal | none | `test_flat_space_in_spherical_coordinates` |
| Schwarzschild | 1916 | static, spherical | none | `test_schwarzschild_is_vacuum`, `test_schwarzschild_kretschmann` |
| de Sitter | 1917 | maximal | $\Lambda$ only | `test_de_sitter_is_lambda_vacuum` |
| FLRW | 1922–37 | homogeneous, isotropic space | perfect fluid | `test_flrw_friedmann_equations`, `test_flrw_acceleration_equation` |
| Kerr | 1963 | stationary, axial | none | not yet (see below) |

Tests live in [`tests/test_solutions.py`](../tests/test_solutions.py).

## Schwarzschild (1916)

```math
ds^2 = -\left(1 - \frac{2M}{r}\right) dt^2 + \left(1 - \frac{2M}{r}\right)^{-1} dr^2 + r^2 d\Omega^2
```

With units restored, $r_s = 2GM/c^2$, about 3 km per solar mass.

- **Birkhoff's theorem.** Any spherically symmetric vacuum solution is (a piece of) Schwarzschild, even if the source is pulsating. A spherical star cannot radiate gravitational waves.
- **$r = 2M$ is not a curvature singularity.** The metric component blows up but
  $K = R_{\rho\sigma\mu\nu}R^{\rho\sigma\mu\nu} = 48 M^2 / r^6$ is finite there. It is a horizon,
  and coordinates such as Kruskal–Szekeres or Eddington–Finkelstein cross it smoothly.
- **$r = 0$ is a real singularity**: $K \to \infty$.
- Landmarks: photon sphere at $r = 3M$, innermost stable circular orbit at $r = 6M$.

## de Sitter (1917)

The vacuum with $\Lambda > 0$, $R_{\mu\nu} = \Lambda g_{\mu\nu}$. In static coordinates

```math
ds^2 = -\left(1 - \frac{\Lambda r^2}{3}\right) dt^2 + \left(1 - \frac{\Lambda r^2}{3}\right)^{-1} dr^2 + r^2 d\Omega^2 ,
```

with a cosmological horizon at $r = \sqrt{3/\Lambda}$. In cosmological
coordinates the same spacetime is $ds^2 = -dt^2 + e^{2Ht} d\vec x^2$ with
$H = \sqrt{\Lambda/3}$: exponential expansion, the far future of our universe
if dark energy is a cosmological constant.

## FLRW: Friedmann–Lemaître–Robertson–Walker

Homogeneous and isotropic space with scale factor $a(t)$ and curvature $k \in \{-1, 0, +1\}$:

```math
ds^2 = -dt^2 + a(t)^2 \left[ \frac{dr^2}{1 - k r^2} + r^2 d\Omega^2 \right] .
```

Its Einstein tensor is

```math
G_{tt} = 3\left( \frac{\dot a^2}{a^2} + \frac{k}{a^2} \right), \qquad
G_{ij} = -\left( \frac{2\ddot a}{a} + \frac{\dot a^2}{a^2} + \frac{k}{a^2} \right) g_{ij} .
```

With a perfect fluid $T_{\mu\nu} = (\rho + p)\, u_\mu u_\nu + p\, g_{\mu\nu}$ at rest
($u_\mu = (-1, 0, 0, 0)$), the $tt$ and $ij$ equations become the **Friedmann equations**.
With $G$ and $c$ restored ($\rho$ a mass density):

```math
H^2 \equiv \left( \frac{\dot a}{a} \right)^2 = \frac{8\pi G}{3}\rho - \frac{k c^2}{a^2} + \frac{\Lambda c^2}{3}
```

```math
\frac{\ddot a}{a} = -\frac{4\pi G}{3} \left( \rho + \frac{3p}{c^2} \right) + \frac{\Lambda c^2}{3}
```

and $\nabla_\mu T^{\mu\nu} = 0$ gives $\dot\rho + 3H(\rho + p/c^2) = 0$.
Acceleration needs $\rho + 3p/c^2 < 0$ or $\Lambda > 0$: this is the
equation behind the 1998 supernova discovery of cosmic acceleration.

## Kerr (1963)

The rotating black hole, mass $M$, angular momentum $J = aM$. In Boyer–Lindquist coordinates

```math
ds^2 = -\left(1 - \frac{2Mr}{\Sigma}\right) dt^2
  - \frac{4 M a r \sin^2\theta}{\Sigma}\, dt\, d\phi
  + \frac{\Sigma}{\Delta} dr^2 + \Sigma\, d\theta^2
  + \left( r^2 + a^2 + \frac{2 M a^2 r \sin^2\theta}{\Sigma} \right) \sin^2\theta\, d\phi^2
```

with $\Sigma = r^2 + a^2 \cos^2\theta$ and $\Delta = r^2 - 2Mr + a^2$. Horizons
at $r_\pm = M \pm \sqrt{M^2 - a^2}$; for $a > M$ there is no horizon (a naked
singularity, forbidden by cosmic censorship if it holds).

**Open contribution:** add `kerr()` to [`einstein/metrics.py`](../einstein/metrics.py)
and a test that its Ricci tensor vanishes. SymPy needs a while to simplify it,
so mark the test with `@pytest.mark.slow`.

## Discussion prompts

- Schwarzschild is "eternal". How does a black hole formed from collapse (Oppenheimer–Snyder) differ?
- Is the FLRW assumption of homogeneity a result or an input? What tests it?
- Kerr has two parameters, mass and spin ("no hair"). What happened to everything else about the star that collapsed?
- de Sitter has a horizon and a temperature, but no black hole. What is radiating?
