# 08 · An exact merger: Kastor–Traschen, and why it is not the answer

> **Benchmark, not the target.** Kastor–Traschen spacetimes solve the
> Einstein–**Maxwell** equations with **Λ > 0**, contain **extremally charged**
> holes, and are asymptotically **de Sitter**. The target of
> [07](07-the-merger-problem.md) is vacuum, Λ = 0, asymptotically flat. Every
> statement below is about the benchmark unless it says otherwise.

Kastor and Traschen (1993) found the first exact solutions in which black holes
coalesce. They are simple enough to check line by line, so they are the right
place to calibrate tools and to see what an exact merger needs.

## 1. The solution

Units $G = c = 1$, Gaussian electromagnetic units (see
[`einstein/fields.py`](../einstein/fields.py)). In the coordinates of Brill,
Horowitz, Kastor and Traschen (1994, "BHKT"),

```math
ds^2 = -\frac{d\tau^2}{U^2} + U^2\left(dx^2 + dy^2 + dz^2\right), \qquad
A = \frac{d\tau}{U}, \qquad
U = H\tau + \sum_i \frac{m_i}{|\vec x - \vec x_i|}, \qquad H = \pm\sqrt{\Lambda/3} .
```

The field equations are $G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi T_{\mu\nu}$ with
$T_{\mu\nu} = \frac{1}{4\pi}\left(F_{\mu\alpha}F_\nu{}^\alpha - \frac14 g_{\mu\nu}F^2\right)$,
and $\nabla_\mu F^{\mu\nu} = 0$.

The original paper uses cosmological time $t$ with $a = e^{Ht}$ and
$\Omega = 1 + \sum m_i/(a r_i)$, $ds^2 = -\Omega^{-2}dt^2 + a^2\Omega^2 d\vec x^2$. The
two forms agree with $U = a\Omega$ and $\tau = a/H$ (so $d\tau = a\,dt$)
(`test_original_cosmological_form_matches`).

**Domain.** The metric is Lorentzian and smooth where $0 < U < \infty$, i.e.
away from the centres $\vec x_i$ and from the hypersurface $U = 0$. For
$H < 0$ and $\tau < 0$, every term in $U$ is positive, so the whole region
$\{\tau < 0\} \times (\mathbb R^3 \setminus \{\vec x_i\})$ is regular. The centres
are not points of spacetime: as in Majumdar–Papapetrou, $r_i \to 0$ is a
horizon at the end of an infinitely long throat, and BHKT extend through it.

## 2. What we verified, and how

All of this is in [`tests/test_kastor_traschen.py`](../tests/test_kastor_traschen.py).

**Reduced system (exact identity).** For an *arbitrary* function
$U(\tau, x, y, z)$, the metric and potential above give

```math
E_{\tau\tau} = -\frac{\Lambda U^3 - 3U^3 U_\tau^2 + 2\Delta U}{U^5}, \qquad
E_{\tau i} = -\frac{2\,\partial_\tau\partial_i U}{U}, \qquad
E_{ij} = -U^2\left(2UU_{\tau\tau} + 3U_\tau^2 - \Lambda\right)\delta_{ij},
```

```math
\nabla_\mu F^{\mu\tau} = \frac{\Delta U}{U^2}, \qquad
\nabla_\mu F^{\mu i} = -\frac{\partial_\tau \partial_i U}{U^2},
```

where $E_{\mu\nu} = G_{\mu\nu} + \Lambda g_{\mu\nu} - 8\pi T_{\mu\nu}$ and
$\Delta$ is the flat Laplacian. (`test_reduced_system_for_arbitrary_U`.)

**Classification within the ansatz.** Solving these:

1. $E_{\tau i} = 0$ gives $U = h(\tau) + V(\vec x)$.
2. $E_{ij} = 0$ gives $2(h + V)h'' + 3h'^2 = \Lambda$ at every point. If $V$ is
   not constant, differentiating in $x$ gives $h'' = 0$, so $h = H\tau$ (a
   constant can be moved into $V$) and $\Lambda = 3H^2$.
3. $E_{\tau\tau} = 0$ then gives $\Delta V = 0$, which is also the Maxwell
   equation.

So **within this ansatz, the Kastor–Traschen family is everything**: $U = H\tau + V$
with $V$ harmonic and $H^2 = \Lambda/3$. With $V = \sum m_i/r_i$ the equations
hold for any number of centres. Two centres are checked twice: through the
reduction (exact), and by evaluating the full, unsimplified residual in exact
rational arithmetic at points where both distances are integers (an
independent diagnostic).

**Corollary (the obstruction we care about).** *With Λ = 0, every solution in
this ansatz is static* (Majumdar–Papapetrou): step 2 forces $h' = 0$.
(`test_lambda_zero_forces_static`.) The motion in Kastor–Traschen comes
entirely from the cosmological constant.

**One centre is Reissner–Nordström–de Sitter with $Q = m$.** The areal radius
$R = rU$ satisfies, exactly,

```math
g^{\mu\nu}\partial_\mu R\, \partial_\nu R = \left(1 - \frac{m}{R}\right)^2 - H^2R^2
= 1 - \frac{2m}{R} + \frac{m^2}{R^2} - \frac{\Lambda}{3}R^2 ,
```

the RNdS horizon function with $Q = m$. The left side is a scalar, so this is a
coordinate-free identification (with the generalised Birkhoff theorem for
Einstein–Maxwell–Λ).

**Mass bound.** The horizons are roots of $(1 - m/R)^2 = H^2R^2$. The
black-hole and cosmological horizons, roots of $|H|R^2 - R + m = 0$, merge at

```math
m = M_{\rm ext} = \frac{1}{4|H|} = \sqrt{\frac{3}{16\Lambda}} .
```

BHKT's coalescence statement assumes total mass $M = \sum m_i < M_{\rm ext}$.

**Scaling.** $(\tau, \vec x, \vec x_i, m_i, H) \to (\lambda\tau, \lambda\vec x, \lambda\vec x_i, \lambda m_i, H/\lambda)$
leaves $U$ unchanged and multiplies $g$ by $\lambda^2$. The physical
parameters are therefore only $Hm_i$ and $H|\vec x_i - \vec x_j|$.

## 3. The contracting branch

With $\tau = a/H$: for $H > 0$, $\tau > 0$ grows with cosmological time and
space expands. For **$H < 0$**, $\tau \in (-\infty, 0)$ still grows toward the
future, but $a = H\tau = |H||\tau|$ shrinks to zero as $\tau \to 0^-$: the
background **contracts**. The centres sit at fixed coordinate positions, so
the physical distance between them shrinks like $a$ while the holes keep
their size. BHKT: "It is only for this case that the spacetime describes black
hole collisions." The expanding branch is the time reverse: white holes
splitting.

Why don't the holes fall together on their own? Each pair has charge equal to
mass, so electric repulsion cancels gravitational attraction exactly, as in
Majumdar–Papapetrou. The holes are "comoving": the merger is driven by the
contraction of the de Sitter background, not by their mutual attraction.

## 4. A simple fact: apparent horizons are constant-mean-curvature surfaces

The slices $\tau = $ const have lapse $N = 1/U$, metric $\gamma_{ij} = U^2\delta_{ij}$
and zero shift, so

```math
K_{ij} = \frac{1}{2N}\partial_\tau\gamma_{ij} = U_\tau\,\gamma_{ij} = H\,\gamma_{ij}, \qquad K = 3H .
```

Every slice is **umbilic** with **constant** mean curvature $3H$ (it is the
flat-slicing value for de Sitter). For any closed surface $S$ in a slice,
$\mathrm{tr}_S K = 2H$, so

```math
\theta_+ = 2H + H_S = 0 \iff H_S = -2H = 2|H| .
```

**MOTSs on Kastor–Traschen slices are exactly the closed surfaces of constant
mean curvature $2|H|$ in the conformally flat metric $U^2\delta$.** For one
centre and coordinate spheres, $\theta_+ = 0$ reproduces $1 - m/R = -HR$: the
RNdS black-hole and cosmological horizons
(`test_single_centre_mots_are_the_rnds_horizons`). We use this in §6 to find
the two-centre apparent horizons by solving a single ODE.

We have not found this observation stated in BHKT or KT (BHKT compute
$\theta$ on centred spheres). It may well be in Nakao–Shiromizu–Hayward (1995),
which we could not access; we claim no novelty for it.

## 5. What is local and what rests on published global analysis

| Statement | Status here |
| --- | --- |
| metric + potential solve Einstein–Maxwell–Λ on $\{U > 0\}$ | **verified** (exact) |
| one centre is RNdS with $Q = m$ | **verified** (scalar identity) |
| horizons merge at $M_{\rm ext} = 1/(4\lvert H\rvert)$ | **verified** (discriminant) |
| MOTSs on slices are CMC surfaces, $H_S = 2\lvert H\rvert$ | **verified** (exact) |
| two-centre MOTSs merge in time | numerical, §6 |
| the coordinates extend through the horizons; extension non-unique; metric $C^2$ at the de Sitter horizon, $C^1$ at inner horizons | BHKT 1994, **not re-derived** |
| for equal masses, the *event* horizon encloses each hole separately at early times and both at late times ("these results prove that the black holes coalesce") | BHKT 1994, analytic bounds on radial null geodesics, **not re-derived** |
| for $M > M_{\rm ext}$: a naked singularity forms to the future of a Cauchy horizon | BHKT 1994, **not re-derived** |

So the coalescence of event horizons in Kastor–Traschen rests on BHKT's
null-geodesic estimates. We reproduce the local field equations and the
apparent-horizon picture.

## 5a. Why Λ → 0 does not give the target

1. **Inside the ansatz, Λ = 0 is static** (corollary in §2). There is nothing
   to take a limit of.
2. **The limit at fixed masses is the static Majumdar–Papapetrou solution.**
   As $H \to 0$, $U \to V$ pointwise. By the scaling in §2, small $\Lambda$ at
   fixed masses is the same spacetime as small masses at fixed $\Lambda$. The
   merger takes a cosmological time of order $|H|^{-1}\ln(d/m)$ for
   separation $d$, which diverges as $H \to 0$.
3. **The matter does not switch off.** Within the ansatz charge equals mass,
   so the vacuum limit $Q \to 0$ is $m \to 0$: no black holes.
4. **The asymptotics are wrong.** The far region is de Sitter in flat slicing,
   whose $\mathscr I^+$ is spacelike; there is no asymptotically flat end.

Points 1 and 3 say the mechanism, not just the solution, fails at Λ = 0. In
[09](09-merger-hypotheses.md) we turn point 1 into a statement that does not
assume the ansatz: zero-shift umbilic slicings of asymptotically flat
spacetimes are static (Lemma A there).

## Discussion prompts

- The holes in Kastor–Traschen never attract each other. Is it still fair to
  call this a "collision"?
- The merger is pushed by a contracting de Sitter background. What would the
  analogous push be in an asymptotically flat spacetime?
- Extremal charge makes the force balance work. Which other "no-force"
  conditions (BPS states in supergravity, for example) could play the same
  role?
