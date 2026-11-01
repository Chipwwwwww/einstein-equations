# 07 · The binary black hole merger problem, stated precisely

Two black holes spiral together, merge and ring down to a single Kerr black
hole. LIGO has heard this happen many times, and numerical relativity
reproduces the waveforms. Yet no exact solution of Einstein's equations
describes it, and no theorem guarantees that such a spacetime exists. This
note pins down what "a solution describing a merger" would mean, so that later
notes ([08](08-kastor-traschen.md), [09](09-merger-hypotheses.md)) can say
exactly what they prove and what they do not.

Units $G = c = 1$, $\Lambda = 0$, conventions of [00](00-conventions.md).
Extrinsic curvature is $K_{ij} = \nabla_i n_j$ for the future unit normal $n$,
so $K > 0$ means the slice is expanding (Wald's sign; MTW use the opposite).

## 1. The equations

The target is a 4-dimensional, asymptotically flat **vacuum** spacetime:

```math
R_{\mu\nu}[g] = 0 .
```

### Initial data and constraints

Einstein's equations are a hyperbolic system only after a gauge is fixed, and
part of them are constraints on the initial slice. Initial data are a
3-manifold $\Sigma$ with a Riemannian metric $\gamma_{ij}$ (the geometry of
space at one instant) and a symmetric tensor $K_{ij}$ (how that geometry is
changing). They must satisfy the Hamiltonian and momentum constraints

```math
{}^{(3)}\!R + K^2 - K_{ij}K^{ij} = 0, \qquad
D_j K^j{}_i - D_i K = 0 ,
```

which are the $G_{\mu\nu}n^\mu n^\nu = 0$ and $G_{\mu i}n^\mu = 0$ components of
the field equations. Physically: the curvature of space is fixed by the energy
density (zero here), and its rate of change by the momentum density (zero).

With lapse $N$ and zero shift, the slices evolve by

```math
\partial_t \gamma_{ij} = 2N K_{ij}, \qquad
\partial_t K_{ij} = D_iD_jN - N\left( {}^{(3)}\!R_{ij} + K K_{ij} - 2K_{ik}K^k{}_j \right)
```

(vacuum, our sign of $K$). At a moment of time symmetry ($K_{ij} = 0$) the
second equation reduces to $\partial_t K_{ij} = D_iD_jN - N\,{}^{(3)}\!R_{ij}$.

### Regularity, topology, asymptotics

- **Asymptotic flatness.** Outside a compact set, $\Sigma$ is a finite union of
  ends diffeomorphic to $\mathbb R^3$ minus a ball, on which
  $\gamma_{ij} - \delta_{ij} = O(r^{-q})$, $\partial\gamma = O(r^{-q-1})$,
  $K_{ij} = O(r^{-q-1})$ with $q > 1/2$, so that the ADM energy and momentum are
  defined (Bartnik 1986). One end is designated physical; the others, if any,
  are hidden inside black holes.
- **Topology for two black holes.** Common choices: $\Sigma = \mathbb R^3$ minus
  two points with three ends (Brill–Lindquist "punctures"), the two-sheeted
  Misner wormhole, or $\mathbb R^3$ minus two balls with inner boundaries
  (excision). The choice is a modelling convention: the black-hole interior
  is not part of what we want to describe.
- **Regularity.** Local existence needs $(\gamma, K) \in H^s_{loc}\times H^{s-1}_{loc}$
  with $s > 5/2$; smooth data suffice for everything below.

### The development

**Theorem** (Choquet-Bruhat 1952; Choquet-Bruhat–Geroch 1969). Smooth vacuum
data satisfying the constraints have a **maximal globally hyperbolic
development** (MGHD), unique up to isometry, in which $\Sigma$ embeds as a
Cauchy surface with induced data $(\gamma, K)$.

This is all that local well-posedness gives. It says **nothing** about how
long the development lasts in any time function, whether it contains a
complete future null infinity, whether its singularities are hidden, or
whether it settles down to Kerr. Penrose's 1965 theorem (null energy
condition + non-compact Cauchy surface + a closed trapped surface) even shows
the MGHD of black-hole data is **future null geodesically incomplete**.

Gauge enters only through how a solution is written down. Statements below
are made about geometric objects (horizons, trapped surfaces, scalars), never
about coordinate components.

## 2. What "two black holes" and "merger" can mean

### Event horizons (global)

Let $\mathscr I^+$ be future null infinity. The **black hole region** is
$\mathcal B = M \setminus J^-(\mathscr I^+)$: the events from which no signal
reaches distant observers. Its boundary $\mathcal H = \partial J^-(\mathscr I^+)$
is the **event horizon**, a null hypersurface (possibly non-smooth) generated
by null geodesics that, once on it, never leave toward the future.

Defining $\mathcal B$ needs **the whole future** of the spacetime, and it is
only meaningful when $\mathscr I^+$ is complete (Christodoulou's notion: the
affine length of its generators is infinite). The horizon is *teleological*:
whether a point is inside depends on what happens later.

**Definition (event-horizon merger).** Let an asymptotically flat vacuum
spacetime with complete $\mathscr I^+$ contain two black holes, i.e. two
families of closed trapped surfaces $T_1(t)$, $T_2(t)$ along a Cauchy foliation
$\Sigma_t$ (each reaching spatial infinity). It *describes their merger* if
there are times $t_1 < t_2$ such that $T_1(t_1)$ and $T_2(t_1)$ lie in
**different** connected components of $\mathcal B \cap \Sigma_{t_1}$, while
$T_1(t_2)$ and $T_2(t_2)$ lie in the **same** component of
$\mathcal B \cap \Sigma_{t_2}$. *Initially distinct* means the components differ
for all $t \le t_1$. (Asking for $\mathcal B \cap \Sigma_{t_2}$ to be connected
would be stronger. Stray components elsewhere would not change the physics.)

The number of components of $\mathcal B \cap \Sigma$ depends on the slice (the
transient toroidal horizons of Emparan–Martínez–Zilhão 2018 are a reminder),
so the definition asks only that *some* pair of slices shows the change. In a
strongly asymptotically predictable spacetime black holes cannot bifurcate
(Hawking & Ellis 1973, §9.2), so two components followed by one is a merger,
never the reverse.

### Apparent horizons (quasi-local, slice-dependent)

On a slice $\Sigma_t$, let $S$ be a closed 2-surface with outward unit normal
$\nu$ in $\Sigma_t$ and mean curvature $H_S = \mathrm{div}_S \nu$. Its future
null normals $\ell_\pm = n \pm \nu$ have expansions

```math
\theta_\pm = \mathrm{tr}_S K \pm H_S .
```

$\theta_+$ is the rate at which the area of an outgoing light front grows.

| Name | Condition | Meaning |
| --- | --- | --- |
| trapped surface | $\theta_+ < 0,\ \theta_- < 0$ | both light fronts shrink |
| outer trapped | $\theta_+ < 0$ | outgoing light front shrinks |
| marginally outer trapped surface (MOTS) | $\theta_+ = 0$ | marginal |
| marginally trapped surface | $\theta_+ = 0,\ \theta_- < 0$ | the usual "apparent horizon" cross-section |

The **apparent horizon** of $\Sigma_t$ is the boundary of its (weakly outer)
trapped region, which is a smooth stable MOTS (Andersson–Metzger 2009).
Finding it needs only the data on one slice: this is what numerical
relativity tracks. But it depends on the slice: Wald and Iyer (1991) showed
that even Schwarzschild has Cauchy slices, reaching arbitrarily close to the
singularity, with no trapped surfaces at all.

**Definition (apparent-horizon merger along a foliation $\Sigma_t$).** The
outermost MOTS of $\Sigma_t$ has two components for $t < t_1$ and a component
enclosing both for $t > t_2$.

### How the two notions are related

- Under null energy and strong asymptotic predictability, closed trapped
  surfaces lie inside $\mathcal B$ (Hawking & Ellis 1973, Prop. 9.2.1; Wald 1984,
  Prop. 12.2.2). So a common trapped surface on $\Sigma_2$ forces the black
  hole region there to contain a connected set enclosing both holes.
- The converse fails: the event horizon can enclose both holes long before any
  common MOTS appears on the chosen slices (in numerical runs it typically
  does), and the outermost MOTS **jumps** discontinuously when a common one
  forms (Andersson–Mars–Metzger–Simon 2009).
- **A common apparent horizon is therefore not a proof of event-horizon
  merger.** Even granting predictability, it shows only that $\mathcal B \cap
  \Sigma_2$ has a component enclosing both holes. To conclude "merger" one
  also needs two components at an earlier time, which needs light rays that
  escape to $\mathscr I^+$ from between the holes: global information again.

## 3. Which claims need which information

| Claim | Information needed |
| --- | --- |
| data satisfy the constraints | the slice only |
| a MOTS / trapped surface exists on a slice | the slice only |
| a development exists for a finite time | local theory + Cauchy stability |
| MGHD exists and is unique | Choquet-Bruhat–Geroch |
| $\mathscr I^+$ is complete (weak cosmic censorship for these data) | global, open in general |
| event horizon exists, number of components | global ($J^-(\mathscr I^+)$) |
| exterior is regular (no naked singularities) | global |
| late-time approach to Kerr | global, plus a stability theorem for the final state |

**Kerr asymptotics is a separate claim.** It needs a final state (parameters
$M_f$, $a_f$) and convergence of the exterior to Kerr in a stated gauge and
norm, at stated rates. Proven results are perturbative: Schwarzschild
(Dafermos–Holzegel–Rodnianski–Taylor 2021, codimension 3), slowly rotating
Kerr (Klainerman–Szeftel 2023; Giorgi–Klainerman–Szeftel 2022, published
2024), Kerr–de Sitter (Hintz–Vasy 2018). None applies to the strongly
nonlinear merger phase itself.

## 4. What would count as a result

| Kind | Example | Status for the vacuum, $\Lambda = 0$ target |
| --- | --- | --- |
| explicit exact metric | a closed-form $g$ with $R_{\mu\nu} = 0$ and a merging horizon | none known |
| rigorous existence theorem | "for data in class X the MGHD has complete $\mathscr I^+$ and $\mathcal B$ merges" | only extreme mass ratio: Hintz 2024 (below) |
| convergent representation with proof | a series or iteration proved to converge to a merger spacetime | none known |
| controlled approximation | expansion with proven error bounds, e.g. in the mass ratio | partial: Hintz's construction is one |
| numerical evidence | convergent numerical relativity | abundant (Pretorius 2005; Campanelli et al., Baker et al. 2006) |
| restricted obstruction | "no merger within ansatz Y" | several, e.g. notes/08 and notes/09 |

The closest rigorous result to the target is **Hintz (2024)**, *Gluing small
black holes along timelike geodesics III*, arXiv:2408.06715. For a small
subextremal Kerr black hole of mass $\epsilon m$ sent along a timelike geodesic
that crosses the event horizon of a slowly rotating unit-mass Kerr black hole,
it constructs vacuum solutions on the region after an initial slice for all
small $\epsilon$. Concatenated with the Kerr stability theorem, the merged hole
settles to a nearby Kerr black hole (Theorem 6.1 for $\Lambda > 0$; Remark 6.2
for $\Lambda = 0$). The paper explicitly leaves open "the construction of the
event horizon, the structure of the black hole region — conjecturally
connected". The black hole there "merges" in the sense of its Definition 4.1:
the small hole's worldline enters the large hole's horizon.

## Discussion prompts

- The event horizon depends on the whole future. Is it the right notion of
  "black hole" for a process we only observe at finite times?
- Apparent horizons jump. Why would a smooth physical process produce a
  discontinuous quasi-local horizon, and what does the jump tell us?
- An exact merger solution would have no symmetry left except possibly an
  axis. Is that by itself a reason to expect no closed form?
