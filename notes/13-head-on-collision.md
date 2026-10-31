# 13 · One head-on collision, constructed

The goal of notes 07–12 was one special but nontrivial example: a
nonspinning, axisymmetric, head-on collision of two black holes in an
asymptotically flat vacuum spacetime, with the merger proved rather than
simulated. This note assembles it. The ingredients are proved in notes 09,
10 and 12; here they are put together, with each one's status.

*Status:* part 1 of the theorem rests on Theorem F′ (note 12 §7). Theorem F′
was verified in session 10 for every $d \in (d_c, 0.767578]$.

No closed-form metric is known for such a collision ([note 07](07-the-merger-problem.md) §4).
So "constructing" the example means three steps:

1. write down initial data;
2. take its vacuum evolution, which exists by the general theory;
3. prove what the horizons do in that evolution.

## The example

**Data.** Equal-mass Brill–Lindquist data with separation $d$:

```math
\Sigma = \mathbb R^3 \setminus \{p_+, p_-\},\quad p_\pm = (0, 0, \pm d/2),\qquad
\gamma = \psi^4\delta,\quad \psi = 1 + \frac{1}{4|x - p_+|} + \frac{1}{4|x - p_-|},\qquad K = 0 .
```

Units are $G = c = 1$, and the ADM mass is $M = 1$ (bare masses $\tfrac12$).

**What makes it the right kind of example.**

- *Vacuum and asymptotically flat.* $\psi$ is harmonic, so the constraints
  hold with $K = 0$. Each puncture is a second asymptotically flat end, the
  inside of a black-hole throat.
- *Nonspinning, head-on.* $K = 0$ is a moment of time symmetry: both holes
  are momentarily at rest, with zero momentum and angular momentum. They
  start to fall toward each other along the line joining them.
- *Axisymmetric.* The data are invariant under rotations about the $z$-axis
  and under the reflection $z \to -z$. So is their evolution, by uniqueness
  of the development.

**Spacetime.** The maximal globally hyperbolic vacuum development $(M, g)$
(Choquet-Bruhat–Geroch). We only use a short stretch of it, the geodesic
slicing $\Sigma_t$: unit lapse and zero shift, starting from $\Sigma_0$. That
stretch is controlled by local well-posedness alone.

**Why "near-critical".** There is a critical separation
$d_c = 0.766197428240 \pm 1.1\times10^{-13}$ (Theorem C, note 09).

- For $d < d_c$ the data already contain a surface enclosing both holes.
- For $d$ just above $d_c$ they do not. The two holes then start apart, and
  only the dynamics can produce a common horizon.

Choosing $d - d_c$ small makes that happen almost at once. This is what lets
a short-time analysis, with no global control, prove that it happens.

## The theorem

> **Theorem H (computer-assisted).** There is $\delta_0 > 0$ such that for every
> separation $d$ with $d_c < d < d_c + \delta_0$ the development of the data
> above has the following properties.
>
> 1. **Two holes at $t = 0$.** Each puncture is enclosed by its own minimal
>    sphere, and the two spheres are disjoint. No closed embedded minimal
>    surface encloses both punctures. On this slice, MOTSs are minimal
>    surfaces, so the apparent horizon of $\Sigma_0$ is exactly the two
>    individual horizons.
> 2. **A common horizon forms.** There is a neighbourhood $\mathcal U$ of the
>    critical surface $S_c$ and a $t_0 > 0$ with the following property. For
>    $0 \le t \le t_0$, the slice $\Sigma_t$ contains a MOTS in $\mathcal U$
>    if and only if $t \ge t_*(d)$, where
>
>    ```math
>    t_*(d) = c\,(d - d_c) + O\big((d - d_c)^2\big), \qquad c \in [2.3768, 2.3936].
>    ```
>
>    That MOTS is an axisymmetric, reflection-symmetric sphere enclosing both
>    punctures.
> 3. **The individual horizons persist.** For $0 \le t \le t_0$ each puncture
>    stays enclosed by its own MOTS on $\Sigma_t$.
>
> So along the geodesic slicing, the apparent horizon has two separate
> components at $t = 0$. By $t = t_*(d)$ it has a component enclosing both
> holes, so the outermost MOTS jumps outward. The birth time of the common
> horizon is linear in $d - d_c$, with a certified coefficient.

### Proof: where each part comes from

| Part | Source | Status |
| --- | --- | --- |
| data are vacuum, AF, symmetric | Brill–Lindquist 1963; harmonic $\psi$ | classical |
| development exists, symmetric, smooth in $(t, d)$ | Choquet-Bruhat 1952, Choquet-Bruhat–Geroch 1969; smooth dependence (note 10, Step 1) | standard, cited |
| fold point, signs of $\partial^2_{z_0}F$, $\partial_dF$, $\partial_tF$ | Theorem C ([note 09](09-merger-hypotheses.md), `research/calc/fold_proof.py`) | computer-assisted (Arb) |
| part 2 (common MOTS iff $t \ge t_*$, in $\mathcal U$) | Proposition B ([note 10](10-onset-theorem-proof.md), Steps 2–5) | written out |
| part 3 (individual MOTSs persist) | note 10, Step 6: barriers for $d \in [0.76, 0.78]$, Andersson–Metzger | written out |
| part 1 (no common surface **anywhere** at $t = 0$) | Theorem F′ ([note 12](12-no-common-horizon.md) §7) for every $d \in (d_c, 0.767578]$ | computer-assisted |

Part 1 is the piece that was missing. Proposition B controls MOTSs only
inside the neighbourhood $\mathcal U$ of $S_c$, and only for $d$ within an
unspecified $\delta_0$ of $d_c$. Theorem F at the single separation
$d = 0.78$ could not be combined with it, because $0.78$ may lie beyond
$\delta_0$. Theorem F′ holds on a whole interval that starts at $d_c$. For
every $d$ near $d_c$ it rules out a common surface anywhere on the slice:
not only near $S_c$, but also the necked surfaces hugging the individual
horizons, and anything farther out. Shrinking $\delta_0$ so that
$d_c + \delta_0 \le 0.767578$ gives the theorem. ∎

## What the example does and does not show

| Question | Answer |
| --- | --- |
| a genuine two-body vacuum collision, AF, axisymmetric, nonspinning? | yes |
| the common apparent horizon is created by the dynamics, not present in the data? | yes (parts 1–2) |
| when it is born? | $t_* = c(d - d_c) + O((d - d_c)^2)$, $c \in [2.3768, 2.3936]$, geodesic slicing; for lapse $\psi^{-2}$, $c \in [10.606, 10.682]$ |
| a specific numerical $d$ for which all this holds? | **no**: $\delta_0$ comes from an implicit-function argument and is not computed |
| two components on a whole interval of times before $t_*$, as note 07's definition asks? | **not written out**. Proved at $t = 0$; for $0 < t < t_*$, part 2 excludes common MOTSs only in $\mathcal U$. A compactness argument for stable MOTSs (Andersson–Metzger) should give two components for $0 \le t < t_1(d)$, but has not been checked |
| event horizons (two components, then one)? | **no**. This needs the global development: complete $\mathscr I^+$ and the final state. Out of reach for these data |
| independent review of the computer-assisted parts? | not done (trust base listed in notes 09 and 12) |

**Comparison with the other examples in this repository.**

- **Kastor–Traschen** (note 08) is an exact merger, but with $\Lambda > 0$
  and charge, and its common horizon is only checked at sample times.
- **Hintz's extreme-mass-ratio solutions** (note 11) give event-horizon
  mergers in vacuum with $\Lambda = 0$. They need a small mass ratio, are
  not explicit, and rest on our reading of Hintz 2026. A small Schwarzschild
  hole falling radially into a large one is also head-on and nonspinning.
  Whether that construction can be chosen axisymmetric was not checked.
- **Theorem H** is the comparable-mass case: equal masses, an explicit
  starting point, an explicit critical separation and onset rate. The price
  is that it controls apparent horizons only, and for $d$ near $d_c$.

## Discussion prompts

- The whole proof lives in a time interval of length $O(d - d_c)$. What
  would it take to make $\delta_0$ explicit? One route is Cauchy–Kovalevskaya
  majorants for the analytic data near $S_c$. What would a numerical value of
  $\delta_0$ say physically?
- Numerical relativity finds that the event horizon of a head-on collision
  encloses both holes long before the common apparent horizon appears. Is
  there a quasi-local quantity on $\Sigma_0$ that would show this early
  merging, without the global development?
- Theorem H uses one slicing family. With lapse $\psi^{-2}$ the onset is
  about 4.5 times later. Which features of the example are independent of
  the slicing, and which are artefacts of it?
