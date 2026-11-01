# 09 · Four hypotheses about exact mergers, and what survived

[07](07-the-merger-problem.md) set the target: a vacuum, asymptotically flat,
$\Lambda = 0$ spacetime in which two black holes merge.
[08](08-kastor-traschen.md) showed how the one exact merger we have
(Kastor–Traschen) works, and why it does not reach the target. This note
records a bounded exploration: four mechanisms, each with a concrete test. Two
were run in depth: H2 (unconventional) and H4 (most promising). H1 was cheap
enough to settle along the way, and H3 was designed but not run.

Conventions: [07](07-the-merger-problem.md) and
[`research/README.md`](../research/README.md). In particular
$K_{ij} = \nabla_i n_j$ and $\theta_\pm = \mathrm{tr}_S K \pm H_S$.

| | Route | Mechanism | Outcome |
| --- | --- | --- | --- |
| H1 | explicit metric | Kastor–Traschen's umbilic, zero-shift slicing at Λ = 0 | **fails**: Lemma A, slices freeze |
| H2 | explicit metric (unconventional) | same, with hyperboloidal (Milne) slices, which escape Lemma A | **fails**: curvature of the slices forbids centres |
| H3 | existence | event horizons of Hintz's extreme-mass-ratio mergers via causal stability | **designed**; two inputs unverifiable this session |
| H4 | existence + controlled expansion | unfolding of the fold at which a common horizon is born | **survives**: Proposition B, with the fold proved by computer (Theorem C) |

## H1 · The Kastor–Traschen mechanism at Λ = 0

1. **Obstacle.** In Kastor–Traschen the motion comes from Λ (note 08, §2).
   Is there an asymptotically flat analogue of its slicing?
2. **Mechanism tested.** Kastor–Traschen slices have zero shift and are
   *umbilic*, $K_{ij} = \frac K3\gamma_{ij}$, with $K$ constant on each slice.
   Any ansatz $\gamma_{ij}(t,x) = \Psi(t,x)^4\mathring\gamma_{ij}(x)$ with zero
   shift has umbilic slices automatically.
3. **Closest work.** The ansatz classification of note 08. Lemma A is
   elementary and probably folklore; we did not find it stated.
4. **Assumptions.** Λ = 0; zero momentum density on the slices (vacuum, or
   a purely electric Maxwell field); slices asymptotically flat at spatial
   infinity.
5. **Test.** One identity (`test_lemma_A_momentum_constraint_for_conformally_flat_zero_shift`):
   for $g = -N^2dt^2 + \Psi^4\delta$ with arbitrary $N$ and $\Psi$, $G_{ti} = -\frac23 N\,\partial_i K$.

**Lemma A.** *Under these assumptions $K_{ij} \equiv 0$ and $\partial_t\gamma_{ij} = 0$:
every slice has the same, time-symmetric geometry.*

*Proof.* The momentum constraint with $K^j{}_i = \frac K3\delta^j_i$ reads
$-\frac23 D_iK = 8\pi j_i = 0$, so $K = K(t)$ on each (connected) slice.
Asymptotic flatness gives $K \to 0$ at infinity, so $K \equiv 0$, hence
$K_{ij} = 0$ and $\partial_t\gamma_{ij} = 2NK_{ij} = 0$. ∎

**Consequence.** Along such a foliation, MOTSs are minimal surfaces of one
fixed Riemannian 3-manifold, so the apparent horizon never changes: **no
apparent-horizon merger along any zero-shift umbilic foliation of an
asymptotically flat Λ = 0 spacetime.** (If, in addition, the normal Killing
initial data $(N, 0)$ are unique up to scale, the spacetime is static and
static black-hole uniqueness leaves only Schwarzschild. We did not prove that
uniqueness, so this stronger statement is a remark, not a result.)

6. **Verdict.** Abandoned. To escape Lemma A one must give up one of its
   assumptions: (a) allow a shift, (b) allow a trace-free $K_{ij}$ (needed
   anyway for gravitational radiation), or (c) use slices that do not reach
   spatial infinity. Option (c) is H2.

## H2 · Hyperboloidal Kastor–Traschen (unconventional)

1. **Obstacle.** Lemma A uses $K \to 0$ at *spatial* infinity.
2. **Mechanism.** Minkowski space has a zero-shift umbilic foliation with
   $K \neq 0$: the Milne wedge, $ds^2 = -dT^2 + T^2h_{-1}$ with $h_{-1}$ the
   hyperbolic metric and $K = 3/T$. Its slices are hyperboloidal (they end on
   $\mathscr I^+$), so Lemma A does not apply. Run backwards in time, the
   contracting Milne wedge (the inside of a past light cone) brings comoving
   observers together at a point, playing the role of contracting de Sitter.
   Hypothesis: extremal holes comoving with Milne, in Einstein–Maxwell with
   $\Lambda = 0$, form an exact merger in **asymptotically flat** space.
3. **Closest work.** Kastor–Traschen (flat slices, Λ > 0); dynamical
   multi-black-hole solutions in FLRW backgrounds need extra fields
   (Maeda–Nozawa and collaborators, not reviewed here).
4. **Assumptions.** Einstein–Maxwell, any Λ; the ansatz below. It keeps the
   merger geometry (comoving holes converging) and asymptotic flatness.
5. **Test.** Ansatz $g = -U^{-2}dt^2 + U^2a(t)^2h_{-1}$, $A = U^{-1}dt$. With
   $W = aU$ and $dt' = a\,dt$ the scale factor drops out:
   $g = -W^{-2}dt'^2 + W^2h_{-1}$, $A = W^{-1}dt'$. For arbitrary $W(t,\vec x)$
   (`test_hyperboloidal_kastor_traschen_reduced_system`):

```math
E_{ti} = -\frac{2\,\partial_t\partial_i W}{W}, \qquad
E_{ij} = -\left(2WW_{tt} + 3W_t^2 - \Lambda - W^{-2}\right) g_{ij} .
```

   $E_{ti} = 0$ forces $W = h(t) + V(\vec x)$. The $x$-derivative of the bracket
   is then $2V_x\,(h'' + W^{-3})$, so wherever $V_x \neq 0$ the function
   $W^{-3} = -h''(t)$ would be independent of $\vec x$, a contradiction. Hence
   **$V$ is constant: no holes**, for any Λ. With no holes, $W = \sqrt{2t}$ is
   exactly Milne (`test_hyperboloidal_background_is_milne`), which confirms
   the ansatz is set up correctly.

6. **Verdict.** Falsified, with an exact obstruction. The curvature of the
   hyperbolic slices adds the term $-W^{-2}$, which depends on position, and
   no time dependence can cancel it. Kastor–Traschen escapes because its
   slices are flat; Lemma A says flat umbilic slices of an asymptotically flat
   Λ = 0 spacetime cannot move. *Assumption that would have to change:* either
   matter whose energy density cancels $W^{-2}$ (a different physical problem),
   or a lapse not tied to the conformal factor by $N = 1/W$ (a much larger
   ansatz, with the Majumdar–Papapetrou force balance no longer available).
   This is a failure *within this ansatz* only.

## H3 · Event horizons of extreme-mass-ratio mergers (executed in session 4: see [11](11-event-horizon-merger.md))

> **Outcome (session 4).** Executed in [note 11](11-event-horizon-merger.md),
> by a different route from the sketch below.
>
> - For Λ > 0 (Hintz Thm 6.1), define the black-hole region relative to the
>   cosmological region. Then the merger holds, the region is connected, and
>   both horizons are localised. **Theorem D** proves this from Hintz's
>   stated estimates alone.
>   **Session 5:** the same holds for Λ = 0, using Hintz's 2026 stability
>   theorem for subextremal Kerr in fixed coordinates.
> - With the textbook $\mathscr I^+$ definition (Λ = 0), **Theorem E** gives
>   the merger exactly as defined in note 07, assuming one far-field
>   property (N) of the late development.
> - (A3) is not needed. (A2) is replaced by a far weaker escape hypothesis,
>   and only at Λ = 0.
> - The key observation: causal statements need only $C^0$ closeness of the
>   metric. Escape routes are timelike curves with margin; traps are time
>   functions.
>
> The text below is the session-1 design, kept as written.


1. **Obstacle.** Hintz (2024, Thm 6.1 / Rem. 6.2) constructs vacuum
   spacetimes in which a small Kerr hole of mass $\epsilon m$ crosses the
   horizon of a unit-mass, slowly rotating Kerr hole, which then settles to
   Kerr. "The construction of the event horizon, the structure of the black
   hole region" is left open there.
2. **Mechanism (conditional proposition).** For all small $\epsilon$, the
   components of $\mathcal B \cap \Sigma_t$ containing the small hole's trapped
   surfaces and the large hole's interior are **different before** the
   crossing and **the same after** it: an event-horizon merger in the sense
   of note 07.
3. **Proof sketch.** (i) Near the geodesic $C$ the rescaled metric is close to
   subextremal Kerr, so closed trapped surfaces exist there (trapping is an
   open condition). They lie in $\mathcal B$ if the spacetime is strongly
   asymptotically predictable (Hawking–Ellis Prop. 9.2.1). (ii) *Before
   crossing*: take a sphere $S_\delta$ of fixed small radius about $C(t_0)$
   with $r > 3M$. In the background, every point of $S_\delta$ emits a null ray
   that moves away from $C$ and has $dr/d\lambda > 0$. In Schwarzschild any
   such ray at $r > 3M$ escapes to infinity. Away from $C$, $g_\epsilon \to g$
   in $C^1$ on compact sets, so the perturbed rays stay close up to the time
   $T + 1$ where the stability theorem takes over. Then $S_\delta \subset J^-(\mathscr I^+)$,
   which separates the small hole's component. (iii) *After crossing*: on
   $\{r < r_H - 3\delta\}$ the background $r$ is a time function with a uniform
   margin. A causal curve can gain at most $O(\delta)$ in $r$ while within
   $\delta$ of $C$, so this region, which contains the small hole, lies in
   $\mathcal B$ and is connected.
4. **Inputs not verified here.** (A2) In the Kerr-stability spacetimes of
   Klainerman–Szeftel and Giorgi–Klainerman–Szeftel, outgoing null geodesics
   leaving a compact set at $r \ge r_H + \delta'$ within a fixed cone of
   directions reach $\mathscr I^+$, uniformly for small perturbations. This
   needs $C^1$ control of the metric along whole null rays, not just decay of
   curvature. (A3) The concatenated spacetime is strongly asymptotically
   predictable.
5. **Smallest calculation.** The background part of (ii) is analytic, as
   above. The decisive inputs (A2) and (A3) are properties of roughly
   1000-page stability proofs that we could not check in a session. That is
   why H3 was not chosen for execution.
6. **Criterion.** Continue if (A2) can be read off the main theorems of
   KS23/GKS22 (or Hintz–Vasy for Λ > 0, where the far region is simpler).
   Abandon if those theorems control the metric only in norms too weak to
   follow null geodesics.

## H4 · How the common horizon is born: a near-critical onset law (survives)

1. **Obstacle.** For vacuum, asymptotically flat, Λ = 0 data, the appearance
   of a common apparent horizon during evolution is known only numerically.
   Andersson–Mars–Metzger–Simon (AMMS 2009, Prop. 5.1) analyse a degenerate
   MOTS in one *given* spacetime. Nobody, as far as we found, has turned that
   into a quantitative statement about a *family* of data.
2. **Mechanism.** Time-symmetric Brill–Lindquist data with separation $d$ have
   a common minimal surface (the common apparent horizon of the $t = 0$ slice)
   only for $d \le d_c$. At $d_c$ an outer (stable) and an inner (unstable)
   common minimal surface merge in a **fold**, leaving a degenerate minimal
   surface $S_c$ whose Jacobi operator $L = -\Delta_S - \mathrm{Ric}(\nu,\nu) - |A|^2$
   has a positive kernel function $\phi_0$. Evolve with any lapse $N > 0$ and
   zero shift. Since $\gamma(t) = \gamma + O(t^2)$ and
   $K(t) = t\dot K + O(t^2)$ with $\dot K_{ij} = D_iD_jN - N R_{ij}$, a surface
   near $S_c$ has

```math
\theta_+ = H_S + t\,p + O(t^2), \qquad p = \mathrm{tr}_S \dot K .
```

   A Lyapunov–Schmidt reduction in the two parameters $(d, t)$ makes the
   common-MOTS condition near $S_c$ the scalar equation
   $\alpha s^2 + \beta(d - d_c) + \gamma t + \dots = 0$ for the amplitude $s$
   along $\phi_0$, with

```math
\beta = \langle \phi_0, \partial_d H\rangle, \qquad
\gamma = \langle \phi_0, p\rangle = -\int_{S_c} N\phi_0\,|A|^2\, dA < 0 .
```

   The second equality is the time-symmetric vacuum case of AMMS's identity.
   It follows from $L\phi_0 = 0$ by integrating by parts:
   $\int\phi_0 p = \int \phi_0(\Delta_S N + N\,\mathrm{Ric}(\nu,\nu)) = -\int N\phi_0|A|^2$.
   So the common MOTS first appears at

```math
\boxed{\; t_*(d) = (d - d_c)\,\frac{\langle \phi_0, \partial_d H\rangle}{\int_{S_c} N\phi_0\,|A|^2\, dA} + O\big((d - d_c)^2\big) \;}
\tag{*}
```

   to the **future**, for every lapse. Only the denominator depends on $N$,
   so the $\phi_0|A|^2$-weighted average of the proper time
   $\int N t_* \phi_0|A|^2 / \int\phi_0|A|^2$ is the same in every slicing, to
   leading order.

3. **Closest work.** AMMS (2009): the sign of $\gamma$ and the tangency of the
   marginally trapped tube to the slice, for a fixed spacetime. Horizon
   pretracking (Schnetter et al.) and Pook-Kolb et al. find the onset
   numerically. Hui & Lin (2025) and earlier authors compute $d_c$. Two
   searches found no statement of (*) or of the linear onset law;
   **potentially new**, and small.
4. **Assumptions.** Vacuum, Λ = 0, asymptotically flat, time-symmetric data:
   the primary target, near its threshold. The fold hypotheses (F) below are
   verified by a computer-assisted proof (Theorem C). The result is local in
   time and about apparent horizons, **not** event horizons.
5. **Calculation** ([`research/calc/bl_fold.py`](../research/calc/bl_fold.py),
   [`tests/test_bl_fold.py`](../tests/test_bl_fold.py)). Brill–Lindquist data
   with bare masses $\frac12$ (ADM mass $M = 1$), punctures at $z = \pm d/2$.
   Axisymmetric shooting for $H_S + tp = 0$.

| Quantity | Value | Check |
| --- | --- | --- |
| $d_c$ | $0.7661974282\,M$ | $= 1.53239486$ in units of the bare mass; Hui & Lin 2025: 1.5323948 (and 1.5323949 from their ref.) |
| area of $S_c$ | $49.101988\,M^2$ | ×4 $= 196.40795$; Hui & Lin: 196.40795 |
| AMMS identity $\langle\phi_0,p\rangle / (-\int N\phi_0\lvert A\rvert^2)$ | $1 + 3\times10^{-7}$ ($N = 1$) | holds for both lapses |
| $\phi_0$ on $S_c$ | between 5.05 and 18.40 | one sign, as a principal eigenfunction must be |
| slope of (*), $N = 1$ | **2.385183** | direct solution: 2.384676, 2.384929, 2.385056 at $d - d_c = 2, 1, 0.5 \times 10^{-4}$; Richardson limit 2.385183 |
| slope of (*), $N = \psi^{-2}$ | **10.644372** | direct: 10.642310, 10.643340, 10.643855; Richardson limit 10.644370 |

   Robustness: $d_c$ moves by $2\times10^{-10}$ when the ODE tolerance is
   loosened 100-fold. The slope moves by $3\times10^{-7}$ when the quadrature
   points are doubled and the finite-difference step tripled. Individual
   horizons: for **every** $d \in [0.76, 0.78]$, interval arithmetic certifies
   $H_S < 0$ on the sphere of radius 0.01 and $H_S > 0$ on the sphere of
   radius 0.3 about each puncture. So each puncture has its own minimal
   surface (Andersson–Metzger Thm 3.1), and the two are disjoint.

   **Is the new horizon trapped, not just marginally outer trapped?** On it
   $\theta_- = 2t_*p$ to leading order. For $N = \psi^{-2}$, $p < 0$ on all of
   $S_c$ (between −0.111 and −0.019), so $\theta_- < 0$: a genuine marginally
   trapped surface. For geodesic slicing ($N = 1$), $p$ ranges over
   $[-0.576, +0.106]$, so $\theta_- > 0$ on part of the surface. That is
   another reminder that quasi-local horizons depend on the slicing.

**The fold, verified by computer.** Let $F(z_0, d, t)$ be the *shooting
function* for lapse $N = 1$. Start the meridian curve of an axisymmetric surface
at the axis point $(0, z_0)$ with horizontal tangent. Solve
$\theta_+ = H_S + t\,p = 0$ until the curve meets $z = 0$, and record
$F = \alpha + \pi/2$ there, where $\alpha$ is the tangent angle. A
reflection-symmetric common MOTS is exactly a zero of $F$.

**Theorem C (computer-assisted).** In the box
$X = (0.6326595528058354, 0.766197428240025) \pm 10^{-12}$, the system
$F = \partial_{z_0}F = 0$ at $t = 0$ has exactly one solution $(z_{0,c}, d_c)$, with

```math
d_c = 0.766197428240 \pm 1.1\times10^{-13}\,M \quad (= 1.532394856480 \pm 2\times10^{-13} \text{ bare masses}).
```

At that point

```math
\partial^2_{z_0}F \in [240.87, 242.70], \qquad
\partial_dF \in [1.2315, 1.2364], \qquad
\partial_tF \in [-0.51816, -0.51655] .
```

Along the critical curve, three things hold. The normal component of
$\partial X/\partial z_0$ is positive. After the axis, $\sin\alpha < 0$ and
$\rho > 0$. And $z > 0$ at the point where the height becomes the independent
variable.

*How.* [`research/calc/fold_proof.py`](../research/calc/fold_proof.py) uses
Arb ball arithmetic throughout, in four stages:

- **Off the axis:** a fixed-point enclosure of the first segment, where
  $\alpha = sA$ and $A = s^{-2}\int_0^s\sigma\Phi$ is a weighted average.
- **Validated Taylor steps:** a Picard a-priori enclosure and a Lagrange
  remainder on each step, using the substitution $\mu = \rho\sin\alpha$. That
  substitution removes the $1/\rho$ singularity exactly, since
  $(\rho\sin\alpha)' = \rho\cos\alpha\,G$.
- **Derivatives:** jets in $(z_0, d, t)$ carry $\partial_{z_0}F$, $\partial_dF$,
  $\partial_tF$, $\partial^2_{z_0}F$ and $\partial_{z_0}\partial_dF$.
- **The fold itself:** a Krawczyk test on $(F, \partial_{z_0}F)$.

The implementation is checked four ways:

1. against ODEs with exact solutions, at several orders;
2. jets against SymPy derivatives;
3. $F$ and its derivatives against the independent floating-point code
   `bl_fold.py`, and against finite differences;
4. against a second, finer step schedule: every enclosure overlaps.

It reproduces Hui & Lin's $d_c$ (1.5323948, and 1.5323949 from their ref.)
as an interval of width $4\times10^{-13}$.

> **History, for trust.** A first "verified" Krawczyk run (session 2) was
> invalid. python-flint truncates power series at length 10 unless told
> otherwise, so Taylor terms and remainders were dropped. The known-solution
> tests exposed it. Session 2 then attributed the next failure to a structural
> barrier at the axis. That diagnosis was **wrong**. The blow-up came from
> low-order Taylor steps once steps grew. The fix was a piecewise schedule:
> order 8 with $h \le s/64$ near the axis, order 18 with
> $h \le \min(s/32, 2^{-9})$ beyond.

**Corollary (the fold hypothesis (F) holds).** (F) asserts four things: a
degenerate common minimal surface $S_c$ exists at $d = d_c$; its Jacobi
operator has a simple zero eigenvalue with positive eigenfunction; the fold is
nondegenerate; and it opens toward $d < d_c$. Indeed, $S_c$ is an embedded,
axisymmetric, reflection-symmetric minimal sphere enclosing both punctures. The field
$\phi_0 = \psi^2\, n\cdot\partial X/\partial z_0$ solves the Jacobi equation on
the half-surface. Because $\partial_{z_0}F = 0$, it meets $z = 0$ with the
symmetric boundary condition, so it extends to a Jacobi field on $S_c$. It is
positive, hence the principal eigenfunction, so $0$ is a simple eigenvalue.
The fold is nondegenerate ($\partial^2_{z_0}F > 0$, $\partial_dF > 0$) and
opens toward $d < d_c$.

**Proposition B.** Take the vacuum development of the Brill–Lindquist data
above, foliated near $S_c$ by slices $\Sigma_t$ with zero shift and a lapse
$N_d > 0$ depending smoothly on $(d, t, x)$. Then there are $\delta_0, t_0 > 0$
and a neighbourhood $\mathcal U$ of $S_c$ such that, for
$d_c < d < d_c + \delta_0$ and $0 \le t \le t_0$:

1. $\Sigma_0$ contains two disjoint minimal surfaces, each enclosing one
   puncture (certified for $d \in [0.76, 0.78]$), and no MOTS in $\mathcal U$;
2. $\Sigma_t$ contains a MOTS in $\mathcal U$ enclosing both punctures **if and
   only if** $t \ge t_*(d)$, where

```math
t_*(d) = c_N\,(d - d_c) + O\big((d - d_c)^2\big), \qquad
c_N = \frac{\langle\phi_0, \partial_dH\rangle}{\int_{S_c} N\phi_0|A|^2\,dA} .
```

Rigorously (Theorem C), $c_1 = \partial_dF/(-\partial_tF) \in [2.3768, 2.3936]$
for geodesic slicing and $c_{\psi^{-2}} \in [10.606, 10.682]$ for $N = \psi^{-2}$
(`certify_slopes()`). Numerically $c_1 = 2.385183$ and $c_{\psi^{-2}} = 10.644372$.
Units are $G = c = M = 1$. The full proof is in
[10](10-onset-theorem-proof.md).

*Proof.* (a) *Smooth dependence.* By Cauchy stability the development and the
slicing depend smoothly on $(d, t)$ near $(d_c, 0)$. So does the shooting
function $\tilde F(z_0, d, t)$ of the true slices. At $t = 0$ it equals $F$,
and $\partial_t\tilde F = \partial_tF$ there, because
$\gamma(t) = \gamma + O(t^2)$ and $K(t) = t\dot K + O(t^2)$. (b) *Fold
analysis* (elementary, given Theorem C). Since $\partial^2_{z_0}\tilde F > 0$,
the implicit function theorem solves $\partial_{z_0}\tilde F = 0$ by
$z_0 = \zeta(d, t)$. Set $m(d, t) = \tilde F(\zeta, d, t)$. Then $m(d_c, 0) = 0$,
$\partial_d m = \partial_d F > 0$ and $\partial_t m = \partial_t F < 0$.
Because $\tilde F$ is convex in $z_0$ near $z_{0,c}$, it has zeros there exactly
when $m \le 0$. That gives item 2 for symmetric MOTSs, with slope
$\partial_dF/(-\partial_tF)$. The same slope also equals the
Lyapunov–Schmidt form (*): both describe the same set. Numerically the two
agree to six digits. (c) *Every MOTS in $\mathcal U$ is symmetric.* By the
Corollary the kernel of the stability operator is one-dimensional and
symmetric. In the Lyapunov–Schmidt reduction (normal graphs
$u = s\phi_0 + w$, with $w$ solved by the implicit function theorem), the
solution for given $(s, d, t)$ is unique. The reflected and rotated solutions
are solutions too, so they coincide with it. (d) Item 1 follows from (b) at
$t = 0$, together with the individual-horizon certificates. ∎

What is rigorous: Theorem C and step (b). What is a standard argument written
as a sketch: steps (a) and (c), i.e. function spaces, smooth dependence of the
zero-shift slicing, and equivariance in the Lyapunov–Schmidt reduction. That
$\Sigma_0$ has **no** common MOTS anywhere, not just in $\mathcal U$, is still
numerical evidence (shooting, and agreement with Hui–Lin's $d_c$).

**What H4 adds to the primary target.** It is an existence statement inside
the target class: vacuum, Λ = 0, asymptotically flat developments in which
the apparent horizon changes from two components to a common one along an
explicit foliation, with a formula for when. It does not address the event
horizon, long-time existence, $\mathscr I^+$ or Kerr asymptotics. Under strong
asymptotic predictability, the marginally trapped common surface of the
$N = \psi^{-2}$ slicing lies in the closure of the black hole region
(Hawking–Ellis §9.2, not re-checked). So the two holes' black-hole regions
are joined by time $t_*$, but whether they were separate *event* horizons at
$t = 0$ is not decided.

6. **Criterion and verdict.** The hypothesis would have been abandoned if the
   two routes to $t_*$ disagreed, or if $\gamma$ had come out positive
   (common horizon appearing to the past). Neither happened, so H4 is the
   route to develop.

## Physical picture

At a moment of time symmetry the holes are momentarily at rest. One might
guess the common horizon appears once they have fallen a little closer. The
separation shrinks like $t^2$, so that guess predicts an onset time of order
$\sqrt{d - d_c}$. The calculation gives $t_* \propto (d - d_c)$, much sooner
when $d$ is close to $d_c$. Trapping is decided by the null expansions, and
those respond at first order to the extrinsic curvature $K \approx t\dot K$,
the rate at which space is changing, before the 3-geometry itself has
changed. Light rays notice the infall before rulers do.

## Novelty check

| Item | Assessment |
| --- | --- |
| Lemma A | elementary; probably folklore; not found stated |
| H2 obstruction | a negative result within one ansatz; not found stated; not important by itself |
| AMMS identity in the time-symmetric case | **known** (AMMS 2009 Prop. 5.1 with $W = \lvert A\rvert^2$); re-derived and checked numerically |
| onset law (*) and Proposition B | **potentially new**; bounded search found no equivalent |
| computer-assisted fold (Theorem C): rigorous $d_c$ to $10^{-13}$ | **potentially new as a rigorous statement**; the value of $d_c$ itself is known numerically |
| BL numbers $d_c$, area | **known**; reproduced to all quoted digits |
| Theorem D (H3, [11](11-event-horizon-merger.md)): event-horizon merger and connected black-hole region in Hintz's Λ > 0 EMR spacetimes | **potentially new**; Hintz lists both as open; elementary given his estimates |

## Discussion prompts

- Proposition B says the common horizon appears at a time linear in
  $d - d_c$. Could a numerical-relativity code, with horizon pretracking,
  test the slope 2.385 in geodesic slicing directly?
- Lemma A and H2 both fail because a constant-mean-curvature umbilic slicing
  is too rigid. What is the least rigid slicing in which an exact merger
  might still be written down?
- H3 needs only that light from outside the small hole escapes. Why was that
  so hard to prove, when it is "obviously" true? (Note 11: it stops being hard
  once timelike curves replace light rays.)
