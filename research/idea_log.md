# Idea log

Hypotheses, the tests run on them, and what happened, including failures,
technical dead ends and ideas dropped before testing. Newest last within each
section. Full write-ups of H1–H4 are in [`notes/09`](../notes/09-merger-hypotheses.md).

## Hypotheses

### H1 · Kastor–Traschen mechanism at Λ = 0 (explicit metric) — failed

- **Idea.** Keep KT's zero-shift umbilic slicing, drop Λ.
- **Test.** Momentum constraint for $g = -N^2dt^2 + \Psi^4\delta$:
  $G_{ti} = -\frac23 N\partial_iK$ (exact identity).
- **Result.** Lemma A: with asymptotic flatness, $K \equiv 0$ and the slices
  freeze. The obstruction is precise: AF forces $K \to 0$, and the momentum
  constraint spreads that to the whole slice.
- **What would have to change:** shift, trace-free $K$, or slices not
  reaching spatial infinity (→ H2).

### H2 · Hyperboloidal (Milne) Kastor–Traschen (explicit metric, unconventional) — failed

- **Idea.** Milne slices are umbilic with $K = 3/T \ne 0$ and end on
  $\mathscr I^+$, so they evade Lemma A. A contracting Milne wedge converges to a
  point, so comoving extremal holes would merge in asymptotically flat space.
- **Test.** Reduced Einstein–Maxwell–Λ system for
  $g = -W^{-2}dt^2 + W^2h_{-1}$, $A = W^{-1}dt$, arbitrary $W(t, \vec x)$.
- **Result.** $E_{ij} \propto 2WW_{tt} + 3W_t^2 - \Lambda - W^{-2}$. The
  curvature term $W^{-2}$ forces $W$ to be independent of $\vec x$: no centres.
  Sanity check: $W = \sqrt{2t}$ is Milne.
- **What would have to change:** extra matter cancelling $W^{-2}$, or a
  lapse not equal to $1/W$, which abandons Majumdar–Papapetrou force balance.

### H3 · Event horizons of Hintz's EMR mergers (existence) — designed in session 1; **executed in session 4 (survives)**

- **Idea.** Causal stability: escaping null rays from a sphere around the
  small hole before crossing; $r$ as a time function inside the large hole
  after crossing.
- **Status (session 1).** Proof sketch written; inputs (A2) (null escape at late times in
  the Kerr-stability spacetimes) and (A3) (strong asymptotic predictability)
  could not be checked that session.
- **Why not chosen then.** The deciding steps sit in 1000-page papers, not in a
  calculation we can run.
- **Session 4: what changed.** Replace escaping *null rays* by escaping
  *timelike curves with margin*. Those survive any $C^0$-small change of the
  metric, and Hintz gives $C^0$ closeness on both scales. Replace "trapped
  surfaces lie in $\mathcal B$" (which needed (A3)) by a direct argument:
  $r$ is a time function inside, and a gain lemma bounds how much $r$ can
  grow near the small hole.
  - For Λ > 0, measure escape against the cosmological region, which lies
    inside Hintz's domain. Then everything needed is in Hintz's Theorem 6.1:
    Theorem D, including connectedness of $\mathcal B$ (notes/11).
  - For Λ = 0, the escape and trap inputs move to the initial slice of the
    Kerr-stability development. KS prove them in their own coordinates
    (Theorem D₀, conditional).
- **Snag found and fixed.** Interpolating Fermi slices and $t_*$-slices with
  an ordinary cut-off fails, because the boost of $C$ makes the error term of
  order 1. A cut-off in $\log|x|$ makes it $O(1/K_0)$ (Lemma 4).

### H4 · Near-critical onset law for the common horizon (existence + controlled) — survives

- **Idea.** Unfold the fold of common minimal surfaces in BL data in the two
  parameters (separation, time). AMMS's sign identity fixes the direction of
  time; Lyapunov–Schmidt gives the onset time to first order.
- **Tests run.** (1) $d_c$ and the critical area against Hui & Lin 2025:
  agreement to all quoted digits. (2) AMMS identity on $S_c$ to $3\times10^{-7}$.
  (3) Slope of (*) against direct solutions of $H_S + tp = 0$, two lapses:
  agreement to 6 significant figures after Richardson extrapolation.
  (4) Resolution and tolerance checks. (5) Individual horizons certified for
  all $d \in [0.76, 0.78]$. (6) Sign of $p$: marginally trapped for
  $N = \psi^{-2}$, not for $N = 1$.
- **Continue/abandon criterion.** Abandon if the two routes disagreed or if
  the horizon appeared to the past. Neither happened.

## Ideas considered and not pursued (this session)

- **Self-similar mergers.** A homothety rescales the ADM mass, so an
  asymptotically flat self-similar vacuum spacetime with a homothety acting
  nontrivially at infinity has zero mass and is flat (positive mass theorem).
  Not pursued: obstruction too cheap to be informative.
- **Boost–rotation-symmetric two-hole spacetimes.** Uniformly accelerated
  holes move on hyperbolae that never meet, and the boost Killing field acts
  as a time translation on the configuration, so there is no merger. A
  symmetry that removes the dynamics; dropped.
- **Kerr–Schild superposition** $g = \eta + \phi_1k_1k_1 + \phi_2k_2k_2$: vacuum
  Kerr–Schild metrics are algebraically special with a single repeated
  principal null direction, which two generic centres break. Dropped.
- **Tidal ($O(\epsilon^2)$) correction to the Emparan–Martínez horizon.** A
  controlled approximation, but a heavy matched-asymptotics calculation.
  Deferred.
- **EMR limit of the Kastor–Traschen event horizon versus Emparan–Martínez
  for extremal Reissner–Nordström.** A nice cross-check of EMR universality
  inside an exact merger, but it needs KT event horizons. Their teleological
  boundary condition lies beyond the coordinate patch (BHKT's extension is
  non-unique). Deferred.

## Technical dead ends (kept so nobody repeats them)

- **Kerr Kretschmann scalar in SymPy**: did not finish in 10 minutes even at
  fixed rational $M$, $a$ (the 4-index contraction is the bottleneck). Ricci
  flatness is checked instead.
- **Two-centre KT residual with `simplify`**: over 90 s per run. Replaced by
  the general-$U$ reduction (a proof) plus exact evaluation at rational
  points (a diagnostic).
- **Shooting for individual KT horizons**: they sit in near-cylindrical
  extremal throats, where axis-to-axis shooting is ill-conditioned. Existence
  is settled by barriers instead.
- **Barriers for the common KT horizon**: coordinate spheres are never fully
  outer-trapped; the halfway surface between inner and outer common MOTSs is
  not either; uniformly scaling the outer MOTS inward leaves untrapped spots,
  because scaling is not a push along the principal eigenfunction. What
  worked: the outer common MOTS of an earlier slice, evaluated on a later
  slice.
- **Fitting MOTS shapes from raw ODE steps**: linear interpolation of solver
  output gives errors that $f''$ amplifies; dense output fixes it.
- **Interval arithmetic with an interval separation $d$**: writing a point as
  $d/2 + R\cos\theta$ and subtracting $d/2$ again loses the cancellation
  (dependency problem) and gives $-\infty$. Displacements must be computed
  relative to each centre directly.
- **python-flint power series are capped at length 10 by default** (`ctx.cap`).
  Longer Taylor series are truncated *silently*. This produced a false
  "Krawczyk verified" fold certificate in session 2. Set `ctx.cap` explicitly,
  and test every validated integrator on ODEs with known solutions at several
  orders (`tests/test_fold_proof.py`).
- **Naive (Moore-type) validated integration near the axis singularity: a
  wrong diagnosis (session 2), corrected in session 3.** Session 2 blamed a
  structural barrier (final width about $C s_{\rm end}^2$, independent of
  $s_1$). Measurements in session 3 refuted that: initial radii up to
  $10^{-22}$ at $s_1$ are not visibly amplified. The blow-up came from the
  per-step remainder of *low-order* steps, jumping each time the step size
  doubled. Fix: a piecewise schedule, order 8 with $h \le s/64$ near the axis
  and order 18 with $h \le \min(s/32, 2^{-9})$ beyond. That made the point
  evaluation accurate to about $10^{-14}$ and the Krawczyk test close
  (Theorem C). Lesson: measure before diagnosing.
- **Parameter boxes wrap badly** (box radius $r = 10^{-10}$ inflates
  $\partial^2_{z_0}F$ to $\pm 480$). The Krawczyk box therefore has to be tiny
  ($10^{-12}$), which only works once the point evaluation is enclosed far
  more tightly than $r$. Taylor models in the parameters would allow bigger
  boxes, but were not needed.
- **Arb's ball `cos`/`sin`** are crude when the radius is comparable to the
  midpoint (radius $10^{-9}$ for $\cos[-3.6\times10^{-6}, 0]$, where $10^{-12}$ is
  possible). Use alternating-series enclosures for small arguments.
