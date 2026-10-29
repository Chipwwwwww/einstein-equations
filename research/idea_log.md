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
- **Session 5.** Hintz's 2026 stability theorem for the whole subextremal
  Kerr family (arXiv:2606.28253) works in fixed coordinates. It states the
  Λ = 0 merger construction as its consequence (3).
  - Theorem D therefore extends to Λ = 0 with "escape to arbitrarily large
    distances". The Klainerman–Szeftel anchoring problem disappears.
  - For the textbook $\mathscr I^+$ version (Theorem E), one gap remains:
    (N), that far points on the late initial slice reach $\mathscr I^+$.
  - Lemma 5 proves (N) where the data are exactly Kerr. Our signals,
    however, travel through the perturbed region.
  - Lemma 6 then derives (N) from the 2026 paper's gauge, in which the
    outgoing null component $h(\partial_0, \partial_0)$ decays faster than
    $r^{-1}$. The asymptotic boost only touches $(dx^1)^2$ (Lemma 3.19).
    Kerr's own outgoing component is $-a^2\sin^2\theta/4r^2$ (checked). This
    rests on our reading of the paper.
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

### Gap 3 attempt (session 6): no common minimal surface at d = 0.78 — partial

- **Reduction (rigorous).** The outermost minimal surface is unique, hence
  symmetric, and its components are spheres.
  - The throat barrier ($H < 0$ on small spheres about each puncture,
    $s \le 0.155$) and the far barrier ($H > 0$ for $|x| \ge 1.05$) are
    certified by interval arithmetic.
  - So a common one is generated by a shooting curve from $(0, z_0)$,
    $z_0 \in [0.545, 1.05]$, meeting $z = 0$ once at a right angle.
- **Landscape (numerical).** There are four regimes: throat-diving; necked
  curves just above the individual horizon's top; plane hits with
  $F \ge 0.017$; and curves leaving the ball.
- **Obstruction.** Just above $z_{\rm top}$ the curves pass a neck of radius
  $\approx 4(z_0 - z_{\rm top})$. This linear law is the signature of a
  passage near a hyperbolic fixed point.
- **Tried and failed, so as not to repeat them:**
  - **A mean-convex level-set function.** The level sets of $\psi$ fail near
    the saddle in the equatorial directions, because $\psi$ is harmonic. A
    custom foliation must keep outward-pushoff mean convexity near the
    individual horizons and still pinch at the saddle. Near those horizons
    no $H < 0$ barrier exists outside them, so this does not avoid the necked
    family either.
  - **Stability of the outermost surface.** Its Jacobi field cannot change
    sign along the profile. Numerically it does change sign, but only just
    after the neck.
  - **The Penrose bound $A \le 16\pi$.** The part of the surface before the
    neck has area of order $2A(S_+) \approx 26 < 16\pi$, so it is not
    violated.
- **Refined picture.** The individual horizon's profile joins the two axis
  fixed points $(0, z_{\rm top})$ and $(0, z_{\rm bot} = 0.2126)$, and read
  backwards it is again a solution. So there is a heteroclinic cycle between
  two saddles, with eigenvalues $\pm1$ in $\sigma$-time.
  - Curves just above $z_{\rm top}$ go down along it and neck near the
    bottom. They come back up outside it ($z_0 - z_{\rm top} = 10^{-4}$:
    the curve rises to $z \approx 0.57$ near the top), and may neck again:
    multi-layer, onion-like embedded spheres.
  - Crossing parity with rays from $p_+$ does not shortcut this: the parity
    changes only after a neck.
- **Route.** In $\sigma$-time ($ds = \rho\,d\sigma$) the system is regular.
  Use isolating blocks at both saddles, validated transport along the two
  connecting orbits, and a validated sweep elsewhere: a
  heteroclinic-cycle (Shilnikov-type) argument.

### Gap 3 attempt (session 7): thin necks are unstable

- **Idea.** Session 6 planned an isolating-block passage argument at the
  saddles. We did not build one. A neck is just as useless to $S^*$ if
  it is *unstable*. Near a waist of radius $\lambda$ the surface is a
  catenoid up to $O(\lambda)$. The catenoid's dilation Jacobi field vanishes
  at $\sigma = \pm1.1997$. So it is enough to certify, for all small
  $\lambda$ at once (λ = 0 included), that the $\lambda$-variation Jacobi
  field changes sign on both sides within $|\sigma| \le 1.5$.
  - Thin necks are then excluded wholesale.
  - Multi-layer necked spheres never have to be followed.
  - The heteroclinic cycle drops out.
- **Approach lemma.** A curve heading for the axis ($u < 0$, $\rho \le \rho_e$)
  either closes on the axis or reaches a waist. The waist radius is at most
  $|v_e| + K\rho_e^2/2$, and its height moves by at most an explicit
  $O(\rho_e)$. This hands the curve to the neck lemma with no further
  integration. It replaces the passage argument.
- **Lohner's method** for the sweep. Plain interval boxes blow up (radius
  $10^{20}$ after 189 steps) on curves that wind around the individual
  horizon. Lohner's parallelepipeds stay at about $10^{-12}$ for point
  starts.
- **Tried and failed, so as not to repeat them:**
  - **Remainder over a wide a-priori box.** In the first Lohner version, the
    Lagrange remainder $c_{p+1}(B)h^{p+1}$ over a wide $B$ dominated: radius
    $10^4$, against a true value of $10^{-8}$. Fix: cap each step at about
    $0.004$ in flat arclength.
  - **Narrow approach section.** Requiring the near-top box to reach
    $\rho \le 0.004$ with $u < 0$ failed. Lohner's enclosure there is about
    300 times the true spread, because errors in the unstable direction grow
    like $e^\sigma$ while the state shrinks like $e^{-\sigma}$. Fix: stop
    earlier (any step with $\rho < 0.02$ will do), at the cost of a wider
    neck rectangle ($\lambda \le 0.00488$).
  - **Wide neck pieces.** A neck piece with $\lambda \in [0.0045, 0.005]$ and a
    $z_w$-width of $0.0025$ gives $\phi(\pm1.5)$ only to $\pm 1$. The
    overestimate grows with the *product* of the two widths. Pieces of
    $2^{-11} \times 2^{-9}$ work, with margin $0.53$.
  - **Not needed:** the $z_0$-Jacobi field (event E4 of the plain-interval
    version). Its sign change comes only after the neck, which N already
    covers.
- **Falsification checks run.**
  - At $d = 0.76 < d_c$, where two common surfaces exist, the sweep code
    leaves the $2^{-14}$ boxes containing them unexcluded and excludes a box
    between them.
  - The neck certificate fails, as it must, when asked for a sign change at
    $\sigma = 1.0 < 1.1997$.

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
