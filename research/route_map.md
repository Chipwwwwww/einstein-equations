# Current research roadmap: an exact vacuum binary-black-hole merger

Updated 2026-11-01. This is the single active roadmap. Progress and evidence
are consolidated in [status.md](status.md). The older route-by-route
exploration is preserved in [the historical roadmap](history_route_map_through_silent_flow.md);
its old “next tasks” and superseded constants are not current instructions.

## Success criterion

Construct or rigorously specify a 4D, asymptotically flat solution of
$R_{\mu\nu}=0$, with no matter, cosmological constant, or supporting
strut, containing two initially separate positive-mass black holes that
merge into one. State which exterior and interior regions the construction
covers, its radiation assumptions, and its future asymptotics. Prefer an
explicit metric; accept a provably convergent effective representation only
with precise domain and error control.

A local vacuum patch, a common MOTS, or an extreme-mass-ratio limiting
horizon is an intermediate result. None alone satisfies this criterion.
A global metric may require an atlas rather than one nonsingular coordinate
formula. That is a coordinate issue, not a reason to weaken the geometry.

## Chosen construction and checkpoints

Use equal-mass, nonspinning, axisymmetric head-on Brill–Lindquist data as the
first finite-ratio laboratory. Retain the EMR gluing literature as the
parallel conceptual route toward a global existence construction; this
roadmap does not assign parallel agents.

| Step | Mathematical task | Current status | Completion test |
| --- | --- | --- | --- |
| 1 | Fix exact data, signs, and merger definition | available | distinguish local metric, MOTS, and event horizon in every theorem |
| 2 | Certify the critical surface and its full stability complement | **completed for the stated BL fold in this session** | positive kernel, exact fold enclosure, all-mode gap $1/(16M^2)$ and inverse bound |
| 3 | Quantify the nonlinear surface problem in the true evolved data | active next target | explicit graph norm, elliptic inverse constants, and pointwise remainder bounds |
| 4 | Exhibit a common trapped/MOTS configuration at explicit $(d,t)$ | open | strict sign/barrier certificate or validated true-evolution bifurcation; verify earlier separate surfaces |
| 5 | Continue a meaningful exterior through merger | open | constraint-preserving horizon-penetrating evolution, a posteriori error bounds, domain-of-dependence coverage |
| 6 | Establish global event-horizon merger and final exterior | open | escape/non-escape barriers tied to future null infinity and late-time control |
| 7 | Seek a simplifying exact metric representation | open throughout | full Ricci check, regular-axis/asymptotic checks, convergence and coverage, comparison against known classes |

Steps 3–4 provide a verifiable geometric checkpoint. Step 5 remains the main
construction bottleneck; step 6 is essential for the ultimate claim.
Success at step 4 must not be reported as completion of steps 5–6.

## Next concrete research calculation

Build a nearby-surface expansion map
$\Theta(u;d,t)$ for normal graphs over $S_c$, using the actual certified
metric and extrinsic curvature. At the critical datum,
$D_u\Theta=L$ and the kernel is one-dimensional.

1. Use the new inverse bound to control the kernel-orthogonal correction,
   while deriving explicit elliptic regularity constants in a norm that
   controls the surface and expansion pointwise.
2. Bound the quadratic remainder on a stated normal-graph ball, including
   parameter and time derivatives of the metric.
3. Prefer a strict outer-trapped witness and an outer untrapped barrier
   if this proves existence with less work than an exact formation-time
   enclosure. Keep the local/global distinction explicit.
4. Choose a concrete $d>d_c$ only after those constants determine a valid
   interval. The heuristic linear formation law is not a certificate.

**Decision rule:** if these constants are dominated by the existing
short-time metric enclosure, switch effort to step 5's hyperbolic
continuation rather than generating more Taylor orders.

## Global construction direction

For comparable masses: a frame or generalized-harmonic formulation with
real Sobolev energy control and horizon-penetrating gauges is the principal
continuation route. A proof must include the constraint subsystem, the
whole domain of dependence, and a restart/error budget. Merely removing
algebraic terms from the equations does not remove propagation or derivative
loss. No commuted-energy certificate is currently available here.

For extreme mass ratios: audit and specialize Hintz's true-solution gluing
construction and its stability attachment, keeping separate (i) existence
of an exact solution, (ii) effective evaluation, and (iii) the event-horizon
topology. The repository's older causal claims remain conditional until
their source-dependent hypotheses are fully verified. Do not treat them as
an explicit metric or silently extend them to equal masses.

## Directions currently deprioritized

- Higher Taylor orders without a new continuation mechanism.
- Fixed-curl-sector boundedness as a convergence test: it is automatic on a
  regular background and gives no uniform summability bound.
- Silent-flow conjugation as a cure for spatial propagation: its principal
  symbol remains similar to the original one.
- Single Kerr–Schild/algebraically special ansätze for the specified BL
  development: its known type-I curvature obstructs that restricted route.
- Finite-order matching or Padé fits labelled “exact.”

These are scoped decisions, not no-go theorems for all merger data or all
mathematical methods. A new route is worth pursuing when it supplies a
specific structure that survives a nontrivial exact test and improves a
construction or global estimate.
