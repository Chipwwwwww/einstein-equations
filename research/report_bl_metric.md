# Report: an evaluable BL metric, and a bounded search for a closed form (session 12)

*2026-10-31. Branch `claude/bl-metric-evaluator`, based on `7c4bff6`.
Mathematics: [`notes/16-bl-metric.md`](../notes/16-bl-metric.md).*

## Metric construction (the primary deliverable)

**What can be evaluated.** The vacuum development of equal-mass, time-symmetric Brill–Lindquist
data, in geodesic slicing. The coordinates are comoving and Cartesian, so the metric is
$g = -dt^2 + \gamma_{ij}(t, x)dx^idx^j$ with $M_{\rm ADM} = 1$ and the punctures at $(0, 0, \pm d/2)$.
`research/calc/bl_metric.evaluate_metric(t, x, y, z, d, ...)` returns certified ball enclosures of:

- every component of $g_{\mu\nu}$ and $g^{\mu\nu}$;
- every first derivative $\partial_\lambda g_{\mu\nu}$;
- every second derivative $\partial_\kappa\partial_\lambda g_{\mu\nu}$, including time derivatives.

**Domain.**

- $(t, x, d) \in (-T_{\min}, T_{\min}) \times U \times I$.
- $T_{\min} = 5.556\times10^{-4}M$.
- $U$ is a compact neighbourhood of the critical surface: the profile boxes inflated by 0.02. It
  contains the $0.02$-tube about $S_c$.
- $I = [d_c - 10^{-13}, 0.767578125]$.

Outside this domain the evaluator refuses, or on request returns an `uncertified` truncated sum.
On the whole certified domain $g$ is Lorentzian.

**Error guarantees** (per evaluation, itemised):

1. coefficient rounding, about $10^{-34}$ at 128 bits;
2. no spatial approximation (the jets are exact);
3. time-truncation bounds for values and for each kind of derivative (notes/16 §3);
4. the radius of $d$, and the imported certificate.

Tolerance requests either succeed or raise an error naming the bound that is too large.

**Demonstration** ($d = 0.7662$, order 5, three points of $U$ including the axis, two patches).

| $t/T_{\min}$ | metric | first derivatives | second derivatives |
| --- | --- | --- | --- |
| 0.1 | $1.1\times10^{-14}$ | $9.6\times10^{-11}$ | $2.9\times10^{-6}$ |
| 0.25 | $6.9\times10^{-10}$ | $9.7\times10^{-7}$ | $1.2\times10^{-2}$ |
| 0.5 | $3.5\times10^{-6}$ | $1.2\times10^{-3}$ | $8$ |

The entries are radii. Each evaluation takes 5–9 s. The low orders are exact and closed form
(notes/16 §4): $\gamma = \psi^4\delta - t^2\mathrm{Ric} + t^4(\tfrac56\mathrm{Ric}\cdot\mathrm{Ric} - \tfrac16\lvert\mathrm{Ric}\rvert^2\gamma_0 - \tfrac1{12}\Delta\mathrm{Ric}) + O(t^6)$.

**What it is.** An exact metric of category B: the complete infinite series converges and solves
$R_{\mu\nu} = 0$ (note 15). The returned values are certified approximations of it. It is a
local development of asymptotically flat data; the patch itself says nothing about
asymptotics.

## Audit and the certified domain

- **Rigor gap repaired.** Scalar bounds in the session-11 certificate used round-to-nearest
  floats: cells, $d$-hull, $\sup\lvert\psi\rvert^{-4}$, the $G$-box radius, $C$ and $T_{\min}$. They
  are now outward rounded. The classic certificate, replayed, moves by $-5.6\times10^{-12}$
  (relative), so $T_{\min} = 2.988\times10^{-4}$ stands.
- **The same solution.** The Taylor coefficients of the certificate's own first-order system
  (built from its `ricci_parts` and `quad_K`, no parity assumed) equal the evaluator's to
  $10^{-28}$.
- **A sharper certificate.** Exact norms of the principal part, $\lVert A(I, \cdot)\rVert = 2$ and
  $\tfrac12\lVert T\rVert = 2.524$, give $T_{\min} = 5.556\times10^{-4}$ (+86%).
- **Limits of this method.**
  - For any weights $C \gtrsim 0.84$, so $T \lesssim 10^{-3}$.
  - Localisation gives no gain.
  - $s_0$ is capped by the punctures' complex singular cones.
  - Recentering cannot beat $s_0/(16C)$.

## B → A experiments

**No category-A expression was found.**

| | Candidate class | Outcome | How decided |
| --- | --- | --- | --- |
| A | universal covariant $t^4$ coefficient | holds (proved); structural only | exact rational arithmetic at 3 points, plus a derivation |
| C1 | conformally flat geodesic slices | **excluded** at $t^2$ | Cotton-tensor enclosure excludes 0; Schwarzschild control |
| C2 | silent flow (ODEs along worldlines; every local-free-fall ansatz) | **excluded** at $t^3$ | $B$-enclosure excludes 0; Schwarzschild control |
| D | Padé in $t^2$ | exploratory: real poles at $t \approx 2.1$–$2.8$, a crunch-type gauge singularity | floating point |

Earlier work (note 14) excluded the algebraically special forms (Petrov type I).

## Horizons (only as far as the physical claims need)

No merger is certified at any explicit $d$. At $d = 0.7662$ the linear law puts the common MOTS at
$t_* \approx 0.011\,T_{\min}$, inside the certified slab. The range $\delta_1$ of note 15's
Corollary J is still not explicit, so this remains a heuristic.

## Rerun versus inherited

**Recomputed this session:**

- both certificates (`ck_bound classic`, `ck_bound sharp`);
- the evaluator demo;
- the structure experiments;
- the tests in `tests/test_bl_metric.py` and `tests/test_ck_bound.py`;
- the full fast suite.

**Inherited, not rerun:** Theorems C, F and F′, Proposition B, Corollary J. The embedding into
the MGHD is cited.

## Gap to a full pre-merger, merger and post-merger spacetime

The certified slab is $5.6\times10^{-4}M$ long. A head-on merger and ringdown take tens of $M$.
The abstract-CK route cannot go much further (ceiling about $10^{-3}$).

## Next decisive calculation

A validated energy-estimate continuation (no derivative loss). Treat the certified series at
$t_1 < T_{\min}$ as an approximate solution. Bound its residual, and apply a certified stability
estimate for the linearised harmonic-gauge system.
