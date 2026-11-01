# Consolidated research progress

Updated 2026-11-01. Current branch: `research/horizon-stability-2026-11-01`.
Base: `3c13047cc665fcc1230854701d710f95dd6af241`, the existing unmerged
`claude/nonlinear-silent-flow` branch. The reviewed main head at the start
was `90f023bb96fe584b5043a9f6ffb4d08d3eef6360`. Both histories are preserved;
this branch does not merge anything into main.

**An explicit, complete vacuum binary-black-hole merger has not been found.**
This session produced a computer-assisted all-mode spectral bound for the
critical common minimal surface. It resolves a specific mathematical estimate
that the earlier roadmap listed as missing. It is not a new merger metric.

The sole current roadmap is [route_map.md](route_map.md). Earlier progress
entries are preserved verbatim in [the historical status](history_status_through_silent_flow.md).
The mathematical proof is [note 20](../notes/20-horizon-spectral-gap.md).

## Strongest result obtained in this session

For the critical common minimal sphere of equal-mass Brill–Lindquist data,
let

```math
L=-\Delta_S-\operatorname{Ric}(\nu,\nu)-|A|^2,
\qquad L\phi_0=0,\quad\phi_0>0.
```

With total ADM mass $M$, the verified conclusions are:

| Quantity | Bound | Status |
| --- | --- | --- |
| Smallest positive eigenvalue, **all angular modes** | $\lambda_{\rm gap}\ge1/(16M^2)$ | analytic proof + interval certificate |
| Inverse on $\phi_0^\perp$ | $\|L^{-1}\|_{L^2\to L^2}\le16M^2$ | consequence of the gap |
| Axisymmetric, equator-even positive spectrum | $\lambda\ge3/(16M^2)$ | same proof in the invariant subspace |
| Nonaxisymmetric Fourier mode $m\ne0$ | $\lambda\ge m^2/(2M^2)$ | centrifugal part of the ground-state energy |
| Whole spectrum, ordered with multiplicity from $k=0$ | $\lambda_k\ge\ell(\ell+1)/(32M^2)$ for $\ell^2\le k<(\ell+1)^2$ | min–max comparison |

These are bounds, not estimated eigenvalues. The zero eigenvalue is retained.
The proof also provides a residual version for a positive approximate mode,
$\lambda_k\ge(a/b)\mu_k+\inf((L\varphi)/\varphi)$; its
application to perturbed/evolved surfaces remains conditional on new bounds.
The result concerns the surface at the initial-data threshold $d_c$, not a
later event-horizon cross-section or spacetime stability.

**Why the proof works.** Set $u=\phi_0v$. The Jacobi quadratic form becomes
$\int\phi_0^2|\nabla v|^2dA$. A radial projection to the unit sphere then
reduces the full spectral problem to two pointwise weights. The validated
profile proves $A\ge24$ and $B\le768$, so the round-sphere eigenvalue 2 gives
$2A/B\ge1/16$. This works for every smooth perturbation and needs no mode
truncation. The critical surface's star-shapedness is verified in the same run.

**Novelty boundary.** Ground-state transforms and min–max/Poincaré comparison
are standard. The explicit bounds and their application to this repository's
critical BL sphere are the contribution here. No claim of publication-level
priority has been established.

## What was actually investigated

- Read the main work on local metric construction, horizon bifurcation,
  spectral gaps, and conditional event-horizon arguments.
- Read the unmerged axisymmetric, curl-hierarchy, and nonlinear silent-flow
  studies. They already explain why fixed-sector regularity is not a
  convergence proof; duplicating that test would not advance the target.
- Rechecked primary literature on extreme-mass-ratio spacetime gluing,
  zero-eigenvalue MOTS evolution, and numerically studied horizon spectra.
- Proposed and proved a ground-state comparison bound, rather than solving
  separate eigenvalue ODEs for all Fourier modes.
- Ran a preliminary floating-point profile experiment only to choose safe
  rational thresholds. Its approximate value 0.07475 was **not** used as
  proof evidence.
- Recomputed the fold root by Krawczyk validation, then checked the new
  inequalities on every interval ODE step and the closed polar collar.

## Verification record

Both final runs passed the fold inclusion, contraction, transversality,
meridian coverage, and all five geometric assertions at **256-bit Arb
precision**, using the final source code.

| Replay | Full-surface boxes including pole collar | Safe diagnostic $\inf A$ | Safe diagnostic $\sup B$ | Exact rational assertions |
| --- | --- | --- | --- | --- |
| [Normal steps](calc/data/horizon_gap.json) | 1,537 | $>25.47$ | $<716.43$ | all passed |
| [Half-sized steps](calc/data/horizon_gap_fine.json) | 3,073 | $>25.53$ | $<699.63$ | all passed |

The theorem uses the conservative **exact integers 24 and 768**, not the
rounded diagnostic values. The two Krawczyk images overlap and lie strictly
inside the original root box; each run independently verifies that same
unique fold. Source hashes in both files match the current implementation.

Commands:

```bash
python -m pytest tests/test_horizon_gap.py tests/test_fold_proof.py -m 'not slow' -q
python -m research.calc.horizon_gap --output research/calc/data/horizon_gap.json
python -m research.calc.horizon_gap --fine --output research/calc/data/horizon_gap_fine.json
```

The fast tests passed: **21 passed, 3 slow tests deselected**. They check the
radial-graph metric, both angular derivative terms in the ground-state
identity (including its approximate-mode residual version), conformal normal-displacement scaling, reparametrization
invariance, invalid-input rejection, the Krawczyk construction on an exact
fold, and the inherited integrator on known exact ODEs.

The full certificate is separate from those unit tests. It replays the
fold calculation, not merely stored log coverage. Source hashes and arithmetic
versions are embedded in each output. Smaller steps are a cross-check of
the same analytic proof and integrator, not an independent formal verification.

## Current state of the ultimate objective

| Requirement | Current evidence | Remaining gap |
| --- | --- | --- |
| Exact vacuum initial data for two centres | Brill–Lindquist formula, established literature | already available; does not determine an event horizon locally |
| Exact local spacetime | inherited convergent-series certificate, corrected in note 17 | only a narrow tube and short time, not a complete exterior |
| Critical common minimal surface | existing fold, replayed here | it is an initial-data threshold |
| All-mode invertibility away from its zero mode | **explicit bound proved this session** | explicit nonlinear/elliptic regularity constants still needed |
| Common MOTS at an explicit $d>d_c$ in the true evolution | qualitative earlier results | pointwise nonlinear residual and true evolved-data bounds |
| First common apparent horizon of the slicing | initial no-common-horizon results are inherited | exclude other common MOTSs at intervening positive times |
| Two event-horizon components merging into one | earlier conditional EMR arguments only | global causal/exterior control; no unconditional new result here |
| Explicit metric over pre-merger, merger and post-merger exterior | not obtained | global PDE construction/continuation and horizon identification |

The inherited local time bounds, $8.683\times10^{-4}M$ on the full tube and
$1.545\times10^{-3}M$ on its waist region, were read from the corrected
repository work and **not recomputed** in this session. They are not extended
by the spectral result. The long global shooting sweeps were also not rerun.

## Precise obstacle and next direction

The immediate unresolved estimate is a pointwise nonlinear bound for the
MOTS expansion on nearby graphs in the **true vacuum-evolved** data.
An $L^2$ spectral inverse does not imply the required $C^{2,\alpha}$ bound.
A quantitative elliptic estimate and a remainder bound must bridge that gap.

The most economical next local target is a certified strictly outer-trapped
surface plus an outer untrapped barrier at one explicit separation and time.
That may avoid locating the exact fold time. It still requires actual
pointwise enclosures; none are asserted here. This target is a checkpoint
on the roadmap, never a substitute for the global merger objective.

For global progress, the roadmap retains horizon-penetrating hyperbolic
continuation and the rigorous EMR gluing route. Adding more early-time
Taylor coefficients or re-testing the silent first sector is deprioritized.

## Sources actually checked

- [Hintz, arXiv:2408.06715, §6 and Remark 6.2](https://arxiv.org/html/2408.06715v1):
  true small-hole vacuum solutions and a merger-plus-stability construction;
  further event-horizon analysis is explicitly left open. This already
  establishes more than local-existence theory, but is not an explicit
  finite-ratio merger metric.
- [Andersson–Mars–Metzger–Simon, arXiv:0811.4721v2](https://arxiv.org/html/0811.4721v2):
  the established framework for zero-principal-eigenvalue MOTS evolution.
- [Pook-Kolb et al., arXiv:2006.03940, §§II.2–III](https://arxiv.org/html/2006.03940v1):
  numerical spectra in head-on mergers; relevant operator and symmetry context.
- [Emparan–Martínez, arXiv:1603.00712](https://arxiv.org/abs/1603.00712):
  abstract checked; the exact horizon construction is in an extreme-mass-ratio
  limit and does not supply the requested global finite-ratio metric.

A search result also identified Xie–Prasad, arXiv:2609.11819, on generic
common-horizon formation. The abstract was available in search, but direct
full-text retrieval failed. It was not used in the proof. The repository's
cited 2026 full-subextremal Kerr stability preprint could not be retrieved
in this session and was likewise not used to upgrade any claim.
