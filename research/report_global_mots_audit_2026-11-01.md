# Research checkpoint: global MOTS audit and literature positioning

**Date:** 2026-11-01. **Branch:** `research/global-mots-audit-2026-11-01`.
**Base:** `claude/hopeful-euler-w9ghbg`, `368fe68d50c9e3c5ea83335c3db443ba123c1082`.
Proofs: [note22](../notes/22-global-mots-audit.md). Literature:
[matrix](literature_matrix.md), [BibTeX](references.bib).
This checkpoint supersedes stronger interpretations of inherited claims125–132.

## Mathematical outcome

No global first-common-horizon bracket or new global vacuum spacetime has
been proved. The main concrete result is a **falsification and replacement**
of the proposed positive-time F' extension. Bounding a variable truncation
residual does not justify sweeping constant expansion levels: an exact
forced-oscillator example disproves that inference. This does not disprove M,
M' or N, whose barriers are actually constructed and use comparison.

Note22 proves the exact evolved-data shooting identity (GS-ODE), including
its regular upper-pole coefficient. It proves that the spatial MOTS stability
operator on an O(2)-invariant smooth sphere has exact drift, so multiplication
conjugates it to `-Delta + Q` (GS-S). The ensuing conditional all-mode gap
(GS-G) and true neck-variation instability test (GS-N) extend the tools
available for positive-time exclusion. These are applications of established
shooting, gauge and spectral techniques, not new general principles.

GS-Q derives the quadratic birth-time coefficient, including the negative
square from motion of the fold minimizer. A contraction/implicit-function
criterion supplies a uniform quadratic remainder once its explicitly listed
bounds are certified. No BL numerical value for that coefficient or uniform
remainder is asserted. The first-order time model misses a possible Ftt
contribution from the second-order metric; extending that model's shooting
coefficients alone is inadequate.

The spectral branch's critical all-mode gap was integrated and independently
replayed: `lambda_gap >= 1/16`, inverse L² norm at most16; nonaxisymmetric
mode bound `m²/2`. This is inherited mathematics, newly reproduced here.

## Branch reconstruction

All remote heads were enumerated before work and rechecked at the checkpoint;
none of these ten changed during the session. Existing branches are preserved.

| Branch | Pinned head | Relation/relevant content |
| --- | --- | --- |
| main | 90f023bb96fe584b5043a9f6ffb4d08d3eef6360 | Consolidated local metric/structural programme; older than M/M'/N branch |
| claude/hopeful-euler-w9ghbg | 368fe68d50c9e3c5ea83335c3db443ba123c1082 | Chosen base: local common-MOTS barriers/foliation, latest session16 claims |
| research/horizon-stability-2026-11-01 | d124cb5b572467a23eb7fdfd3db4c83f7a493d24 | Independent all-mode gap; selected code/tests/data and proof integrated, renumbered note21 |
| claude/nonlinear-silent-flow | 3c13047cc665fcc1230854701d710f95dd6af241 | Independent conditional flow conjugation/curl-sector study; no global continuation. Reviewed, not wholesale merged |
| claude/black-hole-merger-solutions-rpa9hj | c6e7aa8c566d3eb9b83f2a883ddf813c2a8d7088 | Main's merged mathematical history: ansatz obstructions and curl hierarchy |
| claude/bl-bianchi-formulation | dec5d6f1bf128548c4120fa840798f83c4bdd80c | Historical Bianchi formulation work in consolidated ancestry |
| claude/bl-metric-evaluator | bb5c5c68e92f7e8e477d546bb49e3f5e1d4907e8 | Historical certified evaluator |
| claude/bl-metric-extension | 948f0707d86a33f6a15aa593e1c3b85fb0e7885d | Historical spatial/time certificate extensions |
| claude/bl-structure-discovery | 130fb557d3062bb906989a8a8ceeaf2cb08751b6 | Historical structural identities |
| claude/route-j-ck-certificate | cc9e136ee032b696bb07d02b57ceb9f160f0dc17 | Historical CK certificate |

The consistent combined state is the chosen base plus the independent spectral
result, with audit qualifications below. Note10 already proves local symmetry
of nearby MOTSs by equivariant reduction; that is not a discovery here.
The silent-flow results remain linked to their pinned branch rather than
creating contradictory copies of its reports and claim numbering.

## Audit findings and corrections

- **Vacuum constraints:** initial BL harmonic conformal factor and K=0 solve
  the constraints off punctures. The analytic ADM evolution and homogeneous
  Bianchi subsidiary system propagate them. Inherited CK existence/tails
  supply the actual local solution; finite residual tests alone do not.
  This argument does not certify a complete Gaussian slice or exterior.
- **Barriers and topology:** AM09 requires the signs with consistent normals
  on a smooth compact annulus. Smooth reflected closing profiles, pole
  regularity, positive radial-normal projection and nesting are essential.
  Separating homology yields a common component. Outermost on the whole slice
  needs additional exterior hypotheses. No axisymmetric-only existence
  theorem is being substituted for AM09.
- **N and ridges:** the touching-point comparison applies to nonsymmetric
  embedded surfaces bounding the specified interior, wholly in the shell.
  A convex ridge cannot have a smooth inner tangent with two distinct normal
  planes. Neither this nor local symmetry treats shell-straddling candidates.
- **Parameter handoff:** old M passed an exact decimal time to Moore shooting
  but its binary64 rounding to Lohner geometry. These Arb balls are disjoint
  at128bits. The difference is tiny, but no enclosure theorem bridges it.
  Pass the original parameter and outward binary endpoints instead. Pole
  brackets, radial-bin coverage, tail times, gap subtraction and domain bins
  also now use outward enclosures. Conflicting supposed enclosures fail
  loudly instead of silently discarding one. Regression test falsifies the
  old time conversion at both reported times.
- **N endpoints:** printed decimals were rounded outward beyond the tested
  binary64 window. The certificate's stored endpoints define its scope.
- **M' neighbourhood:** its fixed lower endpoint leaves a sliver above dc.
  It does not compute a full delta0-neighbourhood for H; an unquantified
  qualitative neighbourhood cannot be assumed to fill that sliver.
- **Citation:** FR00 §4.4 in the linked arXiv version is Einstein–Dirac;
  frame evolution is in §§4.2–4.3. The mathematical distinction between
  independent connection variables and eliminated coordinate variables matters.

These findings qualify proof provenance; they are not numerical evidence that
the physical M/M'/N conclusions are false. Full corrected end-to-end replays
are still required. The local `[5,6.5]e-6 M` bracket is inherited, audit-qualified,
and must not be promoted to a fully revalidated or global bracket here.

## Reproduction and what actually ran

Environment: Python3.12.14, python-flint0.9.0 (Arb), SymPy1.14.0,
NumPy2.3.5, SciPy1.17.0, pytest9.2.0. Install `requirements.txt` in an isolated
environment if needed. Run from repository root.

```sh
python -m pytest tests/test_bl_barrier.py tests/test_bl_metric.py tests/test_ck_bound.py tests/test_fold_proof.py tests/test_bl_fold.py tests/test_lohner.py tests/test_gap3.py -m 'not slow' -q
python -m pytest tests/test_mots_exact.py tests/test_horizon_gap.py tests/test_bl_barrier.py -m 'not slow' -q
python -m research.calc.horizon_gap --output /tmp/horizon-gap-replay.json
```

| Check | Actual outcome |
| --- | --- |
| First selected inherited fast suite | 73 passed,4 deselected;178s |
| Final affected fast suite | 27 passed,2 deselected;5.02s. Some tests overlap the first suite; do not add counts as distinct tests |
| Full spectral certificate | Passed;256bits,1536 regular boxes plus pole collar;220.78s. New `calc/data/horizon_gap_audit_2026-11-01.json` includes source hashes |
| Corrected near-M replay | Partial only: both endpoint sign changes and6148 expansion boxes per surface passed; inner expansion contained in [-4.162,-3.838]e-8, outer [1.998,2.002]e-6. Runtime session ended before Lohner nesting/domain completion; no final JSON produced. Partial log retained |
| One original-source N fold point | Started, interrupted with no result; not counted as verification |
| Full M at1.2e-5, M', N, F', CK certificates | Inherited; not independently recomputed end-to-end here |

The spectral calculator is byte-for-byte the integrated source, so its replay
hash agrees; its original docstring points to the source branch's note20,
now note21 here. The fine-schedule record is inherited, not independently
replayed. No historical certificate was overwritten.

For the **pending** near-M full corrected replay:

```sh
python - <<'PY'
from research.calc import bl_barrier as b
r = b.main(t='6.5e-6', h_in=-4e-8,
           z_in=.632632142085469, z_out=.6328511544633266,
           out_path='/tmp/bl-barrier-near-corrected.json')
assert r['certified']
PY
```

For M' and N, use `bl_foliate.main_interval` and `bl_foliate.main` with
nonhistorical output paths (inspect signatures before invoking). These are
expensive; they were not used to claim new proof completion. JSON decimal
diagnostics are not standalone rational proof objects: re-execution of the
assertions, source hashes, and analytic reduction are the evidence.

## Literature, novelty and next direction

The mandatory five papers were accessed and their relevant statements/proof
sections inspected; the matrix records versions, hypotheses, limits and the
additional prior art. Hintz24/26 already establish exact merger constructions
in an extreme-mass-ratio setting. Published barrier, local MOTT and shooting
results are used here, not rediscovered. None applies directly to a claimed
global first common apparent horizon of the fixed equal-mass BL slice.

This session produced **genuine but limited progress toward the exact-solution
programme, not merely refinement of an approximation**: it rejects an invalid
proof route, repairs rigorous-arithmetic interfaces, proves exact structural
identities/conditional criteria, and reproduces an all-mode bound. It does
not extend the certified spacetime domain or prove a new merger theorem.
Priority is the precise global lemma in note22 §6: exterior slab/end barriers,
outermost reduction and exhaustive true-ODE classification. For the ultimate
metric, a validated strong-field continuation to a final-state stability slice
is the decisive missing bridge. See [research_gaps.md](research_gaps.md).
