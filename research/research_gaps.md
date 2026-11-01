# Research gaps and decisions — 2026-11-01

The objective remains an explicitly evaluable global vacuum merger spacetime,
including its exterior and causal horizons. Exact solutions need not have
closed formulas: see the [literature matrix](literature_matrix.md), Hintz24/26.
This programme's equal-mass BL data are a different, strong-field problem.

## Executed hypotheses

| Hypothesis | Decisive test and result | Decision |
| --- | --- | --- |
| A constant first-order expansion-level sweep covers all exact MOTSs when the remainder is bounded | Exact forced-oscillator counterexample in note22 §2: all constant forcings miss a root attained by a bounded variable forcing | Rejected. Existing constructed barriers remain legitimate; enumeration of exact surfaces does not. |
| Nonzero K necessarily prevents self-adjoint stability tools | On a smooth O(2)-invariant sphere, the smooth connection form is exact; conjugation cancels drift. Rotational-only counterexample has nonclosed form | Proved with reflection and topology hypotheses, note22 GS-S; standard transformation specialized here. |
| The initial all-mode gap can be used without recertifying infinitely many modes | Weighted energy/mass inequalities plus residual give GS-G by min–max | Conditional proof complete. Evolved geometric inputs absent; do not label the evolved gap certified. |
| Initial thin-neck logic survives on true evolved data | Differentiate the exact MOTS ODE; similarity and strict Dirichlet domain monotonicity imply instability if the variation has two consecutive zeros | GS-N proved conditionally. First-order level variations are not substitutes for true Jacobi fields. |
| More time coefficients alone resolve separation scaling | Implicitly differentiate F=Fz=0; the quadratic term includes movement of the fold minimizer and Ftt of the true metric | GS-Q derived and exact-symbolically tested. Uniform remainder needs F in C³, explicit inverse/contraction bounds and derivative bounds. |

## Next research gates, in order

1. **Repair and finish certificate provenance.** Complete near-M geometry,
   M', and N with preserved parameter intervals. Historical JSON is retained;
   partial output is not a certificate. Audit every scalar-to-interval bridge.
2. **Global finite-time exclusion lemma.** Construct a common geodesic slab
   on a compact punctured exterior and barriers at its ends. Prove that any
   common MOTS yields an enclosing connected stable outermost component;
   establish topology and symmetry, then cover ALL true shooting starts.
   Shell-straddling curves and the neck limit are mandatory. Note22 §6 spells
   out a sufficient implication without presuming the desired conclusion.
3. **True-ODE axis collar.** Implement a regular integral equation using
   GS-ODE's pole coefficient and the certified metric evaluator, with local
   Lipschitz/Taylor bounds. First falsification test: recover conformal initial
   data and reproduce exact expansion zero over a nontrivial collar. Do not
   launch a global sweep until its full domain and collar are certified.
4. **Uniform fold theorem.** Certify the inverse/derivative inequalities in
   note22 (8)–(9), giving an explicit O((d-dc)²) remainder. Determine Ftt from
   the true evolved metric, not the linear-time level equation. This is local
   onset until gate2 is closed.
5. **Metric continuation.** Use a frame/connection symmetric-hyperbolic
   system or a generalized-wave gauge with a posteriori Sobolev energy
   bounds, overlap uniqueness and constraint propagation. The silent-flow
   branch does not remove the principal-symbol obstruction. Control the
   gauge before Gaussian normal coordinates develop caustics.
6. **Future exterior and event horizon.** Certify entry into a final-state
   stability theorem's actual weighted norms, then construct the relevant
   complete null infinity and backward horizon. Local MOTS data alone do not
   supply this. No such transition is proved for these BL data.

## Gluing: useful alternative, precise limitation

ACR21 provides a way to alter far-out data to Kerr while solving the null
constraints and matching charges. With finite propagation and uniqueness,
a modified data set agreeing on a sufficiently large initial region has the
same development on that region's domain of dependence. This is a valid
route to a **different** global spacetime retaining a local merger region;
it cannot silently be called the original BL maximal development. It also
cannot make a strong-field interior small in the norms of a Kerr stability
theorem merely by changing the far exterior. The missing lemma is a finite
strong-field continuation reaching a controlled stability slice.

For a second, extreme-mass-ratio track, making Hintz's gluing/Nash–Moser
construction quantitative and evaluable is more directly grounded in known
global existence than searching for a closed equal-mass formula. It requires
explicit tame constants, inverse construction and error control, not just
writing an iteration symbolically. It is not completed here and does not
replace the BL track.

## Novelty boundary

New to this repository checkpoint: the counterexample to the proposed
constant-residual extension; parameter-handoff repairs; explicit evolved-data
identities and conditional bounds; branch integration and independent gap
replay. No new general relativity principle, global merger theorem, global
first-common-horizon bracket, or uniformly certified birth-time law is claimed.
