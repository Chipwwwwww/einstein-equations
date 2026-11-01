# Primary-literature matrix — 2026-11-01 audit

This is the current comparison, supplementing the historical `literature.md`.
Bibliographic keys resolve in [references.bib](references.bib). Page references
below use printed pages unless labelled PDF. “Read” means the specified sections,
not an independent verification of the entire paper. Search coverage ends on
**1 November 2026**; later November publications cannot yet be audited.
No priority claim is made for this session's analytic propositions.

## Sources and applicable hypotheses

| Key; verified version | Statement and hypotheses inspected | Reading and applicability |
| --- | --- | --- |
| Hintz24; arXiv:2408.06715v1, 13 Aug 2024, 114 pp | Thms 1.1/1.3: smooth globally hyperbolic Einstein background, inextendible timelike geodesic, noncompact Cauchy surface or local absence of Killing fields; true small-hole solutions on prescribed compact sets. Thm 6.1: slowly rotating Kerr–de Sitter, geodesic crosses a horizon in finite time, small subextremal inserted Kerr hole; future relaxation. | Read §§1, 6, including proof pp109–110 and Rem 6.2 (Kerr, Lambda=0). Formal expansion is corrected to a true solution. Small mass ratio, not equal-mass BL; no closed formula. Introduction p4 explicitly defers event-horizon construction and black-region connectivity. |
| AM09; arXiv:0708.4252v3, 23 Oct 2007 | Thms 3.1, 4.1: smooth compact 3D slice data, negative inner/positive outer outgoing expansions with normals toward the exterior end; stable embedded separating MOTS. Thm 7.3 treats the trapped-region boundary. | Read §§2–4 and theorem/proof construction in §7, especially pp11,19,38. No vacuum or symmetry assumption for barrier existence. M/M' must establish a smooth annulus and all boundary signs. It does not identify a shell MOTS as the entire slice's outermost horizon. |
| AMS08; arXiv:0704.2889v2, 22 May 2007 | Thm 9.2: smooth foliation and smooth immersed strictly stable MOTS give a local adapted smooth strictly stable MOTT. §8 concerns inherited symmetry. | Read §§7–9, comparison argument pp19–20 and implicit-function proof pp21–23. Strict stability fails at the critical fold. Neither local uniqueness nor symmetry of a nearby family excludes remote MOTSs. |
| AMMS09; arXiv:0811.4721v2, 6 Apr 2009 | Prop 5.1: connected marginally stable MOTS, NEC and nonzero genericity term W give a locally unique spacelike MOTT tangent to the slice. Thm 5.2 adds outermost, past isolation and (for two future branches) present isolation. | Read §§5–6, C^(2,alpha) Lyapunov–Schmidt proof pp10–13. Qualitative local bifurcation is prior art. Uniform BL parameter bounds and global first-birth exclusion are additional tasks. Do not use past isolation as a premise to prove first birth. |
| FR00; arXiv:gr-qc/0002074v1, 22 Feb 2000 | Hyperbolic reductions, subsidiary constraint systems and local Cauchy theory, including independent frame/connection variables. | Read relevant §§4.2–4.3 equations (4.48), (4.80)–(4.94), PDF pp52,64–67. This is a primary-author survey of formulations, not a large-data merger theorem. In this version §4.4 is Einstein–Dirac; the old ledger's section citation is wrong. Coordinate E,H variables with connection eliminated need not retain frame hyperbolicity. |
| BHP21; arXiv:2104.11343v1, 22 Apr 2021 | General axisymmetric MOTS shooting in meridian geometry, equations (39)–(41); deviation/stability relation (81)–(83). | Full PDF downloaded; §§II–IV inspected, pp3–9. Spatial stability operator and vacuum constraints are explicit. Note22's equation is a convention/evaluator specialization of this known method; the evolved neck argument is not a newly discovered relation. |
| PBH21; arXiv:2104.11344v1, 22 Apr 2021 | Numerical head-on BL evolution, multiple MOTS branches, bifurcations and annihilations. | Read §§II–IV and conclusion; §III gives simulation setup. ArXiv title “Horizons weaving back and forth in time”; journal subtitle “The vanishing of apparent horizons”. These are numerical results, not interval certificates or event-horizon theorems. |
| PK20; arXiv:2006.03940v1, 6 Jun 2020 | Numerical spectra and fluxes in nonspinning merger horizons; self-adjoint simplifications depend on normal connection/deformation choice. | Read §II.B and §III. Note22 uses spatial deformation and an exact-gradient similarity, not an unexamined identification with a null-deformation operator. |
| HuiLin25; published CQG 42, 055008, 11 Feb 2025 | BL common-horizon threshold numerically d=1.5323948 for individual masses one. | Full publisher PDF via INSPIRE; §3.2.1, p13, inspected. Total-mass normalization gives 0.7661974, consistent with repository C. Threshold discovery is not new; rigorous enclosures would be a different contribution. No proof of global nonexistence follows from finder failure. |
| ACR21; arXiv:2107.02456v1, 6 Jul 2021 | Thm 3.1: strongly asymptotically flat ingoing null data (Def 2.12), charges (E,0,L,G), large radius; perturbed inner sphere glued to Kerr sphere along a null hypersurface. | Read §§2.7,3, proof structure §4 and completion §4.4. C² gluing in Thm3.1; higher-order bifurcate gluing Thm3.2. Ten charges must match. Cor3.3 is spacelike gluing. This changes transition-region data, not a continuation theorem for the fixed BL development. |
| Hintz26; arXiv:2606.28253v2, 3 Aug 2026 | Thm1.1: smooth constrained data equal to a subextremal Kerr slice off a compact set, sufficiently small high Sobolev perturbation, imply global exterior relaxation. Thm13.1 permits specified partial polyhomogeneous weighted data. | Read Thm1.1, consequence(3), Thm13.1 and proof opening, PDF pp6–7,322–325. Generalized wave-map gauge, constraint damping, Nash–Moser iteration. This verifies a 2026 extension relevant to Hintz24. No quantitative entry into this small-data regime has been proved for equal-mass BL. |

## Repository theorem comparison

Repository assertions here refer to pinned baseline `368fe68`; audit qualifications
in [note22](../notes/22-global-mots-audit.md) and the report take precedence.
“Specialization” does not assert that every implementation dependency is reverified.

| Repository result | Literature comparison | Established versus remaining gap | Decisive next method |
| --- | --- | --- | --- |
| C (critical fold; note09) | Known BL threshold and saddle-node picture (HuiLin25, BHP21); local bifurcation framework AMS08/AMMS09 | Inherited validated fold and leading onset coefficient; no claim of discovering a new phenomenon. | Certify derivatives of the TRUE reduced map and quantitative implicit-function bounds. |
| F' (note12) | Initial minimal-surface/stability exclusion, unlike time-dependent MOTT existence | Inherited explicit initial-slice near-critical exclusion, with outermost symmetry and neck analysis; not rerun in full here. Positive-time analogue absent. | Exact evolved-data shooting, exterior/end barriers, true Jacobi variation; no constant-residual shortcut. |
| H (note13) | Specializes local MOTT bifurcation theory | Qualitative near-critical onset; uncomputed neighbourhood. Cannot turn M' into coverage of an uncovered sliver. | Uniform inequalities (8)–(9) in note22. |
| M (note20) | Concrete application of AM09 barrier existence | Local common stable MOTS claimed at 1.2e-5 and 6.5e-6. Arithmetic handoff repaired; near-time replay only partial in this session. | Complete corrected nesting/domain replay before upgrading independent-verification status. |
| M' (note20) | Same barrier theorem uniformly in a parameter | Inherited interval [0.76619742824011,0.7662], not all (dc,dc+delta). Complete corrected replay pending. | Preserve decimal/set parameters; certify interval families and cover missing sliver separately. |
| N (note20) | Maximum-principle exclusion, compatible with AMS08 §7 | Shell-only, applies to nonsymmetric embedded candidates wholly in the oriented shell. Rounded display endpoints were wider than tested endpoints. Full corrected replay pending. | Exact endpoints, then classify every outermost common candidate including shell-straddling ones. |
| Local metric (notes15–17) | Constructive local realization of vacuum Cauchy theory (FR00) | Inherited convergent exact series with controlled tails on certified local domains. Neither closed form nor full exterior/EH solution. | Validated finite-slab continuation with independent frame/connection and gauge control. |
| Critical all-mode gap (note21) | Standard ground-state/min–max machinery; spectra studied numerically in PK20 | Independent branch integrated; 1/16 gap and inverse bound16 replayed. This is quantitative certification, not new spectral theory. | Evolved weights/residual in GS-G; quantitative elliptic regularity still needed. |
| This session GS-ODE/S/G/N/Q | Known shooting, gauge transformation, min–max, Sturm comparison and implicit differentiation | Full elementary derivations in note22; conditional evolved bounds, no new global BL existence/exclusion. | Test true ODE axis collar and exterior coverage before an expensive sweep. |

## Discovery trail and verification limits

Backward tracking: AMMS09 cites AMS08; Hintz24 §6 cites final-state stability;
BHP21/PBH21 reference PK20 and AMMS09. Forward/concept searches located
HuiLin25 and Hintz26; INSPIRE's BHP21 citation record was also inspected.
Searches included “MOTS stability merger 2026”, “characteristic gluing Kerr”,
“nonlinear stability full subextremal Kerr”, and BL critical separation.
The 2026 publisher record for Sherif's *Conformal symmetries and MOTS stability*
(DOI 10.1016/j.physletb.2026.140715) was found; full theorem hypotheses were
**not** inspected and it is not used. BL has no established suitable conformal
Killing field here. No MathSciNet/zbMATH full-text access was obtained.
This is a focused prior-art audit, not a proof of bibliographic exhaustiveness.

The five mandatory papers were accessed in full; only the relevant portions
listed above were read. BHP21's web PDF endpoint failed initially, but direct
PDF download succeeded. Abstracts alone are not treated as theorem evidence.
