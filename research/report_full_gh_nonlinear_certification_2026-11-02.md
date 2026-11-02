# Research report — full nonlinear 50-field GH validation (2026-11-02)

**Pinned parent:** PR #9 head 3a7ec7502cefa09522c24e2b853a6bb2cfea0df8. New branch: research/full-gh-nonlinear-certification. Existing PR #7/#8/#9 and main preserved.

**Outcome:** narrow but real **C — CERTIFIED GENERAL 3+1 EINSTEIN GH THEOREM**, for a finite periodic Cauchy slab, undamped gamma1=-1,gamma2=gamma0=0, fixed harmonic gauge H_a=0, fully verified H3 geometric tube and exactly constraint-compatible initial data. The first actual nonflat full-50-field certificate is a manufactured spatially varying approximation of known Kasner. **NOT** a binary merger slab, nonperiodic black-hole exterior or a new theorem relative to mathematical PDE literature.

## P0 independent audit

At audit the remote heads were main 90f023b, PR7 fa5fbce, PR8 ffce85a, PR9 3a7ec75. Theorem25.2 H3 Fourier estimates were checked: H3 C0 embedding 6, gradient embedding 12, H3 matrix product 48, transport commutator 144. The high-order symmetrized energy and Kato first-exit argument are correct under their actual whole-infinite-dimensional tube hypotheses. Tests in PR #9 did NOT prove that the full 50-field nonlinear GH source satisfies those hypotheses, nor did Theorem25.2 itself propagate Einstein constraints. No incorrect Fourier inequality was found in the checked derivation; the missing source-level hypotheses are nontrivial.

## Core mathematics and implementation

**Note28** proves an undamped full-50-field GH a posteriori existence theorem on T3. It reconstructs the ten g_ab, ten Pi_ab and thirty Phi_iab components, the entire gravitational RHS of Lindblom et al. (2006) original equations (35)–(39), the inverse spacetime metric, lapse, shift, connection, and GH principal blocks. GH gauge source is H_a=0. Nonzero damping is explicitly NOT covered.

A rational Banach-algebra intermediate (c,d,L) means for every H3 input in the ball, output deviation <=d and full Fréchet derivative operator norm <=L. Product, reciprocal, lapse inverse square root and Gauss–Jordan inversion all have explicit convergence and rational failure checks. This covers the **whole infinite-dimensional H3 tube**, not sampled points. Known C0 embedding plus a matrix-inertia homotopy gives uniform Lorentzian signature, spatial positivity and GH symmetrizer positivity.

A structural refinement makes the energy growth computable without charging all unknown field derivatives: S=S(g), S A^i=B_i(g), so only g_t=-N Pi+beta Phi and spatial gradients of g contribute to their divergence. The **known full-tube bounds** U_g, X_g, J_S, J_B yield D_adapt=J_S U_g+sum_i(J_Bi) X_g. The usual H3 transport commutator gives Y'<=K Y+(Lambda/sqrt(lambda))*rho, where K=D_adapt/(2lambda)+(Lambda/lambda)(L_F+48sum_i L_Ai V_i+144sum_i M_ci). With E_cert=Lambda rho T/[lambda(1-KT)]<R and exact initial data, the first-exit/Kato continuation is closed. No unbounded unknown-solution norm remains. This is standard symmetric-hyperbolic analysis with explicitly evaluated constants; no publication-level novelty claimed.

**Note29** separately proves an H3 matrix-inverse lemma valid around a *spatially variable* reference metric G, using 48² B R<1 with B>=||G^-1||H3, and 6JR<1 for Lorentzian inertia. It gives ||g^-1||H3 <=B/(1-48²BR) and Lipschitz bound 48²[B/(1-48²BR)]². A nonconstant periodic example and failure tests use exact rational arithmetic. This standard functional calculus still needs physical-domain Sobolev bounds for black holes.

## Full 3+1 nonflat benchmark

Harmonic-time Kasner (2/3,2/3,-1/3):
g_K=diag(-exp(2t), exp(4t/3), exp(4t/3), exp(-2t/3)), with nonzero Riemann R^x_{txt}=2/9. The exact 50-field evolution is independently symbolically checked for all harmonic times. The approximate state adds epsilon*t²*cos(x1) to three distinct GH components psi_01, Pi_22 and Phi_3,22; epsilon=1/10^7. At t=0 the exact ADM Hamiltonian/momentum, GH gauge and reduction/curl constraints are satisfied; the continuum residual is nonzero.

Exact-rational rounded certificate on [0,1/1000], with R=1/10^8:

- Symmetrizer lambda>=489/500, Lambda<=511/500; H3 source bound MF<=10 and whole-tube Frechet Lipschitz LF<=525.
- H3 A_i<=21, derivatives (123,123,119), flux derivative/metric bounds with D_adapt<=66507/200.
- Growth K<750, KT<3/4, continuum residual <=19647/5000000000000 (=3.9294e-9).
- **50-field H3 metric/evolution error <=3346539/206017149996983056 <2e-11 <R**.
- Exact constraint propagation uses published GH homogeneous subsidiary equations and the contracted Bianchi identity. Their full expanded source has not been independently regenerated in code; zero initial constraints and no physical boundary are essential.

This benchmark has an independently known exact metric. It is disclosed as the manufactured comparison and its exact algebraic identity helps compute the residual, but **existence to the certified T follows from the closed energy bootstrap**, not an assumption that the comparator exists. The general GH theorem itself assumes no background symmetry.

## Verification and proof provenance

Proof notes: notes/28-full-gh-h3-source-certificate.md; notes/29-variable-background-inverse-and-obstructions.md.
Source code: research/calc/gh_full_source.py; gh_h3_majorant.py; gh_h3_certificate.py; gh_reference_neumann.py.
Tests: tests/test_gh_full_source_majorant.py (8 tests: full 50-field evolution, independent 4D Ricci and Riemann, generic nonzero spatial-derivative sign/index check, rational majorants and falsification); tests/test_gh_reference_neumann.py (2 tests and failure branches); plus 22 inherited tests from PR #8/#9.

**Actual local result: 32 passed** using Python3.13.5, pytest9.0.2, SymPy1.14. The small rational certificate performs explicit rejection checks even under Python -O. Full inherited repository suite, remote CI and outside proof review were **not** independently completed.

Reproduce from root:

    python -m pytest -q tests/test_gh_full_source_majorant.py tests/test_gh_reference_neumann.py tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py tests/test_full_hs_bootstrap.py tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py tests/test_bl_gh_initial_data.py
    python -m research.calc.gh_h3_certificate
    python -O -m research.calc.gh_h3_certificate
    python -m research.calc.gh_reference_neumann

## Primary literature and limitations

Original formula-level reading: Lindblom–Scheel–Kidder–Owen–Rinne (2006), arXiv gr-qc/0512093 v3, equations 16 and 35–52. Hughes–Kato–Marsden (1977), Arch. Rational Mech. Anal. 63, original sections 2–3, establishes short-time quasilinear well-posedness. Kato (1975), Kato–Ponce (1988) and Friedrich–Rendall (2000) support the first-order energy, commutator and constraint arguments. Other validated numerics and physical-boundary mathematical literature are not exhaustively full-text reviewed: no universal novelty or proof-priority claims. The contribution is an explicit constructive project-level certificate applying established mathematics, not a new universal PDE result.

## Twelve required answers

1. **PR9 H3 sound?** Yes as a conditional full-tube symmetric hyperbolic theorem, subject to the written norms and periodic domain; no counterexample found in the independently checked inequalities.
2. **Full 50-field H3 source/Jacobian?** Yes for the undamped H=0 subfamily and a certified near-Kasner whole H3 ball; no arbitrary gauge/damping certificate.
3. **Closed full-coupled Einstein GH theorem?** Yes on T3 for smooth exact constraint-compatible initial data satisfying the explicit full-source tube bounds, with standard Kato existence, no symmetry reduction, and strict bootstrap.
4. **Gauge/reduction propagation?** Yes analytically via the exact published homogeneous subsidiary system, initial vacuum constraints and boundaryless domain; not independently rederived symbolically in full.
5. **Special spacetime symmetry required by theorem?** No; Kasner benchmark is homogeneous, but the 50-field theorem is not.
6. **Nonflat 3+1 with nonzero error?** Yes; the all-50-variable Kasner certificate with H3 error <2e-11 is a benchmark, not a merger.
7. **Useful time?** Nonzero T=0.001; geometry-adapted K<750 rather than huge generic bound. This is not useful evidence for merger-scale continuation.
8. **Black-hole exterior?** No. Physical-domain Sobolev constants, localized energy and compatible incoming/excision fields remain absent.
9. **BL spacetime enlarged?** No. The earlier exact BL GH initialization has not been extended into a validated positive-time gauge overlap.
10. **Genuinely new result vs PDE literature?** No known novel mathematical principle; source-level majorants and compact constants are new to this repository.
11. **Exact remaining merger obstruction?** Quantitative strong-field H3 calculus on a puncture-free black-hole domain, a localized/no-incoming or constraint-preserving boundary energy theorem, actual numerical-merger continuum GH residuals, controlled BL-to-GH chart overlap and multi-slab propagation through merger and radiative exterior.
12. **Next session?** Establish a finite-domain-of-dependence GH energy theorem with verified Sobolev extension/commutator constants, feed it the spatially variable reference Neumann lemma, and test on horizon-penetrating Schwarzschild before attempting a puncture-excluding BL slab. Do not divert to MOTS-only calculations.

**Classification: C — CERTIFIED GENERAL 3+1 EINSTEIN GH THEOREM in a narrow periodic/undamped harmonic domain.** Not a certified binary-merger slab, full black-hole exterior or global merger spacetime. Draft pending independent verification.
