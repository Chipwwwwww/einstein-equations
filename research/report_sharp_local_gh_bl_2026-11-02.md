# 2026-11-02 — Sharp BL GH continuation: full-second-jet and Bernstein audit

**Parent:** PR #11 head 9ff9dfe1ef1b3ea1c6b3f4ba632c84ac8687051f. **Branch:** research/sharp-local-gh-bl-continuation. **Classification for the specified test:** C — CERTIFIED DYNAMICALLY RESOLVED BL VACUUM EVOLUTION **only on a local, extremely short protected GH slab**. The result is not a merger-scale advance, not geometric enlargement of the previous CK domain, and not novel mathematics relative to Kato/PDE literature.

## 1. Audit finding: PR #11 omitted necessary second time derivatives

The old polynomial u0+t u1+(t²/2)u2 had u2 only in the metric sector; Pi2 and Phi2 were set to zero. For the full 50 GH first-order equations (gamma1=-1, gamma2=gamma0=0, harmonic source H=0), the true second jet has **Pi2_0i=16 p_i(p Delta p+5|grad p|²)/p^13**, Pi2_00=Pi2_ij=0, and Phi2_iab=partial_i g2_ab. Pi1 and g2 include the nonharmonic Delta p terms in the artificial extension. Independent off-axis rational Christoffel variation verifies every Pi2 component, and direct spatial differentiation checks Phi2.

This mistake means the inherited approximation is **not** a complete 50-field second Taylor polynomial; it does NOT by itself invalidate its weaker O(t) residual and local vacuum certificate.

## 2. Structural repair of catastrophic inequalities

1. The old C-infinity flat cutoff used extremely inflated bounds from derivatives of exp(-1/z). A finite C9 polynomial smoothstep suffices for H7 initial GH data (g C9, Phi C8). It agrees exactly with physical BL in the same protected region and is not globally vacuum.
2. Exact rational Bernstein-basis convex-hull derivatives cancel large alternating polynomial coefficients safely. Faà di Bruno is evaluated by derivative-order-specific generating-function convolution, not a worst derivative bound raised to the eighth power.
3. Initial field Hs bounds are separated by order, instead of charging H7 to every source contraction. For s=3,4,5,6,7 their conservative exponents are (22,27,32,37,42).
4. Matrix inverse and lapse use the variable-conformal-reference H3 Neumann/binomial lemma, not a global small-Minkowski assumption.
5. High-order full 50-field GH algebraic F and all A_i are bounded with first and second Frechet derivatives on the entire H3 ball. Pointwise DS and D(SA) use a separate U3+2=24 graph to cover H3-to-C0 embedding. Geometry-dependent energy divergence uses metric g_t and dg only.
6. Full u2=DG(u0)u1 eliminates the linear-in-t residual. H3 remainder is bounded continuously by 10^2856 t².

No new PDE principle relative to Kato, Kato–Ponce or the elementary Bernstein convex-hull theorem was established. The gain is a project-level rigorous certificate using a mathematically valid *different* extension and approximation.

## 3. Exact quantitative comparison

| Bound | PR #11 | Revised |
| --- | ---: | ---: |
| Initial H7 | 10^795 | **10^42** |
| Full GH F Lipschitz H3 | 10^11287 | **10^547** |
| Full GH second Frechet H3 | not constructed | **10^724** |
| Energy K | 10^11302 | **10^560** |
| Continuum GH residual | O(t) | **<=10^2856 t²** |
| Certified T/M | 10^-29045 | **10^-2876** |
| Certified H3 full-state error | <10^-37098 | **<10^-5762** |
| Dynamical metric correction resolved? | No (t² far below error) | **Yes, at one gauge-coordinate point** |
| Extra geometric region beyond prior CK development? | No | **No** |

The new time is **26169 orders of magnitude longer** in the same mass normalization. It remains vastly smaller than 10^-6 M, much smaller than merger timescales, and does not improve the earlier geometric development whose geodesic-slicing local interval is about 8.683e-4 M (different gauge, no validated positive-time overlap).

## 4. Certified full second GH time-jet and closed theorem

Let G(u)=F(u)-A^i(u)partial_i u, u0 the exactly constructed C9 auxiliary BL 50-field initial data on T3. Define u1=G(u0), u2=DG(u0)[u1] and v2=u0+t u1+(t²/2)u2 in **all 50 components**.

Regularity and known-reference bounds:
- u0 in H7, norm <10^42; H3 norm <10^22.
- u1 in H5, with H3/H4/H5 bounds 10^364,10^419,10^474.
- u2 in H4, with H3/H4 bounds 10^906,10^1041.
- Full GH F H3 norm <=10^368; DF <=10^547; D²F <=10^724; all A_i <=10^148, DA_i <=10^318, D²A_i <=10^487.
- Pointwise geometry flux graph (with H3 embedding margin) JS<10^304, JB<10^401, divergence D<10^430; lower-order G<10^551, energy K<10^560.
- Symmetrizer 10^-4 I <=S<=10^4 I on the entire H3 radius R=10^-320 about the known v2; full Lorentzian hyperbolicity is verified by reference-centered Neumann bounds.

From Banach Taylor's theorem, r2(t)=v2_t-G(v2) has exactly the O(t²) expression written in note31 and continuum bound ||r2(t)||H3<=10^2856 t². With T=10^-2876 M, its terminal bound is <=10^-2896. The Kato-local first-exit GH energy bootstrap closes because KT<1/2, v2 remains inside the reference tube, and ||u-v2||_Linf H3 <10^-5762 <R. Existence, uniqueness and full-field error are thus established on the auxiliary torus until T under the stated exact algebra and Sobolev proofs.

The constructed auxiliary data may fail the Einstein constraints outside the physical agreement region; NO global torus vacuum claim is made. Finite speed (<4 for true GH/subsidiary fields throughout the tube) and a Euclidean causal buffer b>=0.1 ensure that the full Einstein/Hamiltonian/momentum, GH gauge/reduction and normal subsidiary constraints propagate exactly on [0,T] times Omega_target. Thus Ric(g)=0 there. This uses the **published GH homogeneous subsidiary system** (not an independent symbolic rederivation of all subsidiary equations), and identifies the interior with the actual BL development up to local diffeomorphism.

## 5. Rigorous dynamical resolution at a concrete point

At x_*=(1/10,0,0), r=sqrt((1/10)^2+(3831/10000)^2), so 39/100<r<2/5. The full GH second jet yields g2_xx=12 p_x²/p^6>1/20 from exact rational inequalities. The coordinate metric correction therefore exceeds T²/40=10^-5752/40, while the certified exact metric-component error is <6*10^-5762. The ratio exceeds 10^10/240>4e7, meeting the requested **resolved dynamical-change** inequality by several orders of magnitude.

This is coordinate-dependent. No pointwise curvature error, gauge-invariant observable or gravitational-wave waveform is certified.

## 6. Independent tests, falsification and provenance

Proof: notes/31-sharp-full-gh-bl-second-jet.md. Code:
- research/calc/bl_sharp_cutoff.py
- research/calc/bl_sharp_certificate.py
- research/calc/bl_full_second_jet.py
- tests/test_bl_sharp_local.py (includes exact Bernstein, C9-not-Cinfinity, nonharmonic auxiliary Laplacian, independently computed all Pi2 connection variations, Phi2 derivative identities, full 50-source graph comparison, strict accelerated-time/space rejection, and dynamical lower bounds)
- research/calc/data/gh_bl_sharp_resolved_certificate_2026-11-02.json.
- Updated claim ledger entries SHBL-01 through SHBL-10.

Local prototype evidence: 7 mathematically focused regression tests passed; two separately executed independent full tensor-jet scripts passed, and 1,000 rational decimal-exponent falsification trials passed. New expanded repository tests were committed; their GitHub CI result must be checked independently before claiming they all passed. Older PR #11 reported 41 local tests. All mathematical certificate guards use explicit exceptions rather than Python assertions and remain active under Python optimized mode.

Run:
    python -m research.calc.bl_sharp_cutoff
    python -m research.calc.bl_sharp_certificate
    python -O -m research.calc.bl_sharp_certificate
    python -m research.calc.bl_full_second_jet
    python -m pytest -q tests/test_bl_sharp_local.py

Literature checked: first-order GH source and constraints in Lindblom et al. 2006, original Kato 1975 quasilinear symmetric-hyperbolic Cauchy work, Kato–Ponce 1988 commutators, Hughes–Kato–Marsden 1977 and Friedrich–Rendall 2000 on constraints/domain dependence. The need for independent mathematical peer review remains.

## 7. All seventeen requested outcome answers

1. PR #11 proof? Its conservative O(t) residual-to-local-vacuum argument is not disproved by the missing jet; its asserted second Taylor completeness is false. Prior global overestimates still demand independent review.
2. Full 50 second derivatives? **Not in PR11**. Now yes, including Pi2_0i and Phi2_iab; off-shell Delta p included.
3. Why H7 and Lipschitz explosion? Flat C-infinity cutoff derivative maxima and raising one global maximum to high derivative powers, then repeating global Banach-algebra product inflation for all full-source factors.
4. Sharper energy? Yes, a **project-level** graded H3/H4/H5 source, Hessian and geometry-only pointwise flux estimate; K<10^560. It uses established PDE methods, not a claimed new theorem.
5. Global artificial-cutoff bottleneck avoided? **Strongly reduced but not eliminated**: still an auxiliary torus with cutoff, now C9 instead of C-infinity.
6. Improved extension? Yes, exact C9 polynomial with Bernstein derivative envelopes and H7<10^42.
7. Higher-order full approximation? Yes, all 50 second time jets, off-shell source checks.
8. Continuum residual? <=10^2856 t²; at terminal time <=10^-2896.
9. Nonlinear K? <10^560.
10. Certified T/M? 10^-2876 (a conservative **resolved-time** witness; not mathematically optimized).
11. Improvement? 29045-2876=26169 decimal orders of magnitude.
12. Time variation exceeds error? Yes at x_*=(0.1,0,0), ratio >10^10/240.
13. Exact vacuum/subsidiary? Yes on protected local causal region by exact physical initial constraints, homogeneous subsidiary propagation, finite speed and closed GH bootstrap; NOT outside cutoff.
14. Enlarged geometric region? No; no validated gauge overlap or positive-time extension of earlier geodesic CK region.
15. Genuinely new estimate vs literature? No original PDE principle established; certified constants and full-jet/cutoff application are new to repository.
16. Merger-scale obstruction? Despite a 26169-order gain, K~10^560 and C2~10^2856 remain enormous and prevent reaching even 10^-6 M. Curvature, high-order continuation and radiative exterior remain uncertified.
17. Next strategy? Prove **direct localized background-relative H3/H4 tame energies** in the physical BL causal lens, with explicit Sobolev extension and commutator constants independent of the artificial cutoff; evaluate physical local p and derivatives by interval spectral models; then validated multislab GH evolution with constraint propagation. This targets the remaining global-algebra constant inflation and actual merger-scale time.

**Classification: C — CERTIFIED DYNAMICALLY RESOLVED BL VACUUM EVOLUTION (extremely local/short).** Do NOT interpret this as B/A merger-scale evolution. Remains draft, pending independent audit and remote CI.

## Independent source graph falsification and repair

An independent traversal of the PRODUCTION full 50-field GH source using the exact analytic majorant nodes showed that the previous hand-transcribed graph was UNDER-COUNTING the H3 source value, Jacobian and Hessian exponent bounds. The independent values exceeded the hand-graph results by 5, 6 and 7 decimal exponents before the final 3-decade margin. We rejected that implementation and now call the actual source in bl_sharp_certificate.graph (geometry override only injects independently certified inverse/lapse/shift bounds). This changes (F,DF,D2F) to (10^368,10^547,10^724), K to <10^560, C2 to 10^2856, the resolved interval to 10^-2876M, and the error to <10^-5762. The earlier 10^-2832 figure is superseded. The reworked test graph independently checks the exact source and rejects an underbound. This incident is recorded explicitly rather than concealed.
