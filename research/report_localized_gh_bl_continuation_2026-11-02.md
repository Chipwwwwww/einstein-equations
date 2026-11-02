# Research report: first quantitative positive-time local BL vacuum development in GH coordinates

**2026-11-02. Starting commit:** PR #10 head f4e56812cba4c057fc89af775641b244e318ec1f. **New branch:** research/localized-gh-bl-continuation. The historical PR #7–#10 and main were not changed. **Classification C — CERTIFIED POSITIVE-TIME BL VACUUM DEVELOPMENT**, **strictly local**, for an astronomically short time and an explicit puncture-excluding region. **No merger, no improvement of the previous geometric coverage, no publishable-new PDE theorem.** Full scientific proof is [note 30](../notes/30-localized-bl-gh-positive-time.md).

## The actual mathematical advance

We constructed smooth auxiliary 50-field GH initial data on a flat 2pi-periodic 3-torus, agreeing exactly with genuine equal-mass, ADM-M=1 BL vacuum initial data on
Omega_phys={|x|<1.1, r_+,r_->0.2}, r_pm=|x-(0,0,+/-0.3831)|. The target is
Omega_*={|x|<=0.9, r_pm>=0.3}, whose Euclidean buffer to artificial data is at least 0.1. The midpoint is included, with p(0)=1+5000/3831>2.3 and gamma_ii(0)>28. The auxiliary extension is positive C-infinity, equals BL throughout Omega_phys, is constant near the torus seams, and generally **violates the Einstein constraints outside Omega_phys**.

Exact derivative recurrences for exp(-1/s) flat cutoffs, reciprocal radii, compositions, powers and spatial derivatives yield conservative global H7 input bounds. A **reference-centered H3 inverse metric/lapse/symmetrizer graph** bypasses PR #10's near-Minkowski smallness assumption. It applies to the whole infinite-dimensional H3 tube, using Banach product 192 valid through H5, Neumann matrix inversion and binomial reciprocal square-root majorants. The graph covers the actual complete 50-field undamped GH source, not only the principal symbol. A physically adapted inverse-metric perturbation bound improves the all-GH/subsidiary characteristic-speed enclosure to <4 rather than the graph's crude 10^4019 bound.

The global approximate GH trajectory is the explicit second time jet
v=u_E,0+t u1+(t²/2)u2 with u1=F(u_E,0)-A^i(u_E,0)∂_i u_E,0, u2^g=-p_E^-2 u1^Pi and the other u2 components 0. Thus the exact source residual vanishes at t=0, and **continuum**, global torus H3 Lipschitz/derivative majorants give a residual bound at every t, not sampled points.

On the genuine BL vacuum region, the metric component evaluator simplifies to:
g00_app=-p^-4-2t²|grad p|²/p^14;
g0i_app=0;
gij_app=p^4 delta_ij + t²(8 p_i p_j-2|grad p|² delta_ij)/p^6.
The complete 50-field approximate state includes Pi_app=t dotPi0 and Phi_app=∂i g0. The evaluator rejects inputs outside the certified rational time and spatial region.

## Quantitative certificate (exact integer/rational inequalities)

All code checks operate under python -O; no implicit float mesh test or Python assert is used as an acceptance condition.

| Quantity | Certified value |
|---|---|
| d, total ADM mass | 3831/5000, M=1 |
| Initial physical region | |x|<1.1 and r_pm>0.2 |
| Protected target | |x|<=0.9, r_pm>=0.3 |
| Protected Euclidean buffer | >=1/10 |
| Entire GH/subsidiary characteristic coordinate speed | <4 |
| Positive certified time | T=10^(-29045) M |
| Initial H7 input component envelope | 10^795 |
| H3 full-tube radius | 10^(-8050) |
| Full 50-field H3 source bound | 10^8871 |
| Full H3 Frechet source Lipschitz bound | 10^11287 |
| Energy rate bound | K<10^11302 |
| Continuum H3 GH residual | rho<=10^(-8060) |
| Nonlinear GH solution error | ||u-v||_{Linf H3}<10^(-37098) |
| Pointwise metric error | <6*10^(-37098) |
| First spatial metric derivative error | <12*10^(-37098) |

The exact global symmetrizable-hyperbolic PDE solution exists and is unique over [0,T] on the auxiliary torus by local Kato existence, full-tube a posteriori H3 energy and a strict first-exit continuation. The auxiliary solution is not necessarily Einstein vacuum outside Omega_phys.

The full GH subsidiary constraints have zero exact initial data **only in Omega_phys**: Hamiltonian and momentum vanish because Delta p=0 and K_ij=0; gauge constraint vanishes from sqrt(-g)g^-1=diag(-p^8,1,1,1); reduction/curl constraints vanish by definition; and the normal subsidiary derivative is directly verified algebraically by dotPi_00=4|grad p|²/p^12, dotPi_ij=4(|grad p|² delta_ij-4p_i p_j)/p^4, dotPi_0i=0. The original homogeneous GH subsidiary symmetric-hyperbolic system of Lindblom et al. (2006), its matched speed bound <4 and local finite propagation then imply constraints vanish exactly throughout [0,T]xOmega_*. Thus Ric(g)=0 there, with no vacuum claim on the remaining torus. Geometric Cauchy uniqueness identifies it with the local development of the asymptotically flat original BL data **up to diffeomorphism**, not through an independently validated coordinate map.

This is rigorously positive time, but cannot be mistaken for merger-scale progress. No geometric region beyond the earlier certified geodesic-slicing interval of approximately 8.683e-4 M is proved, and the two clock values are not directly comparable. H3 g alone does not prove pointwise curvature error; no curvature error is claimed.

## Proof audit and external literature

The earlier PR #10 full source and PR #9 H3 first-exit proof were independently inspected. Beyond prior Kasner checks, we tested generic **off-shell** metrics with nonzero shift and nonzero spatial off-diagonal entries against the independently computed second-order GH wave equation in all 10 symmetric components. No source sign/index contradiction was found. Initial GH normal gauge derivative was verified by independent symbolic contraction. Explicit negative tests reject an uncertified faster slab, out-of-domain metric evaluation, punctures and invalid graph inputs. The selected suite ran locally: **41 passed** (full inherited repository suite and remote CI not yet complete).

Primary mathematical sources: Kato 1975 (symmetrizable hyperbolic local existence), Hughes–Kato–Marsden 1977 (quasilinear local theory), Friedrich–Rendall 2000 (symmetric hyperbolic domains of dependence and Einstein constraints), Lindblom–Scheel–Kidder–Owen–Rinne 2006 (complete GH 50-field first order, homogeneous constraint system and matching propagation speeds). Rinne 2006 and Kreiss–Reula–Sarbach–Winicour 2009 address boundary constructions; **no physical incoming boundary condition is imposed or asserted here**, because a protected finite-speed region avoids that need. The finite-propagation theorem itself is classical; novelty in the current project is its quantitative execution on actual BL data. No exhaustive mathematical-priority claim or independent publication-level theorem is made.

## Reproduction

    python -m research.calc.bl_local_cutoff
    python -m research.calc.bl_torus_majorant
    python -O -m research.calc.bl_torus_majorant
    python -m research.calc.bl_gh_taylor_evaluator
    python -m pytest -q tests/test_bl_local_certification.py tests/test_gh_generic_shift_jet.py
    python -m pytest -q tests/test_gh_full_source_majorant.py tests/test_gh_reference_neumann.py tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py tests/test_full_hs_bootstrap.py tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py tests/test_bl_gh_initial_data.py tests/test_bl_local_certification.py tests/test_gh_generic_shift_jet.py

Proof: notes/30-localized-bl-gh-positive-time.md. Code: research/calc/bl_local_cutoff.py, bl_torus_majorant.py, bl_gh_taylor_evaluator.py, and prior full 50-field gh_full_source.py (independent geometry-majorant override preserves its true source equations). New tests: tests/test_bl_local_certification.py, tests/test_gh_generic_shift_jet.py. Ledger entries appended to research/claim_ledger.md.

## Answers to the required 15 questions

1. PR #10 independently validated? Its source/periodic H3 proof was rechecked with all 10 generic off-shell wave components and existing tests, **not** externally peer-reviewed or globally verified for arbitrary gauge/damping.
2. Quantitative finite-domain theorem? **Yes** — note30, symmetric hyperbolic flux and domain-of-dependence proof with b>CT.
3. Auxiliary nonvacuum outside physical region permissible? **Yes**, exactly because the GH and homogeneous subsidiary characteristic pasts remain in the constraint-compatible region.
4. GH and subsidiary cones controlled? **Yes**, all coordinate speeds <4 throughout certified tube; b=0.1>4T.
5. Explicit C-infinity extension with Sobolev bounds? **Yes**, p_E with flat cutoffs and exact integer/rational derivative recursion through order8, initial H7 bound.
6. Genuine positive-time BL GH approximate evolution? **Yes**, full-50-field time jets and 10-component explicit metric evaluator, not imported from geodesic gauge.
7. Continuum H3 residual certified? **Yes**, rho<=10^-8060 over whole torus/time slab by whole-neighbourhood graph, no mesh samples.
8. Nonlinear a posteriori bootstrap closed? **Yes**, E<10^-37098<R=10^-8050 and K T <1/2.
9. Exact Ric(g)=0 inside protected domain? **Yes**, exact physical initial constraints, algebraically checked normal gauge datum and finite propagation of published GH homogeneous subsidiary system.
10. Metric evaluable with explicit errors? **Yes**, exact rational/symbolic g_app and strict pointwise/first spatial derivative bounds; no pointwise curvature error yet.
11. Physical spacetime and mass-normalized time? [0,10^-29045 M] times Omega_* above.
12. Genuinely beyond previous local metric construction? **No**, and no validated positive-time coordinate overlap map.
13. New mathematical lemma beyond standard literature? **No known novel PDE lemma**; explicit constants and application are repository-level advances.
14. Dominant merger obstacle? Global cutoff H7 bound ~10^795 and full source Jacobian ~10^11287 produce absurd short T; no strong-field useful continuation or radiation exterior.
15. Next strategy? Prove localized, background-relative H3/H4 energies and analytic functional calculus on a causal lens, avoid high-order global cutoff constants, and build validated interval/spectral BL jets; only then study merger/exterior coverage.

**Classification: C — CERTIFIED POSITIVE-TIME BL VACUUM DEVELOPMENT**, but no strong-field binary-merger slab or global metric. Independent review may refine this classification; do not merge unreviewed results into main.

## Stored reproducible certificate

The exact power-of-ten witness is checked into `research/calc/data/gh_bl_local_certificate_2026-11-02.json`; replay `python -O -m research.calc.bl_torus_majorant` and compare output fields, ignoring the provenance-only metadata added to the stored JSON. The stored outcome is not a black-hole-merger certificate. GitHub Actions and independent mathematical review remain outstanding.
