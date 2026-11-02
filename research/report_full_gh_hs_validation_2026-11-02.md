# Full-GH H3 mathematical research checkpoint — 2026-11-02

**Source branch:** PR #8, ffce85ad9943c1da5ebae1a834c5328907623a46.
**New branch:** research/full-gh-hs-validation. PR #7 global-MOTS audit is preserved. Neither main nor PR #8 is modified.

## Principal mathematical result

**VERIFIED GENERAL HYPERBOLIC PDE THEOREM, NOT FULL EINSTEIN.**
Note 25 establishes a quantitative a posteriori existence and uniqueness theorem for arbitrary N-component smooth symmetrizable quasilinear first-order systems on the periodic three-torus. It builds a full H^3 tube around a smooth approximate solution v; assumes a certified continuum H^3 residual and coefficient/symmetrizer/source bounds valid on the ENTIRE tube; and closes the nonlinear H^3 error bootstrap using Kato local existence plus a first-exit contradiction. No true-solution norm enters unbounded.

It uses the normalized Fourier norm with weight (1+|k|^2)^3 and proved elementary constants: H^3->L∞ 6, H^3->gradient L∞ 12, H^3 multiplication 48 and transport commutator 144. The explicit exponential bootstrap is in note25 eq. (25.5). This is a conservative constructive *reformulation of standard PDE energy methods*, not a claimed new PDE principle or original literature theorem.

An exact-rational 3D, two-component nonlinear coupled system with symmetric noncommuting principal matrices yields continuum residual <=393/1250000 and closed H^3 error <=393/124623695 for T=1/100, radius R=1/1000, initial error zero. A time-slab recurrence with finite constants is implemented. Negative tests reject wrong residuals, unmatched initial data, too-wide hyperbolicity tubes and failed bootstrap conditions. The exact independent comparator is used only for regression checking.

## Einstein-specific outcomes

**Note 26:** On all 50 GH first-order fields, the gamma1=-1 first-order principal symbol is symbolically symmetric with the GH symmetrizer. Tests cover gamma2=0 and gamma2=1/4 with the required psi-Pi cross term and Lambda2>gamma2^2. At Schwarzschild Kerr-Schild r=2M, the 50-field symmetrizer has minimum eigenvalue 1/2, N=1/sqrt2, beta_x=1/2. The exact fixed-coordinate GH gauge source H_a=-Gamma_a of Schwarzschild is checked independently. On 3M/2<=r<=5M/2, an **assumed pointwise** metric-error tube <=0.06 certifiably preserves Lorentzian signature and has symmetrizer eigen lower bound 150/359, by exact rational Neumann and Weyl inequalities. This tube does NOT itself follow from a boundary Sobolev theorem on the shell. No full 50-field nonlinear source or nonzero GH evolution residual is verified.

**Note 27:** Time-symmetric equal-mass Brill-Lindquist conformal data with p=1+1/(4r_+)+1/(4r_-), K=0, lapse N=p^-2, zero initial shift/time derivatives, admit explicit 50 GH initial fields g=diag(-p^-4,p^4,p^4,p^4), Pi=0, Phi_i=partial_i g and exact initial harmonic gauge source H=0. Since sqrt(-g)g^-1=diag(-p^8,1,1,1), Gamma_a=0 initially. The Hamiltonian/momentum and first-order reduction constraints hold off punctures. The normal derivative GH subsidiary data, positive-time gauge transition from geodesic slicing and puncture-excluding domains are NOT yet independently certified. This is a new-to-repository exact initialization, not a new relativity theorem.

## Primary mathematical literature and exact hypothesis scope

- Kato (1975), DOI 10.1007/BF00280740: local well-posedness for smooth quasilinear symmetric-hyperbolic Cauchy systems in H^s, s>n/2+1, positive hyperbolicity symmetrizer; original article bibliographic scope checked, not exhaustively full-text proof-read. No ready-made validated numerical Einstein residual constants.
- Hughes, Kato, Marsden (1977), DOI 10.1007/BF00251584: second-order quasilinear hyperbolic solutions in H^{s+1}xH^s, with Einstein case s>n/2 under the specified leading-coefficient hypotheses. Selected original proof and theorem sections inspected. No generic merger validation.
- Kato–Ponce (1988), DOI 10.1002/cpa.3160410704: standard commutator framework. Our integer-H^3 constants are conservative elementary Fourier estimates, not a new harmonic analysis inequality.
- Friedrich–Rendall (2000), arXiv gr-qc/0002074: hyperbolic reductions, geometric constraints and subsidiary propagation; selected original sections consulted. Reduced Einstein equations are not vacuum without exact constraints.
- Lindblom, Scheel, Kidder, Owen, Rinne (2006), arXiv gr-qc/0512093: original eqs (27)–(39), principal 50-field GH and constraint modifications, gamma3=gamma1 gamma2; checked actual principal expressions. Complete GH source and subsidiary equations not yet implemented here.
- Rinne (2006), arXiv gr-qc/0606053, and Kreiss–Reula–Sarbach–Winicour (2009), arXiv 0807.3207: characteristic/constraint-preserving and maximally dissipative boundary theory under exact input hypotheses; the present work does not certify an Einstein moving-excision or outer-radiation boundary.
- Makridakis–Nochetto (2006), Numer. Math. DOI 10.1007/s00211-006-0013-6, and *Numerical existence and uniqueness proof for solutions of nonlinear hyperbolic equations* (JCAM 2001), DOI 10.1016/S0377-0427(00)00563-X: relevant numerical reconstruction/fixed-point theory for specific evolution classes. Only abstracts/theorem scope checked, no generic Einstein transfer established.
- Szeftel (2026), arXiv 2609.40068, and Hintz (2024/26), arXiv 2408.06715: highly specialized nonlinear stability/gluing regimes; their strong-field merger scope cannot be assumed to cover comparable-mass binary simulations. Full long-paper proofs not read in this session.

The task's broad exhaustive mathematical literature demand is **not fully completed**; the above specifies what was actually inspected.

## Reproducible evidence

Code: research/calc/full_hs_bootstrap.py; gh_50_principal.py; gh_ks_source.py; gh_ks_geometry_tube.py; bl_gh_initial_data.py. Proofs: notes/25-full-hs-a-posteriori.md, notes/26-gh-50-principal-and-gap.md, notes/27-bl-gh-initial-data.md.

Run from repo root:

    python -m pytest -q tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py tests/test_full_hs_bootstrap.py tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py tests/test_bl_gh_initial_data.py
    python -m research.calc.full_hs_bootstrap
    python -m research.calc.gh_ks_geometry_tube
    python -m research.calc.bl_gh_initial_data

Actual local combined outcome: **22 passed in 3.25 s** (Python 3.13.5, pytest 9.0.2, SymPy 1.14). Inherited PR #8 subset independently replayed: **9 passed**. Full repository regression and remote CI remain unverified. Tests corroborate algebra, NOT mathematical novelty or a full Einstein vacuum theorem.

## Exact obstacle map

P0: Complete: reconstruct branch provenance and test restricted certificates.
P1: Complete at **general symmetric-hyperbolic** level in H^3 with a closed explicit theorem and test.
P2: Partial: 50-field principal GH symmetrizer, gauge-source identity on Schwarzschild, and exact initial BL constraints; **full nonlinear source**, gauge/reduction subsidiary constraint propagation, and boundary compatibility absent.
P3: Partial geometric only: a nonflat Schwarzschild pointwise hyperbolicity tube; **no full nonzero-error GH 3+1 slab** certified.
P4: Partial exact BL GH initial data off punctures; **no certified overlap or time extension**.
P5/P6: Not achieved; no merger-wide evolution or radiative exterior.

The next mathematically decisive gap is to compute interval-certified H^3 composition/Jacobian bounds for the **complete 50-field GH nonlinear source**, including inverse metric and lapse functional calculus on a whole Sobolev ball, and prove a domain-specific high-order incoming constraint/boundary estimate (or equivalent finite-speed interior theorem). Only then can the note25 bootstrap be applied to true GH vacuum evolution from exact constraints. Long-time nonlinear constants and exterior matching remain open.

## Required final status answers

1. Full-coupled Einstein a posteriori theorem: **NO**.
2. Existence/uniqueness: **YES for general symmetric hyperbolic PDE under explicitly checkable tube hypotheses**; NOT yet a Ricci-flat Einstein proof.
3. Exact spaces/constants: H^3(T^3;R^N), 6, 12, 48, 144, symmetrizer eigen bounds lambda/Lambda and tube constants in note25.
4. Derivative loss: controlled for first-order symmetric-hyperbolic system; actual GH source/physical boundaries not finished.
5. GH reduction/gauge constraints: exact BL initial values checked off punctures; full subsidiary propagation NOT verified.
6. Spherical/cylindrical symmetry: mathematical theorem and 50-field principal algebra require neither.
7. Nonflat 3+1 Einstein benchmark: geometric/principal/gauge checks only, no nonzero-error full evolution validation.
8. BL certified time: **no extension**.
9. Novel mathematics relative to PDE literature: **none established**. Constructive constants and exact certificates are project checkpoints.
10. Remaining merger obstacle: complete GH nonlinear residual/Jacobian and H^3 estimates; chart/boundary/constraint proof, nonlinear time continuation, asymptotically flat exterior.

**Classification: VERIFIED GENERAL HYPERBOLIC PDE THEOREM** (established-method theorem with evaluated constants; not a verified finite-slab full Einstein theorem).


## Addendum: puncture-excluding BL initial hyperbolicity certificate

In the region with both puncture radii $r_\pm\ge1/4$, the exact BL conformal factor satisfies $p\le3$. Thus lapse $N=p^{-2}\ge1/9$ and $\lambda_{\min}(S_{50})\ge1/81$ initially. Under the separately justified uniform *pointwise* metric perturbation bound $\|\delta g\|_{\rm op}\le6/1000$, the Neumann multiplier is $\le243/500$ and $\lambda_{\min}(S_{50})\ge500/40503$ for the perturbed spatial metric. All bounds are exact rational. This does not create a boundary Sobolev inequality or a positive-time GH residual certificate. The explicit no-go test rejects an overly wide radius. Re-run the combined test command above to obtain **22 passed**.
