# 28 — Full 50-field undamped GH Einstein validation on a periodic Kasner slab

**2026-11-02.** Source branch: PR #9 `3a7ec7502cefa09522c24e2b853a6bb2cfea0df8`. **Status: proved conditional general 3+1 GH a posteriori existence theorem and exact-rational closed full-source benchmark on a nonflat Kasner torus.** No publication-level novelty claimed; this is the constructive application of Kato's well-posedness and standard Fourier/Moser algebra. No binary black-hole time extension, puncture/excision/exterior domain, or merger is certified. An independent proof audit is required before merging.

## 1. Independent audit of Theorem 25.2

For Fourier norm `||f||_3^2 = Σ_{k∈Z³}(1+|k|²)^3 |fhat_k|²`, the cube `||k||∞=n` contains `24n²+2` points. Thus (separating the zero mode)

- `C0² ≤ 1+24 ζ(4)+2 ζ(6) <36`, so `||f||∞ ≤6||f||3`.
- `C1² ≤72 ζ(2)+6 ζ(4)<134<144`, so `||∇f||∞≤12||f||3`.
- `w3(p+q) ≤4(w3(p)+w3(q))`, so by Fourier convolution `||fg||3≤48||f||3||g||3`.
- `|w3(p+q)-w3(p)| ≤6|q|(w2(p)+w2(q))`; the corresponding two convolution terms each cost at most `C1||A||3||e||3`, yielding `||[Λ³,A]∂_ie||2≤144||A-C||3||e||3` for any spatially constant matrix `C(t)`.

All estimates apply componentwise and by Frobenius submultiplicativity to arbitrary finite-dimensional matrix/vector fields. Applying them to `e_t+A^i(u)∂_ie = -r+F(u)-F(v)-(A^i(u)-A^i(v))∂_iv`, the energy `Y²=∫(Λ³e)^T S(u)(Λ³e)` has no derivative loss. In particular, **Theorem 25.2 is mathematically sound under the explicitly required whole-infinite-dimensional-tube hypotheses**. The key caveats are: (i) its 3-torus is not a black-hole exterior, (ii) H³ pointwise coefficient bounds must cover every tube element, (iii) the true symmetrizer/time and spatial flux derivatives must be enclosed by the exact PDE on the bootstrap interval, (iv) it is NOT a GH vacuum theorem by itself.

## 2. Exact complete GH convention and source

Metric `(-,+,+,+)`. Fix **the legitimate parameter subfamily** `gamma1=-1, gamma2=gamma0=0, H_a=0`. The 50 fields are `(ψ_ab, Π_ab, Φ_iab)` with each symmetric 4×4 pair, `Π_ab=-n^c∂_cψ_ab`, `Φ_iab=∂_iψ_ab`. Write the inverse spacetime metric `ψ^{ab}`, lapse `N=(-ψ^{00})^-1/2`, normal `n^a=-N ψ^{0a}`, shift `β^i=-ψ^{0i}/ψ^{00}`, inverse spatial metric `h^{ij}=ψ^{ij}+n^in^j`, and `d_0ψ_ab=-NΠ_ab+β^iΦ_iab`, `d_iψ_ab=Φ_iab`. Let `Γ_abc=(d_bψ_ac+d_cψ_ab-d_aψ_bc)/2`.

The complete first-order equations (Lindblom et al. 2006, published/arXiv v3 eqs. 35–39) are

```
∂t ψ_ab = -N Π_ab + β^i Φ_iab,
∂t Π_ab - β^k∂kΠ_ab + N h^{ki}∂kΦ_iab = FΠ_ab,
∂t Φ_iab - β^k∂kΦ_iab + N∂iΠ_ab = FΦ_iab,
```
where

```
FΠ_ab = 2N ψ^{cd}[h^{ij}Φ_icaΦ_jdb - Π_caΠ_db
                    - ψ^{ef} Γ_ace Γ_bdf]
       - (N/2)n^cn^dΠ_cdΠ_ab
       - N n^cΠ_ci h^{ij}Φ_jab,
FΦ_iab = (N/2)n^cn^d Φ_icdΠ_ab
       + N h^{jk} n^c Φ_ijc Φ_kab.
```

All four-index contractions are *actual* four-dimensional contractions; the nonlinear source is not a symmetry-reduced scalar equation. The source code also computes the optional `-2N ∇_(aH_b)` for a fixed prescribed covariant GH source, but this benchmark takes `H=0`. Nonzero damping is **not** certified and is explicitly rejected. The symmetrizer for this subfamily, in 10 independent symmetric-index blocks, is `S=diag(1,1,h^{ij})`. It has positive pointwise lower eigenvalue when `h^{ij}>0`. Original GH paper eqs (27)–(39), (40)–(52) supply the principal and homogeneous subsidiary systems.

## 3. Rigorous full-H³ functional calculus on a nonconstant reference

A scalar intermediate expression is represented by a triple `(c,d,L)` of rational bounds, meaning for **every** input state `u` in the full H³ ball around `v(t)`, `||f(u)-c||_3≤d`, and its Frechet derivative satisfies `||Df(u)[z]||_3≤L||z||_{H³(R^50)}`. The center `c` is spatially constant but the actual approximate field and all tube members may vary arbitrarily in space. Every state in the ball is covered, not a sampled set. For independent 50 input components, set `(c_j,d_j+R,1)` where `d_j` encloses `v_j-c_j` **in the continuum H³ norm**.

For product, with H³ multiplication constant `P=48`, the exact-recursive bounds are

```
(c,d,L)*(c',d',L') =
 (cc', |c|d'+|c'|d+Pdd',
       (|c|+Pd)L'+(|c'|+Pd')L).
```

For inverse `1/f`, require `c≠0` and `Pd<|c|`. The Neumann series `c^-1 Σ(-δ/c)^n`, with `||δ^n||_3≤P^(n-1)d^n`, yields

```
 c_inverse = 1/c,
 d_inverse = d/[|c|²(1-Pd/|c|)],
 L_inverse = (|1/c| + P*d_inverse)^2 L.
```

For `f=1+δ` and `Pd<1`, the binomial series (all absolute coefficients `|(-1/2 choose n)|≤1`) yields an inverse-square-root triple

```
(1, d/(1-Pd), L/(1-Pd)^2).
```

Matrix inversion uses four verified scalar Gauss–Jordan pivots. If a pivot's 48d/|c| reaches 1 the calculation **fails closed**, not with a numerical guess. The full source graph is evaluated with these exact rational rules for **all 50** output components and all three 50-by-50 principal matrices. The validated Fréchet derivative majorant is a whole-tube `H³` Lipschitz constant, by integrating derivatives along segments within the convex ball.

For `J_S` and `J_{SA^i}`, the *same* graph runs in the pointwise algebra with multiplication constant `1` and input radius `6(d_j+R)`, but derivative seed `1` (derivative in finite-dimensional pointwise input). Summing the `5×5` block derivative envelopes bounds the operator norm of the full 50-field derivative, because there are ten identical blocks. The H³ `A^i` matrix bounds use an l¹ sum of norms for each 5×5 block, multiplied by 10, which conservatively dominates its full Frobenius norm. `S` and `SA^i` depend only on the **metric** portion of the 50 fields, not on `Π,Φ`.

**Lorentzian tube:** Relative to constant Minkowski `η`, each `g_ab` component deviation is bounded in H³; `||g-η||op,∞ ≤12 Σ_{a≤b}(d_g_ab+R)` by `C0=6`, symmetry and Frobenius/l¹ norm. If that bound `<1`, the entire line segment `η+θ(g-η)` is invertible and has fixed Lorentzian inertia. The spatial metric is positive, and `g^{00}<0`. Repeated scalar pivots and the lapse binomial expansion are thus both *defined and separately interval-verified*. This is a near-Minkowski / short-Kasner theorem, not a general strong-field matrix inverse theorem with large tube.

## 4. Geometry-adapted high-order energy lemma

The generic note25 constant uses `D_generic=J_S[6M_F+72(V+R)ΣM_A]+12(V+R)ΣJ_B`. This estimates unknown `u_t` in all 50 fields even though the symmetrizer and flux depend only on the ten metric fields. For the **actual GH structure**, define the full-tube known bounds

```
U_g = 6 Σ_{a≤b} ( |center(Fψ_ab)| + deviation(Fψ_ab) ),
X_g = 12 Σ_{a≤b} ( d_g_ab + R ),
D_adapt = J_S*U_g + (Σ_i J_{B_i})*X_g.
```

**Lemma (proved).** On the exact solution's bootstrap interval, `||∂tS(u)+Σ_i∂i(SA^i)(u)||op,∞≤D_adapt`. *Proof:* `S=S(g)` and `SA^i=B_i(g)` algebraically, so `∂tS=DS(g)g_t`, `∂iB_i=DB_i(g)∂ig`; the `g_t` source equation is exactly `Fψ=-NΠ+βΦ` for gamma1=-1, while `||∂ig||∞≤12Σ||g_ab-c_ab||3`. Every coefficient in `U_g`, `X_g`, `J_S`, `J_B` was already certified over the infinite-dimensional ball. QED. No unknown solution norm is added to the exponent. This **sharper, geometry-adapted specialization** is the key difference from note25, not a novel PDE method absent the literature.

The rest of note25's symmetric energy proof now gives

```
Y' ≤ K Y + (Λ/sqrt(λ))*ρ,
K = D_adapt/(2λ) + (Λ/λ)*G,
G = L_F + 48 Σ_i L_Ai V_i + 144 Σ_i M_ci,
```

where `M_ci` is the full-tube H³ deviation of `A^i` from its spatially constant reference `C_i(t)`. The conservative rational comparison `exp(KT)≤1/(1-KT)` (require `KT<1`) gives, for exactly matched initial data,

```
||u-v||_{L∞H³} ≤ Λ ρ T/[λ(1-KT)] < R.
```

First-exit contradiction and Kato continuation prove true existence and uniqueness through `T` (not merely an existing-solution error estimate). All tube data, including λ and Λ, are **computable from v and R**.

## 5. Nonflat 3+1, full-50-field benchmark and exact residual

On `T³` with time `t∈[0,T]`, choose harmonic Kasner

```
g_K = diag(-exp(2t), exp(4t/3), exp(4t/3), exp(-2t/3)),
N=exp(t), beta=0,
Π_K = diag(2exp(t), -4/3 exp(t/3), -4/3 exp(t/3), 2/3 exp(-5t/3)),
Φ_K=0.
```

The Kasner indices `(2/3,2/3,-1/3)` sum to 1 and their squares sum to 1, so the metric is Ricci-flat, and `∂_a(sqrt(-g)g^{ab})=0` exactly. At all t, an independent symbolic identity checks **all 50** nonlinear GH evolution equations `w_t=F(w)` (use `q=exp(t/3)` and `dt=(q/3)dq`). The spacetime is curved: `R^x_{t x t}=(2/3)(1-2/3)=2/9`.

Let `v=w+εt² cos(x1)` in *three distinct 50-field entries*: `ψ_{01}`, `Π_{22}`, and `Φ_{3,22}` (all other entries untouched). At `t=0`, `v(0)=w(0)` and the exact initial Hamiltonian, momentum, harmonic gauge, and reduction constraints are all zero. The approximate solution has spatial variation and nonzero residual; specifically, `r_{ψ01}=2εt cos(x1)` identically, so it is definitely not an exact solution. The approximate `Φ` reduction constraint need not vanish for t>0, but the exact solution's constraints propagate from exact initial data.

With `T=1/1000`, `R=1/10^8`, `ε=1/10^7`, elementary rational exponential bounds `|exp(qt)-1|≤|q|T/(1-|q|T)` provide the H³ input deviations from the constant `t=0` Kasner state. Each active perturbation has norm `||εt²cos(x1)||3=2εt²`, and its time-derivative norm is `4εt`. Thus `||v-w||3≤6εt²<R` and the known compact residual estimate is

```
ρ ≤18εT + 48(Σ M_Ai)(6εT²) + L_F (6εT²).
```

This uses the exact **symbolically checked algebraic Kasner identity** to rewrite the residual, not an assumption about the unknown true solution of the a posteriori theorem. The same source graph is used to enclose `F(v)-F(w)` *for all intermediate states*; no sampled residual substitutes for this continuum bound. The exact comparator is disclosed, and also independently checks the resulting error enclosure.

The **outward-rounded rational proof record** checks the actual full graph below each bound and then computes the energy constants again from the rounded data:

| Certified input | Exact rational upper/lower |
|---|---|
| symmetrizer eigen | `λ ≥489/500`, `Λ≤511/500` |
| H³ source and Jacobian | `M_F≤10`, `L_F≤525` |
| principal H³ bounds | `M_Ai≤21` for each i |
| principal H³ Jacobians | `L_Ai≤(123,123,119)` |
| principal deviations | `M_ci≤(33/500,33/500,29/500)` |
| pointwise derivative | `J_S≤10`, `J_Bi≤13` |
| metric-only derivatives | `U_g≤33`, `X_g≤13/200` |
| adapted energy divergence | `D≤66507/200=332.535` |
| Gronwall rate | `K<750` (computed exact `<747.219`) |
| continuum 50-field residual | `ρ=19647/5000000000000 <4×10^-9` |
| `KT` | `<3/4` |
| final H³ field error | `E_cert=3346539/206017149996983056 <2×10^-11 <R` |

Every comparison uses `Fraction` arithmetic. The short certificate can be independently checked by hand from its displayed rounded constants; the large exact expression DAG is rerun to establish each outward rounding. It is **not** a proof that the same positive-time estimate holds on BL punctured domains or on a strong-field merger slab.

## 6. Exact subsidiary constraints and vacuum recovery

Original Lindblom et al. eqs (40),(43)–(52) define the complete gauge and reduction constraint family `c^A=(C_a,F_a,C_ia,C_iab,C_ijab)` and prove the homogeneous symmetric-hyperbolic subsidiary equations

```
∂tc^A + A^k{}_B^A(u)∂kc^B = F^A{}_B(u,∂u)c^B.
```

The symmetrizer of their eq (52) is positive when the spatial metric is positive and the auxiliary positive `m_ab, Λ²` are selected. On the **boundaryless** 3-torus, zero initial constraints propagate identically by the usual energy uniqueness argument. For a reduced GH solution, the contracted Bianchi identity also gives `□_g C_a+R_a{}^b C_b=0` (for gamma0=0). If `C_a=0` and the vacuum Hamiltonian/momentum constraints `M_a=0` initially, the paper's eq (16) implies the **normal derivative** of C_a vanishes because its remaining terms consist of tangential derivatives of already-zero C_a. Hence `C=0` and `Ric(g)=0` throughout the verified slab. At Kasner time 0, the spatial 3-metric is δ, the extrinsic-curvature eigenvalues (repository positive-K convention) are `(2/3,2/3,-1/3)`, so `R3=0, tr(K)^2=1=tr(K²)` and momentum constraints are zero (homogeneous), while `sqrt(-g)g^{00}=-1` and no spatial dependence gives initial harmonic gauge zero. `Φ=∂g=0` initially, so all reduction/curl constraints vanish. This applies to the full 50-field system, without using a symmetry reduction in the theorem.

These subsidiary equations and their symmetric hyperbolicity are **existing published results** used with verified hypotheses. The 100-plus explicit algebraic source terms of the subsidiary evolution have not been independently regenerated in the code; that is a transparency limitation, not a claim that numerically small constraints imply vacuum.

## 7. Counterexamples and remaining barriers

- **No fictitious periodic binary:** near a BL puncture, `p~const/r`, hence `g_ij~const r^-4` and `∫_{|x|<δ}|g_ij|² dx~∫_0^δ r^-6dr=∞`. The original puncture chart is not even local `L²` across punctures, let alone a global `H³(T³)` metric. Periodization is not a valid workaround.
- **Wide metrics** fail the exact Neumann/inverse-square-root conditions; **long slabs** fail `KT<1`; **unmatched initial data** must explicitly contribute to the bootstrap, and cannot be omitted.
- **Physical boundaries:** no constraint-preserving high-order incoming-field/excision estimates have been verified on a black-hole exterior. Rinne's 2006 Fourier–Laplace boundary-stability argument does not automatically give the missing nonlinear high-order a posteriori certificate there.
- **Strong field:** small coordinate metric deviation from Minkowski is built into the simple exact-rational graph; no generic large-curvature strong-field functional calculus on irregular domains has been established. The current `T=10^-3` is dimensionless harmonic Kasner time, not a binary mass-scaled merger interval.
- **Novelty:** the Banach functional calculus/geometry-adapted energy specialization and rational majorants are newly integrated *in this repository*, but no new mathematical principle relative to PDE literature is established. The known Kasner spacetime was used to make an explicitly checkable manufactured benchmark, not newly discovered.

## 8. Reproduction and provenance

```
python -m pytest -q tests/test_gh_full_source_majorant.py
python -m research.calc.gh_h3_certificate
python -m pytest -q tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py tests/test_full_hs_bootstrap.py tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py tests/test_bl_gh_initial_data.py tests/test_gh_full_source_majorant.py
```

Proof and source locations: `research/calc/gh_full_source.py` (entire nonlinear 50-field system); `research/calc/gh_h3_majorant.py` (all-tube rational analytic Jacobian evaluator); `research/calc/gh_h3_certificate.py` (small rational witness); `tests/test_gh_full_source_majorant.py` (exact full Kasner source, independent arbitrary homogeneous off-diagonal tests, failures). Older PR #7, #8 and #9 stay unmodified. This branch is draft research pending independent mathematical review and external CI.

**Primary sources:** Kato (1975), Arch. Rat. Mech. Anal. 58; Hughes–Kato–Marsden (1977), Arch. Rat. Mech. Anal. 63, sections 2–3; Kato–Ponce (1988), CPAM 41; Friedrich–Rendall (2000), gr-qc/0002074; Lindblom et al. (2006), gr-qc/0512093 v3 (specifically equations 35–39 and 40–52); Rinne (2006), gr-qc/0606053. We did not exhaustively read every referenced mathematics paper's full proof, and no novelty priority search is claimed complete.