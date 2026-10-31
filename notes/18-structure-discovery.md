# 18 · Where the obstruction lives: a bottleneck map, a local/nonlocal split, and foliations

Sessions 11–13 built and repaired a category-B metric. It is the geodesic-slicing development of
equal-mass Brill–Lindquist (BL) data, written as a convergent series in $t$ with certified tails.
Sessions 10–13 excluded several closed-form classes. This note asks a different question. Which
parts of the problem are hard, and does the B-class solution reveal a simplifying structure?

- **§§1–3** fix the baseline, define what would count as essential progress, and map the
  bottlenecks with evidence.
- **§4** audits the scope of every earlier no-go claim.
- **§5** is a portfolio of eight falsifiable hypotheses and the reasons for the choice.
- **§§6–7** are two executed experiments of different kinds:
  - **H1, discovery.** Split the development into a local (silent) part, which is a closed ODE on
    each worldline, and a nonlocal remainder.
  - **H2, gauge.** Can a different slicing, with a bounded maximal-type lapse, have conformally flat
    leaves?
- **§8** documents the tools; **§9** gives the strongest surviving statements; **§10** the
  continuation.

*Status (session 14).*

| Item | Category |
| --- | --- |
| silent system reproduces the development exactly through $t^4$; first deviation at $t^6$ | derived; exact rational checks |
| the deviation through $t^6$ is a single traceless source $(\mathrm{curl}\,H)_3$, with a four-term covariant law | chain derived; law exact at 3 rational points, held out at 3 ball points |
| near the holes the remainder is at most 13 % of each coefficient to $t^{16}$; elsewhere it dominates by $t^{16}$ | numerical (ball midpoints), exploratory |
| the nonlocal remainder has the smallest radius estimate at every point tested; the series' radius is not the local collapse time | numerical, exploratory (fits to $t^{16}$, two runs to $t^{20}$) |
| no slicing with lapse $N = (c_0 + c_+/r_+ + c_-/r_-)/\psi$ (the bounded maximal-type lapses, including geodesic) has leaves conformally flat at $O(t^2)$ | computer-assisted (Arb, box covering), at $d = 0.7662,\ 0.767,\ 1.5$ |
| at each single point such a lapse does exist (two conics, real common zeros); exclusion needs two points | numerical (residual $\le 10^{-16}$) plus a parity argument and rank check |
| category-A metric | **not found** |
| structural breakthrough | **no** |

Conventions are those of note 00 and note 15. Signature $(-,+,+,+)$, MTW/Carroll Riemann sign,
$G = c = M_{\rm ADM} = 1$. The BL bare masses are $\tfrac12$ at $z = \pm d/2$, with
$\psi = 1 + 1/(4r_+) + 1/(4r_-)$. Geodesic slicing means $g = -dt^2 + \gamma_{ij}dx^idx^j$, with the
repo sign $K_{ij} = \tfrac12\partial_t\gamma_{ij}$. Mixed-index matrices are $K = K^i{}_j$ and
$E = E^i{}_j$. For conformally flat $\gamma_0$ these are symmetric in Cartesian coordinates.

## 1 · Baseline: what the evaluator computes, and what is known

**Baseline check (rerun this session).**

- The structure tests pass (5 tests, 28 s).
- `evaluate_metric(2e-4, 0.17, 0.11, 0.5, d = 0.7662, order = 5)` returns *certified*.
  - Metric error $4.8\times10^{-10}$; first derivatives $8.0\times10^{-7}$.
  - It combines the `ext`, `sharp` and `classic` certificates.

No large certificate was rerun.

**What the evaluator computes** (notes/16 §3, notes/17 §4).

- **Output.** Enclosures of $g_{\mu\nu}$, $g^{\mu\nu}$, $\partial g$ and $\partial^2 g$ of the
  geodesic development.
- **Domain.** $(t, x, d) \in (-T, T)\times U\times I$, with
  $I = [d_c - 10^{-13}, 0.767578125]$:
  - $U$ is the 0.02-tube about the critical surface $S_c$, with $T = 8.683\times10^{-4}$;
  - $U_w$ is the waist region, with $T = 1.545\times10^{-3}$.
- **Errors.** Itemised tails from the CK′ schedules.
- **Outside the domain** it refuses.
- **Uncertified mode.** The truncated series can be evaluated anywhere, without guarantees. Fitted
  radius estimates are about 2.1 near the holes and 2.5–4.8 away from them (§6.4).

**Classes excluded before this session, with their exact scope.**

| Claim | Exact scope |
| --- | --- |
| 68–69 | every algebraically special ansatz (single KS, KSK, RT, Kundt, twisting II/D) on an open set meeting the type-I region |
| 87 | conformally flat leaves of the **geodesic** slicing |
| 88 | silent flow ($B = 0$) for the **geodesic** observers |
| 99 | $\gamma(t)$ pointwise algebraic in $(\gamma_0, \mathrm{Ric})$ |
| 100 | closed isotropic tidal laws $E' = f(E)$, and frozen eigenvalue ratios |

**Proposed but never executed.**

- Two-vector Kerr–Schild solvability (route B; task 18).
- The Einstein–Bianchi energy step (notes/17 §3).
- A Geroch-reduction test (route A).
- Alternative data (Misner, Bowen–York).
- The mass-expansion (linear-in-interaction) sector.

**Resting on stored certificates or unreviewed proofs.**

- Theorems C, F and F′, Proposition B and Corollary J: stored logs and certificates, inherited.
- Theorem CK′ (notes/17 §2): proved by us, not externally reviewed.
- The $t^6$ law (claim 98): fitted, not derived by hand.

## 2 · Essential progress

Count a result as *structural* if it does one of the following:

- (a) replaces an infinite family of coefficients by a finite rule;
- (b) identifies a reduction of the field equations, a slicing, or a variable set in which the
  BL development is simpler, and certifies it;
- (c) proves that a whole natural class (not one ansatz) contains no such simplification;
- (d) extends the certified domain by a new mechanism rather than by a better constant.

*Incremental* means tighter constants, more points checked, or a new fitted law at one more order.

By this standard:

- H1's first-deviation law is incremental, an extension of claim 98.
- H1's quantitative split, local dominant near the holes and nonlocal elsewhere, is a structural
  *fact* about where the difficulty lives. It is exploratory evidence, not a theorem.
- H2's exclusion of the whole bounded maximal-type lapse family is a structural no-go (type c),
  of limited scope.

## 3 · Bottleneck map

| Bottleneck | Evidence that it is real | What it blocks |
| --- | --- | --- |
| **coordinate / gauge** | the series is only known in geodesic slicing; geodesic slicing crunches; its leaves are not conformally flat (claim 87); §7: no bounded maximal-type lapse gives conformally flat leaves either | a simpler slicing in which the metric is one scalar |
| **physical coupling (nonlocality)** | the magnetic Weyl tensor is $O(t^3)$, not zero (claim 88); §6: the nonlocal remainder, driven by $\mathrm{curl}\,H$, dominates the coefficients away from the holes by $t^{16}$ and has the smallest fitted radius everywhere | any pointwise (ODE-along-worldlines) closed form, local resummation |
| **initial-data restrictions** | all results use BL (time-symmetric, conformally flat, equal mass, head-on). Conformal flatness kills $H$ at order $t$; time symmetry kills odd orders | statements about generic binaries (spin, momentum, unequal mass) |
| **representation limits** | the representation is a Taylor series in $t$, so it ends at the nearest complex singularity: numerical radius estimates ≈ 2.1 near the holes and larger away from them (claim 72, §6.4), certified radius $\le 1.5\times10^{-3}$ | evaluation at merger times $t_* \sim 2.4(d - d_c)$ only for $d - d_c \lesssim 6\times10^{-4}$; nothing at late times |
| **proof / error-method limits** | abstract CK needs a complex neighbourhood. Its time is capped by the punctures' singular cones and by $C \ge 2\sqrt3\sup\lvert\psi\rvert^{-2}$ (claim 90). The energy method is not built (claim 102) | certified times $\gtrsim 10^{-2}$; an explicit $\delta_1$ |
| **missing global information** | no event horizon, no $\mathscr I^+$, no late-time behaviour from a local series. Theorem H's common MOTS is not shown to be the first apparent horizon (claim 81) | the physical merger statement and radiation |

## 4 · Scope audit of earlier no-go claims

| Earlier statement | What it actually covers | What it does not |
| --- | --- | --- |
| "algebraically special ansätze excluded" (68–69) | the BL development, which is type I on an open dense set | other data whose development might be special somewhere; non-algebraic ansätze |
| "conformally flat slices excluded" (87) | the **geodesic** foliation, at $O(t^2)$ | other foliations. §7 extends this to the bounded maximal-type lapses; arbitrary lapses remain open (§7.5) |
| "silent flow excluded" (88) | the statement $B \equiv 0$ for geodesic observers | the *silent part* as an approximation or organising principle. §6 shows it is exact through $t^4$ and dominant near the holes |
| "pointwise algebraic in Ric excluded" (99) | $\gamma(t) = F(t, \gamma_0, \mathrm{Ric})$ | dependence on $\nabla\mathrm{Ric}$ (it does enter at $t^6$, §6.3); non-pointwise closed forms |
| "isotropic tidal laws excluded" (100) | closed laws $E' = f(E)$ along each observer | laws with source terms. §6 is exactly such a law, with source $\mathrm{curl}\,H$ |
| "series radius ≈ 2.7" (72) | a numerical estimate at one point near $S_c$ | a property of the spacetime. Padé poles are not singularities (claim 89); §6.4 suggests that it is set by the nonlocal remainder and varies with $x$ |
| "CK limits" (90) | one bound formula for the abstract CK method | the spacetime: a method ceiling is not a ceiling of the solution (session 13 corrected the wording) |
| "no category-A metric" (75, 92) | the routes and ansätze listed | category A in general: nothing here proves that no closed form exists |

BL is not all binaries, and geodesic slicing is not all slicings. Every exclusion above is
relative to BL data, and all but claims 68–69 and 99–100 are relative to geodesic slicing.

## 5 · Hypothesis portfolio

| # | Hypothesis (mechanism) | Bottleneck | Tool | Precedent | Novelty vs repo | Smallest decisive calculation | Positive result | Stop condition | Path to a metric |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | **Silent/nonlocal split.** The local tidal ODE (drop $\mathrm{curl}\,H$) carries the development; the nonlocal remainder is small and smoother | coupling, representation | matrix-series ODE plus jet recursion | silent universes [MPS94, BMP95]; local tidal approximation [BJ94]; integrability obstruction [vEULEM97] | first test on vacuum binary data with an exact remainder law | compare the true and silent coefficients to $t^{16}$ at generic points, plus the radius of the remainder | remainder $\ll$ silent with a larger radius, so $g = g_{\rm sil} + $ a small analytic correction | remainder comparable to the coefficients, or no larger radius | semi-explicit B/D-class metric: an ODE per worldline plus a series |
| H2 | **Conformally flat leaves in another slicing.** Some bounded maximal-type lapse makes the leaves conformally flat | gauge | linearised Cotton operator plus box covering of $\mathbb{RP}^2$ | Kerr: no conformally flat slices [GP00, VK04] | first test for BL beyond geodesic slicing | $DC[h_c] = 0$ at points, for all $c$ | a common $c$ at every point, giving the metric $\Phi^4\delta$ to $O(t^2)$ | certified absence of a common zero | one scalar $\Phi(t, x)$ plus shift: a large reduction |
| H3 | **Superposed free falls.** Combine the two single-hole cycloid maps | coupling | exact Schwarzschild free fall | Brill–Lindquist superposition | subsumed by H1: the silent ODE already contains each single hole exactly | — | — | — | — |
| H4 | **Two-vector Kerr–Schild.** $g = \eta + \lambda_1 k_1k_1 + \lambda_2 k_2k_2$ | representation | algebraic solvability at $t^2$ | multi-KS literature | route B, deferred since session 9 | rank of the $t^2$ equations | solvable family | rank deficiency absent | explicit A-class candidate |
| H5 | **Energy step.** Continue past CK with Einstein–Bianchi estimates | proof method | $H^s$ defect bounds | standard hyperbolic theory | planned in session 13 | one certified step | certified $T \sim 10^{-2}$ | constants unattainable | larger B-class domain |
| H6 | **Geroch reduction.** Axisymmetric twist-free vacuum is 2+1 gravity plus a wave map; test the 2+1 sector for simplicity | coupling, gauge | Geroch 1971 reduction | Ernst, Geroch | route A, screened | 2+1 Cotton / wave-map map at $t^2$ | a harmonic-map structure (e.g. totally geodesic) | non-special at $t^2$ | reduce to a scalar wave map |
| H7 | **Alternative data.** Misner or Bowen–York data with simpler developments | data | the same engines | standard | none executed | redo H1/H2 for Misner data | a simpler class | same obstructions | a different A/B family |
| H8 | **Interaction sector.** The $m_+m_-$ cross-terms form a linear hierarchy with a closed form | coupling | mass expansion of the jet recursion | post-Minkowskian-type expansions | new for the repo | closed form of the $m_+m_-$ part at $t^4$, $t^6$ | a recognisable kernel | non-elementary at $t^4$ | perturbative A-class pieces |

**Ranking.** By information per cost: H1 > H2 > H8 > H6 > H4 > H7 > H5.

- H1 and H2 are cheap with the existing engines, are of different kinds (representation versus
  gauge), and each has a sharp positive outcome.
- H1 is the discovery experiment. It looks for a simplifying structure rather than excluding an
  ansatz.
- H5 is the most expensive and is unchanged since session 13.
- H3 is subsumed by H1.

Executed: **H1 and H2**.

## 6 · Experiment H1: the silent part and the nonlocal remainder

### 6.1 Candidate class (fixed before any fit)

Along each comoving geodesic, with $\theta = \mathrm{tr}K$ and $\sigma = K - \tfrac\theta3 I$, the
geodesic-slicing vacuum system in mixed components is

```math
\dot\gamma = 2\gamma K,\qquad \dot K = -E - K^2,\qquad
\dot E = -\theta E + 3(\sigma E)^{\rm TF} + [E, K] + \mathrm{curl}\,H .
```

- The first two equations are exact. They are the evolution equation rewritten with the vacuum
  Gauss equation $E_{ij} = R_{ij} + KK_{ij} - K_{ik}K^k{}_j$.
- The third is the Bianchi identity in a Lie-dragged basis. The commutator converts the Fermi
  derivative, so that $\partial_tE = \nabla_uE + [E, K]$.
- The 1+3 form is that of Ellis and van Elst [EvE98], with $\rho = p = 0$ and no vorticity or
  acceleration.

**The silent system** drops $\mathrm{curl}\,H$, the only term with spatial derivatives. Starting
from $K = 0$ and $E = E_0 = \mathrm{Ric}(\gamma_0)$ (mixed), every matrix stays a polynomial in
$E_0$, so the commutator vanishes. Write $\gamma_{\rm sil} = \gamma_0 M(t)$. The silent metric at a
point is then fixed by $(\gamma_0, E_0)$ there: it is a closed ODE in the three eigenvalues.

**The candidate** is $\gamma = \gamma_{\rm sil} + \rho$, with $\rho$ small and more regular.

### 6.2 What the class must satisfy, and the exact low orders

For conformally flat time-symmetric data, $H = -t\,C(\gamma_0) + O(t^3) = O(t^3)$ (claim 73). So
$\mathrm{curl}\,H$ first enters $\dot E$ at $t^3$.

Write $c_3 = (\mathrm{curl}\,H)_3$ for its $t^3$ coefficient. Propagating it through the system:

- it enters $E$ at $t^4$ as $c_3/4$;
- it enters $K$ at $t^5$ as $-c_3/20$;
- it enters $\gamma_0^{-1}\gamma$ at $t^6$ as $-c_3/60$.

The commutator contributes nothing at these orders, because $E_2$, $K_1$ and $K_3$ are
polynomials in $E_0$. Hence two exact predictions:

```math
\gamma = \gamma_0 M(t) + O(t^6),\qquad
\gamma_0^{-1}\gamma_3 - M_6 = -\tfrac1{15}\,(E_4 - E_4^{\rm sil}) .
```

Both hold **exactly** in rational arithmetic at a BL plane point and at a BL axis point.

- The true coefficients come from the jet recursion of `bl_metric.series_jets`.
- The true $E_4$ comes independently from `bl_structure.electric_series`.
- The $t^6$ coefficients differ.
- At a generic point (256-bit balls) the chain holds to $10^{-68}$.

The deviation is nonzero on the axis too, so this is not a symmetry artefact.

**Silent $t^6$ coefficient** (exact, for any traceless $E_0$):

```math
\gamma_0 M_6 = -\tfrac{13}{180}\,\mathrm{tr}(E_0^2)\,E_0 - \tfrac{41}{540}\,\mathrm{tr}(E_0^3)\,\gamma_0 .
```

Compare the true coefficient (claim 98). The $\mathrm{tr}(\mathrm{Ric}^3)$ coefficient agrees. The
**first failing term** is $\lvert\mathrm{Ric}\rvert^2\mathrm{Ric}$: $-\tfrac{12}{180}$ true versus
$-\tfrac{13}{180}$ silent. Three gradient terms are also missing from the silent coefficient.

### 6.3 The first nonlocal source has a covariant law

Combining the two previous displays:

```math
(\mathrm{curl}\,H)_3 = -\tfrac13\lvert\mathrm{Ric}\rvert^2\mathrm{Ric}
 - \nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric}
 + \tfrac13\lvert\nabla\mathrm{Ric}\rvert^2\gamma_0
 + \tfrac13\,R^{kl}\nabla_k\nabla_lR_{ij} .
```

Here $\nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric} = \nabla_iR_{kl}\nabla_jR^{kl}$. The law is
traceless, as it must be, because $R = 0$ and $\nabla\mathrm{Ric}$ is totally symmetric.

It is checked directly against $4(E_4 - E_4^{\rm sil})$ from the electric-Weyl engine:

- **exactly** at the BL plane and axis points;
- on a Schwarzschild point, where both sides vanish;
- at three **held-out** ball points, residual $\le 5\times10^{-65}$. These include separations
  $d = 0.7$ and $d = 1.5$, which were not used anywhere.

Its evidence category is that of claim 98. It is an identity verified exactly at rational points
and held out at others, not derived by hand.

This answers the question left open in notes/17 (Discussion prompts). In the variables
$(E, K, \gamma)$ the tidal gradients enter as one explicit traceless source. At first order that
source is a cubic polynomial in Ric and its first two covariant derivatives.

### 6.4 High orders: where the silent part dominates

At six generic points (no coordinate zero; $d = 0.7662$, P2 at $d = 0.767$) we computed:

- the true coefficients $\gamma_k$, $k \le 8$ (through $t^{16}$), in 256-bit balls, with the
  largest radius $\le 10^{-55}$;
- the silent ones to $t^{80}$;
- the remainder.

Write $\lvert\cdot\rvert$ for the $\gamma_0$-orthonormal Frobenius norm of mixed matrices.

- $\rho_k$ is the remainder of the coefficient of $t^{2k}$.
- $T_{\rm fit}$ is a radius estimate in $t$ from a least-squares fit of
  $\log\lvert c_k\rvert = a - 2k\log T$ over $4 \le k \le 8$.
- The silent fit uses the same orders.

Everything in the table is computed from ball midpoints, so it is exploratory.

| Point | $E_0$ eigenvalues | $\lvert\rho_3\rvert/\lvert\gamma_3\rvert$ | $\lvert\rho_8\rvert/\lvert\gamma_8\rvert$ | $T_{\rm fit}$ true | $T_{\rm fit}$ silent | silent collapse time | $T_{\rm fit}$ remainder |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P4 $(0.05, 0.02, 0.6)$ | $(-0.585, 0.292, 0.293)$ | 0.010 | 0.131 | 2.137 | 2.111 | 2.054 | 1.755 |
| P1 $(0.17, 0.11, 0.5)$ | $(-0.558, 0.273, 0.286)$ | 0.010 | 0.039 | 2.167 | 2.161 | 2.102 | 1.864 |
| P2 $(0.26, -0.09, 0.12)$ | $(-0.341, 0.118, 0.224)$ | 0.161 | 0.782 | 2.507 | 2.760 | 2.682 | 2.236 |
| P5 $(0.33, 0.21, 0.05)$ | $(-0.180, -0.001, 0.181)$ | 0.514 | 1.007 | 3.042 | 3.710 | 3.589 | 2.939 |
| P3 $(0.6, 0.3, 0.9)$ | $(-0.195, 0.092, 0.103)$ | 0.076 | 0.959 | 3.417 | 3.657 | 3.558 | 2.770 |
| P6 $(1.1, -0.6, 1.4)$ | $(-0.078, 0.038, 0.040)$ | 0.091 | 0.990 | 4.768 | 5.769 | 5.612 | 3.882 |

**Two longer runs, to $t^{20}$** ($N = 10$, about 19 minutes each):

- **P1.** The relative remainder rises to 0.055 and then 0.089 at $k = 9, 10$. The fits are 2.159
  (true) and 1.847 (remainder).
- **P3.** The relative remainder reaches 0.986. The fits are 3.309 (true) and 2.847 (remainder).

At P3 the single ratios $\lvert c_k/c_{k+1}\rvert^{1/2}$ oscillate (3.84, 3.81, 3.58, 2.41, 3.06,
4.13). This is a warning against reading a radius off one ratio. An earlier draft of this section
did exactly that, from the dip at $k = 7$.

**Sign patterns** (a coordinate-dependent heuristic).

- The remainder's $xx$ components keep one sign from $k = 4$ through $t^{16}$ at P1 and P4. At P3
  and P6 they change sign only at the last order computed ($t^{20}$ and $t^{16}$). This suggests a
  singularity near the positive $s = t^2$ axis, i.e. real $t$.
- They alternate at P2 and P5, which suggests the negative axis, i.e. imaginary $t$.

How the collapse times were obtained:

- **Silent collapse time.** Integration of the silent ODE (DOP853, rtol $10^{-12}$) until a scale
  factor reaches $10^{-6}$. Its series to $t^{80}$ agrees (P1: fit 2.117).
- **Type-D control.** For $E_0 = \mu\,\mathrm{diag}(-2, 1, 1)$ the collapse time is
  $\pi/\sqrt{8\mu}$ to $10^{-6}$. This is Schwarzschild free fall, $\pi\sqrt{R_0^3/8M}$
  (`tests/test_bl_silent.py`).

**Reading** (exploratory throughout).

1. **Near the holes** the tidal tensor is nearly of type D, $(-2\mu, \mu, \mu)$. There the silent
   part carries every coefficient through $t^{16}$ up to a remainder of at most 13 % (at most 9 % at
   P1 through $t^{20}$). On the same orders the true and silent fits agree to about 1 %. The local
   collapse is the dominant structure at these orders.
2. **Elsewhere** (the waist P2, the $z \approx 0$ plane P5, and the far zone P3, P6) the remainder is
   comparable to or larger than the silent part by $t^{16}$. The true fits lie 9–18 % below the
   silent ones.
3. **Everywhere** the remainder's fit is the smallest of the three, 13–33 % below the silent fit.
   The remainder's relative size grows with order at every point, near the holes too. So the split
   does not produce a remainder of larger radius.

   If the trend continues, the series' radius is set by the nonlocal remainder, not by the local
   collapse time, at every point tested. Near the holes it would then be about 1.75–1.86, earlier than
   the local collapse time of about 2.1.
4. **The remainder's fit grows with distance from the holes,** from 1.76 to 3.88. At these orders no
   single common radius is visible.

### 6.5 Decision

- **Kept (exact):**
  - the silent system as the exact $t^{\le4}$ structure;
  - the single-source chain at $t^6$;
  - the covariant law for $(\mathrm{curl}\,H)_3$.
- **Rejected, within this scope** (geodesic slicing, BL data, $d \approx 0.766$, coefficients to
  $t^{16}$ at six points and to $t^{20}$ at two):
  - the hypothesis "$g = g_{\rm sil} + $ small, smoother remainder" as a global representation;
  - the hypothesis "the series' radius is the local collapse time".

  The evidence is the remainder's relative size, which grows with order everywhere, and its smaller
  fitted radius. It is numerical and exploratory, not certified.
- **Revised:** the silent part is a good local model near the holes (a D-class approximation there,
  with error measured above). It is not an organising principle for the whole slice.

## 7 · Experiment H2: conformally flat leaves for bounded maximal-type lapses

### 7.1 Candidate class

Take a zero-shift slicing through the time-symmetric slice with lapse $N(t, x)$. Then

```math
\gamma(t) = \gamma_0 + t^2 h + O(t^3),\qquad h = N\,D_iD_jN - N^2\,\mathrm{Ric}(\gamma_0),\quad N = N(0, x).
```

A shift adds $\mathcal L_\beta\gamma_0$ at this order.

**The family.** For $\gamma_0 = \psi^4\delta$, $\Delta_\gamma(\phi/\psi) = \psi^{-5}\Delta\phi$. So
the lapses with $D^2N = 0$ at $t = 0$ (maximal slicing, $D^2N = N\lvert K\rvert^2$) that are bounded
on all three ends are exactly

```math
N_c = \frac{c_0 + c_+/r_+ + c_-/r_-}{\psi},\qquad c \in \mathbb R^3 ,
```

by Bôcher's theorem at the punctures and Liouville's at infinity. The family contains:

- the geodesic slicing, $c = (1, \tfrac14, \tfrac14)$, where $N = 1$;
- the antisymmetric ("odd") lapse $c = (1, -\tfrac14, -\tfrac14)$ of the classic two-hole maximal
  slicings [SCDE76].

### 7.2 Equations

The 3D Cotton tensor $C_{ijk} = \nabla_kS_{ij} - \nabla_jS_{ik}$ vanishes iff the metric is
conformally flat. Since $C[\gamma_0] = 0$, conformally flat leaves require

```math
DC[h_c] = 0 \quad\text{everywhere},
```

where $DC$ is the linearised Cotton operator at $\gamma_0$.

- $DC$ annihilates $\mathcal L_\xi\gamma_0$ and $\phi\gamma_0$.
- With $N\,DDN = \tfrac12\mathcal L_{N\nabla N}\gamma_0 - dN\otimes dN$, this gives the independent
  form $DC[h] = -DC[dN\otimes dN + N^2\mathrm{Ric}]$.
- Pointwise, $DC[h_c](x)$ is a tensor of conics in $c$.

### 7.3 An unexpected pointwise structure

**Rank 2.** At every tested point the nine Cotton components $C_{ijk}$, $j < k$, span only **two**
independent conics. The singular values of the $9\times6$ coefficient matrix are two nonzero and
four zero.

The reason is parity.

- The lapses, $\gamma_0$ and $h_c$ are axisymmetric and invariant under the reflection
  $\varphi \mapsto -\varphi$.
- The Cotton–York tensor is a pseudo-tensor, so only its components odd under that reflection,
  $(\rho\varphi)$ and $(z\varphi)$, can be nonzero.

**Real common zeros.** Two conics in $\mathbb{RP}^2$ meet in up to four points, often real. At
every tested point we found real common zeros numerically (residual $\le 10^{-16}$), for example
$c^\ast = (1, -0.0349, -0.4056)$ and $(1, 0.3377, 1.6039)$ at $(0.17, 0.11, 0.5)$.

**So at each single point there is a lapse in the family for which the leaf is conformally flat
there to $O(t^2)$.** A one-point certificate must fail, and does.

The zeros move with $x$: at $(0.26, -0.09, 0.12)$ one is $(1, -0.110, -0.151)$.

### 7.4 Certificate over two points

Take the conics of two points and cover the faces $\{c_i = 1, \lvert c_j\rvert \le 1\}$. This
suffices because $h_{-c} = h_c$. On each box, some conic from either point must have a ball
enclosure excluding 0 (Arb, 128 bits, enclosures of the conic coefficients at 256 bits).

| Separation | Points | Boxes | Depth | Result |
| --- | --- | --- | --- | --- |
| $d = 0.7662$ | $(0.17, 0.11, 0.5)$, $(0.6, 0.3, 0.9)$ | 247 | 11 | **no common zero** |
| $d = 0.767$ | $(0.26, -0.09, 0.12)$, $(0.05, 0.3, -0.2)$ | 127 | 14 | **no common zero** |
| $d = 1.5$ | $(0.4, 0.2, 1.1)$, $(0.3, 0.5, -0.6)$ | 241 | 12 | **no common zero** |

Controls, all passed (`tests/test_bl_foliation.py`):

- **Pure gauge.** $DC[\mathcal L_\xi\gamma_0]$ and $DC[\phi\gamma_0]$ vanish to $10^{-68}$ for
  polynomial $\xi$ and $\phi$. A non-gauge perturbation gives $1.7\times10^{-2}$.
- **Two independent forms.** The two forms of $h$ agree to $5\times10^{-66}$ in balls, and exactly
  at a rational point.
- **Cross-check with session 12.** The geodesic member $c = (1, \tfrac14, \tfrac14)$ reproduces the
  session-12 Cotton coefficients to 66 digits, for example $C_{001} = -1.0629704725891769\ldots$.
  That value came from a different engine (claim 87).
- **Schwarzschild.** For one hole every lapse $(c_0 + c_1/r)/\psi$ gives $DC = 0$ (spherical
  symmetry). An off-centre $1/\lvert x - q\rvert$ does not.

Treating $d$ as a ball over the whole window $I$ (radius $8\times10^{-4}$) fails. The linearised
Cotton operator cancels large terms, so the enclosures are wider than the values. The statement is
therefore certified at the three listed separations, not on an interval.

### 7.5 Decision and scope

**Rejected (certified):** "some slicing whose initial lapse is bounded and maximal-compatible has
conformally flat leaves." This holds at $d = 0.7662$, $0.767$ and $1.5$, at the first nontrivial
order $t^2$, for any shift.

It does **not** cover:

- lapses outside the family (arbitrary $N_0$);
- slicings not through the time-symmetric slice;
- other data.

**Revised.** For arbitrary axisymmetric $N_0(\rho, z)$ the condition is two fifth-order equations
for one function of two variables. That system is overdetermined but only mildly. Jet counting
gives $2(k - 4)$ equations against $k + 1$ new unknowns at order $k$, so the equations outnumber
the unknowns from $k = 10$ on. Whether it has local solutions is open. It is the natural sequel
(§10).

## 8 · Tools

| Tool | Input → output | Can prove / falsify | Heuristic part | Benchmark with known answer | Binary result |
| --- | --- | --- | --- | --- | --- |
| `bl_silent.silent_series` | $E_0$ (rational or ball) → Taylor coefficients of $(K, E, M)$ of the silent ODE | exact low-order identities; exact silent coefficients | none (exact recursion) | Schwarzschild: silent = true through $t^8$, exactly | silent = true through $t^4$, first failure at $t^6$ |
| `bl_silent.remainder_chain`, `curl_h3_law` | point → $\delta M_6$, $\delta E_4$, the covariant law | exact or ball identities | the law was fitted (via claim 98), then verified | Schwarzschild: both sides 0 | law exact at 2 BL rational points, held out at 3 ball points |
| `bl_silent.high_order`, `silent_crunch` | point → coefficient norms to $t^{16}$, ratio estimates, silent collapse time | nothing (exploratory) | ratio estimates, midpoints | type D: $\pi/\sqrt{8\mu}$ | table of §6.4 |
| `bl_foliation.quadratic_forms` | point, lapse basis → conics $DC[h_{ab}]$ | exact or ball values of the forms | none | gauge controls; Schwarzschild; reproduces claim 87 | rank 2, pointwise zeros |
| `bl_foliation.certify_no_common_zero` | conics at one or more points → covering certificate | absence of a common real zero on $\mathbb{RP}^2$ | none (outward-rounded balls) | fails on single points (where zeros exist), as it must | certified at 3 separations |

Replay:

```bash
python -m research.calc.bl_silent
```

```bash
python -m research.calc.bl_silent long
```

```bash
python -m research.calc.bl_foliation
```

Tests: `tests/test_bl_silent.py` (5 tests) and `tests/test_bl_foliation.py` (3 tests), about 35 s
in total.

## 9 · Strongest surviving statements

1. **Exact (derived, rational checks).**
   - For conformally flat time-symmetric vacuum data in geodesic slicing,
     $\gamma = \gamma_0M_{\rm sil}(t) + O(t^6)$.
   - The $t^6$ deviation is $-\tfrac1{60}\gamma_0(\mathrm{curl}\,H)_3$.
   - Category: identity (R), plus an exact check.
2. **Identity, verified exactly at rational points and held out.** $(\mathrm{curl}\,H)_3$ obeys the
   four-term law of §6.3. Category: as claim 98.
3. **Computer-assisted (certified).** At $d \in \{0.7662, 0.767, 1.5\}$, no slicing with initial
   lapse in the bounded maximal-type family has leaves conformally flat to $O(t^2)$, with any
   shift. At each single point such a lapse exists, numerically.
4. **Exploratory.**
   - Near the holes the local (silent) collapse dominates the development to $t^{16}$–$t^{20}$.
     Away from them the nonlocal remainder dominates.
   - Everywhere the remainder has the smallest fitted radius, and that radius grows with distance
     from the holes.

No A-class candidate survived. The B-class representation gained an organising principle:
local ODE plus a nonlocal source with a known leading law. Its scope and computability are
unchanged.

## 10 · Physical scope and continuation

- **Scope.** Everything is on BL data in a short-time series. Nothing here certifies a merger, a
  horizon or radiation. The silent collapse times (≈ 2.1 near the holes) are local model times,
  not predictions of the true singularity.
- **Continuation note.**
  - Start from §§6.4 and 7.5.
  - Records: `research/calc/data/bl_silent.json`, `bl_silent_N10.json` (the runs to $t^{20}$) and
    `bl_foliation.json`.
- **Next single decisive experiment.** Locate the singular set $t = \tau(x)$ of the development.
  Compute the coefficients to $t^{24}$ ($N = 12$, about an hour per point) at points along a ray
  from a throat to the far zone. Fit the remainder's radius with Domb–Sykes or Padé, and read its
  sign pattern. The decision:
  - Is $\tau$ real (a physical singularity reached before the local collapse, about 1.75–1.86 near
    the holes) or complex?
  - Does it vary smoothly with $x$, as the $t^{16}$ data suggest?

  A real, smoothly varying $\tau(x)$ would make a resummation in $s/\tau(x)^2$ the natural next
  representation. A complex one would require continuation by a nonlocal method (energy estimates,
  or analytic continuation in $x$).

## References

- [BJ94] E. Bertschinger and B. Jain, *Gravitational instability of cold matter*, ApJ 431 (1994) 486.
- [BMP95] M. Bruni, S. Matarrese and O. Pantano, *Dynamics of silent universes*, ApJ 445 (1995) 958.
- [EvE98] G. F. R. Ellis and H. van Elst, *Cosmological models* (Cargèse lectures 1998),
  gr-qc/9812046, §§3–4 (the 1+3 Bianchi equations used in §6.1).
- [GP00] A. Garat and R. H. Price, *Nonexistence of conformally flat slices of the Kerr spacetime*,
  Phys. Rev. D 61 (2000) 124011, gr-qc/0002013.
- [MPS94] S. Matarrese, O. Pantano and D. Sáez, *General relativistic dynamics of irrotational
  dust: cosmological implications*, Phys. Rev. Lett. 72 (1994) 320.
- [SCDE76] L. Smarr, A. Čadež, B. DeWitt and K. Eppley, *Collision of two black holes: theoretical
  framework*, Phys. Rev. D 14 (1976) 2443.
- [vEULEM97] H. van Elst, C. Uggla, W. M. Lesame, G. F. R. Ellis and R. Maartens, *Integrability of
  irrotational silent cosmological models*, Class. Quantum Grav. 14 (1997) 1151, gr-qc/9611002.
  The integrability conditions of nonlinear silent models are not automatically satisfied, and
  inhomogeneous Petrov type-I silent solutions are conjectured not to exist. This is consistent
  with claim 88 and with §6.
- [VK04] J. A. Valiente Kroon, *On the nonexistence of conformally flat slices in the Kerr and other
  stationary spacetimes*, Phys. Rev. Lett. 92 (2004) 041101, gr-qc/0310048.

**Novelty**, classified conservatively:

- The silent system is standard [MPS94, BMP95, EvE98].
- Its use as an exact low-order skeleton of a vacuum binary development, and the explicit first
  source law, were not found in the sources we checked. That was a bounded search, so no
  originality is claimed.
- H2 is a BL analogue of the Kerr nonexistence results [GP00, VK04], restricted to one lapse family.

## Discussion prompts

- The silent ODE depends only on the three tidal eigenvalues at a point. Near the holes it
  reproduces the development to a few percent through $t^{16}$. Is there a matched description:
  silent near the holes, linearised far away, glued in between?
- Each point admits a lapse that makes the leaf conformally flat *there*. Can the pointwise zeros
  $c^\ast(x)$ be promoted to a field-dependent lapse $N_0(x)$ solving $DC[h] = 0$? It would no
  longer be maximal. What are the compatibility conditions of that overdetermined system?
- If the far-zone radius of the series is set by the near-hole singularity, what is the shape of
  the complex singular set $t = \tau(x)$? Does it relate to the formation of the common apparent
  horizon?
