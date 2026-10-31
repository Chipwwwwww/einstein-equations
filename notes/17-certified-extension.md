# 17 · A corrected and extended certificate, better derivative bounds, and the series' geometry

This note does three things.

- **Repairs the foundation** of notes 15–16. The constant of the abstract Cauchy–Kovalevskaya (CK)
  theorem had been misread, so every certified time there was too large (§1).
- **Proves and applies a generalised theorem**, CK′ (§2), to the same system. The certified
  interval becomes longer than the corrected baseline: $T = 8.683\times10^{-4}$ on the whole tube and
  $T = 1.545\times10^{-3}$ on the waist region (§3).
- **Improves the evaluator's derivative bounds** (§4), and runs two bounded searches for closed-form
  structure (§§5–6): a covariant $t^6$ law and geometric (Weyl) variables.

It ends with a single-hole control that fixes the meaning of the series' convergence time (§7).

*Status (session 13).*

| Item | Category |
| --- | --- |
| BG constant misread in notes 15–16; corrected times $2.885\times10^{-4}$, $5.115\times10^{-4}$ | audit result (primary source re-read on the page) |
| Theorem CK′ | proved here (§2), a generalisation of Baouendi–Goulaouic's proof |
| $T = 8.683\times10^{-4}$ on $U$ (tube), $T = 1.545\times10^{-3}$ on $U_w$ (waist), every $d \in I$ | computer-assisted (Arb), candidate |
| evaluator: multi-certificate tails, extra derivative bounds | software; tests |
| $\gamma_2$ algebraic in Ric; $\gamma_3$ covariant law (5 terms) for conformally flat data | $\gamma_2$ derived; $\gamma_3$ fitted, held-out and exactly checked |
| $\gamma$ not pointwise algebraic in Ric at $t^6$ | certified (ball minor $\ne 0$) |
| $E(t) = E_0 - t^2(\tfrac32E_0^2 - \tfrac12\mathrm{tr}E_0^2\,I) + O(t^4)$; eigenframe rotates at $t^4$; speciality drifts at $t^2$ | derived and exact; certified exclusions |
| category-A metric | **not found** |
| merger at an explicit $d$ | **not certified** |

Conventions are unchanged (note 00, note 15):

- signature $(-,+,+,+)$, $G = c = M_{\rm ADM} = 1$;
- geodesic slicing, $g = -dt^2 + \gamma_{ij}dx^idx^j$;
- $K_{ij} = \tfrac12\partial_t\gamma_{ij}$;
- Brill–Lindquist (BL) data $\gamma_0 = \psi^4\delta$, $\psi = 1 + \tfrac1{4r_+} + \tfrac1{4r_-}$, punctures at $(0, 0, \pm d/2)$.

## 1 · Audit: what was wrong and what was repaired

**The abstract CK constant (consequential).**

- Session 11 reconstructed the time condition of Baouendi–Goulaouic (BG, *Journées EDP* 1977) from
  a garbled text extraction, as $a < \min(1/(16C), R/(8M))$ (normalised $s_0 = 1$).
- Read on the page, their condition (12) is

  ```math
  a < \inf\Big(\frac T2,\ \frac{R}{16RC + 8M}\Big),
  ```

  with the fixed point taken in the closure of a ball of $E_{2a}$. Their Lemma 3 needs one argument
  in $E_{2a}$; the extraction had lost that hypothesis.
- With the same certified $C$ and $K_0$, the times become
  - $2.885\times10^{-4}$ for note 15, instead of $2.988\times10^{-4}$;
  - $5.115\times10^{-4}$ for note 16, instead of $5.557\times10^{-4}$.
- The bound $\sup\lvert u\rvert_w \le 4aM$ survives, because $8aC < \tfrac12$ under (12).
- **Every downstream number** that used the old times is superseded: the evaluator's tails, the
  demonstration, and the heuristic $\delta_1$. That includes the "baseline $T_0 = 5.556\times10^{-4}$"
  this iteration was asked to improve on. Corrected figures appear below.

**Silent NaNs in ball powers (repaired).**

- python-flint returns `nan` for `x ** 2` when the midpoint of `x` is exactly 0, while `x * x` is
  fine.
- `frob_upper` squared absolute values that can have zero midpoint (purely real or purely
  imaginary entries). Python's `max` silently drops a `nan`.
- So a cell with a `nan` bound *could* have been skipped. All such powers are now products, and any
  non-finite cell bound is a failed cell (`ok: false`), never skipped.
- All certificates were recomputed under the guard. The tube result at the session-12 parameters
  is unchanged ($T_{\rm bg} = 5.103\times10^{-4}$ at cell 0.01), so no skipped cell had mattered
  there.

**Checked and unchanged.**

- The tail-bound summation indices (re-derived term by term, notes/16 §3).
- The relation between jet weight and derivative order: values through $g_N$, derivatives through
  $g_{N-1}$, Ricci exact to weight $W - 2$.
- Inverse-metric enclosures (`arb_mat.inv`), domain membership for points and boxes.
- Compatibility of the coefficient engine with the certified first-order system (re-run:
  $10^{-28}$).

**Smaller fixes.**

- The demonstration formed $t$ at 53-bit default precision (harmless balls of radius about
  $10^{-20}$); now at working precision.
- The $s$-grid for spatial Cauchy estimates could be empty for $q > 63/64$; it now always contains
  a valid point.

**Critical re-reading of note 16's "ceiling" and "recentering" claims.**

| Claim (note 16 §2) | What it actually is |
| --- | --- |
| "$C \ge 2\sqrt3\sup\lvert\psi\rvert^{-2}$, so $T \lesssim 10^{-3}$" | A theorem about one bound formula (block sup norms, Cauchy estimate, zero-order terms charged as $s_0\times$Lipschitz) at $s_0 = 0.014$. It is **not** a ceiling: the time also depends on $s_0$, the theorem's constants and the zero-order bookkeeping, all changed in §2. |
| "larger $s_0$ is worse" | A numerical parameter search with one covering. It depends on box overestimation, reduced in §3. |
| "recentering cannot beat $s_0/(16C)$" | Conditional on $C_1 \ge C$ after the restart. Even so, a continuous restart gains at most $\int_0^{s_0}ds/(16C(s))$, which is close to $s_0/(16C)$ here. Not a universal obstruction. |
| "complex near-zeros of $\psi$ cap $s_0$" (session-12 draft) | **wrong**: $\inf\lvert\psi\rvert$ only falls from 2.02 to 1.75 between $s_0 = 0.014$ and $0.08$. The cap comes from the punctures' complex singular cones $q_\pm = 0$, reached by the inflated cells and the corners of the imaginary boxes. |

## 2 · Theorem CK′ (proved here)

Let $(X_s)_{0<s\le s_0}$ be a decreasing scale of Banach spaces. Each $X_s$ consists of bounded
holomorphic maps on $\Omega_s$, normed by a maximum over blocks $A$. Let $F$ be autonomous and
holomorphic from the $R$-ball of $X_s$ into $X_{s'}$ for $s' < s$. Assume, for $0 < s' < s \le s_0$,
$x, y$ in the $R$-ball of $X_s$, and every block $A$:

```math
\text{(H1}''\text{)}\quad \lVert F_A(x) - F_A(y)\rVert_{s'} \le \frac{C_A\lVert x - y\rVert_s}{s - s'} + L_A\lVert x - y\rVert_{s'},
\qquad \text{(H2)}\quad \lVert F(0)\rVert_s \le M .
```

$C_A$ measures the derivative loss. $L_A$ is a zero-order Lipschitz constant at the *output* level.

> **Theorem CK′.** Let $0 < \beta < 1$ and $k = \max_A\big(4C_A/(\beta s_0) + L_A/(1 - \beta)\big)$. If
>
> ```math
> 2ak < 1 \qquad\text{and}\qquad 3aM < R\,(1 - \beta)(1 - 2ak),
> ```
>
> then $u' = F(u)$, $u(0) = 0$ has a solution $u(t) \in X_s$, holomorphic in $t$ for
> $\lvert t\rvert < a(1 - s/s_0)$, with
>
> ```math
> \sup\lVert u(t)\rVert_s \le \frac{aM}{(1 - \beta)(1 - ak)} \;(< R/3).
> ```
>
> So $a_{\max}(\beta) = (1 - \beta)/\big(3M + 2(1 - \beta)k\big)$ with $R = 1$, and every
> $a < a_{\max}$ works.

*Proof.* Normalise $s_0 = 1$ (replace $C_A$ by $C_A/s_0$). Let $E_b$ be the holomorphic maps $v$ with
$v(t) \in X_s$ for $\lvert t\rvert < b(1 - s)$, normed by

```math
|||v|||_b = \sup_{s,\ \lvert t\rvert < b(1-s)} (1 - s)\Big(1 - \frac{\lvert t\rvert}{b(1 - s)}\Big)^\beta\lVert v(t)\rVert_s .
```

This is BG's (7) with exponent $\beta$ in place of $\tfrac12$. It is a Banach space, and uniform
limits of holomorphic maps are holomorphic.

- *Lemma 1′.* $\lVert\int_0^tv\rVert_s \le b|||v|||_b/(1 - \beta)$, because
  $\int_0^{A}(1 - \sigma/A)^{-\beta}d\sigma = A/(1 - \beta)$ with $A = b(1 - s)$.
- *Lemma 2′.* With $s(\sigma) = \tfrac12(1 + s - \lvert\sigma\rvert/b)$, the weight at $(s(\sigma), \sigma)$
  equals $\tfrac{A + \sigma}{2b}\big(\tfrac{A - \sigma}{A + \sigma}\big)^\beta$, and $s(\sigma) - s = (A - \sigma)/(2b)$. So

  ```math
  \int_0^{\lvert t\rvert}\frac{\lVert v(\sigma)\rVert_{s(\sigma)}}{s(\sigma) - s}\,d\sigma \le 4b^2|||v|||_b\,A^{\beta-1}\int_0^{\lvert t\rvert}(A - \sigma)^{-1-\beta}d\sigma
  \le \frac{4b}\beta\,\frac{|||v|||_b}{(1 - s)(1 - \lvert t\rvert/A)^\beta}.
  ```

  For $\beta = \tfrac12$ this is BG's Lemma 2.
- *Lemma 3′* (BG's Lemma 3, with (H1″)). Let $u \in E_a$ and $v \in E_{2a}$ with
  $a|||u|||_a/(1-\beta) + 2a|||v|||_{2a}/(1-\beta) < R$. Then for $\lvert t\rvert < a(1 - s)$ and every
  block,

  ```math
  \lVert F_A(\textstyle\int_0^tu) - F_A(\int_0^tv)\rVert_s \le C_A\int_0^t\frac{\lVert u - v\rVert_{s(\sigma)}}{s(\sigma) - s}d\sigma + L_A\int_0^t\lVert u - v\rVert_sd\sigma .
  ```

  The proof is BG's telescoping (9) over $w_j = \int_0^{t_j}u + \int_{t_j}^tv$. Membership of $w_j$
  in the $R$-ball of $X_{s_j}$ follows from Lemma 1′; the second argument's being in $E_{2a}$ is
  what makes $\int_{t_j}^tv$ defined at level $s_j$. Each step uses (H1″) between the levels
  $s_j$ and $s$. The zero-order pieces add up to $L_A\int_0^t\lVert u - v\rVert_s$.

Write $G(v)(t) = F(\int_0^tv)$. The lemmas give, row by row and then taking the maximum:

- (10′) for $u \in E_b$ with $b|||u|||_b/(1-\beta) < R$: $|||G(u)|||_b \le bk|||u|||_b + M$;
- (11′) for $u \in E_a$, $v \in E_{2a}$ under the Lemma 3′ condition:
  $|||G(u) - G(v)|||_a \le ak|||u - v|||_a$.

Now fix the ball. Let $\rho = M/(1 - 2ak)$ and $\mathfrak B = \{v \in E_{2a} : |||v|||_{2a} \le \rho\}$.

- By (10′) with $b = 2a$ (allowed because $2a\rho/(1-\beta) < 3a\rho/(1-\beta) < R$),
  $G(\mathfrak B) \subset \mathfrak B$.
- On $\mathfrak B$ we have $|||\cdot|||_a \le |||\cdot|||_{2a}$, and the Lemma 3′ condition holds since
  $3a\rho/(1-\beta) < R$. So (11′) makes $G$ a contraction in $|||\cdot|||_a$, with factor $ak < \tfrac12$.
- $G$ extends to the closure $\bar{\mathfrak B}$ in $E_a$, mapping it into itself. For $u \in \bar{\mathfrak B}$,
  approximate by $u_n \in \mathfrak B$ and use (11′) with the pair $(u, u_n)$. The contraction property
  passes to the closure the same way.
- Banach's theorem gives a fixed point $v^*$. Then $u = \int_0^tv^*$ solves the problem.
- (10′) at $b = a$ gives $|||v^*|||_a \le M/(1 - ak)$, and Lemma 1′ gives the stated bound.

Holomorphy in $t$ is BG's Remark: take complex $t$ and integrate along rays. ∎

Two remarks on the hypotheses and the comparison with BG:

- With $\beta = \tfrac12$ and $L = 0$, the conditions read $16aC < 1$ and $6aM < R(1 - 16aC)$, that is
  $a < R/(16RC + 6M)$. That is slightly better than BG's printed $8M$, because the ball's radius is
  optimised.
- **Checking (H1″) for the BL system.** Note 15's bounds already split into the two kinds.
  - *Loss* terms come from Cauchy estimates: $2\nabla u_K$ (constant $2\sqrt3$); $A(G, \nabla\Delta g)$
    (the exact norms of note 16 §2); and $A(\Delta G, \nabla u_g)$ (constant $\tfrac12\lVert T\rVert\Gamma^2\sqrt3R_g$).
  - *Zero-order* terms are pointwise Lipschitz bounds at each $w \in \Omega_{s'}$: $2u_K$;
    $B(G, g)$ and $Q(G, K)$ through their Jacobians; and $A(\Delta G, \nabla g_0)$ (constant
    $\tfrac12\lVert T\rVert\Gamma^2G_1$).
  - The weighted row constants are $C_g = (R_K/R_g)2\sqrt3$,
    $C_K = (R_\gamma/R_K)c^{\rm loss}_{K\gamma} + (R_g/R_K)c^{\rm loss}_{Kg}$,
    $L_\gamma = 2R_K/R_\gamma$, and
    $L_K = (R_\gamma/R_K)c^0_{K\gamma} + (R_g/R_K)J_{Bg} + J_{QK}$. The constant
    $M = \sup\lvert\mathrm{Ric}(\gamma_0)\rvert_F/R_K$.

## 3 · The extension: certificates

**Covering improvements** (`ck_bound.certify(..., im_split=m)`).

- *Symmetry.* Complex conjugation and the reflection $y \to -y$ preserve every bound: they are
  norms of tensors built from real-analytic, reflection-symmetric data. So the imaginary parts can
  be restricted to $\operatorname{Im}x, \operatorname{Im}y \ge 0$, a quarter of the previous cube.
- *Sub-boxes.* The quadrant $[0, s_0]^2 \times [-s_0, s_0]$ is split into $m \times m \times 2m$ boxes.
  Boxes with no point of norm $< s_0$ are dropped.
- *Parallelism.* The cells run in parallel (4 processes). Any non-finite bound fails the run.

**Parameters** (searches recorded in the session-13 log; any choice is valid, the search only
selects).

- *Whole tube.* $s_0 = 0.017$, $(R_\gamma, R_g, R_K) = (0.25, 22, 10)$, $m = 2$.
- *Waist region.* $U_w$ is the profile boxes with $z_{\rm hi} \le 0.15$, inflated by $0.02$, rotated
  and reflected: the part of the tube around the symmetry plane. Its parameters are $s_0 = 0.03$,
  $(0.2, 18, 7)$, $m = 3$.
- Larger $s_0$ fails for the whole tube by $s_0 \approx 0.03$. Cells inflated by $w + s_0$ approach
  the singular cones of the punctures, and $K_0$ grows from 200 to 1,600 to 9,760 at
  $s_0 = 0.014, 0.03, 0.04$.

**Results** (`python -m research.calc.ck_bound {classic|sharp|ext|waist}`, cell side 0.005,
128-bit balls, outward rounded).

| mode | region | $s_0$ | weights | imag. covering | cells | $C$ | $K_0$ | BG (12) time | CK′ time | $B$ at that time | wall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `classic` | tube | 0.012 | (1.0, 20.0, 10.0) | cube | 2815 | 2.5100 | 148.77 | 2.885e-4 | **4.918e-4** (β = 0.7) | 0.0454 | 86 s |
| `sharp` | tube | 0.014 | (0.25, 22.0, 10.0) | cube | 3032 | 1.5746 | 194.46 | 5.114e-4 | **7.053e-4** (β = 0.8) | 0.1137 | 87 s |
| `ext` | tube | 0.017 | (0.25, 22.0, 10.0) | 2 | 3249 | 1.7034 | 135.39 | 5.842e-4 | **8.683e-4** (β = 0.75) | 0.0824 | 379 s |
| `waist` | waist | 0.03 | (0.2, 18.0, 7.0) | 3 | 1137 | 2.3567 | 44.58 | 7.646e-4 | **15.45e-4** (β = 0.6) | 0.0458 | 435 s |

All runs use 4 processes and peak at 44–47 MB per process. Times are rounded down; the stored values are the certified floats. The `classic` and `sharp` rows are the session 11–12 constants with the corrected time; the old claims were $2.988\times10^{-4}$ and $5.557\times10^{-4}$.

So, **on the whole tube** the certified time rises from the corrected baseline
$5.115\times10^{-4}$ to $8.683\times10^{-4}$ (×1.70). **On the waist region** it rises to $1.545\times10^{-3}$
(×3.02), passing the stretch target $2T_0$.

These are certified extensions of one and the same analytic solution: the CK solutions on $U$ and
on $U_w$ agree on the common domain, by uniqueness of analytic solutions. A longer interval
**does not certify a merger at any explicit $d$** (§8).

*The energy-method alternative, inspected and not executed.* An a posteriori stability estimate
for a hyperbolic formulation would avoid the loss of a derivative per Cauchy estimate. It is not
available with the present representation, for four reasons:

1. *Formulation.* The fixed-lapse, zero-shift ADM system the series solves is only weakly
   hyperbolic. An energy estimate needs a different reduction. Options are the generalised
   harmonic gauge (which needs a coordinate change from the geodesic chart, with a Jacobian to
   control), or the Einstein–Bianchi system in Gauss gauge, which is symmetric hyperbolic in this
   same chart.
2. *Norms.* Pointwise bounds on $\partial^2g$ need $H^s$ with $s > 7/2$ on shrinking
   domain-of-dependence regions, with explicit Sobolev and commutator constants.
3. *Residual.* The defect of the truncated series must be bounded in that norm over the whole
   region, not at points.
4. *Closure.* A continuation (bootstrap) argument is needed.

None of these is in place, and a sketch is not an extension. Here the analytic route is extended
instead, because its hypotheses are fully verified.

## 4 · The evaluator: derivative bounds before and after

Changes in `research/calc/bl_metric.py`:

- *Several certificates.* Every certificate whose region contains the point (or box) and whose time
  covers $t$ contributes; each tail bound takes the minimum over them and over every $(a, B)$ pair
  of the CK′ schedules.
- *A second $\partial_t\partial_x\gamma$ bound.* A Cauchy estimate *in $t$* on the $g$-block,
  intersected with the Cauchy-in-$x$ bound on the $K$-block.
- *A valid $s$-grid* at every $q < 1$.
- *Cached spatial jets.* They do not depend on $t$, so a time scan reuses them.
- *The field-equation enclosure of $\partial_t^2\gamma$* (note 16) is kept.

**Comparison** (`python -m research.calc.bl_metric compare`, record `bl_metric_compare.json`).
$d = 0.7662$. Times are fractions of the session-12 $T_0 = 5.557\times10^{-4}$. Entries are the
largest radius over the components of $g$ (metric), of $\partial g$ (first) and of $\partial^2g$
(second).

| point | $t/T_0$ | quantity | before, claimed (invalid) | before, corrected | after, order 5 | after, order 7 |
| --- | --- | --- | --- | --- | --- | --- |
| generic $(0.17, 0.11, 0.5)$ | 0.1 | $g$ | 1.09e-14 | 2.7e-14 | 9.8e-17 | 1.6e-21 |
| generic $(0.17, 0.11, 0.5)$ | 0.1 | $\partial g$ | 9.61e-11 | 2.0e-10 | 2.1e-12 | 3.5e-17 |
| generic $(0.17, 0.11, 0.5)$ | 0.1 | $\partial^2 g$ | 2.88e-06 | 5.6e-06 | 8.1e-08 | 1.9e-12 |
| generic $(0.17, 0.11, 0.5)$ | 0.25 | $g$ | 6.87e-10 | 1.7e-09 | 5.9e-12 | 3.9e-15 |
| generic $(0.17, 0.11, 0.5)$ | 0.25 | $\partial g$ | 9.67e-07 | 2.1e-06 | 2.0e-08 | 1.3e-11 |
| generic $(0.17, 0.11, 0.5)$ | 0.25 | $\partial^2 g$ | 1.18e-02 | 2.3e-02 | 3.2e-04 | 2.9e-07 |
| generic $(0.17, 0.11, 0.5)$ | 0.5 | $g$ | 3.52e-06 | 9.3e-06 | 2.6e-08 | 2.8e-10 |
| generic $(0.17, 0.11, 0.5)$ | 0.5 | $\partial g$ | 1.24e-03 | 2.8e-03 | 2.3e-05 | 2.4e-07 |
| generic $(0.17, 0.11, 0.5)$ | 0.5 | $\partial^2 g$ | 8.25e+00 | 1.7e+01 | 1.8e-01 | 2.6e-03 |
| generic $(0.17, 0.11, 0.5)$ | 1.0 | $g$ | — | refused | 1.6e-04 | 2.8e-05 |
| generic $(0.17, 0.11, 0.5)$ | 1.0 | $\partial g$ | — | refused | 3.5e-02 | 5.9e-03 |
| generic $(0.17, 0.11, 0.5)$ | 1.0 | $\partial^2 g$ | — | refused | 1.6e+02 | 3.5e+01 |
| axis $(0, 0, 0.62)$ | 0.1 | $g$ | 1.09e-14 | 2.7e-14 | 9.8e-17 | 1.6e-21 |
| axis $(0, 0, 0.62)$ | 0.1 | $\partial g$ | 9.61e-11 | 2.0e-10 | 2.1e-12 | 3.5e-17 |
| axis $(0, 0, 0.62)$ | 0.1 | $\partial^2 g$ | 2.88e-06 | 5.6e-06 | 8.1e-08 | 1.9e-12 |
| axis $(0, 0, 0.62)$ | 0.25 | $g$ | 6.87e-10 | 1.7e-09 | 5.9e-12 | 3.9e-15 |
| axis $(0, 0, 0.62)$ | 0.25 | $\partial g$ | 9.67e-07 | 2.1e-06 | 2.0e-08 | 1.3e-11 |
| axis $(0, 0, 0.62)$ | 0.25 | $\partial^2 g$ | 1.18e-02 | 2.3e-02 | 3.2e-04 | 2.9e-07 |
| axis $(0, 0, 0.62)$ | 0.5 | $g$ | 3.52e-06 | 9.3e-06 | 2.6e-08 | 2.8e-10 |
| axis $(0, 0, 0.62)$ | 0.5 | $\partial g$ | 1.24e-03 | 2.8e-03 | 2.3e-05 | 2.4e-07 |
| axis $(0, 0, 0.62)$ | 0.5 | $\partial^2 g$ | 8.25e+00 | 1.7e+01 | 1.8e-01 | 2.6e-03 |
| axis $(0, 0, 0.62)$ | 1.0 | $g$ | — | refused | 1.6e-04 | 2.8e-05 |
| axis $(0, 0, 0.62)$ | 1.0 | $\partial g$ | — | refused | 3.5e-02 | 5.9e-03 |
| axis $(0, 0, 0.62)$ | 1.0 | $\partial^2 g$ | — | refused | 1.6e+02 | 3.5e+01 |
| waist $(0.26, 0, 0.05)$ | 0.1 | $g$ | — | 2.7e-14 | 4.3e-20 | 7.2e-26 |
| waist $(0.26, 0, 0.05)$ | 0.1 | $\partial g$ | — | 2.0e-10 | 3.0e-15 | 5.0e-21 |
| waist $(0.26, 0, 0.05)$ | 0.1 | $\partial^2 g$ | — | 5.6e-06 | 9.9e-11 | 2.3e-16 |
| waist $(0.26, 0, 0.05)$ | 0.25 | $g$ | — | 1.7e-09 | 2.6e-15 | 1.7e-19 |
| waist $(0.26, 0, 0.05)$ | 0.25 | $\partial g$ | — | 2.1e-06 | 2.9e-11 | 1.9e-15 |
| waist $(0.26, 0, 0.05)$ | 0.25 | $\partial^2 g$ | — | 2.3e-02 | 3.8e-07 | 3.5e-11 |
| waist $(0.26, 0, 0.05)$ | 0.5 | $g$ | — | 9.3e-06 | 1.1e-11 | 1.1e-14 |
| waist $(0.26, 0, 0.05)$ | 0.5 | $\partial g$ | — | 2.8e-03 | 3.0e-08 | 3.1e-11 |
| waist $(0.26, 0, 0.05)$ | 0.5 | $\partial^2 g$ | — | 1.7e+01 | 2.0e-04 | 2.9e-07 |
| waist $(0.26, 0, 0.05)$ | 1.0 | $g$ | — | refused | 4.9e-08 | 8.2e-10 |
| waist $(0.26, 0, 0.05)$ | 1.0 | $\partial g$ | — | refused | 3.4e-05 | 5.7e-07 |
| waist $(0.26, 0, 0.05)$ | 1.0 | $\partial^2 g$ | — | refused | 1.2e-01 | 2.7e-03 |
| waist $(0.26, 0, 0.05)$ | 2.0 | $g$ | — | refused | 3.6e-04 | 9.7e-05 |
| waist $(0.26, 0, 0.05)$ | 2.0 | $\partial g$ | — | refused | 6.3e-02 | 1.7e-02 |
| waist $(0.26, 0, 0.05)$ | 2.0 | $\partial^2 g$ | — | refused | 1.3e+02 | 4.6e+01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.1 | $g$ | 4.21e-02 | 4.2e-02 | 4.2e-02 | 4.2e-02 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.1 | $\partial g$ | 4.58e-01 | 4.6e-01 | 4.6e-01 | 4.6e-01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.1 | $\partial^2 g$ | 1.37e+01 | 1.4e+01 | 1.4e+01 | 1.4e+01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.25 | $g$ | 4.21e-02 | 4.2e-02 | 4.2e-02 | 4.2e-02 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.25 | $\partial g$ | 4.58e-01 | 4.6e-01 | 4.6e-01 | 4.6e-01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.25 | $\partial^2 g$ | 1.37e+01 | 1.4e+01 | 1.4e+01 | 1.4e+01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.5 | $g$ | 4.21e-02 | 4.2e-02 | 4.2e-02 | 4.2e-02 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.5 | $\partial g$ | 4.59e-01 | 4.6e-01 | 4.6e-01 | 4.6e-01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 0.5 | $\partial^2 g$ | 1.86e+01 | 2.5e+01 | 1.4e+01 | 1.4e+01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 1.0 | $g$ | — | refused | 4.2e-02 | 4.2e-02 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 1.0 | $\partial g$ | — | refused | 4.9e-01 | 4.6e-01 |
| box $(0.17, 0.11, 0.5)\pm10^{-4}$ | 1.0 | $\partial^2 g$ | — | refused | 1.6e+02 | 4.0e+01 |

Over the $\pm10^{-4}$ box at $t = T_0/4$, the enclosure radius versus the half-range over 27 sample points (`bl_metric_boxvariation.json`) is: $g_{xx}$ 0.04211 vs 0.04207 (×1.001); $\partial_xg_{xx}$ 0.458 vs 0.454 (×1.007); $\partial_x^2g_{xx}$ 13.7 vs 5.2 (×2.6, no mean-value form at that order). So the box radii are genuine variation, apart from a factor 2.6 on second derivatives.

**Reading the table.**

- The time-truncation bounds dominate every radius. Coefficient rounding is about $10^{-34}$, and
  there is no spatial approximation.
- At a box, "metric" includes the genuine variation of $\gamma$ over the box
  (about $\lvert\nabla\gamma\rvert \times$ half-width). The interval overestimation beyond that is
  small, thanks to the mean-value form (note 16).
- "before (claimed)" is session 12's output with the invalid constant, kept for reference only. It
  is in `bl_metric_baseline_session12.json`.
- **Second derivatives at $T_0/2$:** generic and axis points $1.75\times10^{1}$ (corrected baseline, order 5) → 1.81e-01 (order 5) → 2.59e-03 (order 7), a factor $6737$; the session-12 claim was 8.25. At the waist point, 2.9e-07 (order 7). Metric: 9.3e-06 → 2.8e-10; first derivatives: 2.8e-03 → 2.4e-07.
- **Cost.** the jets do not depend on $t$. The first evaluation at a point costs about 7 s at order 5 and 80–100 s at order 7; a box costs twice that (centre plus box). Further times at the same point reuse the cached jets and take about 0.02 s. Peak memory of the comparison run: 74 MB. Arithmetic precision is 128 bits throughout; raising it changes nothing, since the time tails dominate.
- The evaluator still refuses unsupported tolerances and out-of-domain certified requests (tests).

## 5 · Coefficient structure at $t^6$ (bounded closure test)

For conformally flat, time-symmetric vacuum data ($R = 0$ and vanishing Cotton tensor), the
gradient $\nabla_kR_{ij}$ is totally symmetric. Commuting derivatives and using the 3D
decomposition of Riemann into Ricci gives

```math
\Delta\mathrm{Ric} = 3\,\mathrm{Ric}\cdot\mathrm{Ric} - \lvert\mathrm{Ric}\rvert^2\gamma_0,
\qquad
\gamma_2 = \tfrac1{12}\big(7\,\mathrm{Ric}\cdot\mathrm{Ric} - \lvert\mathrm{Ric}\rvert^2\gamma_0\big).
```

So through $t^4$ the BL development is *pointwise algebraic in the tidal tensor*. The general
formula of note 16 reduces to this one.

**Candidate class (fixed before fitting).** $\gamma_3$ is a linear combination, with universal
rational coefficients, of the 14 natural symmetric 2-tensors of the right order:

- $\Delta^2\mathrm{Ric}$;
- eleven quadratic terms in $\mathrm{Ric}$, $\nabla\mathrm{Ric}$ and $\nabla\nabla\mathrm{Ric}$;
- $\lvert\mathrm{Ric}\rvert^2\mathrm{Ric}$ and $\mathrm{tr}(\mathrm{Ric}^3)\gamma_0$.

$\mathrm{Ric}^3$ is removed by Cayley–Hamilton for traceless $\mathrm{Ric}$.

**Results** (`bl_structure.curvature_values`, `t6_basis`, `cf_identities`).

- *Rank.* At 300 bits the 14 tensors span only a 5-dimensional space on BL data: the singular
  values drop from 0.39 to $10^{-90}$. The nine relations come out as simple rationals, e.g.
  $\Delta^2\mathrm{Ric} = 3\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} + 6\,\nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric} - 2\lvert\nabla\mathrm{Ric}\rvert^2\gamma_0$.
  They follow from conformal flatness.
- *Fit.* Five generic points, $d \in \{0.70, 0.75, 0.7662, 0.80\}$, residual $7\times10^{-90}$:

  ```math
  \gamma_3 = -\tfrac1{15}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\,\gamma_0
  + \tfrac1{60}\nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric} - \tfrac1{180}\lvert\nabla\mathrm{Ric}\rvert^2\gamma_0
  - \tfrac1{180}R^{kl}\nabla_k\nabla_lR_{ij}.
  ```

- *Falsification tests passed.* The formula holds:
  - on two held-out generic points (balls, $\sim10^{-69}$, plus two more at $10^{-89}$);
  - **exactly** in rational arithmetic, at a rational-distance point on the plane and one on the
    axis (BL, $d = 3831/5000$), and at a Schwarzschild point.
- *Status.* An identity for conformally flat time-symmetric vacuum data, **conjectured** in general.
  It is checked exactly at three points and on held-out balls, but not derived by hand.
- *Closure test.* The class "$\gamma(t, x)$ is a pointwise algebraic function of $(\gamma_0, \mathrm{Ric})$"
  holds through $t^4$ and is **excluded at $t^6$**. At two generic points a $3\times3$ minor of
  $[\lvert\mathrm{Ric}\rvert^2\mathrm{Ric},\ \mathrm{tr}(\mathrm{Ric}^3)\gamma_0,\ \gamma_3]$ has a ball
  excluding 0 ($1.2284\ldots$ and $-0.0093153\ldots$). Gradients of the tidal field enter at $t^6$,
  through two quadratic $\nabla\mathrm{Ric}$ terms and one $R^{kl}\nabla_k\nabla_l$ term.
- *Recurrence.* No recurrence in the order was tested beyond this: two exactly known orders and one
  fitted order are too few to falsify a recurrence class honestly.

## 6 · Geometric variables (spacetime invariants along the comoving observers)

**Variables.**

- $E_{ij}$ and $B_{ij}$, the electric and magnetic Weyl tensors relative to the comoving geodesic
  observers $n = \partial_t$. These depend on the observer field and the slicing.
- The complex invariants $I = \tfrac12\mathrm{tr}Q^2$ and $J = \tfrac16\mathrm{tr}Q^3$, with $Q = E + iB$.
  These are spacetime invariants.
- The speciality index $S = 27J^2/I^3$, invariant, with $S = 1$ for Petrov types D and II.
- The eigenframe of $E$, relative to the observers.

**Results** (`bl_structure.electric_series`, `weyl_speciality`, `eigenframe_tests`; $B = O(t^3)$
does not enter below $t^3$).

- *A closed law at $t^2$* (derived from §5, exact at the rational plane point, $10^{-71}$ at generic
  points):

  ```math
  E^i{}_j(t) = E_0 - t^2\big(\tfrac32E_0^2 - \tfrac12\mathrm{tr}(E_0^2)\,I\big) + O(t^4).
  ```

  So at this order the tidal tensor evolves by an isotropic matrix law, as in "silent" dust models.
- *The speciality drifts.* $dS/d(t^2) \ne 0$ at generic points: $0.0039223\ldots$,
  $0.346039\ldots$, and $0.412815\ldots$ on the plane, all balls excluding 0. Schwarzschild has
  $S \equiv 1$, and the BL axis has $S = 1$ by symmetry. This excludes the class "tidal eigenvalue
  *ratios* frozen along each worldline" (isospectral up to scaling), which contains every
  type-D-like free fall.
- *The eigenframe rotates at $t^4$.* $|[E_0, E_2]|^2 = 9.81\times10^{-6}$ and $1.18\times10^{-6}$
  at two generic points (balls excluding 0); it is 0 for Schwarzschild. Any closed, pointwise,
  isotropic evolution law $E' = f(E)$ along the observers keeps $[E_0, E_n] = 0$. **So every such
  law is excluded at $t^4$**: gradients of the tidal field enter the dynamics there. This matches
  §5's $t^6$ result for $\gamma$ ($E$ involves two more derivatives).
- *Invariants as coordinates.* $\operatorname{Re}I$ and $\operatorname{Re}J$ are even in $t$;
  $\operatorname{Im}I$ and $\operatorname{Im}J$ are $O(t^3)$. So $(t, x) \mapsto$ (invariants) has zero
  $t$-derivative on $\Sigma_0$ and is not locally invertible there. On the axis, $S \equiv 1$ as well.
  Invariant-based coordinates would have to exclude a neighbourhood of the time-symmetric slice
  and of the axis.
- *Reconstruction.* Even a closed scalar or tidal law would not give the metric: one would still
  need the 1+3 constraint and propagation equations (for the shear and the expansion, with their
  spatial derivatives). Here the closed laws themselves fail first.

**Outcome.** These variables simplify the first non-trivial order (an algebraic law). They expose
no finite closed system: tidal-field gradients enter at $t^4$ for $E$ and at $t^6$ for $\gamma$.
**No category-A candidate emerged.** Only the classes stated are excluded.

## 7 · What the convergence time means: the single-hole control

Take Schwarzschild in geodesic slicing from rest at areal radius $R_0$. Then $R = R_0F(\tau)$ with
$F = \cos^2(\eta/2)$, $\tau = (\eta + \sin\eta)/\sqrt8$ and $\tau = t\sqrt{M/R_0^3}$
(`bl_structure.schwarzschild_crunch`, tested).

- $\gamma_{rr} = \psi^4H^2$ with $H = F - \tfrac32\tau F' \ge 1$ on $(0, \pi)$. The comoving
  congruence has **no caustic** before the end.
- $d\tau/d\eta = 0$ only at $\eta = \pi$ (mod $2\pi$). So the Taylor radius in $t$ at that comoving
  point is exactly $t_c = \pi\sqrt{R_0^3/8M}$.
- At $t_c$ the observer reaches $R = 0$, and the Kretschmann scalar $48M^2/R^6$ diverges ($10^{45}$
  at $1 - 10^{-4}$ of the way).

So for one hole the series stops converging where the slicing meets the **$r = 0$ curvature
singularity**. Notes 14 and 16 called this a gauge breakdown; that is now corrected. Whether a
slicing reaches the singularity, and when, is gauge-dependent: maximal slicing avoids it. The event
itself is physical.

For BL, the numerical radius (ratio and Padé estimates $\approx 2.1$–$2.8$) is not interpreted.
Padé poles cannot distinguish a curvature singularity, a caustic and a complex singularity, and no
certified bound exists there.

## 8 · Physical scope

- **What is certified.** A **local** development of BL data: $U \times (-T, T)$ and
  $U_w \times (-T_w, T_w)$, for every $d \in I$. It embeds in the maximal development of
  asymptotically flat data, on the domain of dependence (note 15 §6; coordinate light speed
  $\le 0.25$). The patch itself says nothing about asymptotics.
- **No merger at an explicit $d$.** At $d = 0.7662$ the linear onset law predicts $t_*$ at about
  $0.01$ of the certified time, but the range of validity $\delta_1$ is still not explicit (note 15
  §7.3). So $d = 0.7662$ is not a certified merger example. A longer certified time does not
  change this.
- **Not addressed:** the first common apparent horizon anywhere on a slice, event horizons, and a
  full inspiral–merger–ringdown spacetime.
- **Equal-mass head-on BL only.** Nothing here applies to spins, other mass ratios or orbits.

## 9 · Status and next step

| Claim | Category |
| --- | --- |
| BG constant: notes 15–16 times too large; corrected $2.885\times10^{-4}$, $5.115\times10^{-4}$ | audit (page re-read) |
| Theorem CK′ | proved (§2) |
| $T = 8.683\times10^{-4}$ on $U$, $T = 1.545\times10^{-3}$ on $U_w$ | computer-assisted (rerun in session 13); candidate |
| derivative-error comparison | computed (§4) |
| $\Delta\mathrm{Ric}$ identity and algebraic $\gamma_2$ for conformally flat data | derived; exact checks |
| $\gamma_3$ five-term law | fitted; held-out balls and exact checks at 3 points; not derived |
| algebraic-in-Ric closure excluded at $t^6$; isotropic tidal laws excluded at $t^4$; speciality not frozen | certified (balls excluding 0), Schwarzschild controls |
| single-hole series radius = time to the curvature singularity | exact (classical cycloid) |
| category-A metric; merger at explicit $d$ | **no**; **no** |

**Next decisive calculation.** The Einstein–Bianchi system in Gauss gauge is symmetric hyperbolic
in the present chart, so no coordinate change is needed. Use the certified series at a time
$t_1 < T$ as an approximate solution, and carry out one validated energy-estimate step for that
system. That needs two things:

- a certified $H^4$ bound on the series' defect over a domain-of-dependence region;
- explicit commutator and Sobolev constants.

Unlike any refinement of CK, this route has no derivative-loss ceiling.

## Discussion prompts

- The constant error survived two sessions because the text extraction lost the hypothesis
  $v \in E_{2a}$. Which other hypotheses in a computer-assisted proof are read only from extracted
  text, and how should a project guard against that?
- Through $t^4$ the BL development is pointwise algebraic in the tidal tensor, and the tidal tensor
  obeys an isotropic law at $t^2$. Is there a reformulation, for example in terms of $(E, B)$ and
  the shear, in which the tidal gradients that enter at $t^4$ appear as a small, explicit source?
- For one hole, the series at a comoving point converges exactly up to the moment that observer
  hits $r = 0$. What would it take to decide whether the BL series' numerical radius of about 2.5
  marks a curvature singularity, a caustic or a complex singularity?
