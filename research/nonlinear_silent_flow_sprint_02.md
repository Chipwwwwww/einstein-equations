# Nonlinear silent-flow conjugation and a checkable convergence condition (sprint 02)

*Branch `claude/nonlinear-silent-flow`, from `c6e7aa8`. That commit contains the axisymmetric sprint
and curl-hierarchy sprint 01; `main` is still `90f023b`. Checker:
`research/calc/silent_flow_checks.py`.*

**Labels.**

- **P**: proved here (analytic derivation).
- **C**: conditional on stated hypotheses.
- **E**: computation run this sprint.
- **I**: inherited, read and not rerun.

**Scope.** Conformally flat, time-symmetric vacuum data in Gauss gauge. Everything is analytic except
§6 (one BL check) and §5's numbers (inherited certificate constants).

## 0 · Answer in one paragraph

Composing with the exact pointwise silent flow, $\Phi = F_t(V)$, removes the whole algebraic part $A$
from the equations. The transformed system is purely quasi-linear first order,
$\partial_tV = \varepsilon\sum_i\tilde B^i(t, V)\partial_iV$, with **no zero-order term**. Its source at
the initial state is $J_t^{-1}L(\Phi_{\rm sil})$, which is $O(t^2)$ and nonzero.

Its principal matrices are similar to the original ones ($\tilde B^i = J^{-1}B^iJ$). So the
characteristic structure is unchanged, and the loss constant can only grow by $\kappa(J)$.

The abstract CK theorem with $\varepsilon$ as a holomorphic parameter turns this into a full nonlinear
**conditional theorem** (§4):

- the sector series converges at $\varepsilon = 1$, with explicit tail bounds, for
  $t < T_{\rm conj} = \rho_0R(1 - s)/(16RC_1 + 8M_Q)$;
- the reconstruction goes through $\Phi = F_t(V)$.

$T_{\rm conj}$ has **the same CK scale** as before: analytic width over loss constant. The
conjugation's genuine contributions are:

- an exact removal of the zero-order constraints;
- an $O(t^2)$ source;
- a rigorous way to bound *all* sector convolutions at once, by holomorphy in $\varepsilon$.

It does not change the mechanism. The decisive constant is the loss bound $B_*$ of the E–H principal
part on complex cells; the complex-ODE quantities ($J$, $D^2F$) are not decisive.

## 1 · Scope corrections to the previous sprints

1. **Proposition 4 (sprint 01) controls an auxiliary linear hierarchy** (coefficients frozen at the
   silent background, curl only on the $(E, H)$ sectors). It does not control the full nonlinear
   hierarchy.

   The full linearisation at a silent background ($\bar H = 0$) is given in §2. It contains a
   *metric-variation loss term* $-D_\gamma(\mathrm{curl}\,\bar E)[Y_\gamma]$. The phrase "replace $C_L$
   by $C_L + C_N$" is withdrawn; §4 replaces it.
2. **The ×1.55 / ×1.80 ceilings (sprint 01, claim 126)** are an arithmetic upper bound obtained by
   fixing the CK′ norms, weights, analytic width, loss rows and sufficient-condition formula, and then
   idealising away all zero-order and source constraints. They are not:
   - the intrinsic maximal gain of the curl-degree method;
   - a bound on all resummations;
   - a convergence or physical limit of BL.

   Removing $A$ does not make the curl forcing vanish: §3, item 5.
3. **Complex domains.** Sprint 01 wrote $t_{\rm sil}(z) \le T$ for complex $z$. That is not defined and
   is withdrawn. It is replaced by a *common regularity domain*:
   - a complex neighbourhood $D_\rho$ on which, for all real $0 \le t \le T$, the background (silent
     solution), the inverse metric and the needed variations are holomorphic with uniform bounds;
   - with the principal branch of $\sqrt{\det\gamma}$ (§3.4).

   The claim that the admissible width "shrinks to 0 at a definite rate as $T \to$ the real collapse
   time" is a **conjecture**. Different collapse times at a few real points do not prove it.
4. **Axisymmetric sprint, scope.** KFA excluded the tested Killing-norm-gradient frame, and claim 100
   excludes alignment with the tidal frame. That does **not** exclude every other symmetry-adapted
   frame. The meridian orientation angle is a coordinate needed to *describe* $S$; it is not an
   independent physical degree of freedom.

## 2 · Full linearisation of the spatial operator at a silent background

Write $c(\gamma; A)_{ij} = \mathrm{sym}_{ij}\,\epsilon_i{}^{kl}D^\gamma_kA_{lj}$, with
$\epsilon_i{}^{kl} = \gamma_{im}\epsilon^{mkl}$ and $\epsilon^{mkl} = [mkl]/\sqrt{\det\gamma}$. Then
$L(\Phi) = (0, 0, c(\gamma; H), -c(\gamma; E))$. At $\bar\Phi$ with $\bar H = 0$:

```math
DL(\bar\Phi)[Y] = \bigl(0,\ 0,\ c(\bar\gamma; Y_H),\ -c(\bar\gamma; Y_E) - D_\gamma c(\bar\gamma; \bar E)[Y_\gamma]\bigr)
```

The term $D_\gamma c(\bar\gamma; \bar H)[Y_\gamma]$ is absent because $\bar H = 0$. The metric variation,
with $h = Y_\gamma$ (P), is

```math
D_\gamma c(\bar\gamma; \bar E)[h]_{ij} = \mathrm{sym}_{ij}\Bigl(
 \underbrace{h_{im}\bar\epsilon^{mkl}\bar D_k\bar E_{lj}}_{\text{lowering}}
 \underbrace{-\tfrac12(\bar\gamma^{pq}h_{pq})\,\bar\epsilon_i{}^{kl}\bar D_k\bar E_{lj}}_{\text{volume form}}
 \underbrace{-\,\bar\epsilon_i{}^{kl}\,\delta\Gamma^p{}_{kj}[h]\,\bar E_{lp}}_{\text{connection}}\Bigr),
\qquad \delta\Gamma^p{}_{kj} = \tfrac12\bar\gamma^{pq}\bigl(\bar D_kh_{qj} + \bar D_jh_{qk} - \bar D_qh_{kj}\bigr).
```

The variation $\delta\Gamma^p{}_{kl}$ in $D_kE_{lj}$ cancels against the antisymmetric
$\epsilon^{kl}$. The inverse metric enters through $\delta\Gamma$, linearly here. For the *nonlinear*
operator at $\gamma = \bar\gamma + h$ the full expansions are needed:

- $\gamma^{-1} = \sum_n(-\bar\gamma^{-1}h)^n\bar\gamma^{-1}$;
- $\det(\bar\gamma + h)^{-1/2} = \det\bar\gamma^{-1/2}\det(I + \bar\gamma^{-1}h)^{-1/2}$.

Both converge when $\lvert\bar\gamma^{-1}h\rvert \le r$ with $r < 2^{1/3} - 1$, which keeps
$\det(I + X)$ in the disc $\lvert w - 1\rvert < 1$ (principal branch).

**Two consequences.**

- *At $t = 0$* the lowering and volume terms vanish. $\bar\epsilon_i{}^{kl}\bar D_k\bar E_{lj} = 0$
  entirely: its symmetric part is Cotton–York $= 0$, and its antisymmetric part is proportional to
  $\epsilon_{ijm}(\mathrm{div}\,E)^m = 0$ by the vacuum constraint at time symmetry.
- The **connection term is a first-order (loss) term in $h$ with coefficient $\bar E \ne 0$.** It is
  the $H$-row ← $\partial\gamma$-column entry responsible for the Jordan block of the coordinate symbol
  (sprint 01 §2.4).

**Earliest order where it acts** (derived before checking).

- By the lag lemma, $\gamma - \gamma_{\rm sil} = O(\varepsilon^2t^6)$.
- So $\varepsilon D_\gamma c(\cdot; \bar E)[\gamma - \gamma_{\rm sil}] = O(\varepsilon^3t^6)$ in $\dot H$.
- Hence it first changes $H$ at $\varepsilon^3t^7$, $E$ at $\varepsilon^4t^8$ (through
  $\mathrm{curl}\,H$), and $\gamma$ at $t^{10}$.

Checked in §6 (E).

## 3 · Nonlinear silent-flow conjugation

**Setting.** $\partial_t\Phi = A(\Phi) + \varepsilon L(\Phi, \partial\Phi)$, with $\Phi = (\gamma, K, E, H)$:
24 complex functions, full symmetric matrices, no reduction. Let $F_t$ be the flow of the pointwise ODE
$\dot\Phi = A(\Phi)$, and $J_t(V) = D_VF_t(V)$ its state Jacobian (24 × 24).

**Proposition 1 (conjugation; P).** Assume:

- (a) $A$ is holomorphic on an open set $\mathcal U \subset \mathbb C^{24}$ where $\gamma$ is invertible;
- (b) for the states $V$ considered and $s \in [0, t]$, the flow $F_s(V) \in \mathcal U$ exists;
- (c) $V(t, x)$ is $C^1$ in $t$ and holomorphic in $x$.

Then $\Phi(t, x) := F_t(V(t, x))$ solves the $\varepsilon$-system iff

```math
\partial_tV = \varepsilon\,J_t(V)^{-1}\,L\bigl(F_t(V),\ J_t(V)\,\partial_xV\bigr),\qquad V(0) = \Phi(0).
```

*Proof.*

- $\partial_t\Phi = (\partial_tF_t)(V) + J_t(V)\partial_tV = A(F_t(V)) + J_t\partial_tV$.
- $\partial_x\Phi = J_t(V)\partial_xV$, because $F_t$ has no explicit $x$-dependence: the silent ODE
  is autonomous and pointwise.
- $J_t$ is invertible, since $\det J_t = \exp\int_0^t\mathrm{tr}\,DA(F_s(V))\,ds$ (Liouville). $\square$

**Principal matrices.** $L$ is exactly quasi-linear first order with no zero-order part:
$L(\Phi, \partial\Phi) = \sum_iB^i(\Phi)\partial_i\Phi$. Every term of $c(\gamma; A)$ contains exactly
one derivative, of $A$ or of $\gamma$ (Christoffels), and $(\det\gamma)^{-1/2}$ contains none. Hence

```math
\partial_tV = \varepsilon\,Q(t, V, \partial V),\qquad Q = \sum_i\tilde B^i(t, V)\,\partial_iV,\qquad
\tilde B^i(t, V) = J_t(V)^{-1}B^i(F_t(V))\,J_t(V).
```

**Answers to the six questions** (P).

1. **Absorbed:** all of $A$. That is the pointwise algebraic dynamics: the $K^2$, $\theta E$ and
   $\sigma E$ terms, and the silent collapse.
2. **Remaining:** all of $L$, including the metric-variation terms of §2.
3. **Symbol:** $\tilde B(\xi) = J^{-1}B(\xi; F_t(V))J$ is *similar* to the original symbol at the
   physical state. Eigenvalues and Jordan structure are identical, so the weak hyperbolicity of the
   coordinate form is **not** repaired.
4. **Norms:** $\lVert\tilde B^i\rVert \le \kappa(J)\lVert B^i\rVert$ with $\kappa = \lVert J\rVert\lVert J^{-1}\rVert \ge 1$. Lipschitz constants acquire
   $D_V\tilde B = -J^{-1}(DJ)J^{-1}BJ + J^{-1}(D_\Phi B\cdot J)J + J^{-1}B(DJ)$, which needs the
   state Hessian $D^2F_t$.
5. **Source at the initial state:** $Q(t, V_{\rm in}, \partial V_{\rm in}) = J_t^{-1}L(\Phi_{\rm sil}(t), \partial\Phi_{\rm sil}(t)) = J_t^{-1}(0, 0, 0, -\mathrm{curl}_{\gamma_{\rm sil}}E_{\rm sil})$.
   It is nonzero for $t > 0$ and of order $t^2$:
   - $\mathrm{curl}\,E_0 = 0$;
   - the $t^1$ term vanishes because the silent flow is even;
   - the $t^2$ term is $-3H_3$, whose curl is CGR (59).
6. **Background curl source:** preserved exactly, as the $O(t^2)$ forcing.

The transformation is a change of unknowns. It is not a new exact solution.

## 4 · Analytic-scale estimate and the conditional theorem

**Scale and state ball.**

- $D \subset \mathbb R^3$ bounded open (for example the tube $U$), $D_\rho = \{z \in \mathbb C^3 : \mathrm{dist}(z, D) < \rho\}$.
- $X_\rho$: bounded holomorphic $\mathbb C^{24}$-valued maps on $D_\rho$, with
  $\lVert\cdot\rVert_\rho = \sup_{D_\rho}\lvert\cdot\rvert_w$ and the weighted norm
  $\lvert\Phi\rvert_w = \max_X\lvert\Phi_X\rvert/R_X$.
- Ball: $\mathfrak B(R) = \{V \in X_\rho : \lVert V - V_{\rm in}\rVert_\rho \le R\}$, where
  $V_{\rm in} = \Phi_{\rm init}$ is holomorphic on $D_{\rho_0}$.
- State set: $\mathcal K = \{v : \lvert v - \Phi_{\rm init}(z)\rvert_w \le R \text{ for some } z \in D_{\rho_0}\}$, a compact subset of $\mathcal U$.

**Background bounds, all over $t \in [0, T]$ and $v \in \mathcal K$** (each is a statement about the
pointwise ODE or about explicit algebraic coefficients):

| Quantity | Meaning |
| --- | --- |
| (B1) | $F_t(v)$ exists in $\mathcal U$ |
| (B2) | $J_* = \sup\lvert J_t\rvert$, $J^{-1}_* = \sup\lvert J_t^{-1}\rvert$, $J'_* = \sup\lvert D_vJ_t\rvert$ |
| (B3) | for $\Phi \in F(\mathcal K)$: $\lvert\gamma^{-1}\rvert \le 1/g_*$, $\lvert\det\gamma\rvert \ge \delta_*$, and $\det\gamma$ in a sector avoiding $(-\infty, 0]$, so the principal $\sqrt{\ }$ is holomorphic |
| (B4) | $B_* = \sum_i\sup\lvert B^i\rvert$, $B'_* = \sum_i\sup\lvert D_\Phi B^i\rvert$ over $F(\mathcal K)$; these are explicit rational functions of $\gamma$, $\gamma^{-1}$, $(\det\gamma)^{-1/2}$, $E$, $H$, controlled by (B3) |

**Lemma 2 (loss estimate; P).** Let $\kappa_J = J_*J^{-1}_*$ and
$\Lambda = J^{-1\,2}_*J'_*B_*J_* + J^{-1}_*B'_*J_*^2 + J^{-1}_*B_*J'_*$. For
$V, W \in \mathfrak B(R) \subset X_\rho$ and $\rho' < \rho \le \rho_0$:

```math
\lVert Q(t, V) - Q(t, W)\rVert_{\rho'} \le \frac{C_1}{\rho - \rho'}\lVert V - W\rVert_\rho,\qquad
C_1 = \kappa_JB_* + \Lambda\bigl(\lVert V_{\rm in}\rVert_{\rho_0} + R\bigr),
```

```math
\lVert Q(t, V_{\rm in})\rVert_{\rho'} \le \frac{M_Q}{\rho_0 - \rho'},\qquad M_Q \le \kappa_JB_*\lVert V_{\rm in}\rVert_{\rho_0}.
```

There is no zero-order constant ($C_0 = 0$).

*Proof.* Write $Q(V) - Q(W) = \sum\tilde B^i(V)\partial_i(V - W) + \sum(\tilde B^i(V) - \tilde B^i(W))\partial_iW$.

- First sum: Cauchy on a one-variable disc of radius $\rho - \rho'$.
- Second sum: the mean value theorem on the convex ball, together with
  $\lVert\partial_iW\rVert_{\rho'} \le (\lVert V_{\rm in}\rVert + R)/(\rho - \rho')$ and
  $\lVert V - W\rVert_{\rho'} \le \lVert V - W\rVert_\rho$.
- Source: $Q(V_{\rm in}) = \sum\tilde B^i\partial_iV_{\rm in}$ and Cauchy. $\square$

Refinement of the source: since it is $O(t^2)$, $M_Q$ may be replaced by
$\tfrac{T^2}{2}\sup\lvert\partial_t^2(J_t^{-1}L(\Phi_{\rm sil}))\rvert$ times the same Cauchy factor.

**All sector convolutions at once.** No sector is estimated individually. With $\varepsilon$ a complex
parameter, $\varepsilon Q$ satisfies the hypotheses below with $C = \lvert\varepsilon\rvert C_1$ and
$M = \lvert\varepsilon\rvert M_Q$.

- The Picard iterates are holomorphic in $\varepsilon$ and converge uniformly, so the solution
  $V(t; \varepsilon)$ is holomorphic in $\varepsilon$.
- $\Phi(t; \varepsilon) = F_t(V(t; \varepsilon))$ is holomorphic too, and solves the
  $\varepsilon$-system (Proposition 1).
- By uniqueness of the order-by-order recursion, its $\varepsilon$-Taylor coefficients **are** the
  sectors $\Phi_m$.
- So the $V$-series and the $\Phi$-series differ, and the metric sectors are read off from $\Phi$, not
  from $V$.

**Abstract theorem used** (I: Baouendi–Goulaouic, corrected form, notes/17 §1). Scale $X_s$,
$0 < s \le 1$, with $\rho = \rho_0s$. Hypotheses, mapped item by item:

| BG hypothesis | Here |
| --- | --- |
| $F(t, u)$ continuous in $t$, maps $\{\lVert u\rVert_s < R\}$ into $X_{s'}$ | $F(t, u) = \varepsilon Q(t, V_{\rm in} + u)$; holomorphic by (B1)–(B4) |
| (H1) $\lVert F(t, u) - F(t, v)\rVert_{s'} \le C\lVert u - v\rVert_s/(s - s')$ | Lemma 2, $C = \lvert\varepsilon\rvert C_1/\rho_0$ |
| (H2) $\lVert F(t, 0)\rVert_s \le M/(1 - s)$ | Lemma 2, $M = \lvert\varepsilon\rvert M_Q/\rho_0$ |

Conclusion: a unique solution for $\lvert t\rvert < a(1 - s)$ whenever $a < R/(16RC + 8M)$.

**Theorem 3 (conditional; C on (B1)–(B4)).** Let $0 < s < 1$, $0 \le t \le T$, and

```math
R_\varepsilon(t) = \frac{\rho_0R(1 - s)}{t\,(16RC_1 + 8M_Q)},\qquad
T_{\rm conj} = \min\Bigl(T,\ \frac{\rho_0R(1 - s)}{16RC_1 + 8M_Q}\Bigr).
```

1. The sector series $\sum_m\varepsilon^m\Phi_m(t)$ converges in $X_{\rho_0s}$ for $\lvert\varepsilon\rvert < R_\varepsilon(t)$.
   In particular it converges at $\varepsilon = 1$ for $t < T_{\rm conj}$.
2. **Tails:** for $1 < r < R_\varepsilon(t)$,
   $\lVert\sum_{m>N}\Phi_m(t)\rVert \le S\,r^{-N-1}/(1 - r^{-1})$, with $S = \sup_{\mathcal K}\lvert F_t\rvert$ (Cauchy in $\varepsilon$).
3. **Derivatives:** spatial derivatives converge on smaller strips (Cauchy); $\partial_t$ follows from
   the equation.
4. **Nondegeneracy:** $\gamma(F_t(V))$ is invertible, by (B3).
5. **Physics at $\varepsilon = 1$.** $\Phi$ solves the coordinate E–H system with Einstein data. On the
   common domain with the ADM CK solution it coincides with that solution: the ADM solution generates
   an E–H solution through the Gauss and Bianchi identities, and the analytic E–H solution is unique.
   There it is therefore a vacuum metric.
   - Beyond the ADM domain, physicality needs constraint propagation **for the coordinate E–H form in
     the analytic class**. That was not proved. Friedrich–Rendall cover the frame form. This is
     gap G1.

Theorem 3 is **not closed for BL**. (B1)–(B4) were not certified, and $\rho_0$, $R$, the weights were
not chosen.

## 5 · Is it worth certifying? Comparison with CK′ and sprint 01

| Construction | Loss constant | Zero-order and source | Width | Time scale |
| --- | --- | --- | --- | --- |
| CK′ (ADM, certified) | weighted ADM loss rows, $C_{\rm loss} = 1.57$ (tube) | $M$, $L_A$ present (they cost ≤ ×1.55 / ×1.80) | $s_0 = 0.017$ / 0.03 | $8.7\times10^{-4}$ / $1.5\times10^{-3}$ (I) |
| sprint 01, Prop. 4 | $C_L$ of frozen curl | — (auxiliary linear hierarchy) | $s_0$ | not a full-system statement |
| this sprint, Thm 3 | $C_1 = \kappa_JB_* + \Lambda(\cdot)$, E–H principal part | none; $M_Q = O(T^2)$ | $\rho_0$ | $\rho_0R/(16RC_1 + 8M_Q)$ |

1. **Is there a mechanism that could enlarge the region?** Only $C_1 < C_{\rm loss}$. That would be a
   property of the E–H *formulation*, whose coordinate speeds are $\le\lvert\psi\rvert^{-2}$ times the
   curl symbol norm 1. Conjugation itself can only *raise* the loss constant, by $\kappa_J \ge 1$.
2. **Which constant decides?** $B_*$, the principal coefficient bound on the complex cells.
   - Its curl part is $\gtrsim\lvert\psi\rvert^{-2} \approx 0.243$ on the tube (I: `psi_abs_lower` 2.027).
     It cannot be reduced by weights: the $E$–$H$ pair couples symmetrically.
   - Its connection part $\propto\lvert\bar E\rvert R_\gamma/R_H$ can be reduced by weights, since
     there is no zero-order coupling.
   - At $T \sim 10^{-2}$, $J = I + O(T)$ and $\Lambda = O(T)$. So $J$ and $D^2F$ are **not**
     decisive.
3. **How much improvement would matter?** Idealise: $C_1 = 0.243$, same widths, the CK′-type form
   $\rho_0/(8C_1)$.
   - The ceiling is $8.7\times10^{-3}$ (tube) and $1.5\times10^{-2}$ (waist), about ×10 the certified
     times. This is arithmetic, not a bound on the method.
   - Even that leaves the result local. Via $t_* \approx 2.4(d - d_c)$ it moves the near-critical
     window from $d - d_c \lesssim 6\times10^{-4}$ to $\lesssim 6\times10^{-3}$, with $\delta_1$ still
     not explicit. It is three orders of magnitude short of the Taylor radius (≈ 2).
4. **Worth a complex-ODE certification?** **No.** The complex-ODE quantities ($J$, $J^{-1}$, $D^2F$)
   are trivial at these times. At most, a cheap interval bound of $B_*$ (with weights) on the existing
   complex cells would decide whether the ×≤10 ceiling is approached. It would not change the category
   of the result: it is another CK-scale local existence estimate.

## 6 · Check run this sprint (E)

`python -m research.calc.silent_flow_checks`: generic BL point $(0.17, 0.11, 0.5)$, $d = 0.7662$,
256-bit balls, jets of weight 8.

- The full system is compared with the frozen-metric system (curl taken in $\gamma_{\rm sil}$). The
  difference is exactly the metric-variation part.

| Field | First nonzero difference | Earlier orders | $\varepsilon = 2$ / $\varepsilon = 1$ ratio |
| --- | --- | --- | --- |
| $H$ | $t^7$: $9.4\times10^{-5}$ | $\le 10^{-74}$ | 8.000000 (pure $\varepsilon^3$) |
| $E$ | $t^8$: $3.3\times10^{-5}$ | $\le 10^{-73}$ | 16.000000 (pure $\varepsilon^4$) |
| $\gamma$ | none through $t^8$ (as predicted: first possible $t^{10}$) | — | — |

All three predictions of §2 hold. These are ball-midpoint differences, with the earlier orders at
roundoff level.

## 7 · Status

| Statement | Label |
| --- | --- |
| full linearisation of $L$ at a silent background (§2) | P |
| conjugation identity and principal matrices; similarity of symbols; source $O(t^2)$ | P |
| Lemma 2 (loss estimate, $C_0 = 0$) | P (given B1–B4) |
| Theorem 3 (nonlinear sector convergence, tails, reconstruction) | C on (B1)–(B4); physics beyond the ADM domain also needs gap G1 |
| metric-variation term first acts at $\varepsilon^3t^7$ ($H$), $\varepsilon^4t^8$ ($E$) | P (derivation) + E (BL check) |
| ceilings ×1.55 / 1.80 (sprint 01) and ×≈10 (here) | idealised arithmetic under fixed frameworks |
| certified metric region | **unchanged** |

## 8 · Single next direction (from the obstruction)

The obstruction is the analytic-scale **loss of one derivative per unit time against a fixed analytic
width**. It is untouched by any pointwise conjugation.

The one direction that changes the growth mechanism is to replace the analytic scale by **energy
(Sobolev) estimates for the frame form**, which is symmetric hyperbolic. There the controlling ratio
becomes causal width over propagation speed, i.e. tube size over $\lvert\psi\rvert^{-2}$, instead of
analytic width over loss constant.

The operator that changes is the principal $(E, H)$ curl block, now estimated through the
Bel–Robinson-type energy of notes/19 §5 instead of Cauchy inequalities. This is not started here.
