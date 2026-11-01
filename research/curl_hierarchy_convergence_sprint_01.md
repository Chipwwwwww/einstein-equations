# Curl-degree hierarchy: convergence mechanism, analyticity obstruction, controllable reconstruction (sprint 01)

*2026-11-01. Branch `claude/black-hole-merger-solutions-rpa9hj`, start commit `7a13442`. That commit is
the reviewed baseline `90f023b` plus the axisymmetric sprint, which does not touch this question.
Checker: `research/calc/curl_hierarchy_checks.py`.*

> **Scope corrections (sprint 02, `research/nonlinear_silent_flow_sprint_02.md` §1).**
>
> 1. Proposition 4 controls an *auxiliary linear hierarchy* with frozen curl coefficients. It does not
>    control the full nonlinear hierarchy, whose linearisation also contains the metric-variation loss
>    term $-D_\gamma(\mathrm{curl}\,\bar E)[Y_\gamma]$ (sprint 02 §2). The "$C_L + C_N$" sentence below is
>    withdrawn.
> 2. The ×1.55 / ×1.80 ceilings are arithmetic bounds within fixed CK′ norms, weights, width, loss rows
>    and sufficient-condition formula, with zero-order and source constraints idealised away. They are
>    not the intrinsic gain of the method, not a bound on all resummations, and not a limit of BL.
>    Removing $A$ does not remove the curl forcing.
> 3. "$t_{\rm sil}(z) \le T$" for complex $z$ is undefined and is withdrawn. Read $s_{\max}$ as the width
>    of a common regularity domain on which the background, the inverse metric and the needed variations
>    are holomorphic with uniform bounds for $0 \le t \le T$. That this width shrinks to 0 at a definite
>    rate is a **conjecture**, not a consequence of the different real collapse times.

## 1 · Question, scope, labels

**Question.** Suppose each curl-degree sector $\Phi_m$ of the session-15 hierarchy (notes/19 §4) is
computed as a full function of $t$. Does the representation $\sum_m\Phi_m$ enlarge the controllable
region of the metric? Or is the difficulty (spatial analyticity, propagation, nonlinearity) only moved
into the infinite sum over $m$?

**Scope.**

| Statement | Applies to |
| --- | --- |
| §§3–5 (lemmas, propositions) | any conformally flat, time-symmetric vacuum data, in Gauss gauge, on regions where the silent background is regular |
| numbers in §2.5 and §6 | equal-mass BL data, stored certificates (read, not rerun) |
| the transport model of §4 | calibration only; not a statement about BL |

**Labels.**

- **P**: proved here.
- **C**: conditional on named hypotheses.
- **E**: exact computation.
- **N**: exploratory numerics or heuristics.
- **I**: inherited.

## 2 · Corrections to session 15, and the literature

### 2.1 Truncation order (notes/19 §7 was wrong)

**Lemma 2 (P).** For conformally flat data, $G_j = O(t^{2j+4})$ for $j \ge 1$. Hence

```math
\gamma - \sum_{j=0}^{J} G_j = \sum_{j > J} G_j = O\bigl(t^{2J+6}\bigr).
```

So the even coefficients are guaranteed **through $t^{2J+4}$**. The $t^{2J+6}$ coefficient is *not*
guaranteed: it contains the leading term of $G_{J+1}$.

*Check with $J = 0$ (E, inherited claim 103).* The silent solution reproduces $\gamma$ through $t^4$
and fails at $t^6$; the deviation is $-\tfrac1{60}\gamma_0(\mathrm{curl}\,H)_3 \ne 0$.

Notes/19 §7 said "exact through $t^{2j+6}$". It is corrected to "through $t^{2J+4}$, error
$O(t^{2J+6})$".

### 2.2 Sector ratios at a fixed power of $t$

At $t^{2k}$ the lag lemma allows only $j \le k - 2$, so finitely many sectors appear.

- A decrease over these finitely many numbers is **not evidence** that the series over $j$ converges.
- It is also **not necessary** for convergence.
- "$j = 1$ is about 3 % of the silent part" (notes/19 §4.4) is a ratio of **Taylor coefficients at
  one point (P1) and one power ($t^{12}$)**. It is not a metric error.

Notes/19 §4.4 and claim 120 are corrected.

### 2.3 Boundedness of the first sector is automatic

By Lemma 3 below, on any interval where the silent background and the spatial derivatives it needs
are regular, every fixed sector $\Phi_m$ is bounded. It solves a linear forced ODE with bounded
coefficients. So "$G_1$ stays bounded until $t_{\rm sil}$" decides **nothing** about:

- uniformity in $m$;
- convergence at $\varepsilon = 1$;
- whether the representation goes past the Taylor radius.

The "next decisive task" of notes/19 §8 and report_session15 is withdrawn (§8 replaces it).

### 2.4 Principal-symbol quantifiers (E)

Coordinate form ($\gamma$, $K$, $E$, $H$, Christoffels from $\gamma$), exact rationals, at $\gamma = I$,
$E = \mathrm{diag}(-1, -1, 2)$, $H = 0$:

| direction $\xi$ | algebraic mult. of 0 | geometric mult. | rank $P$ / rank $P^2$ |
| --- | --- | --- | --- |
| $e_3$ | 16 | 16 | 8 / 8 (diagonalisable) |
| $(3/5, 0, 4/5)$ | 16 | 15 | 9 / 8 (Jordan block) |

**Correct statement.** At this state there is a direction with a Jordan block. Strong hyperbolicity
needs uniform diagonalisability in *all* directions, so the coordinate form is not strongly hyperbolic
there.

Session 15 wrote "at any state with $E \ne 0$, Jordan block". That was based on one random float
sample, so it is replaced: neither "every direction" nor "every state" is claimed. Conversely, the
diagonalisable direction $e_3$ does not make the system hyperbolic. Claim 115 and notes/19 §2.3 are
corrected.

### 2.5 Clifton–Gregoris–Rosquist 2017, eq. (59)

**Their assumptions.** Eqs. (46)–(47): $\Theta = \sigma = H = 0$ and $\dot\Theta = \dot H = 0$ at
$t = 0$, hence $\mathrm{div}\,E = 0 = \mathrm{curl}\,E$. Their slices are conformally flat.

**Their conventions.**

- $\mathrm{curl}\,S_{ab} = \epsilon_{cd(a}D^cS_{b)}{}^d$, as in our (15).
- $H = \mathrm{curl}\,\sigma$ (their eq. 22); ours is $H = \mathrm{curl}\,K$, the same since
  $\mathrm{curl}(\theta\gamma) = 0$.
- Their (19)–(20) match our equations (notes/19 §2.1).

This is the same setting and the same sign conventions as ours.

**Their result.** $\partial_t^3(\mathrm{curl}\,H)\vert_0 = 4\,\mathrm{curl}\,\mathrm{curl}(E^2)$, with
$(E^2)_{ab} = E_a{}^cE_{bc}$. The repo works with Taylor coefficients,
$(\mathrm{curl}\,H)_3 = \tfrac1{3!}\partial_t^3(\mathrm{curl}\,H)\vert_0$, so CGR (59) reads

```math
(\mathrm{curl}\,H)_3 = \tfrac23\,\mathrm{curl}\,\mathrm{curl}\bigl(\mathrm{Ric}^2\bigr).
```

**Checks performed.** `python -m research.calc.curl_hierarchy_checks`, exact rational arithmetic.

1. The repository's exact chain $(\mathrm{curl}\,H)_3 = 4\gamma_0(E_4 - E_4^{\rm sil})$ equals
   $\tfrac23\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ at the BL plane and axis rational points (both
   nonzero).
2. The four-term law of claim 105 equals $\tfrac23\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ at all
   63 weight-6 harmonic-jet points.

   By the reduction of notes/19 §3.3, the two expressions are therefore identical for **every**
   conformally flat time-symmetric vacuum data set. This is P, computer-assisted, exact.

**Classification.**

| Item | Class |
| --- | --- |
| $(\mathrm{curl}\,H)_3 = \tfrac23\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ | **known**: CGR 2017, eq. (59). Claim 105 is its transcription to Taylor coefficients plus its conformally flat expansion |
| the four-term expansion of $\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ via $\nabla\mathrm{Ric}$ totally symmetric and $\Delta\mathrm{Ric} = 3\mathrm{Ric}^2 - \lvert\mathrm{Ric}\rvert^2\gamma_0$ | known result rewritten in this gauge |
| the **metric** coefficient $\gamma_3 = -\tfrac{13}{180}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\gamma_0 - \tfrac1{90}\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ (Theorem 1 restated) and its universal exact verification | formula and verification added by this project; **novelty not established**. CGR use curl H on LRS curves, not the metric coefficient, but a wider search is needed |
| the Einstein–Bianchi system in Gauss gauge, symmetric hyperbolicity of its frame form | known [Friedrich–Rendall §4.4; Ellis–van Elst]; nothing new |

## 3 · The hierarchy and a finite-sector regularity lemma

**Unknowns.** $\Phi = (\gamma_{ij}, K_{ij}, E_{ij}, H_{ij})$: lower comoving components, symmetric,
24 functions.

**System** (notes/19 §2.1):

```math
\partial_t\Phi = A(\Phi) + \varepsilon\,L(\Phi, \partial\Phi),\qquad
A(\Phi) = \bigl(2K,\ -E + K\!\cdot\!K,\ -2\theta E + \tfrac52\{K, E\} - \gamma\,\mathrm{tr}(KE),\ -2\theta H + \tfrac52\{K, H\} - \gamma\,\mathrm{tr}(KH)\bigr),
```

```math
L(\Phi, \partial\Phi) = \bigl(0,\ 0,\ \mathrm{curl}_\gamma H,\ -\mathrm{curl}_\gamma E\bigr),
```

- **$A$** is rational in $\Phi$ (through $\gamma^{-1}$) and contains no derivatives.
- **$L$** is linear in $(\partial E, \partial H)$, and bilinear in $(\partial\gamma, (E, H))$ through the
  Christoffel symbols and $\det^{-1/2}$.
- **$\varepsilon$** is a formal parameter multiplying every curl.
- **Data:** $\Phi(0) = (\gamma_0, 0, \mathrm{Ric}(\gamma_0), 0)$, independent of $\varepsilon$.
- **Parity:** $\Phi(\cdot\,; -\varepsilon) = P\,\Phi(\cdot\,; \varepsilon)$ with $P = \mathrm{diag}(1, 1, 1, -1)$.
  So $\gamma$, $K$, $E$ are even in $\varepsilon$ and $H$ is odd.
- **$\varepsilon = 1$** is the vacuum Einstein system. The constraints propagate only there.

For $\varepsilon \ne 1$ the deformed system is not Einstein. The grading belongs to *this* formulation:
putting $\varepsilon$ on the curls of a different formulation (frame form, ADM) defines a different
deformation, and nothing here identifies the two.

**Sectors.** With $\Phi = \sum_{m\ge0}\varepsilon^m\Phi_m$:

- $\partial_t\Phi_0 = A(\Phi_0)$ is the silent system, with $H_0 = 0$.
- For $m \ge 1$:

```math
\partial_t\Phi_m = DA(\Phi_0)\Phi_m + S_m,\quad S_m = \mathcal N_m(\Phi_1, \dots, \Phi_{m-1}) + \bigl[L(\Phi, \partial\Phi)\bigr]_{m-1},\quad \Phi_m(0) = 0,
```

```math
\Phi_m(t, x) = \int_0^t U_0(t, s; x)\,S_m(s, x)\,ds,\qquad \partial_tU_0 = DA(\Phi_0(t, x))\,U_0,\ U_0(s, s) = I .
```

Here $\mathcal N_m$ collects the higher derivatives of $A$ at $\Phi_0$, applied to products of lower
sectors. It contains no spatial derivatives. $[L]_{m-1}$ contains exactly one spatial derivative of
the lower sectors.

**Lemma 3 (finite-sector regularity; P).** Let $D \subset \mathbb R^3$ be open, $T > 0$ and
$k \ge 1$. Assume:

- (a) $\Phi_0 \in C^k$ in $x$, uniformly on $[0, T]\times D$;
- (b) $\gamma_0(t, x) \ge g_*\,I$ with $g_* > 0$ there. This holds for
  $T < \inf_{\bar D}t_{\rm sil}$, with a margin.

Then for each $1 \le m \le k$, $\Phi_m$ exists on $[0, T]\times D$, is unique and $C^{k-m}$ in $x$, and

```math
\sup_{[0,T]\times D}\lvert\partial_x^\alpha\Phi_m\rvert \le C_{m,\alpha}\bigl(T,\ g_*,\ \lVert\Phi_0\rVert_{C^{m+\lvert\alpha\rvert}},\ \sup\lVert\partial_x^{\le m+\lvert\alpha\rvert-1}U_0\rVert\bigr).
```

*Proof.* Induction on $m$. $S_m$ is a polynomial in $\Phi_{<m}$ and in *one* derivative of
$\Phi_{m-1}$ and $\gamma_{<m}$, with coefficients rational in $\Phi_0$ and denominators bounded below
by (b). So $S_m \in C^{k-m}$. $U_0$ inherits the $C^k$ regularity of $\Phi_0$ by smooth dependence on
parameters. Grönwall on $[0, T]$ bounds the integral. $\square$

**What is automatic and what is not.**

- Fixed-$m$ boundedness is automatic. The required regularity costs one derivative per $\varepsilon$.
- The constant $C_{m,0}$ grows with $m$ through $\lVert\Phi_0\rVert_{C^m}$, i.e. $m$-th spatial
  derivatives of the silent background.
- Lemma 3 says nothing about $\sum_mC_{m,0}$. All of the difficulty sits in that growth.

## 4 · A solvable model, and what transfers to BL

**Transport model.** $\partial_tu = \varepsilon c\,\partial_xu$, $u(0) = f$. Then $u = f(x + \varepsilon ct)$
and $u_m = (ct)^mf^{(m)}(x)/m!$.

- Every sector is a polynomial in $t$, defined for all $t$.
- The sum is the Taylor series of $\varepsilon \mapsto f(x + \varepsilon ct)$. Its radius in
  $\varepsilon$ is $\mathrm{dist}(x, \mathrm{sing}\,f)/(ct)$, where the distance is to the complex
  singularities of $f$.
- So the sum converges at $\varepsilon = 1$ iff $ct < \mathrm{dist}_{\mathbb C}(x, \mathrm{sing}\,f)$.
  For $f = 1/(1 + x^2)$: $ct < \lvert x \pm i\rvert$.
- The $m!$ from the $m$ derivatives cancels the $1/m!$ from the $m$ time integrations, which leaves
  geometric growth.

**General mechanism (P).** Fix $(t, x)$. The sector sum is a power series in $\varepsilon$.

- If it converges at $\varepsilon = 1$, it converges for $\lvert\varepsilon\rvert < 1$, where it is
  holomorphic in $\varepsilon$ (Abel).
- If the Taylor coefficients in $\varepsilon$ are the $\Phi_m$ of an $\varepsilon$-holomorphic solution
  of the deformed system, then the radius in $\varepsilon$ equals the distance from 0 to that
  solution's nearest singularity in $\varepsilon$ (Cauchy–Hadamard).
- For complex $\varepsilon$ the characteristic speeds $\varepsilon c$ become complex, and the deformed
  solution moves information *into complex $x$*.

So **convergence at $\varepsilon = 1$ up to time $t$ is, in this mechanism, a statement of spatial
analyticity over a complex width of about $c\,t$.** That is exactly the kind of condition that limits
the Taylor series in $t$ (CK). The model calibrates this; it does not prove anything about BL.

**What BL changes, and what it does not.**

- *Trivial linear sector* (claim 117). For conformally flat data the linearisation around flat space
  gives $E_{\rm lin} = \mathrm{const}$ and $H_{\rm lin} = 0$. Data singularities are therefore not
  propagated at linear order, and the model's data-singularity mechanism does not operate on BL
  directly.

  A heuristic check: the model's scale at P1, $\min r_\pm\,\psi^2/\sqrt2 \approx 0.91$ (complex distance
  to the punctures' cones $q_\pm = 0$ times $\psi^2$, the coordinate speed$^{-1}$), lies well below the
  observed Taylor radius $\approx 2.15$ (claim 106, N).
- *Delay* (Lemma 2). Sector $j$ starts at $t^{2j+4}$. In a delayed-source model
  ($\partial_tE = \varepsilon DH + N(E)$, $\partial_tH = -\varepsilon DE$, $DE_0 = 0$, $DN(E_0) \ne 0$), the
  $2j$ derivatives in sector $j$ are again balanced by the $2j + 4$ time integrations. The ratio
  $(2j)!/(2j+4)!$ is polynomial in $j$, so **the geometric rate $(\kappa t/s)^{2j}$ is unchanged**. The
  delay only shifts powers. A better rate would need cancellations in $\mathrm{curl}^{2j}$ acting on
  silent-generated sources, and none is known.
- *The term that matters.* In the actual hierarchy the operator iterated $m$ times is
  $\Psi \mapsto \int_0^tU_0(t, s; x)\,\mathrm{curl}_{\gamma_0}\Psi\,ds$. It acts on **nonlinear silent
  quantities**: the first source is $\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ (CGR 59).

  Its $m$-fold growth is controlled by the complex analyticity width of the **silent background**
  $\Phi_0(t, \cdot)$, not of the data. That width is
  $s_{\max}(t) = \mathrm{dist}_{\mathbb C}\bigl(D, \{z : t_{\rm sil}(z) \le t\}\bigr)$, and it **decreases
  to 0** as $t \to t_{\rm sil}$, because the local collapse time varies in $x$.

  This is the BL-relevant obstruction. It is absent if $t_{\rm sil}$ is constant (homogeneous
  collapse) and present whenever $\nabla t_{\rm sil} \ne 0$: for BL, $t_{\rm sil}$ ranges over 2.05–5.6
  across the session-14 points.

## 5 · Main proposition (route A): a uniform-in-$m$ bound and its exact limitation

*Why route A.* §4 identifies the obstruction as uniform analytic control of $m$-fold derivative
losses. Route B (an a-posteriori residual estimate at $\varepsilon = 1$) needs the frame form, since
the coordinate form is not strongly hyperbolic (§2.4), together with commuted energies. Notes/19 §5
has only the $L^2$ identity. Route A answers the question asked directly.

**Scale.** $X_s$ is the space of bounded holomorphic functions on
$D_s = \{z \in \mathbb C^3 : \mathrm{dist}(z, D) < s\}$, with the sup norm, for $0 < s \le s_0$. For
$s' < s$, Cauchy gives $\lVert\partial f\rVert_{s'} \le \lVert f\rVert_s/(s - s')$.

**Proposition 4 (linear principal part; P).** Assume on $[0, T]$:

- **(H1)** $\Phi_0(t, \cdot)$ is holomorphic on $D_{s_0}$, with $\lvert\det\gamma_0\rvert \ge \delta_* > 0$
  there;
- **(H2)** $\sup_{0\le\sigma\le t\le T}\lVert U_0(t, \sigma; \cdot)\rVert_{s_0} \le B$, the sup norm of the
  matrix function;
- **(H3)** the principal operator $L_0\Psi := (0, 0, \mathrm{curl}_{\gamma_0}\Psi_H, -\mathrm{curl}_{\gamma_0}\Psi_E)$,
  with coefficients frozen at $\Phi_0$, satisfies $\lVert L_0\Psi\rVert_{s'} \le C_L\lVert\Psi\rVert_s/(s - s')$.

Define the linear hierarchy $\Psi_0 = \Phi_0$, $\Psi_m = \int_0^tU_0\,L_0\Psi_{m-1}$. Then for
$0 < s < s_0$:

```math
\lVert\Psi_m(t)\rVert_{s} \le \lVert\Phi_0\rVert_{s_0}\left(\frac{e\,B\,C_L\,t}{s_0 - s}\right)^{m},
```

so $\sum_m\varepsilon^m\Psi_m$ converges in $X_s$ for $\lvert\varepsilon\rvert < (s_0 - s)/(eBC_Lt)$.
In particular, $t < (s_0 - s)/(eBC_L)$ is **sufficient** for convergence at $\varepsilon = 1$. It is
not necessary. On smaller
strips the sum can be differentiated term by term in $x$ (Cauchy) and in $t$.

*Proof.* Use nested widths $s_j = s + (s_0 - s)(m - j)/m$. Each step costs
$BC_L/(s_{j-1} - s_j) = BC_Lm/(s_0 - s)$ and one time integration. Hence
$\lVert\Psi_m(t)\rVert_s \le \lVert\Phi_0\rVert_{s_0}(BC_Lm/(s_0 - s))^mt^m/m!$, and
$m^m/m! \le e^m$. $\square$

**Conditional extension (C).** The full hierarchy adds the nonlinear mixing $\mathcal N_m$ and the
dependence of the curl coefficients on $\gamma_{\ge1}$. A Cauchy–Kovalevskaya majorant for these terms
would give a bound of the same form, with $C_L$ replaced by $C_L + C_N(\lVert\Phi_0\rVert_{s_0}, \delta_*)$ *(withdrawn: see the scope corrections above)*.
This is not proved here; it is the standard majorant step, but its constants were not computed.

*If the conditional extension holds:*

- the sum at $\varepsilon = 1$ solves the Einstein–Bianchi system;
- it equals the unique analytic solution, so it is the CK series of the evaluator and the constraints
  hold;
- $\gamma$ stays nondegenerate as long as the tail $\lVert\sum_{m\ge1}\rVert_s$ stays below the
  nondegeneracy margin of $\gamma_0$.

**The limitation, stated exactly (P, given H1–H3).**

1. (H1) must hold over the whole interval. The silent background at complex $z$ exists only while
   $t < t_{\rm sil}(z)$, so the admissible width is at most $s_{\max}(T)$ of §4, and $s_{\max}(T) \to 0$
   as $T \to \inf t_{\rm sil}$.
2. The sufficient condition $eBC_LT < s_{\max}(T) - s$ therefore stops **strictly before** the local
   collapse time. Its location is set by the spatial variation of $t_{\rm sil}$ and by the complex
   propagator bound $B$.
3. Compared with CK, the sector structure replaces the zero-order constants (the $M$, $L_A$ rows of
   CK′) by the exact silent propagator. It **does not touch the loss term**, which carries the factor
   $1/(s - s')$.

**Quantitative consequence for the existing certificates** (I + arithmetic). In the CK′ framework
used for the certificates (same norms, weights and loss rows), set all zero-order constants to zero:
$L_A = 0$, $M = 0$, $\beta \to 1$. The best time is then $s_0/(8C_{\rm loss})$.

| Region | $s_0$ | $C_{\rm loss}$ | ceiling $s_0/(8C_{\rm loss})$ | certified now | ratio |
| --- | --- | --- | --- | --- | --- |
| tube | 0.017 | 1.5746 | $1.35\times10^{-3}$ | $8.68\times10^{-4}$ | at most 1.55 |
| waist | 0.03 | 1.3472 | $2.78\times10^{-3}$ | $1.545\times10^{-3}$ | at most 1.80 |

So with the certified loss constants, *exact treatment of the local (silent) dynamics* can at best
raise the certified time by a factor below 2. Reaching the Taylor radius ($\approx 2$, i.e. three
orders of magnitude more) would need a far smaller loss constant or a far larger analyticity width.
The E–H formulation's own loss constant (coordinate speeds $\lesssim\psi^{-2}$ times the curl symbol
norm 1, plus Christoffel terms) was **not computed**. A smaller value would be a gain of the
*formulation*, not of the sector decomposition.

**Heuristic illustration (N, not evidence).** Set $B = 1$, $C_L \approx \psi^{-2}(P1) \approx 0.18$, and
$\lvert\nabla t_{\rm sil}\rvert \approx 0.27$ (from the session-14 values at P1 and P4). The condition
$eBC_LT = s_{\max}(T)$, with $s_{\max} \approx (t_{\rm sil} - T)/\lvert\nabla t_{\rm sil}\rvert$, gives
$T \approx 1.9$. That is the same order as the remainder radius $\approx 1.85$ at P1 (claim 106). None
of the inputs is certified.

## 6 · Status of the statements

| Statement | Label |
| --- | --- |
| truncation $\gamma - \sum_{\le J}G_j = O(t^{2J+6})$, coefficients guaranteed through $t^{2J+4}$ | P |
| fixed-order sector ratios are not evidence of convergence | P (logic) |
| Lemma 3: finite-sector regularity, constants not uniform in $m$ | P |
| exact symbol multiplicities at $E = \mathrm{diag}(-1, -1, 2)$: (16, 16) along $e_3$, (16, 15) along $(3/5, 0, 4/5)$ | E |
| claim 105 equals CGR (59): at BL points and for all conformally flat vacuum data | E / P (computer-assisted); the result itself is **known** |
| $\gamma_3 = -\tfrac{13}{180}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\gamma_0 - \tfrac1{90}\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$ | P (Theorem 1 + CGR check); novelty unverified |
| convergence at $\varepsilon = 1$ means spatial analyticity of complex width about $ct$ (model and Abel/Cauchy–Hadamard) | P (mechanism) |
| Proposition 4 (linear principal part, geometric bound) | P |
| full nonlinear hierarchy bound | C (majorant step not done) |
| ceilings 1.35e-3 / 2.78e-3 for exact zero-order treatment in CK′ with certified loss rows | I + arithmetic |
| $T \approx 1.9$ illustration | N |
| certified region | **unchanged**: tube $8.683\times10^{-4}$, waist $1.545\times10^{-3}$ |

## 7 · What this adds, and what it rules out

- **Answer to the question.** Computing each sector exactly in $t$ moves the difficulty into the sum
  over $m$. There the $m$-fold derivative loss meets the shrinking complex analyticity width of the
  silent background.
- **What sector completeness buys.** Only the removal of the zero-order (non-derivative) constraints.
  With the certified loss constants that is a factor below 2 in certifiable time, not an extension to
  the collapse time.
- **Inferences ruled out** (each previously at risk in this project):
  - first-sector boundedness ⇒ viability;
  - fixed-order decrease ⇒ convergence;
  - real-time regularity of every sector ⇒ a larger complex Taylor radius;
  - one random symbol sample ⇒ a universal multiplicity statement.
- **The obstruction is not the weak hyperbolicity of the coordinate form.** It appears equally in the
  frame form. It is the analytic-scale loss of one derivative per order of coupling.

## 8 · The single next mathematical task

**Certify (H1)–(H2) for the BL silent ODE on complex boxes around the tube.** Concretely:

- a lower bound on the complex analyticity width $s_{\max}(T)$ of the silent solution, i.e. the
  distance from the real tube to the complex set $\{t_{\rm sil}(z) \le T\}$;
- an upper bound on the complex propagator $B(T)$.

Both are statements about an explicit 9-dimensional ODE (eigen-scale factors, $K$, $E$) with analytic
dependence on $z$, so they can be certified with the existing interval tools. Together with
Proposition 4 and the majorant step they yield a certified sector-convergence time. Comparing it with
$s_0/(8C_{\rm loss})$ then decides whether the hierarchy offers any certified gain at all.

If $s_{\max}$ falls linearly with slope $\lvert\nabla t_{\rm sil}\rvert$ as §4 predicts, the answer is
already bounded by §5. Then the route should be modified: change the background (a resummation around
a better-than-silent local model), not refine the sectors.

## Reproduce

```bash
python -m research.calc.curl_hierarchy_checks
```

It runs the CGR (59) checks (exact at BL points; 63 harmonic-jet points) and the exact symbol
multiplicities.
