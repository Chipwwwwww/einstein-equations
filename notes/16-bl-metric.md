# 16 · The near-critical BL spacetime as an evaluable metric

Note 15 proved that near the critical surface $S_c$ the Brill–Lindquist development *equals*
its convergent time series, for a certified time. This note turns that into a usable object:
a **metric evaluator**. Given $(t, x, y, z, d)$, it returns every component of $g_{\mu\nu}$, its
inverse, and its first and second derivatives, each with a certified enclosure and an itemised
error budget. It refuses points outside the certified domain.

It also reports three side results:

- an audit of the certificate, with one rigor gap repaired (§1);
- an 86% longer certified time (§2);
- four bounded experiments that use the series to look for a closed form (§5). **None was
  found.**

*Status (session 12).*

| Object | Category |
| --- | --- |
| the metric on $U \times (-T, T)$, $d \in I$, as the complete infinite series | **B**, exact: the series converges and its sum solves $R_{\mu\nu} = 0$ |
| the truncated series returned by the evaluator, with tail bounds | certified approximation of that metric (computer-assisted) |
| a closed-form (category-A) metric | **not found**; four candidate classes tested (§5) |
| a common apparent horizon at a specific numerical $d$ | **not certified** (§7) |

Conventions are those of note 00 and note 15: signature $(-,+,+,+)$, $G = c = 1$, $M_{\rm ADM} = 1$,
$K_{ij} = \frac1{2N}\partial_t\gamma_{ij}$. The objects are recovered exactly from the code:

- **The data.** The Brill–Lindquist data are $\gamma_0 = \psi^4\delta$ with
  $\psi = 1 + \frac1{4r_+} + \frac1{4r_-}$ and $r_\pm = \lvert x - (0, 0, \pm d/2)\rvert$ (bare masses $\tfrac12$),
  and $K = 0$.
- **The gauge.** Geodesic slicing: lapse 1, shift 0, spatial coordinates comoving with the normal
  geodesics of $\Sigma_0$. So

  ```math
  g = -dt^2 + \gamma_{ij}(t, x)\,dx^idx^j,\qquad g_{00} = -1,\ g_{0i} = 0 \ \text{exactly}.
  ```

- **The coordinates.** $t$ is proper time along the normal geodesics, and $(x, y, z)$ are their
  Cartesian labels on $\Sigma_0$. This chart is regular on the axis and is the chart of the
  Cauchy–Kovalevskaya proof. The cylindrical components of note 14 are used only as an
  independent check, off the axis.

## 1 · Audit of the certificate the evaluator relies on

I re-read note 15's proof and `research/calc/ck_bound.py` against what the evaluator uses.

| Item | Finding |
| --- | --- |
| first-order reduction, constraint propagation, time reflection | as in note 15 §2; re-checked |
| abstract CK constants | re-derived in note 15 §3; unchanged |
| norm conversions, inverse-metric bound | $\lvert G - \psi^{-4}I\rvert_{\rm op} \le \sup\lvert\psi\rvert^{-4}q/(1-q)$, entries bounded by the operator norm; correct |
| **conversion of ball results to stored scalars** | **gap, repaired.** The session-11 code formed several bounds in floating point with round-to-nearest. These were the cell midpoints and radii, the hull of the $d$-interval, $\sup\lvert\psi\rvert^{-4}$, the radius of the $G$-box, the coefficients $c_{AB}$, $C$ and $T_{\min}$. Adjacent cells could leave gaps of about $10^{-19}$, and the stored $T_{\min}$ was not rounded down. Now every bound is formed in ball arithmetic and rounded outward (`up`/`down`). Cells are exact intervals $[i\,h, (i+1)h]$, with index ranges computed with outward rounding, and $d$ is the exact hull. |
| effect of the repair | the replayed classic certificate covers the same 2,815 cells and gives $T_{\min} = 2.98804567489\times10^{-4}$, $5.6\times10^{-12}$ (relative) below the stored value. Note 15's rounded-down $2.988\times10^{-4}$ stands |
| complex-domain coverage | the cells cover the componentwise inflation, by $0.02 + s_0$, of the 361 profile boxes. That is exactly the complex $s_0$-neighbourhood of the compact set $U$ below, after rotation and reflection. Note 15 stated the result on the tube of radius $0.02$, which $U$ contains |
| Cartesian proof vs coefficient engine | the new engine is Cartesian, as in the proof. It reproduces the Taylor coefficients of the *first-order system the certificate is proved for*, built with the certificate's own `ricci_parts` and `quad_K` and with no parity assumed, to $10^{-31}$ (§6). So the coefficients and the tail theorem refer to the same solution |
| a fact corrected in passing | $S_c$ meets the plane $z = 0$ at $\rho = 0.2523$, not at its widest point $\rho = 0.2821$. Points such as $(0.28, 0, 0)$ lie outside the tube |

**The certified spatial domain.** Each of the 361 Lohner step boxes $B_i$ in the meridian
half-plane contains a point of the profile of $S_c$. Let

```math
U = \{x : (\rho(x), \lvert z(x)\rvert) \in \textstyle\bigcup_i (B_i + [-0.02, 0.02]^2)\} \ \supset\ \{x : \operatorname{dist}(x, S_c) \le 0.02\}.
```

$U$ is compact, axisymmetric and reflection-symmetric, and membership is decidable in ball
arithmetic (`bl_metric.in_domain`).

## 2 · A sharper certificate: exact norms of the principal part

The dominant constant in note 15 was $c_{Kg} = 9\,\Gamma_{\rm op} + s_0J_{Bg}$. It came from the
bound $\lvert A(G, D)\rvert_F \le 3\lvert G\rvert_F\lvert D\rvert_F$ for the principal part
$A_{ij} = \tfrac12G^{kl}T_{klij}(D)$ of the Ricci tensor. Two exact computations sharpen it.

- **The principal part at $G = I$.** $\lVert A(I, \cdot)\rVert = 2$: the integer matrix $4LL^T$ has
  characteristic polynomial $(x - 16)(x - 10)^5(x - 6)^3$.
- **The general bilinear bound.** $\lvert A(G, D)\rvert_F \le \tfrac12\lVert T\rVert\,\lvert G\rvert_F\lvert D\rvert_F$
  with $\lVert T\rVert^2 = 14 + 2\sqrt{33}$, since $T^TT$ has characteristic polynomial
  $x^{54}(x - 8)^{15}(x - 16)^6(x^2 - 28x + 64)^3$. So $\tfrac12\lVert T\rVert = 2.5244$, not 3.

Both polynomials are recomputed in `tests/test_ck_bound.py`, from $L$ and $T$ read off the
certificate's own `ricci_parts`.

Write $G = \psi^{-4}I + E$ with $\lvert E\rvert_F \le \sqrt3\,\varepsilon_{\rm op}$. Then

```math
c_{Kg} = \sqrt3\big(2\sup\lvert\psi\rvert^{-4} + \tfrac12\lVert T\rVert\sqrt3\,\varepsilon_{\rm op}\big) + s_0J_{Bg},
```

and $3 \to \tfrac12\lVert T\rVert$ in $c_{K\gamma}$. A small search over $s_0$ and the weights
(record: `research/calc/data/ck_bound_sharp.json`) gives:

> **Proposition T′ (computer-assisted).** Proposition T of note 15 holds with
>
> ```math
> T_{\min} = 5.556\times10^{-4}\,M\quad(5.55699634\ldots\times10^{-4}, \text{rounded down}),
> ```
>
> on the domain $U$ above, uniformly for $d \in I = [d_c - 10^{-13}, 0.767578125]$. The parameters are
> $s_0 = 0.014$ and weights $(R_\gamma, R_g, R_K) = (0.25, 22, 10)$.
>
> The certified constants are:
>
> - $C = 1.5746$ and $K_0 = 194.46$ (3,032 cells, 128-bit balls, outward rounded);
> - $B = 4T_{\min}K_0/R_K = 0.04322$;
> - $\inf\lvert\psi\rvert \ge 2.0269$ on $\Omega_{s_0}(U)$.

The time is 86% longer than in note 15, and $C$ is 37% smaller.

**Replay.**

```bash
python -m research.calc.ck_bound sharp
```

It takes about 2.5 minutes. The classic run, `python -m research.calc.ck_bound classic`,
reproduces note 15's certificate with the rounding repair.

**Two tested ways to go further, and what blocks them.**

1. **Localisation does not help.** The theorem applies to any compact $U$, so I certified the waist
   region alone (56 profile boxes, $z \le 0.15$). The result is again $T = 5.557\times10^{-4}$.
   The binding row of $C$ is the $g$-row, $R_K/R_g\cdot2\sqrt3$, which does not depend on the
   domain. Larger $s_0$ is worse, ending in a collapse: $T = 4.65\times10^{-4}$ at $s_0 = 0.03$,
   $2.7\times10^{-5}$ at $0.05$, and $10^{-14}$ at $0.08$. $\inf\lvert\psi\rvert$ only drops
   from 2.02 to 1.75. The collapse comes from the punctures' complex singular cones
   $q_\pm = \sum_i(w_i - p_{\pm,i})^2 = 0$: inflating by $s_0$ brings cells nearer the punctures,
   and the imaginary boxes $[-s_0, s_0]^3$ reach $\sqrt3\,s_0$ at their corners.
2. **The ceiling of this method.** For any weights,

   ```math
   C \ge \max(\text{row } g, \text{row } K) \ge \sqrt{2\sqrt3\,c_{Kg}} \ge 2\sqrt3\,\sup\lvert\psi\rvert^{-2} \approx 0.84,
   ```

   so $T \le s_0/(16 \times 0.84) \approx 1.0\times10^{-3}$ at $s_0 = 0.014$. Larger $s_0$ makes
   $K_0$ and $c_{K\gamma}$ grow. This is a ceiling of the method, not of the solution: the
   numerical convergence radius is about 2.7.
3. **Recentering cannot beat it.** At time $t_1$ the abstract-CK solution is known to be
   holomorphic only on $\Omega_s(U)$ with $s < s_0(1 - t_1/a)$. Restarting from $t_1$ with that
   width and a Lipschitz constant $C_1 \ge C$ reaches at most

   ```math
   t_1 + \frac{s_0(1 - t_1/a)}{16C} = a .
   ```

   The analytic width is used up at rate $16C$ per unit time, whatever the restart schedule.
   Going further needs a method that does not lose a derivative per step. Energy estimates for
   the hyperbolic system, with validated numerics, are the natural candidate (§8).

## 3 · Error bounds for values and derivatives

Write $u = (\gamma - \gamma_0, \partial\gamma - \partial\gamma_0, K)$ and $T_s = T_{\min}(1 - s/s_0)$.
For $0 < s < s_0$, Proposition T′ with Cauchy's estimate in $t$ gives the coefficient bound

```math
\sup_{\Omega_s(U)}\lvert u_n\rvert_w \le B\,T_s^{-n},\qquad
\lvert u\rvert_w = \max\big(\lvert u_\gamma\rvert_F/R_\gamma,\ \lvert u_g\rvert_F/R_g,\ \lvert u_K\rvert_F/R_K\big).
```

The evaluator at order $N$ sums $g_0, \ldots, g_N$ for values, and $g_0, \ldots, g_{N-1}$ for
derivatives. Let $q = \lvert t\rvert/T_{\min}$ and $q_s = \lvert t\rvert/T_s$. The left-out
tails are bounded in Frobenius norm, so every component inherits the same bound:

| quantity | tail bound | source |
| --- | --- | --- |
| $\gamma_{ij}$ | $R_\gamma B\,q^{2N+2}/(1 - q^2)$ | $\gamma$-block, even powers |
| $\partial_k\gamma_{ij}$ | $R_gB\,q^{2N}/(1 - q^2)$ | the $g$-block *is* $\partial\gamma$ (no Cauchy loss) |
| $\partial_t\gamma_{ij} = 2K_{ij}$ | $2R_KB\,q^{2N+1}/(1 - q^2)$ | $K$-block, odd powers |
| $\partial_k\partial_l\gamma_{ij}$ | $\sqrt3R_gB\,s^{-1}q_s^{2N}/(1 - q_s^2)$, minimised over $s$ | Cauchy in $x$ on the $g$-block (spatial loss $s$) |
| $\partial_t\partial_k\gamma_{ij}$ | $2\sqrt3R_KB\,s^{-1}q_s^{2N-1}/(1 - q_s^2)$, minimised over $s$ | Cauchy in $x$ on the $K$-block |
| $\partial_t^2\gamma_{ij}$ | (a) $2R_KB\,T_{\min}^{-1}\sum_{i\ge N}(2i+1)q^{2i}$; (b) the field equation | Cauchy in $t$ (a); intersected with (b) |

Option (b) is valid because the solution satisfies
$\partial_t^2\gamma = 2(-\mathrm{Ric}(\gamma) + 2K\gamma^{-1}K - \mathrm{tr}(\gamma^{-1}K)K)$
exactly. Evaluating the right side on the certified enclosures of $\gamma$, $\partial\gamma$,
$\partial\partial\gamma$ and $K$ gives a second enclosure; the evaluator keeps the intersection.
Disjoint enclosures would mean a bug, and they raise an error. At $q = 1/4$ this tightens
$\partial_t^2\gamma_{xx}$ from $\pm5.9$ to $\pm0.54$.

The error budget returned with each evaluation keeps the four sources separate:

1. **Coefficient arithmetic.** The jets compute the Taylor coefficients at the point *exactly*,
   up to ball rounding: about $10^{-34}$ at 128 bits.
2. **Spatial approximation.** None.
3. **Time truncation.** The table above.
4. **Parameters and imported certificate.** The radius of $d$ if $d$ is a ball. The certificate
   file, with its replay command.

## 4 · The metric

**Low orders in closed form.** For *any* time-symmetric vacuum data in geodesic slicing,

```math
\gamma(t) = \gamma_0 - t^2\,\mathrm{Ric} + t^4\Big(\tfrac56\,\mathrm{Ric}\cdot\mathrm{Ric} - \tfrac16\lvert\mathrm{Ric}\rvert^2\gamma_0 - \tfrac1{12}\Delta\mathrm{Ric}\Big) + O(t^6).
```

Here $\mathrm{Ric}$ is the Ricci tensor of $\gamma_0$, $(\mathrm{Ric}\cdot\mathrm{Ric})_{ij} = R_{ik}\gamma_0^{kl}R_{lj}$,
and $\Delta = \nabla^k\nabla_k$ is the rough Laplacian of $\gamma_0$.

*Proof.* The recursion gives $\gamma_2 = \frac16\big(D\mathrm{Ric}[\mathrm{Ric}] + 2\,\mathrm{Ric}\cdot\mathrm{Ric}\big)$,
using $R(\gamma_0) = 0$ (Hamiltonian constraint). Take
$D\mathrm{Ric}[h]_{ij} = \tfrac12(\nabla_k\nabla_ih_j{}^k + \nabla_k\nabla_jh_i{}^k - \Delta h_{ij} - \nabla_i\nabla_j\mathrm{tr}\,h)$
with $h = \mathrm{Ric}$:

- $\mathrm{tr}\,h = 0$ and $\nabla^kR_{jk} = \tfrac12\nabla_jR = 0$;
- commuting the derivatives and using the 3D decomposition of Riemann into Ricci (with $R = 0$)
  gives $\nabla_k\nabla_iR_j{}^k = 3(\mathrm{Ric}\cdot\mathrm{Ric})_{ij} - \lvert\mathrm{Ric}\rvert^2\gamma_{ij}$.

So $D\mathrm{Ric}[\mathrm{Ric}] = 3\,\mathrm{Ric}\cdot\mathrm{Ric} - \lvert\mathrm{Ric}\rvert^2\gamma_0 - \tfrac12\Delta\mathrm{Ric}$. ∎

The identity is also checked *exactly*, in rational arithmetic, at a rational-distance point on
the plane and one on the axis (BL, $d = 3831/5000$), and at a Schwarzschild point. At a generic
point the ball residual is $5\times10^{-70}$ (Experiment A, §5). It is presumably known; no
novelty is claimed.

For BL data, with $\psi$ harmonic,

```math
R_{ij} = -\frac2\psi\,\partial_i\partial_j\psi + \frac6{\psi^2}\,\partial_i\psi\,\partial_j\psi - \frac2{\psi^2}\lvert\nabla\psi\rvert^2\delta_{ij},\qquad
\partial_i\psi = -\sum_\pm\frac{(x - p_\pm)_i}{4r_\pm^3},\quad
\partial_i\partial_j\psi = \sum_\pm\frac{3(x - p_\pm)_i(x - p_\pm)_j - r_\pm^2\delta_{ij}}{4r_\pm^5}.
```

This is an explicit $t^2$ coefficient. The $t^4$ coefficient involves derivatives of $\psi$ up to
fourth order.

*Exact values* at the rational point $x = (205597/800000, 0, 0)$, with $r_\pm = 369053/800000$ and
$d = 3831/5000$ (the exact rationals are in `bl_metric_demo.json`):

- $\gamma = 18.856892380\,\delta$
  $+\ t^2\,\mathrm{diag}(-1.915055271, -4.099021914, 6.014077185)$
  $+\ t^4\,\mathrm{diag}(-0.136848658, 0.269465233, 0.868583067) + O(t^6)$.

**Higher orders.** The general term comes from the recursion

```math
(2k+1)\kappa_k = [-\mathrm{Ric}(\gamma)]_k + [2\kappa\gamma^{-1}\kappa - \mathrm{tr}(\gamma^{-1}\kappa)\kappa]_{k-1},\qquad
g_{k+1} = \kappa_k/(k+1).
```

It is evaluated exactly in the jet algebra of `research/calc/bl_metric.py`. The cost is 0.35 s
at order 3, 1.5 s at order 4 and 6 s at order 5, per point.

**The evaluator.**
`bl_metric.evaluate_metric(t, x, y, z, d, tolerance=None, precision=128, order=None, ...)`.

- **Inputs.** Each argument may be a number, an exact decimal string, or a ball. A ball means
  "all values in it", so the same call encloses a spatial patch or a time interval.
- **Returns.**
  - `g`, `g_inv`, `dg[l][m][n]` $= \partial_lg_{mn}$ and `ddg[k][l][m][n]`, all as balls.
  - `status`, `lorentzian`, `error_budget`, `max_radius_by` (metric, first derivatives, all).
  - `order`, `precision`, `domain`, `certificate`, `seconds`.
- **Domain checks** (`DomainError` otherwise):
  - $d \in I$;
  - $x \in U$ (decided in ball arithmetic);
  - $\lvert t\rvert < T_{\min}$.

  With `allow_uncertified=True`, it instead returns the truncated sum with `status = "uncertified"`
  and no tail bound.
- **Tolerance.** With `tolerance` and no order, the order rises until the radii of the chosen
  quantities (`tolerance_on`: metric, first, or all) are below it. Otherwise it raises
  `DomainError`, naming the bound that is too large. It never reports unsupported digits.

**Lorentzian on the whole certified domain.** On $U$,
$\lvert\gamma(t) - \gamma_0\rvert_F \le R_\gamma B = 0.0108$, while $\gamma_0 \ge 2.0269^4 = 16.88$.
So $\gamma(t)$ is positive definite and $g$ is Lorentzian wherever the evaluator certifies. The
evaluator also checks Sylvester's criterion on each returned enclosure.

**Demonstration** (`python -m research.calc.bl_metric`, record `research/calc/data/bl_metric_demo.json`).

- $d = 0.7662$, so $d - d_c = 2.57\times10^{-6}$.
- Order 5: the metric through $t^{10}$, derivatives through $t^8$.
- Three points of $U$: the axis $(0, 0, 0.62)$ near the top of $S_c$, the plane $(0.26, 0, 0)$ near
  the waist, and a generic point $(0.17, 0.11, 0.5)$.
- Three times, $t/T_{\min} = 0.1, 0.25, 0.5$.
- Two patches: the boxes $(0.17, 0.11, 0.5) \pm 0.001$ and $\pm 0.0001$ (each coordinate), with $\lvert t\rvert \le T_{\min}/2$.

The numbers are in the table of §4.1.

### 4.1 Demonstration results

| point | $t/T_{\min}$ | $\gamma_{xx}$ | $\gamma_{xz}$ | $\partial_t\gamma_{xx}$ | radius: metric / first / second derivatives | time |
| --- | --- | --- | --- | --- | --- | --- |
| axis, near the top of S_c (0, 0, 0.62) | 0.1 | `[28.2049728851492 +/- 3.99e-14]` | `[+/- 1.10e-14]` | `[-0.0009118190 +/- 3.58e-11]` | 1.1e-14 / 9.6e-11 / 2.9e-06 | 7.6 s |
| axis, near the top of S_c (0, 0, 0.62) | 0.25 | `[28.204972752 +/- 8.29e-10]` | `[+/- 6.88e-10]` | `[-0.002280 +/- 6.73e-7]` | 6.9e-10 / 9.7e-07 / 1.2e-02 | 7.9 s |
| axis, near the top of S_c (0, 0, 0.62) | 0.5 | `[28.20497 +/- 5.80e-6]` | `[+/- 3.52e-6]` | `[+/- 5.13e-3]` | 3.5e-06 / 1.2e-03 / 8.3e+00 | 7.3 s |
| plane z = 0, near the waist of S_c (0.26, 0, 0) | 0.1 | `[18.7149346996347 +/- 2.17e-14]` | `[+/- 1.10e-14]` | `[-0.0002053492 +/- 1.59e-11]` | 1.1e-14 / 9.6e-11 / 2.9e-06 | 7.7 s |
| plane z = 0, near the waist of S_c (0.26, 0, 0) | 0.25 | `[18.71493467 +/- 1.01e-9]` | `[+/- 6.88e-10]` | `[-0.000513 +/- 5.93e-7]` | 6.9e-10 / 9.7e-07 / 1.2e-02 | 7.7 s |
| plane z = 0, near the waist of S_c (0.26, 0, 0) | 0.5 | `[18.71493 +/- 8.09e-6]` | `[+/- 3.52e-6]` | `[-0.001 +/- 5.90e-4]` | 3.5e-06 / 1.2e-03 / 8.3e+00 | 8.2 s |
| generic point near S_c (0.17, 0.11, 0.5) | 0.1 | `[30.2491910537313 +/- 2.31e-14]` | `[2.72278e-8 +/- 4.34e-14]` | `[0.0005965091 +/- 4.06e-11]` | 1.1e-14 / 9.6e-11 / 2.9e-06 | 9.9 s |
| generic point near S_c (0.17, 0.11, 0.5) | 0.25 | `[30.249191141 +/- 9.43e-10]` | `[1.70e-7 +/- 8.61e-10]` | `[0.001491 +/- 4.93e-7]` | 6.9e-10 / 9.7e-07 / 1.2e-02 | 9.5 s |
| generic point near S_c (0.17, 0.11, 0.5) | 0.5 | `[30.24919 +/- 4.97e-6]` | `[+/- 4.20e-6]` | `[0.003 +/- 5.81e-4]` | 3.5e-06 / 1.2e-03 / 8.3e+00 | 8.8 s |
| **patch** $(0.17, 0.11, 0.5) \pm 0.001$, $\lvert t\rvert \le T_{\min}/2$ | $[-0.5, 0.5]$ | `[3.0e+1 +/- 0.679]` | `[+/- 5.55e-6]` | — | 4.2e-01 / 4.9e+00 / 1.5e+02 | 8.8 s |
| **patch** $(0.17, 0.11, 0.5) \pm 0.0001$, $\lvert t\rvert \le T_{\min}/2$ | $[-0.5, 0.5]$ | `[30.2 +/- 0.0914]` | `[+/- 4.22e-6]` | — | 4.2e-02 / 4.6e-01 / 1.9e+01 | 9.3 s |

- Every row is `status = certified`, with `lorentzian = True`. The balls are shown as arb
  prints them: the decimal display rounds the radius outward. The radius columns are the exact
  stored radii.
- The radii are dominated by the time-truncation bounds, which are uniform over $U$, so they
  coincide at the three points. The coefficient-arithmetic part is about $10^{-34}$.
- At $t = T_{\min}/10$ the metric is certified to $10^{-14}$, and its second derivatives to
  $3\times10^{-6}$. At $T_{\min}/2$ the metric is certified to $4\times10^{-6}$; the second
  derivatives are no longer useful there ($\pm 8$), and higher order or smaller $t$ is needed.
- The patch rows enclose all components over a box and a time interval at once.
- Their widths are genuine spatial variation, not overestimation: $\lvert\nabla\gamma_{xx}\rvert \approx 4\psi^3\lvert\nabla\psi\rvert \approx 235$
  there. Shrinking the box tenfold shrinks the radius tenfold.

## 5 · From B to A: four bounded experiments

Each experiment fixed its candidate class before computing. Each was decided by an exact or
ball-arithmetic computation (`python -m research.calc.bl_structure`, record
`research/calc/data/bl_structure.json`).

| Experiment | Candidate class | Decisive quantity | Outcome |
| --- | --- | --- | --- |
| A · coefficient structure | the $t^4$ coefficient is a universal covariant expression $a\,\Delta\mathrm{Ric} + b\,\mathrm{Ric}\cdot\mathrm{Ric} + c\,\lvert\mathrm{Ric}\rvert^2\gamma_0$ | exact equality at rational points | **holds**, with $(a, b, c) = (-\tfrac1{12}, \tfrac56, -\tfrac16)$, proved (§4). It is a structural fact about the coefficients, not a closed form for the metric |
| C1 · finite ansatz | the geodesic slices stay conformally flat, $\gamma(t) = \Phi^4\cdot$(flat); this includes Schwarzschild | $t^2$ coefficient of the Cotton tensor $C_{ijk}$ of $\gamma(t)$ | **excluded**. At $t^0$, $C = 0$ (conformally flat data). At $t^2$, $C_{xxy} = -1.06297\ldots \pm 8\times10^{-67} \ne 0$ at $(0.17, 0.11, 0.5)$, $d = 0.7662$. Schwarzschild control: $\le 5\times10^{-70}$ |
| C2 · finite ansatz | silent flow: magnetic Weyl $B_{ij} = 0$ relative to the comoving observers. The equations then reduce to ODEs along each worldline. This includes Schwarzschild free fall and every "local free fall" ansatz | $t^3$ coefficient of $B$ | **excluded**. $B = O(t^3)$ with $B^{(3)}_{xx} = -0.150407\ldots \pm 5\times10^{-68} \ne 0$. Schwarzschild control: $\le 2\times10^{-71}$ |
| D · resummation (exploratory) | Padé $[4/4]$ in $s = t^2$ of the coefficients through $t^{16}$ | pole positions | **no closed form; no error control.** The smallest positive real pole gives $t \approx 2.6$–$2.8$ for $\gamma_{yy}$, $\gamma_{zz}$ at $(0.26, 0, 0)$, and $t \approx 2.1$–$2.2$ for all components at $(0.05, 0, 0.6)$. This matches the ratio estimates (2.1–2.9) and suggests a real-time, crunch-type singularity of geodesic slicing at $t$ of order $M$. Not stable across components: $\gamma_{xx}$ at $(0.26, 0, 0)$ has none. Schwarzschild control: pole at $t = 4.18$, exact crunch $3.75$ |

**What C1 and C2 exclude.**

- **C1.** The *geodesic* slices of the BL development are not conformally flat for small
  $t \ne 0$. So no metric $-dt^2 + \Phi^4\delta$, in any spatial chart, describes this foliation.
  Other slicings are not excluded.
- **C2.** No silent (ODE-along-worldlines) reduction reproduces it. That includes any
  "superposition of free falls" with a time-independent eigenframe.
- **Neither result excludes A-class metrics in general.** They exclude these two classes for this
  spacetime.
- Together with Proposition P of note 14 (Petrov type I, which rules out algebraically special
  ansätze), three structurally different ways to obtain a closed form are now closed for the BL
  head-on spacetime.

**Experiment B** (a variable change guided by invariants) was not executed beyond what C2 uses
(the magnetic Weyl tensor). It is listed in §9.

## 6 · Independent validation

All checks are in `tests/test_bl_metric.py`, `tests/test_ck_bound.py` and `tests/test_bl_series.py`.

1. **Exact Schwarzschild free fall in every Cartesian component.** Geodesic slicing from rest is
   exactly
   $\gamma_{ij} = \psi^4[F^2(\delta_{ij} - n_in_j) + (F - \tfrac32\tau F')^2n_in_j]$, with
   $F'' = -1/F^2$, $\tau = t\sqrt{M/R_0^3}$ and $R_0 = r\psi^2$. This is a closed-form (category A)
   metric for one hole, the Novikov-type comoving form. The engine reproduces it through $t^6$ to
   $10^{-25}$.
2. **The independent cylindrical float code** (`bl_series`) agrees through $t^8$ to $10^{-9}$, and
   the $\rho\varphi$, $z\varphi$ components vanish.
3. **The first-order system of the certificate** (its own Ricci reduction, no parity assumed) gives
   the same coefficients to $10^{-28}$. The odd ones vanish.
4. **Constraints.** The Hamiltonian and momentum constraints vanish order by order to $10^{-25}$.
   This tests the code; the exact claim rests on constraint propagation (note 15 §2).
5. **Symmetry.** On the axis the metric is diagonal with $\gamma_{xx} = \gamma_{yy}$. Rotated points
   give rotated tensors. The metric is even in $t$ and $\partial_t\gamma$ is odd.
6. **Overlapping enclosures.** Certified enclosures computed at 128 and 192 bits, at orders 2 and 3,
   and for a patch versus points inside it all overlap or contain each other.
7. **Refusals.** $\lvert t\rvert \ge T_{\min}$, a point far from $S_c$, and $d = 0.78$ are
   rejected. An unreachable tolerance fails with the offending bound named.

These residual checks test the implementation. That the metric *exactly* solves
$R_{\mu\nu} = 0$ rests on the analytic construction (note 15 §§2–3), not on them.

## 7 · What the metric does and does not show physically

- **A certified local development of asymptotically flat data.**
  - The data are asymptotically flat, and the local region embeds in their maximal globally
    hyperbolic development (note 15 §6, Choquet-Bruhat–Geroch, cited).
  - The patch $U \times (-T_{\min}, T_{\min})$ itself says nothing about asymptotics.
- **No certified merger at $d = 0.7662$.**
  - The linear law predicts the birth of the common MOTS at
    $t_* \approx 2.385\,(d - d_c) = 6.1\times10^{-6} = 0.011\,T_{\min}$, inside the certified slab.
  - But the range $\delta_1$ in note 15's Corollary J is not explicit, and the quadratic
    remainder is not bounded.
  - So $d = 0.7662$ is a **certified local BL development**, not a certified merger example.
- **Gauge vs physics.** Geodesic slicing crunches at $t \approx \pi M$ (Schwarzschild throat).
  That is a gauge breakdown, not a curvature singularity. Experiment D's real positive Padé poles
  are consistent with such a crunch scale. The certified slab is four orders of magnitude shorter
  than either.

## 8 · Literature and novelty

| Work | What it gives | Difference from this note |
| --- | --- | --- |
| Cauchy–Kovalevskaya; Nirenberg 1972; Nishida 1977; Baouendi–Goulaouic 1977 | analytic local existence; abstract theorem with constants | the method used here, with constants certified for BL |
| Choquet-Bruhat 1952 and its descendants | local existence in Sobolev spaces (energy estimates) | no derivative loss, but no explicit evaluable formula |
| Hintz 2024/2026 (note 11) | true solutions for EMR mergers by gluing and nonlinear iteration | an existence proof (C), not an evaluable formula (B) |
| Datta–Reiterer–Trubowitz, arXiv:1203.3766 | computer-assisted construction of the Choptuik spacetime, with rigorous error bounds | closest methodological precedent; spherical symmetry, self-similar |
| numerical relativity of head-on BL/Misner collisions (e.g. Anninos et al. 1993) | full dynamics, including the merger | no rigorous error bounds |

A bounded search found no prior rigorous, evaluable representation of a BL (or any
two-black-hole) development with certified errors.

**Novelty (candidate).** The evaluator and its certificate are a validated local representation
of a two-black-hole vacuum spacetime with explicit error bounds. That is a software and
verification result about a known spacetime, not a new solution. Two items are candidates:

- the sharp certificate;
- the exclusion of conformally flat and silent geodesic foliations for BL.

The $t^4$ formula is presumably known.

## 9 · Status

| Claim | Category |
| --- | --- |
| Proposition T′: $T_{\min} = 5.556\times10^{-4}$ on $U$, $d \in I$ | computer-assisted (Arb); candidate |
| rounding repair of note 15's certificate | done; $T_{\min} = 2.988\times10^{-4}$ confirmed |
| evaluator: certified enclosures of $g$, $g^{-1}$, $\partial g$, $\partial^2g$ on $U \times (-T_{\min}, T_{\min})$ | computer-assisted; implementation tested (§6) |
| exact sum solves $R_{\mu\nu} = 0$ | analytic (note 15) |
| universal $t^4$ coefficient | proved; exact checks |
| C1, C2 exclusions for BL | computer-assisted (ball arithmetic, nonzero enclosures) |
| category-A metric | **not found** |
| merger at a specific $d$ | **not certified** |
| method ceiling $T \lesssim 10^{-3}$; recentering gives no gain | analysis |

**Next decisive calculation.** A validated *energy-estimate* step that does not lose
derivatives. Start from the certified series data at $t_1 < T_{\min}$ as an analytic
approximate solution. Bound its residual and use a certified stability estimate for the
linearised harmonic-gauge system. This would extend the certified time toward the merger at
an explicit $d$ ($T \gtrsim 0.0033$ for the whole Theorem F′ interval) without the
$s_0/(16C)$ ceiling.

## Discussion prompts

- The certified slab is $5.6\times10^{-4}M$ thick, yet the common apparent horizon at
  $d = 0.7662$ is predicted to form inside it. How much of a "merger" can happen in so short a
  time, and what does that say about apparent horizons?
- The geodesic slices fail to stay conformally flat at order $t^2$, and silence fails at order
  $t^3$. Is there a slicing in which the BL development *does* stay conformally flat to higher
  order, for example maximal slicing?
- The single-hole development has the closed form $\psi^4[F^2P_\perp + (F - \tfrac32\tau F')^2P_\parallel]$.
  Which features of it could survive in a closed form for two holes, given that C1 and C2 fail?
