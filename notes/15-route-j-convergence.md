# 15 · Route J: a certified convergence time near the critical surface

Note 14 §5 wrote the Brill–Lindquist (BL) development in geodesic slicing as a power series in
$t$, and §6 observed that for $d$ close to $d_c$ the apparent-horizon merger of Theorem H happens
inside that series' domain of convergence, *if* the domain has a lower bound near the critical
surface $S_c$. The bound there was a heuristic, $T \approx 2.7$, from ratio tests. This note
replaces it by a **certified** bound.

*Status (session 11).*

- **Proposition T** (computer-assisted): an explicit convergence time
  $T_{\min} = 2.988\times10^{-4}\,M$, uniform over an explicit compact neighbourhood $U$ of $S_c$
  and over an explicit interval of separations.
- **Corollary J**: the common-MOTS formation of Proposition B takes place inside that certified
  region. This is only *qualitative*: the range of $d$ for which it holds is not computed.
- The estimate that would make it quantitative is stated in §7.3. It is the blocking inequality.

Conventions are those of [note 00](00-conventions.md) and `research/README.md`:

- signature $(-,+,+,+)$;
- $G = c = 1$ and $M_{\rm ADM} = 1$, so the bare masses are $\tfrac12$;
- $K_{ij} = \frac1{2N}\partial_t\gamma_{ij}$;
- $\theta_\pm = \mathrm{tr}_SK \pm H_S$.

Lengths and times are in units of $M$. All coordinates are the comoving Cartesian
coordinates $x = (x^1, x^2, x^3)$ of the initial slice, with the punctures at $(0, 0, \pm d/2)$.

## 0 · What kind of object each statement is about

Each statement below concerns exactly one of these objects; the column says which statements.

| Object | Here |
| --- | --- |
| initial data $(\gamma_0, K_0) = (\psi^4\delta, 0)$ | closed form |
| formal series $\sum_n u_n t^n$ (the recursion of note 14 §5) | explicit coefficients |
| **convergent** series: the formal series converges and its sum solves the vacuum equations | Proposition T, on $U \times \{\lvert t\rvert < T_{\min}\}$ |
| numerical approximation (truncated series, floating point) | `bl_series.py`; not used in any proof here |
| closed-form spacetime metric | none |
| local common-MOTS bifurcation (Proposition B, near $S_c$) | Corollary J, qualitative in $d$ |
| globally first common apparent horizon on $\Sigma_t$ | **not proved** for $0 < t < t_*$ (§7.4) |
| event-horizon merger | **not addressed** for BL data (note 13) |

## 1 · Statement

Let $S_c$ be the critical minimal surface of Theorem C ([note 09](09-merger-hypotheses.md)), at
$d = d_c$. It is an embedded, axisymmetric, reflection-symmetric sphere enclosing both
punctures, with meridian profile in $\rho \in [0, 0.2821]$, $\lvert z\rvert \le 0.63266$. Define

```math
U = \{x \in \mathbb R^3 : \operatorname{dist}(x, S_c) \le 0.02\},\qquad
I = [d_c - 10^{-13},\ 0.767578125] \supset (d_c, 0.767578125],
```

with $d_c - 10^{-13} = 0.7661974282399$. The upper end is the end of the interval of Theorem F′.

> **Proposition T (computer-assisted).** For every $d \in I$, the vacuum Einstein equations in
> geodesic slicing (lapse 1, shift 0), with equal-mass time-symmetric BL data at separation $d$,
> have a unique real-analytic solution $(\gamma(t), K(t))$ on $U \times (-T_{\min}, T_{\min})$, where
>
> ```math
> T_{\min} = 2.988\times10^{-4}\,M \qquad (\text{precisely } 2.98804567\ldots\times10^{-4},\ \text{rounded down}).
> ```
>
> 1. *Convergence.* The solution extends holomorphically to $\Omega_{s}(U) \times \{\lvert t\rvert < a(1 - s/s_0)\}$,
>    for every $a < T_{\min}$ and $0 < s \le s_0 = 0.012$. Here
>    $\Omega_s(U) = \{w \in \mathbb C^3 : \operatorname{dist}(w, U) < s\}$ (Hermitian distance).
>    Its Taylor series in $t$ is the series of note 14 §5. It converges on $U$ for
>    $\lvert t\rvert < T_{\min}$, uniformly for $x \in U$, $d \in I$ and $\lvert t\rvert \le \tau$, for every $\tau < T_{\min}$.
> 2. *Bounds.* Write $u = (\gamma - \gamma_0, \partial\gamma - \partial\gamma_0, K)$ and
>    $\lvert u\rvert_w = \max(\lvert u_\gamma\rvert_F/1, \lvert u_g\rvert_F/20, \lvert u_K\rvert_F/10)$
>    (Frobenius norms of the coordinate components). Then
>    $\sup_{U}\lvert u(t)\rvert_w \le 4T_{\min}K_0/10 = 0.01778$ for $\lvert t\rvert < T_{\min}$, with $K_0 = 148.77$.
>    The Taylor coefficients satisfy $\sup_U \lvert u_n\rvert_w \le 0.01778\,T_{\min}^{-n}$, so for
>    $\lvert t\rvert \le \tau < T_{\min}$
>
>    ```math
>    \sup_U\Big|u(t) - \sum_{n<N}u_n t^n\Big|_w \le 0.01778\,\frac{(\tau/T_{\min})^N}{1 - \tau/T_{\min}} .
>    ```
> 3. *Vacuum.* The solution satisfies the Hamiltonian and momentum constraints. So
>    $g = -dt^2 + \gamma_{ij}(t, x)\,dx^idx^j$ is a real-analytic vacuum metric on
>    $\mathring U \times (-T_{\min}, T_{\min})$. In it, $\gamma$ is even in $t$ and $K$ is odd.
>    The metric is axisymmetric and $z$-reflection symmetric.
> 4. *Relation to the maximal development.* The domain of dependence of $\mathring U$ in this
>    spacetime contains $\{(t, x) : \lvert t\rvert < T_{\min},\ \operatorname{dist}(x, \partial U) > \lvert t\rvert/4\}$.
>    It embeds isometrically into the maximal globally hyperbolic development (MGHD) of the BL
>    data, and $t$ is proper time along the normal geodesics of $\Sigma_0$ there.

*Evidence.*

- Parts 1–2 rest on two things: the abstract Cauchy–Kovalevskaya (CK) theorem with explicit
  constants (§3), and bounds certified in ball arithmetic (§§4–5).
- Part 3 is the standard constraint-propagation argument (§2), with the vacuum constraints also
  checked order by order.
- Part 4 uses the bound on the coordinate light speed (§6) and the standard uniqueness of
  developments, which is **cited, not re-checked**.
- The certificate: `research/calc/ck_bound.py`, record `research/calc/data/ck_bound.json`,
  replay `python -m research.calc.ck_bound` (about 2 minutes).

**What is new.** The series representation itself is standard Cauchy–Kovalevskaya
(Proposition S). New here is an *explicit, certified, uniform* time of validity on an explicit
neighbourhood of the critical surface, for an explicit interval of $d$. We found no such
computation in the literature we checked; this is a **candidate** new result.

## 2 · The system and what its solution is

**First-order reduction.** The unknowns are $\gamma_{ij}$, $g_{kij}$ (standing for
$\partial_k\gamma_{ij}$) and $K_{ij}$. Write $G = \gamma^{-1}$. The system is

```math
\partial_t\gamma_{ij} = 2K_{ij},\qquad \partial_t g_{kij} = 2\partial_kK_{ij},\qquad
\partial_t K_{ij} = -\big(A_{ij} + B_{ij}\big) + 2K_{ik}G^{kl}K_{lj} - (G^{kl}K_{kl})K_{ij}.
```

Here $A + B$ is the Ricci tensor written in terms of $(\gamma, g, \partial g)$.

- $A$ is linear in $\partial g$:

  ```math
  A_{ij} = \tfrac12G^{kl}\big(\partial_kg_{ijl} + \partial_kg_{jil} - \partial_kg_{lij} - \partial_jg_{ikl} - \partial_jg_{kil} + \partial_jg_{lik}\big).
  ```

- $B$ is algebraic in $(G, g)$: the terms $\partial_kG^{kl}\Gamma_{lij} - \partial_jG^{kl}\Gamma_{lik} + \Gamma\Gamma - \Gamma\Gamma$,
  with $\Gamma_{lij} = \frac12(g_{ijl} + g_{jil} - g_{lij})$ and
  $\partial_mG^{kl} = -G^{ka}g_{mab}G^{bl}$.

When $g = \partial\gamma$, $A + B$ is exactly $R_{ij}(\gamma)$. This is checked against the
closed-form Ricci tensor of $\psi^4\delta$ at real and complex points to $10^{-25}$
(`tests/test_ck_bound.py`).

The system is of the form $\partial_tu = F(u)$, with $F$ first order in $x$. That is the setting
of the abstract CK theorem (§3).

**Its solution solves the ADM equations.**

- $\partial_t(g_{kij} - \partial_k\gamma_{ij}) = 2\partial_kK_{ij} - \partial_k(2K_{ij}) = 0$, and
  $g = \partial\gamma$ at $t = 0$. So $g = \partial\gamma$ for all $t$, and $A + B = \mathrm{Ric}(\gamma)$.
- Hence $(\gamma, K)$ solves the vacuum ADM evolution equations with $N = 1$, $\beta = 0$:
  $\partial_tK_{ij} = -R_{ij} + 2K_{ik}K^k{}_j - KK_{ij}$.

**Symmetric unknowns.** Apply the theorem to the symmetrised system, with $F_K$ replaced by its
symmetric part. Orthogonal projection does not increase Frobenius norms, so the constants of §4
hold for it unchanged. Along its solution, $\gamma$ and $K$ are symmetric and
$A + B = \mathrm{Ric}(\gamma)$, so the projection does nothing. Hence the symmetrised solution
solves the unsymmetrised system, and by uniqueness it is *the* solution.

**Time reflection.** $(\gamma(-t), g(-t), -K(-t))$ solves the same system with the same data,
because the quadratic term is even in $K$. By uniqueness, $\gamma$ and $g$ are even in $t$ and
$K$ is odd.

**Reality and symmetry.** The data are real on real points and $F$ has real coefficients. So
$w \mapsto \overline{u(\bar t, \bar w)}$ is again a solution, and uniqueness makes $u$ real for
real $(t, x)$. The covering domain and the data are invariant under rotations about the $z$-axis
and under $z \to -z$, so by uniqueness the solution is too.

**Constraints.** The constraints are

```math
\mathcal C = R(\gamma) + K^2 - K_{ij}K^{ij},\qquad \mathcal C_i = D^jK_{ij} - D_iK .
```

At $t = 0$ both vanish: $R(\psi^4\delta) = -8\psi^{-5}\Delta\psi = 0$, and $K = 0$.

- For any solution of the ADM evolution equations, $(\mathcal C, \mathcal C_i)$ obey a linear,
  homogeneous, first-order system: Frittelli, *Phys. Rev. D* 55 (1997) 5992, eqs. (12)–(13),
  with the evolution residuals $E_{ab} = 0$ (**read**). The system is solved for the time
  derivatives, and its coefficients are analytic in $(t, x)$.
- So every $\partial_t^n(\mathcal C, \mathcal C_i)$ at $t = 0$ is a linear combination of spatial
  derivatives of lower time derivatives. By induction, all of them vanish on $\Omega_s(U)$.
- $(\mathcal C, \mathcal C_i)$ are holomorphic in $t$ on a disc, because they are built from the
  analytic solution. So they vanish identically.

Independent check: the Hamiltonian and momentum constraints of the computed series vanish order
by order to $10^{-13}$ and $10^{-15}$ respectively. This holds for BL and for Schwarzschild
(`tests/test_bl_series.py`).

**Coordinates.** The theorem is applied in the Cartesian coordinates $x$, which are regular on
the axis. The $(\rho, z)$ form of `bl_series.py` is used only off the axis, to evaluate
coefficients. It plays no role in Proposition T.

## 3 · The abstract Cauchy–Kovalevskaya theorem, with constants

Let $(X_s)_{0<s\le s_0}$ be a decreasing scale of Banach spaces:
$X_s \subset X_{s'}$ and $\lVert\cdot\rVert_{s'} \le \lVert\cdot\rVert_s$ for $s' < s$. Let
$F$ satisfy, for all $0 < s' < s \le s_0$ and $\lVert u\rVert_s, \lVert v\rVert_s < R$,

```math
\text{(H1)}\quad \lVert F(u) - F(v)\rVert_{s'} \le \frac{C\,\lVert u - v\rVert_s}{s - s'},\qquad
\text{(H2)}\quad \lVert F(0)\rVert_s \le \frac{M}{s_0 - s}.
```

**Theorem (Nishida 1977; proof of Baouendi–Goulaouic 1977).** For every
$a < s_0\min\!\big(\tfrac1{16C}, \tfrac{R}{8M}\big)$, the problem $u' = F(u)$, $u(0) = 0$ has a
unique solution $u(t) \in X_s$ for $\lvert t\rvert < a(1 - s/s_0)$, with
$\lVert u(t)\rVert_s \le 4aM/s_0 < R/2$. If $F$ is holomorphic, so is $u$ in complex $t$ (their
Remark).

*Constants.* We read the proof in M. S. Baouendi, C. Goulaouic, *Le théorème de Nishida pour le
problème de Cauchy abstrait par une méthode de point fixe*, Journées Équations aux dérivées
partielles (1977), 1–8 (numdam JEDP_1977____1_0), normalised to $s_0 = 1$. The
extraction of their condition (12) was garbled, so it is re-derived here from the readable
lemmas and inequalities. The norm is

```math
|||v|||_a = \sup_{0<s<1,\ |t|<a(1-s)}(1 - s)\sqrt{1 - \frac{|t|}{a(1-s)}}\ \lVert v(t)\rVert_s .
```

1. Lemma 1: $\lVert\int_0^tv\rVert_s \le 2a|||v|||_a$.
2. Inequalities (10)–(11) for $G(v)(t) = F(\int_0^tv)$:
   - $|||G(v)|||_a \le 8aC|||v|||_a + M$;
   - $|||G(u) - G(v)|||_a \le 8aC|||u - v|||_a$.

   These hold for $|||u|||_a, |||v|||_a < R/(4a)$; the constraint comes from Lemma 3.
3. Take the closed ball $|||v|||_a \le 2M$. Under $a < R/(8M)$ it lies inside
   $\{|||v|||_a < R/(4a)\}$.
4. Under $a < 1/(16C)$, $G$ maps the ball into itself, since $8aC\cdot2M + M < 2M$. It is a
   contraction with factor $8aC < \tfrac12$. These two conditions are (12).
5. Then $\lVert u(t)\rVert_s = \lVert\int_0^tv\rVert_s \le 2a\cdot2M = 4aM$.

Uniqueness (their §5) needs $a < 1/(8C)$, which is weaker. Rescaling $s \to s/s_0$ turns
$C, M$ into $C/s_0, M/s_0$. This gives the condition stated above.

**Application.**

- Take $R = 1$ in the weighted norm $\lvert\cdot\rvert_w$ of §1.
- Take $M = s_0K_0/R_K$, where $K_0 \ge \sup_{\Omega_{s_0}}\lvert\mathrm{Ric}(\gamma_0)\rvert_F$.
  This works because $F(0) = (0, 0, -\mathrm{Ric}(\gamma_0))$, and $s_0 - s \le s_0$.
- Then

```math
T_{\min} = \min\Big(\frac{s_0}{16\,C},\ \frac{R_K}{8K_0}\Big),\qquad
\sup\lvert u(t)\rvert_w \le \frac{4T_{\min}K_0}{R_K}.
```

Every $\lvert t\rvert < T_{\min}$ is reached: choose $a \in (\lvert t\rvert, T_{\min})$ and $s$
small. **Cauchy's estimate in $t$** on the disc of radius $a(1 - s/s_0)$, with
$a \to T_{\min}$ and $s \to 0$, gives the coefficient bound in Proposition T(2) and hence the
truncation bound.

## 4 · The scale and the bounds

**Scale.** $X_s$ is the space of bounded holomorphic maps $\Omega_s(U) \to \mathbb C^{9+27+9}$
with the sup of $\lvert\cdot\rvert_w$. The weights are $(R_\gamma, R_g, R_K) = (1, 20, 10)$.

**Cauchy estimates.** For $w \in \Omega_{s'}$, the Hermitian ball of radius $s - s'$ about $w$
lies in $\Omega_s$. So $\lvert\partial_kf(w)\rvert \le \lVert f\rVert_s/(s - s')$ for
vector-valued $f$, and

```math
\lvert\nabla f\rvert_F \le \sqrt3\,\lVert f\rVert_s/(s - s').
```

**Block coefficients.** Write $\lvert F_A(u) - F_A(v)\rvert \le \sum_Bc_{AB}\lVert u_B - v_B\rVert_s/(s - s')$
in unweighted Frobenius norms. Each bound below is proved by the triangle inequality, the mean
value theorem along a segment, and $s - s' \le s_0$.

| $c_{AB}$ | bound | why |
| --- | --- | --- |
| $c_{\gamma K}$ | $2s_0$ | $F_\gamma = 2u_K$ loses no derivative |
| $c_{gK}$ | $2\sqrt3$ | $F_g = 2\nabla u_K$, Cauchy |
| $c_{Kg}$ | $3\sqrt3\,\Gamma_F + s_0J_{Bg}$ | $\lvert A(G, D)\rvert_F \le 3\lvert G\rvert_F\lvert D\rvert_F$ (six index permutations of $D$, Cauchy–Schwarz) plus Cauchy; $B$ by the mean value theorem |
| $c_{K\gamma}$ | $s_0(J_{BG} + J_{QG})\Gamma_{\rm op}^2 + 3\Gamma_{\rm op}^2(s_0G_1 + \sqrt3R_g)$ | $G_u - G_v = G_u(\gamma_v - \gamma_u)G_v$, so $\lvert\Delta G\rvert_F \le \Gamma_{\rm op}^2\lvert\Delta\gamma\rvert_F$; the $A$-term uses $\lvert\nabla g_u\rvert_F \le G_1 + \sqrt3R_g/(s - s')$ |
| $c_{KK}$ | $s_0J_{QK}$ | quadratic term, mean value theorem |

The quantities in the table:

- $\Gamma_{\rm op} \ge \sup\lvert G\rvert_{\rm op}$, with $\Gamma_F = \sqrt3\,\Gamma_{\rm op}$.
  - Write $G = \psi^{-4}(1 + E)^{-1}$ with $E = \psi^{-4}u_\gamma$.
  - Then $\lvert E\rvert_{\rm op} \le q = R_\gamma\sup\lvert\psi\rvert^{-4} < 1$ (checked), and
    $\Gamma_{\rm op} = \sup\lvert\psi\rvert^{-4}/(1 - q)$.
  - The entries of $G - \psi^{-4}I$ are bounded by $\sup\lvert\psi\rvert^{-4}q/(1 - q)$.
- $G_1 \ge \sup\lvert\partial\partial\gamma_0\rvert_F$.
- The $J$'s are Frobenius norms of the Jacobians of $B(G, g)$ and of $Q(G, K) = 2KGK - \mathrm{tr}(GK)K$:
  - with respect to the entries of $G$ ($J_{BG}$, $J_{QG}$), of $g$ ($J_{Bg}$) and of $K$ ($J_{QK}$);
  - enclosed over the convex sets $G \in \psi^{-4}I + [\text{entry box}]$,
    $g \in g_0 + [\text{entry box of radius } R_g]$, $K \in [\text{entry box of radius } R_K]$.
  - These contain every segment the mean value theorem uses.

For the weighted norm, the Lipschitz constant is the row maximum

```math
C = \max\Big(\frac{R_K}{R_\gamma}c_{\gamma K},\ \frac{R_K}{R_g}c_{gK},\ \frac{R_\gamma}{R_K}c_{K\gamma} + \frac{R_g}{R_K}c_{Kg} + c_{KK}\Big).
```

**Analyticity of the data on $\Omega_{s_0}(U)$.**

- The data are $\psi = 1 + \sum_\pm\frac14q_\pm^{-1/2}$, with $q_\pm = \sum_i(w_i - p_{\pm,i})^2$, on
  the principal branch.
- This is holomorphic and equals the real $\psi$ on real points wherever $\operatorname{Re}q_\pm > 0$.
  That is checked on every cell.
- $\psi$ is bounded away from 0 ($\lvert\psi\rvert \ge 2.0317$), also checked.

## 5 · The computation and its checks

**Covering.**

- The profile of $S_c$ in $z \ge 0$ is enclosed by the 361 step boxes of a validated Lohner run of
  the shooting system. The run starts from Theorem C's box ($z_0 \in [0.6326595528055, 0.6326595528065]$,
  $d \in [d_c \pm 10^{-13}]$) and continues until the whole set has crossed the plane.
- Every $w \in \Omega_{s_0}(U)$ can be moved, by a real rotation about the axis and possibly
  $z \to -z$, to a point with $\operatorname{Re}w = (\rho, 0, z)$. Then $(\rho, z)$ is within
  $0.02 + s_0$ of the profile (the meridian distance is at most the 3D distance), and
  $\lvert\operatorname{Im}w\rvert < s_0$.
- All the norms above are invariant under these maps: the data are, and the formulas are
  $O(3)$-covariant with unitary actions on the entries.
- So the suprema are taken over 2,815 cells
  $\{\operatorname{Re}x \in [\rho\text{-cell}], \operatorname{Re}y = 0, \operatorname{Re}z \in [z\text{-cell}], \operatorname{Im}w \in [-s_0, s_0]^3\}$
  of side 0.005, with $d$ the interval $I$.

**Arithmetic.**

- python-flint `acb` and `arb` balls at 128 bits. The precision is set inside a context manager
  and restored.
- Jacobians are computed by forward-mode duals over balls.
- A cell on which any check fails (branch, $\lvert\psi\rvert > 0$, $q < 1$) makes the whole run
  return `ok: false`. There is no partial result.

**Result** (`research/calc/data/ck_bound.json`, 63 s for the cells):

| quantity | value |
| --- | --- |
| $c_{\gamma K}, c_{gK}, c_{K\gamma}, c_{Kg}, c_{KK}$ | 0.024, 3.4641, 4.8875, 0.8190, 0.3832 |
| $C$ (weighted) | 2.5100 (row $K$: $0.4888 + 1.6381 + 0.3832$) |
| $K_0 = \sup\lvert\mathrm{Ric}(\gamma_0)\rvert_F$ | 148.77 (cell $\rho \in [0.125, 0.13]$, $z \in [0.54, 0.545]$) |
| $s_0/(16C)$, $R_K/(8K_0)$ | $2.988\times10^{-4}$, $8.40\times10^{-3}$ |
| $T_{\min}$ | $2.988\times10^{-4}$ (the Lipschitz term binds) |
| $\sup\lvert u\rvert_w$ | $\le 0.01778$ |

The choice $s_0 = 0.012$, with weights $(1, 20, 10)$, came from a small search over
$s_0 \in \{0.012, 0.015, 0.018\}$ and six weight vectors. Any choice gives a valid
certificate; the search only selects the best one.

**Checks** (`tests/test_ck_bound.py`, about 40 s).

1. $A + B$ equals the closed-form Ricci tensor at a real and at a complex point ($< 10^{-25}$).
2. The algebraic Lipschitz bounds hold at random complex points of a cell. The test uses an
   independent NumPy implementation of $B$ and $Q$. A mutation that shrinks the bounds 100-fold
   is caught.
3. Refusals: $R_\gamma = 40$ violates $q < 1$, and a cell at a puncture violates the branch
   condition. Both return `ok: false`.
4. The flint precision is restored after a run.
5. The certified profile boxes contain the independently integrated floating-point critical
   curve.
6. **Schwarzschild benchmark** (one centre, tube about the isotropic throat sphere):
   - $K_0$ lies between the exact $\max\sqrt6\,\psi^4/R^3$ on the tube and three times it;
   - $T_{\min}$ is below the exact convergence radius $\pi\sqrt{R_{\rm in}^3/8}$, which is set
     by the crunch of geodesic slicing.

## 6 · Consequences used later

- **Lorentzian.** On $U$, $\lvert\gamma(t) - \gamma_0\rvert_F \le 0.0178$, while
  $\gamma_0 = \psi^4\delta$ has $\psi^4 \ge 17$. So $\gamma(t)$ is positive definite.
- **Coordinate light speed.**
  - A causal curve has $\lvert dx/dt\rvert^2 \le \lvert G\rvert_{\rm op} \le \Gamma_{\rm op} = 0.0623$,
    so its coordinate speed is at most $0.25$.
  - This gives the domain-of-dependence statement of Proposition T(4).
  - For example, every point within $0.0199$ of $S_c$ is in it for all $\lvert t\rvert < T_{\min}$.
- **Embedding in the MGHD.**
  - The domain of dependence of $\mathring U$ in the CK spacetime is a globally hyperbolic
    vacuum development of $(\mathring U, \gamma_0, 0)$. By the uniqueness part of
    Choquet-Bruhat–Geroch, it embeds into the MGHD of the full data, compatibly with $\Sigma_0$.
    *Cited, not re-checked.*
  - The geodesic-slicing coordinates are determined by $\Sigma_0$, so they correspond.
- **Spatial derivatives.** At $s = s_0/2$, the bounds hold on $\Omega_{s_0/2}(U)$ for
  $\lvert t\rvert < T_{\min}/2$. Cauchy's estimate then bounds
  $\lvert\partial_x^\alpha u\rvert \le \alpha!\,(2/s_0)^{\lvert\alpha\rvert}\times(\text{weights})\times0.0178$
  on $U$. These are the inputs the bifurcation estimate (§7.3) needs.

## 7 · Connection to the MOTS bifurcation

### 7.1 Provenance of the slope

Theorem C (`research/calc/fold_proof.py`) certifies, at the fold point $(z_{0,c}, d_c)$, the
derivatives

```math
\partial_dF \in [1.2315, 1.2364],\qquad \partial_tF \in [-0.51816, -0.51655],\qquad
\partial^2_{z_0}F \in [240.87, 242.70],
```

so $c = \partial_dF/(-\partial_tF) \in [2.3768, 2.3936]$.

- Its shooting function uses $\theta_+ = H + t\,p$ with $p = \mathrm{tr}_S\kappa_0$ and
  $\kappa_0 = -\mathrm{Ric}(\gamma_0)$. That is the true expansion to first order in $t$:
  $K(t) = t\kappa_0 + O(t^3)$ and $\gamma(t) = \gamma_0 + O(t^2)$.
- The slope involves only first derivatives at $t = 0$, so $c$ is the slope of the **true**
  bifurcation function.
- In Proposition B (note 10), the $O((d - d_c)^2)$ remainder in $t_*(d)$ comes from the implicit
  function theorem. Its constant is not computed.

### 7.2 What is proved

> **Corollary J (qualitative in $d$).** There is $\delta_1 > 0$ such that for every
> $d \in (d_c, d_c + \delta_1)$:
>
> - the common MOTS of Proposition B forms at $t_*(d) < T_{\min}$;
> - for $0 \le t \le t_*(d)$, all MOTSs of Proposition B's neighbourhood lie in the region
>   $\{\operatorname{dist}(x, S_c) < 0.0199\}$ of Proposition T;
> - so the whole local formation event — no MOTS near $S_c$ for $t < t_*$, the MOTS born at
>   $t_*$ — takes place where the spacetime metric equals its convergent series, with the
>   explicit error bound of Proposition T(2).
>
> $\delta_1$ is **not explicit**.

*Proof.*

- Proposition B holds for $d \in (d_c, d_c + \delta_0)$ and $0 \le t \le t_0$, with MOTSs that are
  $C^{2,\alpha}$-small normal graphs over $S_c$. Shrinking its neighbourhood keeps it valid, so we
  may take the graphs of height $< 0.0199$ and $t_0 < T_{\min}$.
- Since $t_*(d) \to 0$ as $d \downarrow d_c$, choose $\delta_1 \le \delta_0$ with $t_* < t_0$ on
  $(d_c, d_c + \delta_1)$.
- On $\{\operatorname{dist}(x, S_c) < 0.0199\} \times [0, T_{\min})$, the geodesic-slicing
  development of the MGHD is the CK solution (§6). ∎

*Heuristic only.* The linear law gives $t_* < T_{\min}$ for
$d - d_c < T_{\min}/2.3936 = 1.25\times10^{-4}$. This is the expected size of $\delta_1$. It is
*not* a proof: the quadratic remainder is not bounded. On the full interval of Theorem F′,
$t_*$ reaches about $0.0033$, eleven times $T_{\min}$. So Proposition T does not cover the
merger for all of that interval.

### 7.3 The blocking estimate

An explicit $\delta$ with $\lvert t_*(d) - c(d - d_c)\rvert \le C_2(d - d_c)^2$ on
$(d_c, d_c + \delta)$ needs quantitative versions of two steps of Proposition B.

**(i) From the full MOTS equation to the shooting function.**

- Proposition B reduces MOTSs near $S_c$ to zeros of an axisymmetric shooting function. It does
  so by an equivariant Lyapunov–Schmidt reduction: the kernel of the Jacobi operator $L_{S_c}$ is
  spanned by the positive $\phi_0$.
- A quantitative version needs a certified lower bound $\lambda_2(L_{S_c}) \ge \underline\lambda > 0$
  on the second eigenvalue, over *all* modes including the non-axisymmetric ones.
- It also needs explicit $C^{2,\alpha}$ bounds on the nonlinear part of the expansion map over
  graphs of size $\le\eta$, for the data $(\gamma(t), K(t))$ of Proposition T.
- Neither is computed. This step decides the size of the neighbourhood in which "if and only if"
  holds.

**(ii) Second derivatives of the true bifurcation function.**

- Let $\Phi(z_0, d, t)$ be the shooting function of the true slice $\Sigma_t$, and $\Phi_1$ the
  linearised one of Theorem C. Proposition T(2) and §6 bound $\gamma(t) - \gamma_0$, $K(t) - t\kappa_0$
  and their spatial derivatives by $O(t^2)$, with explicit constants.
- What is missing is a *validated shooting solver for non-conformally-flat axisymmetric data with
  $K \ne 0$*. It would turn those bounds into certified bounds on

  ```math
  \sup_B\lvert\partial_t^2\Phi\rvert,\quad \sup_B\lvert\partial_t\partial_d\Phi\rvert,\quad \sup_B\lvert\partial_t\partial_{z_0}\Phi\rvert,\quad \sup_B\lvert\partial_d^2\Phi\rvert,\quad \sup_B\lvert\partial_d\partial_{z_0}\Phi\rvert,\quad \inf_B\partial_{z_0}^2\Phi > 0
  ```

  on a box $B = [z_{0,c} \pm r] \times [d_c, d_c + \delta] \times [0, T_{\min}]$.
- With those bounds, the quantitative fold argument gives $\delta$ and $C_2$. The $t = 0$ face
  of $B$ is already certified, by `gap3_fold.py`: $\partial^2_{z_0}F \ge 125.1$ and
  $\partial_dF \ge 0.177$ on the fold window.

The precise blocking inequality for a statement at a specific $d$ is

```math
\lvert t_*(d) - c(d - d_c)\rvert \le C_2(d - d_c)^2 \ \text{ with } \ c(d - d_c) + C_2(d - d_c)^2 < T_{\min},
```

where $C_2$ is built from the six bounds above and $\lambda_2(L_{S_c})$ from (i). None of these is
certified yet.

### 7.4 The global gap

- Corollary J is about MOTSs **near $S_c$** only.
- For $0 < t < t_*$, no statement is made about common MOTSs elsewhere on $\Sigma_t$. Theorem F′
  excludes them at $t = 0$ only.
- So the common MOTS of Corollary J is **not shown to be the first common apparent horizon** of
  the slicing.
- Nothing here concerns event horizons. The BL apparent-horizon analysis gives no event-horizon
  statement (note 13).

## 8 · Audit of route E / Cauchy–Kovalevskaya

| Item | Finding |
| --- | --- |
| Kovalevskaya form | First-order system $\partial_tu = F(u, \partial_xu)$, analytic in its arguments. The note 14 form ($\partial_t^2\gamma = \ldots$) is equivalent, because $g = \partial\gamma$ propagates exactly (§2). |
| evolution equations | Signs as in `tests/test_validation.py::test_adm_evolution_equation_sign`. The recursion of note 14 §5 is the Taylor expansion of this system. |
| time symmetry | $\gamma$ even, $K$ odd, by uniqueness (§2). |
| constraint propagation | Frittelli 1997 eqs. (12)–(13), **read**: linear and homogeneous in $(\mathcal C, \mathcal C_i)$. Vanishing follows by analyticity in $t$. Checked numerically, both constraints. |
| noncharacteristic condition | $t = 0$ is noncharacteristic for the first-order system: it is solved for $\partial_t$. |
| analytic domain | $\Omega_{s_0}(U)$ with $\operatorname{Re}q_\pm > 0$ checked on every cell; the punctures' complex singular cones are excluded. |
| patches and axis | One Cartesian patch, regular on the axis. Symmetry is used only to reduce the covering. |
| uniformity in $d$ | $d$ is an interval in every bound, so the result is uniform on $I$. Note 14's earlier "$T$ is analytic in $d$, hence bounded below" was not a proof; it is now replaced by the uniform certificate. |
| existence vs effective bound | Classical CK gives existence only. The abstract theorem with the constants of §3 gives the effective bound. |
| relation to the MGHD | §6: domain of dependence plus Choquet-Bruhat–Geroch uniqueness (cited). |
| $T \approx 2.7$ (note 14) | It remains a numerical ratio estimate. The certified bound is about 9,000 times smaller. The gap comes from the worst-case loss of one derivative per Cauchy estimate and from the constant 16 of the abstract theorem. |

## 9 · Status of the claims in this note

| Claim | Category |
| --- | --- |
| Proposition T ($T_{\min} = 2.988\times10^{-4}$ on $U$, $d \in I$) | computer-assisted (Arb balls); candidate new |
| constants of the abstract CK theorem | re-derived from Baouendi–Goulaouic (lemmas read; their (12) re-derived) |
| constraint propagation | literature (Frittelli 1997, read) plus a numerical check |
| embedding of the CK region in the MGHD | standard; cited, not re-checked |
| Corollary J (local MOTS formation inside the certified region) | analytic consequence of Proposition B and T; $\delta_1$ not explicit |
| $\delta_1 \approx 1.25\times10^{-4}$ | heuristic (linear law) |
| explicit $\delta$, $C_2$ | **not achieved**; blocking estimate in §7.3 |
| first common apparent horizon; event horizons | **not addressed** |

## Discussion prompts

- The certified time is set by $s_0/(16C)$. That is a worst-case loss of one derivative,
  not the geodesic crunch that sets $T \approx 2.7$. Would a scale adapted to the hyperbolic
  structure close most of the factor of 9,000? One example is analytic norms weighted along
  characteristic cones, as in Ovsyannikov–Treves.
- Corollary J needs the spectral gap $\lambda_2(L_{S_c})$ over non-axisymmetric modes. Each
  Fourier mode $e^{im\varphi}$ gives an ODE eigenvalue problem on the profile. Is a validated
  bound for $m \ge 1$ easier than the $m = 0$ fold computation?
- The common MOTS of Corollary J is born at a time of order $10^{-4}M$ after a moment of time
  symmetry. Does anything physical, such as the area or the Hawking mass of the new MOTS, depend
  on the slicing at leading order in $d - d_c$?
