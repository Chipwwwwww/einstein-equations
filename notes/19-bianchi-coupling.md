# 19 · The spatial coupling in Einstein–Bianchi form: a verified formulation, a proved $t^6$ law, and a curl-degree hierarchy

Session 14 (notes/18) found two things about the BL development:

- its local, *silent* part (a closed ODE on each worldline) is exact through $t^4$;
- the nonlocal remainder carries the rest, and it sets the series' radius.

Its single nonlocal source is $\mathrm{curl}\,H$. This note takes that coupling as the object of
study.

- **§2** derives the Einstein–Bianchi system in Gauss gauge with the repository's conventions. It
  checks the system against the primary sources and against the existing evaluator, and corrects the
  note-17 statement about its hyperbolicity.
- **§3** proves the $t^6$ law of note 17 for **all** conformally flat time-symmetric vacuum data.
- **§4** organises the whole coefficient hierarchy by curl degree, proves a lag lemma, and runs the
  first decisive experiment of the main route.
- **§5** proves the energy identity for the alternative route.
- **§§6–8** give the portfolio, the role of each lemma, and the continuation.

*Status (session 15).*

| Item | Category |
| --- | --- |
| E–H system in Gauss gauge, repo conventions (§2.1) | derived; agrees term by term with Ellis–van Elst (40)–(41) and, after converting the signature, with Friedrich–Rendall's $P_{ij}$, $Q_{kl}$; reproduces the ADM coefficients **exactly** through $t^8$ |
| frame form symmetric hyperbolic, speeds $0, \pm\tfrac12, \pm1$ | inherited (Friedrich–Rendall §4.4), symbol recomputed here |
| coordinate form ($\gamma$, $K$, $E$, $H$, Christoffels from $\gamma$) **not** strongly hyperbolic | computed (exact): at $E = \mathrm{diag}(-1,-1,2)$, $H = 0$, a Jordan block at speed 0 along $(3/5,0,4/5)$, none along $e_3$. *(Session 15 wrote "whenever $E \ne 0$" from one float sample; corrected in the curl-hierarchy sprint, §2.4.)* |
| $\Delta\mathrm{Ric}$, $t^4$ and $t^6$ laws for all conformally flat time-symmetric vacuum data | **proved** (reduction §3.3 plus 72 exact evaluations) |
| $(\mathrm{curl}\,H)_3$ law (claim 105) | proved (follows from §3 and claims 103–104) |
| lag lemma: the curl-degree-$2j$ sector of $\gamma$ is $O(t^{2j+4})$ for conformally flat data | proved (induction, §4.1); exact checks |
| curl-degree sectors at three generic points to $t^{12}$: lag lemma confirmed; past $j = 1$ the sectors decrease geometrically (ratios about 0.2–0.4, one exception) | numerical, exploratory (§4.4) |
| energy identity and flux bound | proved (§5); exact jet check |
| category-A metric; certified merger | **no**; **no** |

Conventions are those of note 00 and note 15:

- signature $(-,+,+,+)$, MTW/Carroll Riemann sign, $G = c = M_{\rm ADM} = 1$;
- geodesic slicing $g = -dt^2 + \gamma_{ij}dx^idx^j$ in comoving coordinates;
- $K_{ij} = \tfrac12\partial_t\gamma_{ij} = \nabla_iu_j$, with $u = \partial_t$;
- $K = \sigma + \tfrac\theta3\gamma$, $\theta = \mathrm{tr}K$.

Sources quoted in another signature are converted.

## 1 · Baseline

- **Base.** PR #5's head `bbc2368` (`claude/bl-structure-discovery`). It is newer than PR #4's
  `c26cb9f` and contains the silent/nonlocal split this note builds on.
- **Baseline rerun:**
  - `tests/test_bl_silent.py` passes;
  - the session-14 full fast suite was green at that head (151 passed).
- **Corrections inherited and respected:**
  - the CK times of notes 15–16 are superseded (notes/17 §1);
  - the certified times remain $8.683\times10^{-4}$ (tube) and $1.545\times10^{-3}$ (waist);
  - the $t^6$ law was "fitted" (claim 98) until §3 below.

## 2 · The Einstein–Bianchi system in Gauss gauge

### 2.1 Derivation in the repository's conventions

**Observer, frame and gauge.**

- $u = \partial_t$ is geodesic and hypersurface-orthogonal: $a = 0$, $\omega = 0$.
- Coordinates are comoving (the evaluator's chart).
- $E_{ab} = C_{acbd}u^cu^d$.
- $H_{ab} = (\mathrm{curl}\,K)_{ab}$, with $(\mathrm{curl}\,A)_{ij} = \mathrm{sym}(\epsilon_i{}^{kl}D_kA_{lj})$
  and $\epsilon^{mkl} = [mkl]/\sqrt{\det\gamma}$. This is Ellis–van Elst's constraint (C3), eq. 34, with
  $\omega = 0$. It fixes the sign of $H$ by the same orientation that the curl uses.

**Frame versus coordinate derivatives.** Ellis–van Elst's $\dot E_{\langle ab\rangle}$ is the
projected covariant derivative along $u$. For geodesic $u$, $\nabla_uE$ is already spatial and
trace-free. In comoving coordinates $\nabla_u\partial_i = K^k{}_i\partial_k$, so

```math
(\nabla_uE)_{ij} = \partial_tE_{ij} - K_i{}^kE_{kj} - K_j{}^kE_{ik}.
```

Ellis–van Elst (40)–(41) [EvE98], with $\mu = p = q = \pi = 0$ and $a = \omega = 0$, read
$\nabla_uE = \mathrm{curl}\,H - \theta E + 3\,\mathrm{STF}(\sigma E)$, and the same with
$(E, H) \to (H, -E)$. In lower-index comoving components:

```math
\begin{aligned}
\partial_t\gamma_{ij} &= 2K_{ij}, \\
\partial_tK_{ij} &= -E_{ij} + K_{ik}K^k{}_j, \\
\partial_tE_{ij} &= \phantom{-}\varepsilon\,(\mathrm{curl}\,H)_{ij} - 2\theta E_{ij} + \tfrac52(K\!\cdot\!E + E\!\cdot\!K)_{ij} - \gamma_{ij}\,\mathrm{tr}(K\!\cdot\!E), \\
\partial_tH_{ij} &= -\varepsilon\,(\mathrm{curl}\,E)_{ij} - 2\theta H_{ij} + \tfrac52(K\!\cdot\!H + H\!\cdot\!K)_{ij} - \gamma_{ij}\,\mathrm{tr}(K\!\cdot\!H),
\end{aligned}
```

with $(K\cdot E)_{ij} = K_i{}^kE_{kj}$ and $\varepsilon = 1$. The parameter $\varepsilon$ multiplies every
curl and is used in §4.

**Where each equation comes from.**

- The second equation is the ADM evolution equation rewritten with the vacuum Gauss equation
  $E_{ij} = R_{ij} + KK_{ij} - K_{ik}K^k{}_j$.
- The last two are the once-contracted Bianchi identities.

**Primary-source check.** Friedrich–Rendall [FR00, §4.4, p. 52] give the Gauss-gauge reduction for
pressure-free matter, which includes vacuum. Their gauge conditions are $\Gamma_0{}^j{}_k = 0$ and
$e^\mu{}_0 = \delta^\mu{}_0$. Their unknowns are $(e^\mu{}_a, \Gamma_b{}^i{}_j, E_{ij}, B_{ij})$, and the
reduced system is $T_0{}^j{}_a = 0$, $\Delta_{ij0a} = 0$, $P_{ij} = 0$, $Q_{kl} = 0$.

Their signature is $(+,-,-,-)$. The spatial metric $h = -\gamma$ is negative definite; lower-index $E$,
$B$, $\chi$ change sign while mixed ones do not; and $\epsilon\epsilon$ products equal
$-\det[h_{\cdot\cdot}]$. With these conversions their $P_{ij}$ and $Q_{kl}$ (read on the page) become
exactly the two Bianchi equations above. A transcription *without* the conversion gives
$\partial_tH = -\mathrm{curl}\,E + \tfrac12(K\cdot H + H\cdot K) + \gamma\,\mathrm{tr}(K\cdot H)$. That
version is wrong, and §2.6 shows the jet comparison detects it at $t^8$.

### 2.2 Constraints and their propagation

The constraints are:

- Gauss: $E = \mathrm{Ric} + \theta K - K\cdot K$;
- (C3): $H = \mathrm{curl}\,K$;
- momentum: $D^jK_{ij} = D_i\theta$;
- $\mathrm{div}\,E$: $D^bE_{ab} = \epsilon_{abc}\sigma^b{}_dH^{cd}$;
- $\mathrm{div}\,H$: $D^bH_{ab} = -\epsilon_{abc}\sigma^b{}_dE^{cd}$ [EvE98, (45)–(46)];
- in the frame form, also torsion-freeness and the curvature (structure) equations.

Friedrich–Rendall derive the subsidiary system for the zero quantities in this gauge and state its
characteristic polynomial: the factors $\xi_0^K(\xi_0^2 + \tfrac14h^{ab}\xi_a\xi_b)^L$, i.e. speeds 0
and $\tfrac12$. It is symmetric hyperbolic with right-hand sides $N(z)$ linear in the zero
quantities. So constraints satisfied initially remain satisfied on the domain of dependence. They
note that boundaries are not covered.

*Check here:* along the BL E–H jets (rational point), $H = \mathrm{curl}\,K$ through $t^6$ and the
Gauss constraint through $t^5$ hold **exactly**.

### 2.3 Principal symbol: which form is hyperbolic

| Form | Unknowns | Principal part | Result |
| --- | --- | --- | --- |
| frame [FR00] | $e_a{}^i$, $\Gamma$, $K$, $E$, $H$ | only $(E, H)$: $P(\xi) = \begin{pmatrix}0 & C_\xi\\ -C_\xi & 0\end{pmatrix}$ | $C_\xi$ antisymmetric on STF, so symmetric hyperbolic; speeds $0, 0, \pm\tfrac12, \pm\tfrac12, \pm1, \pm1$ (computed) |
| coordinate (this note's jet engine) | $\gamma$, $K$, $E$, $H$ | as above **plus** $\partial\gamma$ entering through the Christoffels in $\mathrm{curl}$ | exact example $E = \mathrm{diag}(-1,-1,2)$, $H = 0$: along $\xi = (3/5,0,4/5)$ eigenvalue 0 has algebraic multiplicity 16, geometric 15 (Jordan block), so **not strongly hyperbolic**; along $\xi = e_3$ the symbol is diagonalisable (16/16). No statement for all states or all directions is made. *(Corrected in research/curl_hierarchy_convergence_sprint_01.md §2.4.)* |

`python -m research.calc.bianchi_symbol` computes the table.

**Correction to note 17.** "Symmetric hyperbolic in this same chart" is true for Friedrich's frame
unknowns. There the frame components $e_a{}^i$ live in the comoving chart, and the connection
coefficients are *independent* unknowns transported by ODEs:
$e_0(\Gamma^a{}_{bc}) = R^a{}_{b0c} - K_c{}^d\Gamma^a{}_{bd}$, whose source is the Weyl tensor. It
is false for the natural coordinate implementation, in which the Christoffels are derivatives of
$\gamma$.

The jet engine below uses the coordinate form. That is harmless for formal power series, where
Cauchy–Kovalevskaya applies, but an energy method must use the frame form.

### 2.4 Initial values from BL data

- $\gamma_0 = \psi^4\delta$, $K = 0$, $E = \mathrm{Ric}(\gamma_0)$, $H = 0$. Here $H(0) = \mathrm{curl}\,K(0) = 0$,
  and $\mathrm{curl}\,E(0) = \mathrm{curl}\,\mathrm{Ric} = $ Cotton–York $= 0$ (conformal flatness).
- Frame form: $e_a{}^i = \psi^{-2}\delta_a^i$, with $\Gamma$ from $\psi$.

### 2.5 Reconstruction of the metric

$\gamma^{ij} = \delta^{ab}e_a{}^ie_b{}^j$, in **the same comoving chart and the same proper time**
$t$ as the evaluator. No coordinate change is involved, so events correspond identically.

### 2.6 Comparison with the evaluator's recursion

`research/calc/bl_bianchi.py` solves the system above as $(x, y, z, t)$-jets. At the rational BL
plane point ($d = 3831/5000$) its $\gamma$ coefficients equal those of `bl_metric.series_jets`
**exactly** through $t^8$, and all odd coefficients vanish. The unconverted transcription of §2.1
agrees through $t^6$ and fails at $t^8$. This is the first order where $K\cdot H$ reaches $\gamma$
($H_3K_1$ in $\dot H$ at $t^4$, then $E$ at $t^6$, $K$ at $t^7$, $\gamma$ at $t^8$).

### 2.7 The four questions for this formulation

1. **What becomes easier.**
   - In the frame form, every metric, frame, connection and $K$ variable obeys an ODE along the
     worldlines.
   - The only spatial coupling is a linear-principal-part curl system for $(E, H)$.
   - It is symmetric hyperbolic with speeds $\le 1$, so energy estimates lose no derivatives
     (§5), and domains of dependence are explicit.
2. **Obstruction removed.** The derivative loss of the second-order ADM form in Gauss gauge, which
   is only weakly hyperbolic and limits the method to Cauchy–Kovalevskaya and analytic data.
3. **How the metric is recovered.** From the frame, in the same chart (§2.5).
4. **What remains as hard as before.**
   - The nonlinear algebraic terms (products with $K$) and the geodesic crunch of the gauge itself.
   - Constraint propagation needs a domain of dependence without boundary.
   - Energy estimates still need $L^\infty$ control of $K$ and $\Gamma$, plus higher-order commuted
     energies.

   Nothing here is a solution formula.

## 3 · The $t^6$ law is a theorem

### 3.1 Statement

**Theorem 1.** Let $\gamma_0$ be a conformally flat Riemannian 3-metric with $R = 0$, i.e.
time-symmetric vacuum data. Let $\gamma(t) = \sum_k\gamma_kt^{2k}$ be its geodesic-slicing
development. Then

```math
\begin{aligned}
\Delta\mathrm{Ric} &= 3\,\mathrm{Ric}\cdot\mathrm{Ric} - \lvert\mathrm{Ric}\rvert^2\gamma_0, \qquad
\gamma_2 = \tfrac1{12}\bigl(7\,\mathrm{Ric}\cdot\mathrm{Ric} - \lvert\mathrm{Ric}\rvert^2\gamma_0\bigr),\\
\gamma_3 &= -\tfrac1{15}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\gamma_0
 + \tfrac1{60}\nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric} - \tfrac1{180}\lvert\nabla\mathrm{Ric}\rvert^2\gamma_0
 - \tfrac1{180}R^{kl}\nabla_k\nabla_lR_{ij}.
\end{aligned}
```

Every coefficient is universal: nothing depends on the BL family.

### 3.2 Where the terms come from (structural derivation)

By claims 103–104 and §2.1:

```math
\gamma_3 = \gamma_0M_6^{\rm sil} - \tfrac1{60}(\mathrm{curl}\,H)_3,\qquad
\gamma_0M_6^{\rm sil} = -\tfrac{13}{180}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\gamma_0 .
```

Since $H_1 = -\mathrm{curl}\,E_0 = 0$ and every algebraic term of $\dot H$ contains $H$:

```math
H_3 = -\tfrac13\,[\mathrm{curl}_{\gamma(t)}E(t)]_{t^2}
    = -\tfrac13\bigl(\mathrm{curl}_{\gamma_0}E_2 + (\delta_{\gamma_1}\mathrm{curl})(\mathrm{Ric})\bigr),\qquad
(\mathrm{curl}\,H)_3 = \mathrm{curl}_{\gamma_0}H_3,
```

with $E_2 = -\tfrac52\mathrm{Ric}^2 + \tfrac12\lvert\mathrm{Ric}\rvert^2\gamma_0$ and $\gamma_1 = -\mathrm{Ric}$.

So the entire non-silent content of $\gamma_3$ is a double curl, one of which is varied, of
expressions quadratic in Ric. That is why gradients enter only as $\nabla\mathrm{Ric}\otimes\nabla\mathrm{Ric}$
and $\mathrm{Ric}\cdot\nabla\nabla\mathrm{Ric}$. It is also why the $\mathrm{tr}(\mathrm{Ric}^3)\gamma_0$
coefficient is purely silent: a curl produces no pure-trace term.

The reduction below makes this rigorous without hand-expanding $\mathrm{curl}\,\mathrm{curl}$.

### 3.3 Proof

1. **Normal form.** Conformally flat means locally $\gamma_0 = \phi^4\delta$. Then
   $R = -8\phi^{-5}\Delta\phi = 0$ means $\phi$ is harmonic.
2. **Naturality.** Both sides are natural tensors. The Gauss-gauge development commutes with
   isometries of the data, and the Taylor coefficient $\gamma_k(x)$ depends only on the $2k$-jet
   of the data at $x$. So it suffices to check the identities at $x = 0$ after any isometry.
   - Translations, rotations and special conformal maps $F_b(x) = (x + b\lvert x\rvert^2)/(1 + 2b\cdot x + \lvert b\rvert^2\lvert x\rvert^2)$
     act on $\phi$ by $\phi \mapsto (\phi\circ F)\lvert DF\rvert^{1/2}$. This preserves the form
     $\phi^4\delta$ and harmonicity: the Kelvin weight in dimension 3 is $\tfrac12$.
   - $\partial(\phi\circ F_b\,\lvert DF_b\rvert^{1/2})(0) = \nabla\phi(0) - \phi(0)\,b$, so
     $b = \nabla\phi(0)/\phi(0)$ removes the gradient.
   - A constant rescaling multiplies both sides by the same power, since both have weight 6.

   So we may assume $\phi(0) = 1$, $\nabla\phi(0) = 0$, and a diagonal Hessian.
3. **Polynomial reduction.** Write $\phi = 1 + q_2 + q_3 + \dots$, with $q_k$ harmonic and
   homogeneous of degree $k$, and $q_2 = a(x^2 - z^2) + b(y^2 - z^2)$.
   - The recursion only differentiates and multiplies. Hence every residual of derivative order $w$
     ($w = 4$ for the first two identities, $w = 6$ for $\gamma_3$) is, at $x = 0$, a polynomial in the
     coefficients that is weighted-homogeneous of weight $w$. The coefficients of $q_k$ have weight $k$.
   - So, for $w = 6$, $P = P_6(q_6) + P_{24}(q_2, q_4) + P_{33}(q_3) + P_{222}(a, b)$, with no $q_5$
     term because no weight-1 coefficient is left.
   - For $w = 4$, $P = P_4(q_4) + P_{22}(a, b)$.
4. **Finite check.** $P \equiv 0$ if and only if it vanishes on the following:
   - each basis vector of $q_4$ and of $q_6$ (linear parts);
   - $q_2 = (1, 0), (0, 1)$ combined with each basis vector of $q_4$ (bilinear part);
   - $q_3 = e_i$ and $e_i + e_j$ (polarisation);
   - $(a, b) \in \{(1,0), (0,1), (1,1)\}$ for the binary quadratic;
   - $(a, b) \in \{(1,0), (0,1), (1,1), (1,-1)\}$ for the binary cubic, because a nonzero binary
     cubic has at most three projective zeros.

   That is 72 evaluations. Each one is the exact rational value of the recursion's coefficient and
   of the right-hand side at $0$, computed by `research/calc/cf_universal.py` with the existing engine
   and exact rational harmonic bases. All 72 vanish.

**Controls.**

- Non-harmonic perturbations break the identities: $\phi = 1 + x^4$ breaks the weight-4 identities,
  $1 + x^3$ breaks the $\gamma_3$ law, and $1 + x^2$ breaks all three.
- $\gamma_3 \ne 0$ at 45 of the 72 evaluation points.

$\square$

The evidence category is a proof whose finite part is an exact rational computation, so it is
computer-assisted. It contains no floating point and no sampling. This promotes **claim 98**. With
claims 103–104 it also promotes **claim 105**: the $(\mathrm{curl}\,H)_3$ law.

*Literature (curl-hierarchy sprint §2.5).* Claim 105 is the **known** result of Clifton, Gregoris and
Rosquist (Gen. Rel. Grav. 2017, eq. 59), $\partial_t^3(\mathrm{curl}\,H)\vert_0 = 4\,\mathrm{curl}\,\mathrm{curl}(E^2)$,
under the same assumptions. In Taylor coefficients it reads
$(\mathrm{curl}\,H)_3 = \tfrac23\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$, and it is expanded here with the
conformally flat identities. The equality is checked exactly, for all conformally flat vacuum data.

So $\gamma_3 = -\tfrac{13}{180}\lvert\mathrm{Ric}\rvert^2\mathrm{Ric} - \tfrac{41}{540}\mathrm{tr}(\mathrm{Ric}^3)\gamma_0 - \tfrac1{90}\mathrm{curl}\,\mathrm{curl}(\mathrm{Ric}^2)$.
This metric-level form and its universal verification are this project's; their novelty is not
established.

### 3.4 What it reveals about higher orders

- **No term linear in curvature survives beyond $t^2$.** At linear order around flat space,
  $\mathrm{Ric}_{\rm lin} = -2\,\partial\partial\phi_{\rm lin}$ and $\Delta^m\mathrm{Ric}_{\rm lin} = 0$ for
  $m \ge 1$ ($\phi_{\rm lin}$ is harmonic). Hence $E_{\rm lin}(t) = \mathrm{Ric}_{\rm lin}$,
  $H_{\rm lin} = 0$ and $\gamma_{\rm lin} = \gamma_0 - t^2\mathrm{Ric}$ exactly. Every higher coefficient is
  at least quadratic in curvature. This is a forced cancellation, proved.
- **Conformal flatness delays the coupling by two orders.** $\mathrm{curl}\,E_0 = 0$ makes $H = O(t^3)$
  instead of $O(t)$. This is the base case of the lag lemma (§4.1).
- **Natural gradings.**
  - Derivative weight: $\gamma_k$ has weight $2k$.
  - Curvature degree: between 2 and $k$ for $k \ge 2$.
  - Curl degree (§4).

  Only the last tracks the spatial coupling.

## 4 · The curl-degree hierarchy (main route)

### 4.1 Grading and lag lemma

Give $\varepsilon$ in §2.1 the role of a bookkeeping parameter. The solution is a formal series in
$t$ whose coefficients are polynomials in $\varepsilon$. Changing $\varepsilon \to -\varepsilon$
together with $H \to -H$ leaves the system invariant. Hence $\gamma$, $K$, $E$ are even in
$\varepsilon$ and $H$ is odd.

Write $G_{k,j}$ for the $\varepsilon^{2j}$ part of the $t^{2k}$ coefficient of $\gamma$. This is the
part produced by exactly $2j$ applications of curl.

**Lemma 2 (lag).** For conformally flat time-symmetric data, the $\varepsilon^m$ part ($m \ge 1$) of
the $t^n$ coefficient vanishes unless $n \ge m + L_F$, where $L_E = L_H = 2$, $L_K = 3$, $L_\gamma = 4$.
In particular $G_{k,j} = 0$ for $j > k - 2$, i.e. $\gamma_k$ has degree at most $k - 2$ in
$\varepsilon^2$.

*Proof (induction on $n$).*

- Each $\varepsilon$ comes with one curl. A curl in $\dot E$ raises the $t$-order by one and the
  $\varepsilon$-degree by one. A product of two fields adds the $\varepsilon$-degrees and the
  $t$-orders, plus one for the integration.
- The lags are preserved by every term once they hold for the inputs, because $K$ starts at $t^1$.
- The only place the base case matters is $\varepsilon\,\mathrm{curl}\,E$ applied to the silent
  ($\varepsilon^0$) part of $E$, in $\dot H$. There one needs $[\mathrm{curl}_{\gamma(t)}E(t)]^{(0)}_{t^0} = 0$
  and $[\cdot]_{t^1} = 0$.
- The first is $\mathrm{curl}\,\mathrm{Ric} = $ Cotton–York $= 0$, which is where conformal flatness
  enters. The second holds because the silent solution is even in $t$.
- Then $H$ has lag 2, $E$ lag 2 (through $\varepsilon\,\mathrm{curl}\,H$), $K$ lag 3, and $\gamma$
  lag 4. $\square$

The bound is sharp at the first orders. Exactly at the rational point: the $\varepsilon^2$ sector
of $\gamma$ first appears at $t^6$ and the $\varepsilon^4$ sector at $t^8$, and $\varepsilon = 0$
reproduces the silent ODE through $t^8$.

### 4.2 Each sector is a linear transport problem

Write $\Phi = (\gamma, K, E, H)$, $\dot\Phi = A(\Phi) + \varepsilon L(\Phi)$, where $A$ is algebraic
and $L = (0, 0, \mathrm{curl}\,H, -\mathrm{curl}\,E)$. Expand $\Phi = \sum_m\varepsilon^m\Phi_m$:

```math
\dot\Phi_0 = A(\Phi_0)\ \ (\text{silent ODE}),\qquad
\dot\Phi_m = DA(\Phi_0)\,\Phi_m + \mathcal N_m(\Phi_1, \dots, \Phi_{m-1}) + [L(\Phi)]_{m-1}.
```

Each sector is the solution of the **linearised silent ODE along each worldline**, driven by a
source made of products and curls of lower sectors. This is the transport equation that organises
the hierarchy.

Every sector can therefore be computed as a function of $t$, by ODE integration, up to the silent
collapse time, with no Taylor truncation in $t$. The cost: sector $m$ needs about $m$ spatial
derivatives of the silent solution, i.e. its variational equations.

*If the $\varepsilon$-series converged at $\varepsilon = 1$ beyond the Taylor radius, this would be
a new constructive representation (category B).* §4.4 tests the first necessary condition.

### 4.3 Answers to the four questions for the hierarchy

1. **Easier:** each sector is an ODE in $t$ (exact in $t$).
2. **Obstruction addressed:** the Taylor truncation in $t$.
3. **Metric:** $\gamma = \sum_j G_j(t, x)$, in the same chart.
4. **As hard as before:** convergence in $j$. The $j$-th source needs about $2j$ spatial
   derivatives, so convergence needs analyticity in $x$, just as Cauchy–Kovalevskaya does.

### 4.4 First decisive experiment: sector sizes at generic BL points

**Run.** `python -m research.calc.bl_bianchi sectors` (about 14 minutes per point, three points in
parallel, record `research/calc/data/bl_bianchi_sectors.json`).

- The E–H jets are run to $t^{12}$ at $\varepsilon = 0, 1, \dots, 5$, plus $\varepsilon = 6$ as a check, in
  256-bit balls.
- $G_{k,j}$ is obtained by exact-in-arithmetic Vandermonde interpolation in $\varepsilon^2$.
- The extrapolation to $\varepsilon = 6$ agrees to $\le 2\times10^{-53}$, so the degree bound of Lemma 2
  holds at these generic points.
- Every entry the lemma forces to vanish is $\le 10^{-69}$.

The table gives Frobenius norms (mixed, ball midpoints, so exploratory) of $G_{6,j}$, the $t^{12}$
coefficient split by curl degree $2j$.

| Point | $j = 0$ (silent) | $j = 1$ | $j = 2$ | $j = 3$ | $j = 4$ |
| --- | --- | --- | --- | --- | --- |
| P1 $(0.17, 0.11, 0.5)$, near a hole | $2.69\times10^{-4}$ | $8.65\times10^{-6}$ | $2.35\times10^{-6}$ | $5.47\times10^{-7}$ | $1.71\times10^{-7}$ |
| P2 $(0.26, -0.09, 0.12)$, waist | $1.40\times10^{-5}$ | $1.62\times10^{-5}$ | $1.17\times10^{-5}$ | $3.96\times10^{-6}$ | $9.68\times10^{-7}$ |
| P3 $(0.6, 0.3, 0.9)$, far | $4.86\times10^{-7}$ | $1.59\times10^{-7}$ | $6.57\times10^{-8}$ | $6.32\times10^{-8}$ | $1.49\times10^{-8}$ |

**Reading.**

1. **Ratios at fixed order.** From $j = 1$ on, successive ratios are 0.27, 0.23, 0.31 (P1) and
   0.72, 0.34, 0.24 (P2). At P3 they are 0.41, 0.96, 0.24, so one step is not monotone.
   - At the coefficient level, the finitely many sectors present at $t^{12}$ decrease past $j = 1$.
     By the lag lemma only $j \le k - 2$ can occur at a fixed order. So this decrease is *not*
     evidence that the series over $j$ converges, and it is not necessary for convergence
     (curl-hierarchy sprint §2.2).
   - At P1, the coefficient of $t^{12}$ in $j = 1$ is 3.2 % of that in the silent part, and
     $j \ge 2$ adds 1.1 %. These are coefficient ratios at one point and one power of $t$, **not
     metric errors**.
   - At the waist the first three sectors are comparable.
2. **Radius of each sector in $t$.** Each sector's coefficient sequence in $k$ starts at
   $k = j + 2$. Its single-step radius estimates *increase* with $k$, for example P3, $j = 1$: 2.74,
   2.84, 3.06. This is onset growth, and with three or four terms per sector it does not determine
   the sector's radius.
3. **Decision for M1.**
   - *Not falsified*: the first necessary condition, decrease in $j$ at fixed order, holds at the
     points tested, with the exceptions noted.
   - *Not confirmed*: whether the sector functions remain analytic up to the local collapse time
     $t_{\rm sil}(x)$, while their sum has the smaller radius of notes/18, cannot be decided from
     Taylor coefficients. It needs the sector **functions** (§4.2), computed as ODE solutions in $t$.

   That is the next decisive task (§8).

## 5 · Alternative route: the energy identity

**Lemma 3 (curvature energy).** For solutions of §2.1 (with $\varepsilon = 1$), with
$w = \lvert E\rvert^2 + \lvert H\rvert^2$ and $F^k = \epsilon_i{}^{kl}E^{ij}H_{lj}$:

```math
\partial_t\bigl(\sqrt{\gamma}\,w\bigr) - 2\,\partial_k\bigl(\sqrt{\gamma}\,F^k\bigr)
 = \sqrt{\gamma}\,\bigl(-\theta\,w + 6\,\mathrm{tr}(\sigma(E^2 + H^2))\bigr),
\qquad \lvert 2F\cdot n\rvert \le w\ \ (\lvert n\rvert_\gamma = 1).
```

*Proof.*

- $\partial_t\lvert E\rvert^2 = -4\,\mathrm{tr}(KE^2) + 2E^{ij}\partial_tE_{ij} = 2E\cdot\mathrm{curl}\,H - 4\theta\lvert E\rvert^2 + 6\,\mathrm{tr}(KE^2)$,
  and likewise for $H$.
- $\partial_t\sqrt\gamma = \theta\sqrt\gamma$.
- $E\cdot\mathrm{curl}\,H - H\cdot\mathrm{curl}\,E = D_k(\epsilon_i{}^{kl}E^{ij}H_{lj})$, by swapping indices
  in $\epsilon$.
- $K = \sigma + \tfrac\theta3\gamma$.
- The flux bound: $2F\cdot n = 2\langle E, C_nH\rangle$ with $\lVert C_n\rVert = 1$, the maximal
  speed, so $\lvert 2F\cdot n\rvert \le 2\lvert E\rvert\lvert H\rvert \le w$. $\square$

*Check:* the identity holds exactly for the BL jets at the rational point, $t^0$–$t^5$. The flux
ratio sampled over $2\times10^4$ random states is at most 0.946.

**Consequence.** On a region shrinking at unit speed,
$\frac{d}{dt}\int w \le \int(-\theta + 6\sigma_{\max})\,w$. That is $L^2$ control of the curvature by
Grönwall, given $L^\infty$ bounds on $K$. This is a standard Bel–Robinson-type estimate, re-derived
in our conventions.

*If it is extended* to commuted energies of order $s > 5/2$, with ODE bounds for the frame and
connection, *it gives a continuation criterion and a lossless a-posteriori estimate.* That is the
step notes/17 §3 lacked. Not done here:

- the commuted energies;
- the size of the domain of dependence of the tube (width 0.02, coordinate speed $\le\psi^{-2}$);
- the restart error from the truncated series.

## 6 · Portfolio and selection

| | Hypothesis (precise statement sought) | Why it might hold | Precedent | New lemma or tool | First decisive calculation | Falsified if | Link to the metric |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **M1 (main)** | The curl-degree expansion $\gamma = \sum_j G_j$ converges at $\varepsilon = 1$ on a time interval set by spatial analyticity, beyond the $t$-Taylor radius | near the holes the silent part dominates (notes/18); each sector is an exact ODE in $t$ | gradient and long-wavelength expansions in cosmology; silent models [MPS94, BMP95] | lag lemma (proved), transport form (§4.2) | sector sizes $G_{k,j}$ at generic points (§4.4) | the higher sectors are not smaller than the lower ones at fixed $k$, or a sector's $t$-radius is not above the full series' | a category-B representation evaluable as ODE solutions in $t$ |
| **M2 (alternative)** | A lossless energy estimate for the frame E–H system closes on the tube for times beyond the CK certificate | symmetric hyperbolic, speeds $\le 1$ [FR00] | Bel–Robinson / Christodoulou–Klainerman | Lemma 3 (proved); commuted energies (open) | the $L^2$ identity (done); next, the $H^s$ commutator constants | the domain of dependence of the tube is too small, or the constants explode | certified continuation of the existing B-class metric |
| M3 | A closed form for each sector at fixed $j$ | each sector is a linear ODE with explicit silent coefficients | variation of constants | — | does $G_{\cdot,1}$ have a closed form in $t$? | non-elementary at $j = 1$ | piecewise A-class building blocks |
| M4 | Different data (Misner, unequal mass) simplify the sectors | Cotton $= 0$ is what delays coupling; other conformally flat data keep it | — | — | — | — | — |

M1 and M2 were selected: one constructive, one analytic. M1's first decisive calculation was run
(§4.4). M2's base lemma was proved (§5).

## 7 · Role of each result

- **Theorem 1** (proved). It enables certified covariant control of the first three nonzero
  orders for any conformally flat vacuum data. It gives an exact local initial layer, valid for all
  such data, against which approximations and restart data can be checked.
- **Lemma 2** (proved). It enables a sharp bookkeeping of where spatial coupling enters. Truncating
  the curl-degree expansion at $J$ is exact through $t^{2J+4}$, with error $O(t^{2J+6})$. It gives a
  hierarchy of transport problems with known error order. *(Session 15 wrote "through $t^{2j+6}$",
  which is wrong. With $J = 0$, the silent solution fails at $t^6$.)*
- **The frame/coordinate hyperbolicity distinction** (computed). An energy method must use frame
  and connection unknowns; it rules out a coordinate implementation.
- **Lemma 3** (proved). It is the base energy of the alternative route; with commuted versions it
  gives a continuation criterion.

## 8 · Physical scope and continuation

- **Scope.** Everything concerns BL data and short times in the geodesic gauge. No merger is
  certified, no horizon statement is new, and no radiation is computed.
- **Before a certified merger:**
  - an explicit $\delta_1$ (notes/15);
  - certified MOTS formation inside a certified metric domain;
  - global (event-horizon) information.

  None was advanced by this note.
- **Continuation note.**
  - Start from §4.4 and §5.
  - Records: `research/calc/data/cf_universal.json`, `bl_bianchi_sectors.json`.
  - Tools: `bl_bianchi.py` (E–H jets with $\varepsilon$, energy check), `cf_universal.py`
    (harmonic-jet proofs), `bianchi_symbol.py` (symbols).

*(Withdrawn in the curl-hierarchy sprint, §2.3. Boundedness of a fixed sector on an interval where
the silent background is regular is automatic: it solves a linear forced ODE with bounded
coefficients. So the test below cannot decide uniform control, convergence at $\varepsilon = 1$, or
passage beyond the Taylor radius. The replacement task is in
research/curl_hierarchy_convergence_sprint_01.md §8.)*

**Next single decisive task (superseded).** Integrate the first curl-degree sector ($j = 1$) as a function of
$t$ at a near-hole point. Use the transport equation of §4.2: solve the silent ODE together with its
spatial 2-jet (variational equations), feed $\mathrm{curl}^2$ of it as the source, and integrate
numerically up to the silent collapse time.

*Exact system for that task.* Truncate at $\varepsilon^2$ and write
$(\gamma, K, E) = (\gamma^{(0)}, K^{(0)}, E^{(0)}) + \varepsilon^2(\gamma^{(2)}, K^{(2)}, E^{(2)})$ and
$H = \varepsilon H^{(1)}$. Then:

```math
\begin{aligned}
&\text{silent:}\ \ \dot\gamma^{(0)} = 2K^{(0)},\ \dot K^{(0)} = -E^{(0)} + K^{(0)}\!\cdot K^{(0)},\
\dot E^{(0)} = -2\theta^{(0)}E^{(0)} + \tfrac52\{K^{(0)}, E^{(0)}\} - \gamma^{(0)}\mathrm{tr}(K^{(0)}E^{(0)}),\\
&\dot H^{(1)} = -\mathrm{curl}_{\gamma^{(0)}}E^{(0)} - 2\theta^{(0)}H^{(1)} + \tfrac52\{K^{(0)}, H^{(1)}\} - \gamma^{(0)}\mathrm{tr}(K^{(0)}H^{(1)}),\\
&\dot\Phi^{(2)} = DA(\Phi^{(0)})\,\Phi^{(2)} + \bigl(0,\ 0,\ \mathrm{curl}_{\gamma^{(0)}}H^{(1)}\bigr),\qquad \Phi^{(2)}(0) = 0,\ H^{(1)}(0) = 0 .
\end{aligned}
```

Here $\{A, B\} = A\cdot B + B\cdot A$, the dot products use $\gamma^{(0)}$, and $DA$ is the
linearisation of the algebraic right-hand side.

- The silent system is pointwise. Carried with spatial jets of order 2 it loses no order.
- $H^{(1)}$ then needs the 2-jet of $E^{(0)}$ and is exact to order 1, so the source of $\Phi^{(2)}$
  is exact at the point.
- So the $j = 1$ sector at a point is an ODE system of fixed finite size in $t$. It can be
  integrated numerically, or expanded in $t$ to high order with spatial jets of weight 2 only.

The test: does the $j = 1$ sector stay bounded until $t_{\rm sil}(x) \approx 2.1$?

- **If yes**, the sectors carry only local singularities, the smaller radius of the full series is
  produced by the sum over $j$, and M1 becomes a convergence question for a hierarchy of explicit
  transport problems.
- **If no**, nonlocal singularities already appear at first order in the coupling, and M1 is
  rejected as a route to a longer domain.

## References

- [EvE98] G. F. R. Ellis and H. van Elst, *Cosmological models*, Cargèse lectures 1998,
  gr-qc/9812046, eqs. (34), (40)–(46). Read in the PDF.
- [FR00] H. Friedrich and A. D. Rendall, *The Cauchy problem for the Einstein equations*, Lect. Notes
  Phys. 540 (2000) 127, gr-qc/0002074, §4.4, pp. 51–57. Read on the rendered page; signature
  $(+,-,-,-)$, converted here.
- [MPS94], [BMP95]: as in notes/18.

## Discussion prompts

- The lag lemma uses only $\mathrm{curl}\,\mathrm{Ric} = 0$. For data with nonzero Cotton–York
  tensor (spin, boosts), the coupling starts at $t^3$ in $H$. Does the hierarchy then lose its
  separation of scales?
- Each curl-degree sector is a linear transport problem along the silent flow. Is there a norm in
  which the $j$-th source is bounded by $C^j\,j!^{\alpha}$ with $\alpha < 1$, giving Gevrey-type
  rather than analytic requirements?
- The frame form is symmetric hyperbolic and the coordinate form is not, although they describe
  the same solution. Which variable in the frame form plays the role of the missing eigenvector?
