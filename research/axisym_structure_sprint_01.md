# Axisymmetric structure sprint 01

*2026-11-01. Stages 0 → 1 → 2 → A1 → A2 → A3 only. Baseline: `main` at `90f023b`, after merging PRs
#1–#6. Branch `claude/black-hole-merger-solutions-rpa9hj`. Checker:
`research/calc/axisym_sprint01.py`.*

**Verdict.** The one selected hypothesis (Killing-frame alignment of the first nonlocal source) is
**falsified within its stated scope**. The evidence is a certified ball excluding zero at three generic
off-axis, off-plane points, at two separations. Two other candidates were rejected at selection, because
results already in the repository exclude them.

## 1 · Baseline and files used

- **Commit.** `90f023b` (`main`). It contains sessions 1–15.
- **Files read:**
  - `research/report_session13.md`, `research/report_session14.md`;
  - `notes/17-certified-extension.md`, `notes/18-structure-discovery.md`, `notes/19-bianchi-coupling.md`;
  - `research/claim_ledger.md` (claims 68–121).
- **Code used:**
  - `research/calc/bl_silent.py` (`remainder_chain`, `curl_h3_law`);
  - `research/calc/bl_structure.py` (`curvature_values`, `t6_basis`);
  - `research/calc/bl_metric.py` (jets, `working_precision`).

## 2 · Target (stage 0) and physical case (stage 1)

**Target of the programme.**

- A metric of 4D vacuum GR, with $\Lambda = 0$ and an asymptotically flat exterior.
- Uncharged, nonspinning black holes in a head-on, axisymmetric configuration, with finite positive
  mass ratio.
- It must cover a meaningful exterior region before, during and after merger, with stated
  initial-data and radiation assumptions.

**Outcome classes.**

| Class | Meaning |
| --- | --- |
| A | closed form or a genuinely simplifying analytic representation |
| B | convergent, effectively computable representation with rigorous error control over the required region |
| C | local existence or a local convergent expansion |
| D | finite-order approximation |
| E | horizon or geometric theorem without a metric |

The repository's best result is **C**: a certified local series ($T \le 1.5\times10^{-3}M$), plus
E-type horizon theorems. This sprint is a structural test, not a construction.

**Case tested.** Equal-mass Brill–Lindquist data.

- **Data.**
  - $\gamma_0 = \psi^4\delta$, with $\psi = 1 + \tfrac1{4r_+} + \tfrac1{4r_-}$ and bare masses
    $\tfrac12$ at $z = \pm d/2$.
  - $K_{ij} = 0$ (time symmetric).
  - Units $G = c = M_{\rm ADM} = 1$, with $\psi \to 1$ at infinity.
- **Constraints.** $R(\gamma_0) = 0$ (since $\Delta\psi = 0$), and the momentum constraint holds
  trivially.
- **Gauge.**
  - Geodesic slicing: lapse 1, shift 0, comoving isotropic Cartesian chart.
  - $g = -dt^2 + \gamma_{ij}(t, x)dx^idx^j$, with $K_{ij} = \tfrac12\partial_t\gamma_{ij}$.
  - Signature $(-,+,+,+)$, MTW Riemann sign.
- **Symmetry.**
  - Axisymmetry about $z$.
  - Twist-free: $\partial_\phi$ is hypersurface-orthogonal, which is preserved from time-symmetric
    data.
  - Reflection $z \to -z$.
- **Radiation fixed by the data.** Conformal flatness plus time symmetry fix the initial
  gravitational-wave content to that of the conformally flat choice, the "junk" content. At $t = 0$:
  $H = 0$, $E = \mathrm{Ric}(\gamma_0)$, and $\mathrm{curl}\,E = 0$.

Masses and $d$ alone do not determine the spacetime; the data family does.

## 3 · Evidence and literature audit (stage 2)

### Repository results

| Claim | Assumptions | Evidence | Status | Relevance |
| --- | --- | --- | --- | --- |
| 68–69 Petrov type I on an open dense set, so algebraically special ansätze are excluded | BL, any $d > 0$ | closed-form eigenvalue differences on $z = 0$; exact check | proved (small) | excludes Kerr–Schild-type routes for *this* spacetime only |
| 87 geodesic leaves not conformally flat | BL, geodesic slicing, $O(t^2)$ | ball excluding 0 | certified | gauge-specific |
| 108 no bounded maximal-type lapse gives conformally flat leaves | BL, $d \in \{0.7662, 0.767, 1.5\}$, $O(t^2)$ | box covering | certified | still a family, not all slicings |
| 88, 103 silent flow exact through $t^4$; first deviation at $t^6$ from the single source $(\mathrm{curl}\,H)_3$ | conformally flat time-symmetric data, geodesic slicing | derivation + exact rational checks | proved | the nonlocal coupling to be organised |
| 105 / 116 $(\mathrm{curl}\,H)_3$ law; $t^6$ law (Theorem 1) | all conformally flat time-symmetric vacuum data | reduction + 72 exact evaluations | proved (computer-assisted, exact) | used in A3 as an independent check |
| 99 "not pointwise algebraic in Ric" | polynomial (constant-coefficient) basis | ball minors | certified, **narrow scope** | see the scope note below |
| 100 tidal eigenframe rotates at $t^4$, $[E_0, E_2] \ne 0$ | BL, generic points | balls | certified | — |
| 106, 120 silent dominance near the holes; curl-degree sector sizes | BL, finitely many points and orders | ball midpoints | exploratory numerics | — |
| 95 certified domain | $d \in I$, tube $U$, $\lvert t\rvert < 8.683\times10^{-4}$ (waist $1.545\times10^{-3}$) | CK′ certificate | computer-assisted, inherited (not rerun) | all local; no merger time reached |
| session 14 axis result | BL axis point $(0, 0, 31/50)$ | exact rationals: $E_4 \ne E_4^{\rm sil}$ | exact | used in the screening below |

**Scope note produced in this sprint (exact algebra plus inherited certificate).** Claim 99 excludes
$\gamma_3$ only from a constant-coefficient basis. It can be strengthened:

- Since $\gamma_0^{-1}\gamma_3 = M_6^{\rm sil} - c_3/60$ and $E_2 = E_2^{\rm sil} + c_3/4$, with the
  silent parts commuting with $E_0 = \mathrm{Ric}$, we get
  $[\mathrm{Ric}, \gamma_0^{-1}\gamma_3] = -\tfrac1{15}[E_0, E_2]$.
- This is nonzero by claim 100.
- So $\gamma_3 \notin \mathrm{span}\{\gamma_0, \mathrm{Ric}, \mathrm{Ric}^2\}$ even with **arbitrary
  invariant-dependent scalar coefficients**, at those points.

### Literature (bounded search; read in the PDF text)

- **Brink–Zimmerman–Hinderer**, arXiv:1303.1919, eqs. (2.3)–(2.10). Geroch reduction with
  $\bar h = g - \lambda^{-1}\xi\xi$.
  - Twist-free vacuum: $\bar R_{\alpha\beta} = \bar D_\alpha\bar D_\beta\lambda/2\lambda - \bar D_\alpha\lambda\bar D_\beta\lambda/4\lambda^2$
    and $\bar D^2\lambda = \lvert\bar D\lambda\rvert^2/2\lambda$.
  - With $\Lambda = \sqrt\lambda$: $\mathrm{Ric}(\bar h) = \Lambda^{-1}\mathrm{Hess}\,\Lambda$ and $\Box_{\bar h}\Lambda = 0$.
  - They also build a null-triad hierarchy. A single scalar ($\lambda$) carries the dynamics.
- **Gudapati**, arXiv:1702.07950. The reduced system is 2+1 Einstein–wave-map. For the reduction of
  asymptotically flat 3+1 spacetimes, the 2+1 mass and the wave-map energy **diverge**, so finite-energy
  2+1 theorems cannot be imported without changes. This was checked; none was imported.
- **Clifton–Gregoris–Rosquist**, arXiv:1607.00775, eqs. (17)–(24) and §4.2–4.3. The 1+3 system has the
  same form as ours.
  - On a reflection-symmetric surface the axial-tensor and scalar parts of $H$ vanish and only its
    vector part survives. $\mathrm{curl}\,H$ (polar) is not forced to vanish.
  - On intersection curves $H = 0$, but $\mathrm{curl}\,H$ need not vanish.
  - In twist-free axisymmetry *every* meridian plane is such a surface. Hence $H = e_\phi \odot w$, and
    $S := \mathrm{curl}\,H$ has only $(\rho\rho, \rho z, zz, \phi\phi)$ components. This is known
    structure.

"Not found in the sources checked" is the only novelty claim made below.

## 4 · Axisymmetric formulation (A1)

**Standard (Geroch; BZH).**

- Killing field $\xi = \partial_\phi$, norm $\lambda = \Lambda^2 = g(\xi, \xi)$, twist zero.
- $g = g_3 + \Lambda^2d\phi^2$, where $g_3$ is the orbit-space metric. In our gauge,
  $g_3 = -dt^2 + q_{ab}dx^adx^b$ with $x^a = (\rho, z)$, so geodesic slicing restricts to Gauss gauge
  on the orbit space.
- **Reduced variables:** $U = (q_{ab}, \Lambda)$.
- **Equations:** $\Box_{g_3}\Lambda = 0$ and $\mathrm{Ric}(g_3) = \Lambda^{-1}\mathrm{Hess}_{g_3}\Lambda$.
  The $t$-components give evolution equations; the $(tt)$ and $(ta)$ components are the constraints.
- **Reconstruction:** $g = g_3 + \Lambda^2d\phi^2$.
- **Axis regularity:** $\Lambda = \psi^2\rho\,(1 + O(\rho^2))$ and elementary flatness
  $\lvert\nabla\Lambda\rvert \to 1$ on the axis.
- **Asymptotics:** $\Lambda \sim \rho$, $q \to \delta$.
- **Remaining gauge freedom:** none (the comoving geodesic chart is fixed by the data).
- Twist-free does **not** mean $H = 0$, nor no radiation (Control 2).

**Derived here, used in A2** (a direct computation in the frame $\{e_\rho, e_z, e_\phi\}$, not
independently cross-checked):

- $\ddot\Lambda = -E^\phi{}_\phi\,\Lambda$ exactly: $\xi$ is a Jacobi field along the comoving
  geodesics.
- $(\mathrm{curl}\,H)_{\phi\phi} = \Lambda^2\epsilon_{(2)}^{bc}\partial_b(W_c/\Lambda)$, with $W = H(\cdot, \xi)$.

**Controls** (`python -m research.calc.axisym_sprint01`, SymPy, exact).

1. **Schwarzschild** (isotropic, static, meridian polar). $g_3 = -N^2dt^2 + \psi^4(dr^2 + r^2d\theta^2)$
   and $\Lambda = \psi^2r\sin\theta$ give $\Box\Lambda = 0$ and $\mathrm{Ric} - \Lambda^{-1}\mathrm{Hess}\,\Lambda = 0$
   identically.
2. **Genuine radiation.** Teukolsky's even-parity $m = 0$ linearised quadrupole wave on Minkowski
   (PRD 26, 745 (1982)), with $F(x) = x^7$:
   - linearised Ricci $= 0$ exactly, which also confirms the recalled formulas;
   - no $d\phi$ cross terms, so it is twist-free;
   - $\delta\Lambda/\Lambda \ne 0$;
   - transverse curvature $r\,R^{\hat\theta}{}_{\hat t\hat\theta\hat t} \to 1890\,u\,\sin^2\theta = \tfrac38F^{(6)}(u)\sin^2\theta$
     at fixed $u = t - r$, i.e. $O(1/r)$ with the quadrupole pattern.

   So the twist-free reduction retains a genuine radiative degree of freedom. It is not pure gauge
   (Riemann $\ne 0$) and not a mass perturbation ($\ell = 2$).

## 5 · Screening and the frozen hypothesis (A2)

**Screened and rejected** (already excluded by repository results):

1. **Axial-Jacobi silence.** The axial norm follows the silent ODE, equivalently
   $(\mathrm{curl}\,H)_{\phi\phi} = 0$ (by the derivation above, $W/\Lambda$ closed).
   - On the axis, $S_3 = (\mathrm{curl}\,H)_3$ is LRS, $\mathrm{diag}(a, a, -2a)$, and the session-14
     axis point gives $a \ne 0$ exactly.
   - Regularity then forces $S_{\phi\phi} \to a \ne 0$ near the axis. Excluded.
2. **Nonlocal correction carried only by $\lambda$, or meridian-isotropic.**
   - $S_3$ is trace-free, so "meridian part zero" forces $S_3 = 0$, contradicting claims 103–105.
   - The isotropic form $c(e_\phi e_\phi - \tfrac12q)$ requires $c \to 0$ on the axis, contradicting
     $a \ne 0$. Excluded.

**Frozen hypothesis card (KFA, Killing-frame alignment).** Frozen before the A3 computation.

1. **Statement.** Let $\Lambda_0 = \psi^2\rho$ be the initial axial norm, $n = \nabla\Lambda_0/\lvert\nabla\Lambda_0\rvert$
   (meridian), $m = e_\phi \times n$, and $R(t) := \gamma_0^{-1}(\gamma(t) - \gamma_{\rm sil}(t))$, the
   nonlocal remainder (mixed, comoving). Then $C[U] \equiv 0$, with $C$ as in the formula after this
   card.
2. **Class.** Geodesic development of equal-mass BL data, any $d > 0$, all off-axis points.
3. **Covariance.** Geometric given the Gauss gauge, since the frame comes from the Killing field.
   At $t^6$ the time-$t$ and $t = 0$ frames agree.
4. **Jet orders.** $R$ starts at $t^6$; we use the stored Taylor coefficients $a_k$ of $t^{2k}$.
   $\partial_t^6C\vert_0 = 6!\,a_6$ with $a_6 = -\tfrac1{60}S_3$, which needs the ADM jets to $N = 3$
   (weight 6), already exact in the repo. Lower orders vanish identically (claim 103).
5. **Motivation.**
   - The only radiative scalar of the reduction is $\Lambda$, and the Killing 2-form
     $F = d\xi^\flat = 2\Lambda\,d\Lambda\wedge d\phi$ is a source-free Maxwell field in vacuum. Its
     principal frame at $t = 0$ is $\{n, m, e_\phi\}$.
   - $C_{abcd}\xi^d = \nabla F$ (Killing identity).
   - So one might expect the first nonlocal source to be aligned with $F$.
6. **Simplification if true.** $S_3$, and the first nonlocal correction to the metric, would be two
   scalars in an explicitly known frame instead of a three-component field. That is a candidate
   building block for an effective representation.
7. **Nearest known result.** CGR parity restrictions ($e_\phi$ is an eigenvector — known and
   automatic). KFA adds the second eigenvector $n$; not found in BZH, Gudapati or CGR.
8. **Cheapest test.** $S_3(n, m)$ at a generic point (off the axis, off $z = 0$) in 256-bit balls.
9. **Failure.** A ball for $S_3(n, m)$ that excludes 0.

The statement of item 1, as a formula:

```math
C[U](t) := R(t)(n, m) \equiv 0\quad\text{for all }t \qquad(\text{equivalently, } n\text{ and }e_\phi\text{ are eigenvectors of }R(t)).
```

## 6 · Calculation performed (A3)

**Method.**

- $S_3$ is computed by the exact chain $S_3 = 4\gamma_0(E_4 - E_4^{\rm sil})$ (`bl_silent.remainder_chain`).
- It is cross-checked against the proved covariant law (`curl_h3_law`); they agree to $\le 10^{-65}$.
- $n$ comes from the closed-form $\nabla\psi$.
- Components are taken in the $\gamma_0$-orthonormal frame.
- Arb, 256 bits. The values are midpoints with true ball radii.

| Point $(x, y, z)$ | $d$ | $S_3(n, m)$ | $\lvert S_3\rvert$ | $S_3(n, e_\phi)$ (parity check) | verdict |
| --- | --- | --- | --- | --- | --- |
| $(0.17, 0.11, 0.5)$ | 0.7662 | $0.009403057987$ (rad $6\times10^{-67}$) | 0.01672 | $< 5\times10^{-67}$ | **fails** |
| $(0.31, -0.2, 0.27)$ | 0.7662 | $-0.008589318752$ (rad $9\times10^{-69}$) | 0.04715 | $< 3\times10^{-68}$ | **fails** |
| $(0.4, 0.2, 1.1)$ | 1.5 | $-0.004553665314$ (rad $3\times10^{-68}$) | 0.00719 | $< 2\times10^{-68}$ | **fails** |
| $(0.26, -0.09, 0)$, plane $z = 0$ | 0.767 | $0$ (rad $2\times10^{-70}$) | 0.07857 | 0 | passes by symmetry only |

The first failing order is $t^6$, the first order at which the remainder exists:
$\partial_t^6C\vert_0 = 720\cdot(-\tfrac1{60})S_3(n, m) = -12\,S_3(n, m) \ne 0$.

**The responsible terms.** The $(n, m)$ components of the four terms of the proved law, at
$(0.17, 0.11, 0.5)$:

| Term | $(n, m)$ component |
| --- | --- |
| $-\tfrac13\lvert\mathrm{Ric}\rvert^2\mathrm{Ric}$ | $-0.01342$ |
| $-\nabla_i\mathrm{Ric}\cdot\nabla_j\mathrm{Ric}$ | $+0.00211$ |
| $+\tfrac13\lvert\nabla\mathrm{Ric}\rvert^2\gamma_0$ | $0$ (pure trace) |
| $+\tfrac13R^{kl}\nabla_k\nabla_lR_{ij}$ | $+0.02071$ |

- Already $\mathrm{Ric}(n, m) = 0.0861 \ne 0$: the Killing-norm gradient is not a principal tidal
  direction, even at $t = 0$.
- The algebraic term inherits that misalignment, and the two gradient terms add their own.
- There is no cancellation. The failure is generic, not a fine-tuned residue.
- These are the known terms of the $t^6$ law, so no new coupling term is claimed.

## 7 · Verdict and scope

- **Falsified within the stated scope.** KFA fails at $t^6$ for the geodesic development of
  equal-mass BL data at $d = 0.7662$ and $1.5$, at three generic off-axis, off-plane points. This is
  certified (Arb balls excluding 0).
- It holds trivially on $z = 0$ and on the axis, by symmetry, which is exactly the kind of hiding the
  generic points avoid.
- This does not exclude other frames, other data, other gauges, or the exact-merger programme.
- Candidates 1–2 were **rejected at selection**: they are inconsistent with existing exact results,
  through axis regularity and tracelessness.

## 8 · What changed in our understanding

- In twist-free axisymmetry the first nonlocal source of the BL development has **no
  symmetry-adapted eigenframe beyond $e_\phi$**:
  - not the tidal frame (claim 100);
  - not the Killing-norm frame (this sprint).

  Its meridian orientation is a genuine third function.
- On the axis it is nonzero and LRS. That excludes, at once, every candidate that removes the
  coupling from the $\lambda$ sector or from the meridian sector.
- Claim 99's exclusion extends to arbitrary invariant-dependent coefficients (§3).
- The twist-free reduction retains radiation (Control 2). Restricting to it is not the reason the
  coupling is hard.

## 9 · One recommended next step

Verify, exactly on the BL jets, the identity that the magnetic Weyl 1-form is determined by the
axial norm:

```math
H(\cdot, e_\phi) = \Lambda^{-1}\star_{(2)}\bigl(d\dot\Lambda - k\cdot d\Lambda\bigr),
```

with $k$ the meridian extrinsic curvature. This is derived from $H = \mathrm{curl}\,K$ and the momentum
constraint; the derivation is sketched, not checked.

If it holds, $\mathrm{curl}\,H$, the entire nonlocal source, is a second-order operator acting on the
reduced radiative scalar $\Lambda$ (and $k$). The spatial coupling of the 4D problem would then be
expressed exactly through the reduced wave equation $\Box_{g_3}\Lambda = 0$: an explicit link between
the session 14–15 hierarchy and the Geroch reduction.

## Reproduce

```bash
python -m research.calc.axisym_sprint01
```

This takes about 2 minutes. It prints the controls and the A3 table.
