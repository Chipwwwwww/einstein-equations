# 14 · Toward an explicit dynamical binary: screening the routes

The question of this note: can one write down, explicitly, a dynamical spacetime containing two
black holes that attract, radiate and merge? The setting is 4D vacuum, $R_{\mu\nu} = 0$,
$\Lambda = 0$, asymptotically flat (AF), $G = c = 1$. The first target is the simplest
nontrivial case: nonspinning holes, axisymmetric, head-on.

**Answer of this session.**

- No explicit (closed-form) binary metric was found. None is known in the literature checked.
- The session gives two obstructions, each with an exact scope. One is proved here; the other is
  read from the literature.
- It gives one exact but local construction: a convergent series.
- It identifies one route, near-critical Brill–Lindquist (BL), along which an exact convergent
  representation of a spacetime *containing an apparent-horizon merger* appears within reach.

The route map with all ten routes, A–J, is
[`research/route_map.md`](../research/route_map.md).

## 1 · What counts as a result

| Category | Meaning |
| --- | --- |
| A explicit exact metric | known functions or quadratures, all field equations verified on a stated domain |
| B convergent constructive representation | a series or iteration with convergence, regularity and the full equations established |
| C existence / geometric theorem | proof without an evaluable metric |
| D controlled approximation | truncation with a stated domain and error assessment |
| E formal or numerical candidate | no convergence or proof |

Levels of what a spacetime covers:

1. exact data on one slice;
2. a local dynamical spacetime;
3. a substantial pre-merger interval;
4. a rigorously established merger, with the horizon notion stated;
5. broad binary parameters.

The repository so far:

| Level | Result | Category |
| --- | --- | --- |
| 1 | BL data (classical) | A |
| 4, apparent horizons, near-critical BL | Theorem H, note 13 | C, computer-assisted, with a non-explicit $\delta_0$ |
| 4, event horizons, extreme mass ratio | Theorems D/E, note 11 | C, conditional on our reading of Hintz 2026 |

## 2 · Literature baseline for explicit binaries

| Known result | What it gives | How it was checked |
| --- | --- | --- |
| Kastor–Traschen (1993) | exact merging black holes, but $\Lambda > 0$ and charged | reproduced, note 08 |
| Bach–Weyl (1922), Israel–Khan (1964) | static two-hole vacuum metrics | need a strut (conical singularity) on the axis; secondary sources |
| Neugebauer–Hennig (2009), Hennig–Neugebauer (2011) | no stationary axisymmetric two-hole vacuum equilibrium | abstracts (literature.md) |
| Bičák–Schmidt (1984, 1989); Pravda–Pravdová review (2000) | axisymmetric, locally AF, radiative vacuum with a second Killing field ⇒ boost-rotation symmetric; acceleration needs nodal singularities, or negative masses | review section read (gr-qc/0003067); Bičák review §11 read (gr-qc/0004016); 1984 paper not read |
| Kerr–Schild (1965), Gürses–Gürsey (1975), Xanthopoulos | vacuum KS: $k$ geodesic and shear-free, field equations linear, algebraically special | abstracts and citing papers |
| superposed KS data | not vacuum; used as approximate data in numerical relativity | arXiv:2405.03289, gr-qc/0305071, gr-qc/0205060 (abstracts) |
| Wessels (1998) | exact analytic characteristic *data* for axisymmetric non-rotating vacuum | abstract, CQG 15 2509 |
| Emparan–Martínez (2016), Hintz (2024–26) | exact EMR event horizon; EMR gluing (existence) | read in sessions 4–9 |

A broad search found no claimed exact, regular, AF, dynamical binary-black-hole vacuum solution.
It used extended web search plus the citations above.

## 3 · Result P: the head-on BL spacetime is not algebraically special

> **Proposition P (exact).** For every separation $d > 0$, the vacuum development of equal-mass
> time-symmetric BL data is of Petrov type I on an open dense subset of the initial slice,
> and hence on an open spacetime neighbourhood of it. On the symmetry axis it is type D.

*Proof.* At a moment of time symmetry the magnetic Weyl part vanishes, and the electric part is
the slice's Ricci tensor, $E_{ij} = R_{ij}(\gamma)$. So the Weyl tensor is algebraically special
exactly where the real, symmetric, traceless $E$ has a repeated eigenvalue.

On the plane $z = 0$, with $r$ the distance to either puncture (bare masses $\tfrac12$), $E$ is
diagonal in the radial, azimuthal and axial directions. The eigenvalue differences are

```math
\lambda_\rho - \lambda_\varphi = -\frac{3(2r - d)(2r + d)}{r^3(2r+1)^2},\qquad
\lambda_\rho - \lambda_z = -\frac{3(8r^3 - 4d^2 r - d^2)}{2r^4(2r+1)^2},\qquad
\lambda_\varphi - \lambda_z = \frac{3d^2}{2r^4(2r+1)} .
```

These are coordinate components; the eigenvalues with respect to $\gamma$ are these divided by
$\psi^4$. They are distinct except at the origin and on one circle. The discriminant of $E$ is
real-analytic away from the punctures and not identically zero, so its zero set has empty
interior. Type I at a point is an open condition. At $d = 0.78$ and $(x, y, z) = (0.8, 0, 0)$ the
eigenvalues are distinct rationals; exact arithmetic is in `research/calc/petrov_bl.py` and
`tests/test_petrov_bl.py`. ∎

**Check.** At $d = 0$ (one Schwarzschild hole), $\lambda_\varphi = \lambda_z = -\lambda_\rho/2$:
type D, as it must be.

**Consequence.** Every exact-solution family built on a shear-free geodesic null congruence is
algebraically special (Goldberg–Sachs). That includes:

- Kerr–Schild with one null vector, and Kerr–Schild–Kundt;
- Robinson–Trautman and Kundt;
- twisting type II/D.

No such family reproduces the head-on BL spacetime on any open set meeting the type-I region.
This is a no-go for *that* spacetime, not for every head-on binary. But it shows that the most
popular exact-solution technology cannot describe the simplest binary whose data are known
exactly.

## 4 · Result G: integrable methods need a symmetry a regular binary does not have

Ernst reduction, inverse scattering (Belinski–Zakharov) and Bäcklund transformations need two
commuting Killing fields. By §2, an axisymmetric, locally AF, radiative vacuum spacetime with a
second Killing field forming a 2D Lie algebra with the axial one is boost-rotation symmetric.
Its sources are uniformly accelerated: by nodal singularities (strings or struts), by negative
masses, or by external fields that are not AF.

A head-on collision of two positive masses with a regular axis therefore has exactly one
continuous symmetry. In that sector the vacuum equations reduce to 2+1 gravity coupled to one
wave equation (Geroch reduction), with no known integrable structure. This reasoning comes from
the literature (category: established), with the caveat that we read the 1984 theorem only
through the reviews.

## 5 · Result S: the BL head-on spacetime as a convergent series

Take geodesic slicing (unit lapse, zero shift) and coordinates comoving with the observers at
rest at $t = 0$. With $K_{ij} = \tfrac12\partial_t\gamma_{ij}$, the vacuum equations are

```math
\partial_t\gamma_{ij} = 2K_{ij},\qquad
\partial_t K_{ij} = -R_{ij}(\gamma) + 2K_{ik}K^k{}_j - K\,K_{ij},
```

plus the constraints. Time symmetry makes $\gamma$ even in $t$ and $K$ odd. Write
$\gamma = \sum_k g_k\,t^{2k}$ and $K = t\sum_k \kappa_k\,t^{2k}$. Then

```math
g_{k+1} = \frac{\kappa_k}{k+1},\qquad
(2k+1)\,\kappa_k = -[\mathrm{Ric}(\gamma)]_k + \big[2\kappa\gamma^{-1}\kappa - \mathrm{tr}(\gamma^{-1}\kappa)\,\kappa\big]_{k-1},
```

starting from $g_0 = \psi^4\delta$, so that $g_1 = -\mathrm{Ric}(g_0)$. Every coefficient is an
explicit differential expression in $\psi = 1 + \frac{1}{4r_+} + \frac{1}{4r_-}$.

> **Proposition S (standard theory applied; category B, level 2).** Near every point of the
> initial slice the series converges for $|t| < T(x)$, with $T(x) > 0$. Its sum solves the full
> vacuum Einstein equations there.

*Proof.* Written second order in time,
$\partial_t^2\gamma = F(\gamma, \partial\gamma, \partial^2\gamma, \partial_t\gamma)$ is in
Kovalevskaya form with analytic $F$. The BL data are real-analytic away from the punctures. By
Cauchy–Kovalevskaya there is a unique analytic solution near each point, and its Taylor series
is the recursion above. The constraint quantities obey a linear homogeneous Kovalevskaya-type
system (contracted Bianchi identities) and vanish at $t = 0$, so they vanish identically. ∎

**Validation** (`research/calc/bl_series.py`, `tests/test_bl_series.py`). The code is an
independent jet engine: truncated power series in $(\rho, z, t^2)$.

- For one centre it reproduces the exact Schwarzschild free fall from rest,
  $\ddot r = -M/r^2$, through $t^{10}$, to $10^{-12}$.
- The Hamiltonian constraint holds order by order to rounding ($\le 10^{-13}$).
- $g_1$ agrees with the independent Cartesian formula of Proposition P.

**Numbers** (floating point, numerical evidence only;
[`research/calc/data/bl_series_dc.json`](../research/calc/data/bl_series_dc.json)). At $d = d_c$
and $(\rho, z) = (0.25, 0)$, near where the critical surface crosses the plane, the coefficients
through $t^{14}$ give ratios $|g_k/g_{k+1}| \approx 7.5$ for the $zz$ component. That suggests a
convergence radius $T \approx 2.7$ there. For comparison, Schwarzschild at the throat has
$T = \pi M$, set by the crunch of geodesic slicing.

**A structural fact.** The magnetic part of the Weyl tensor is

```math
B_{ij} = -t\,\mathrm{curl}\,R_{ij}(\gamma_0) + O(t^3) = -t\,C_{ij}(\gamma_0) + O(t^3),
```

where $C$ is the Cotton–York tensor. That vanishes for conformally flat $\gamma_0$, so
$B = O(t^3)$. The $t^3$ coefficient is nonzero for BL, while for Schwarzschild $B \equiv 0$
(a control). This was checked numerically; the identity at order $t$ is exact.

**What S is and is not.**

- It is an exact representation, local in time: it covers $|t| < T$, and the geodesic slicing
  ends at its crunch.
- It is not a closed form.
- It covers neither a substantial pre-merger interval nor the late-time radiation at
  $\mathscr I^+$.

## 6 · Route J: an explicit representation of a merging region

Combine S with Theorem H ([note 13](13-head-on-collision.md)). For $d$ just above $d_c$, a
common apparent horizon appears at $t_*(d) = c\,(d - d_c) + O((d - d_c)^2)$, with
$c \in [2.3768, 2.3936]$, near the critical surface $S_c$. The radius $T$ depends analytically on
$d$, so it is bounded below near $S_c$ for $d$ near $d_c$. Hence, for $d - d_c$ small:

> the whole apparent-horizon merger region ($0 \le t \le t_*$, near $S_c$) is covered by the
> convergent series of S.

That is a category-B representation of a spacetime region containing a level-4 (apparent-horizon)
merger, with a category-C proof of the merger. On the interval of Theorem F′, $t_* \le 0.0033$,
against $T \approx 2.7$: a factor of about 800.

**Missing, in order of difficulty.**

1. A rigorous lower bound on $T$ near $S_c$. The plan: interval (Arb) jets for the recursion, plus
   a Cauchy–Kovalevskaya majorant for the tail.
2. Certified bounds on $\partial_t^2 F$ from the series, which would make Proposition B's
   $\delta_0$ explicit and give a theorem for a specific $d$.
3. Common MOTSs away from $S_c$ for $0 < t < t_*$.
4. Event horizons, which are out of reach for these data.

## 7 · Radiation: what can and cannot be said

- BL data are time symmetric. Whatever radiation they contain is a standing mixture of incoming
  and outgoing waves at $t = 0$, not a clean outgoing signal. It is not specified independently
  of the masses and separation (it is fixed by conformal flatness).
- Proposition S is local in time and in space near each point. It reaches neither null
  infinity nor Bondi news. No waveform is predicted.
- To obtain one, one would need either the global development (with $\mathscr I^+$), or
  characteristic data on an outgoing null cone, where the radiation content is free data
  (route F).
- The early-time fact $B = O(t^3)$ is a local statement about the gravitomagnetic field near the
  slice. It is not a radiation measurement.

## 8 · Status of the claims in this note

| Claim | Category |
| --- | --- |
| P: Petrov type I on an open dense set, every $d > 0$ | analytically proved; exact arithmetic |
| P ⇒ no algebraically special ansatz for the BL development | proved (with Goldberg–Sachs) |
| G: two-Killing integrable methods give accelerated sources (struts or negative mass) | literature (reviews read; 1984 theorem not read) |
| S: convergent series, full vacuum equations | standard theorem applied; code validated |
| $T \approx 2.7$ near $S_c$ | numerical evidence |
| $B = O(t^3)$ for conformally flat time-symmetric data | exact at order $t$; nonzero $t^3$ coefficient numerical |
| J: the merger region lies inside the convergence domain for $d \to d_c$ | proved qualitatively (analytic dependence); explicit bounds missing |
| explicit binary metric (category A) | **not found** |

## Discussion prompts

- Proposition P rules out algebraically special forms for the BL collision. Is there a natural
  *non*-algebraically-special exact-solution technique left that does not need a second Killing
  field?
- The convergence radius of the time series is set by the crunch of geodesic slicing, which is a
  gauge effect. Would maximal slicing give a series that reaches much further in time? What
  would replace Cauchy–Kovalevskaya there, given that the lapse solves an elliptic equation?
- For conformally flat time-symmetric data, the gravitomagnetic field starts at order $t^3$. Is
  there a physical reading of the $t^3$ coefficient, for example as the onset of the radiation
  the data were hiding?
