# Report: bottlenecks, a local/nonlocal split, and foliations, session 14

*2026-10-31. Branch `claude/bl-structure-discovery`, based on `c26cb9f` (PR #4's head). Mathematics:
[`notes/18-structure-discovery.md`](../notes/18-structure-discovery.md).*

**There is no structural breakthrough, and no category-A candidate survived.** The session produced
a map of where the difficulty lives, two executed experiments of different kinds, and one certified
no-go of limited scope.

## What new ideas were tested

**H1 (discovery). Is the BL development a local ODE plus a small nonlocal correction?**

The silent system drops $\mathrm{curl}\,H$ from the exact 1+3 equations. It is then a closed ODE
along each comoving geodesic, fixed by the three tidal eigenvalues there. The candidate class,
fixed in advance, was $\gamma = \gamma_{\rm sil} + \rho$, with $\rho$ small and of larger radius.

**H2 (gauge). Is there a slicing with conformally flat leaves?**

Session 12 excluded this only for the geodesic slicing. H2 tests every bounded maximal-type initial
lapse, $N = (c_0 + c_+/r_+ + c_-/r_-)/\psi$. That is the complete set of bounded lapses with
$D^2N = 0$ at $t = 0$, and it includes the geodesic one.

## What calculations decided them

**H1.**

- Exact rational checks at BL plane and axis points:
  - silent = true through $t^4$;
  - the $t^6$ deviation is $-\tfrac1{60}\gamma_0(\mathrm{curl}\,H)_3$, cross-checked against the
    independent electric-Weyl engine.
- The decision itself: true versus silent coefficients through $t^{16}$ at six generic points and
  through $t^{20}$ at two, with the remainder's size and fitted radius (notes/18 §6.4).
  - Near the holes the remainder is ≤ 13 %.
  - Elsewhere it dominates by $t^{16}$.
  - Everywhere its fitted radius is the smallest, and it grows relative to the silent part with
    order.

**H2.**

- The linearised Cotton tensor of the $t^2$ leaf correction is a tensor of conics in $c$.
- In axisymmetry only two of them are independent, and at every single point they have real common
  zeros.
- So the certificate covers $\mathbb{RP}^2$ with Arb boxes for the conics of **two** points.
- Result: no common zero at $d = 0.7662$, $0.767$ and $1.5$.

## Did it change the structure of the problem?

Partly, and only as knowledge of where the difficulty is, not as a new representation.

- **H1 locates the obstruction.**
  - The local (silent) dynamics is exact through $t^4$ and dominant near the holes.
  - The nonlocal coupling ($\mathrm{curl}\,H$, the magnetic Weyl tensor) carries everything else.
  - It has the smallest radius estimate everywhere, so the series' radius is not the local collapse
    time (exploratory).
- **H2 removes a whole natural gauge family** as a route to a one-scalar metric.

Neither improves an existing implementation. The evaluator is unchanged.

## What new metric information can be obtained

- **A pointwise local model with measured error.** Near the holes, $\gamma \approx \gamma_0M_{\rm sil}(t)$,
  where $M_{\rm sil}$ solves a 3-eigenvalue ODE. The error is $O(t^6)$, with known leading term
  $-\tfrac1{60}\gamma_0(\mathrm{curl}\,H)_3t^6$ and a closed covariant law. Through $t^{16}$ the
  relative remainder is ≤ 13 % at the near-hole points tested.
- **Diagnostics.** The silent collapse time at any point (≈ 2.05–2.10 near the holes), and fitted
  radius estimates of the nonlocal remainder: 1.76–3.88, growing away from the holes (exploratory).

## A-class candidates

None survived.

## B-class representation

- **Organising principle:** gained. The development is a local ODE plus a nonlocal source with a
  known first law.
- **Scope and computability:** unchanged. The certified domain and the evaluator are as in session
  13.

## Which hypotheses failed, and within what exact scope

| Hypothesis | Failed within |
| --- | --- |
| H1 global form ("silent + small, smoother remainder"; "radius = local collapse time") | BL, geodesic slicing, $d \approx 0.766$, coefficients to $t^{16}$ at six generic points and to $t^{20}$ at two; numerical and exploratory |
| H2 | the family $N_c$, first order $t^2$, any shift, at $d \in \{0.7662, 0.767, 1.5\}$; certified |

Not covered by these failures:

- H1: other data, other slicings, higher orders.
- H2: arbitrary lapses, slicings not through the time-symmetric slice, $d$-intervals (claim 110).

## Which claims remain conjectural or unreviewed

- The $(\mathrm{curl}\,H)_3$ law (claim 105) and the $t^6$ law it rests on (claim 98): identities
  verified exactly and held out, not derived by hand.
- The high-order comparison (claim 106): exploratory.
- Theorem CK′ (session 13): proved here, not externally reviewed.
- Inherited and not rerun: the stored certificates and Theorems C, F, F′ and H.

## Next single decisive experiment

Locate the singular set $t = \tau(x)$ of the development. Compute the coefficients to $t^{24}$
($N = 12$) at points along a ray from a throat to the far zone. Fit the remainder's radius and read
its sign pattern.

- **Real, smoothly varying $\tau$:** a singularity reached before the local collapse. The
  $t^{16}$ data suggest about 1.75–1.86 near the holes. A resummation in $s/\tau(x)^2$ is then the
  natural next representation.
- **Complex $\tau$:** continuation must be nonlocal (an energy method, or analytic continuation in
  $x$).

**A correction made during the session.** An earlier draft read a far-zone radius of 2.4 off one
coefficient ratio. The $t^{20}$ run at P3 showed that the ratios oscillate there. The note now uses
least-squares fits over several orders.

## Rerun versus inherited

- **Rerun:**
  - the structure tests and one certified evaluator call (baseline);
  - all new calculations;
  - the new tests (8) and the full fast suite.
- **Inherited:** all certificates; Theorems C, F, F′ and H; Corollary J.

## Reproduce

```bash
python -m research.calc.bl_silent
```

```bash
python -m research.calc.bl_foliation
```

```bash
python -m pytest tests/test_bl_silent.py tests/test_bl_foliation.py
```
