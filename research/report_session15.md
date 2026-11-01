# Report: the spatial coupling in Einstein–Bianchi form, session 15

*2026-11-01. Branch `claude/bl-bianchi-formulation`, based on `bbc2368` (PR #5's head). Mathematics:
[`notes/19-bianchi-coupling.md`](../notes/19-bianchi-coupling.md).*

**No category-A metric and no certified merger.** The session produced:

- a verified first-order formulation of the coupling;
- a corrected hyperbolicity statement;
- a theorem replacing a fitted law;
- a proved structural lemma for all orders;
- an exploratory first test of the chosen constructive route.

## Obstacle investigated

The spatial coupling of the gravitational field. In the geodesic development of BL data, all
nonlocality enters through $\mathrm{curl}\,H$ (notes/18). The question was whether that coupling can be
put into a form that is better organised than the local Taylor framework.

## What was obtained

1. **The Einstein–Bianchi system in Gauss gauge, in repo conventions.**
   - It is derived from Ellis–van Elst (40)–(41).
   - It matches Friedrich–Rendall's $P_{ij}$, $Q_{kl}$ after converting their $(+,-,-,-)$ signature.
     Both sources were read on the page.
   - It reproduces the evaluator's recursion **exactly** through $t^8$ at a rational BL point.
   - The constraints ($H = \mathrm{curl}\,K$, Gauss) hold along the solution.
2. **Which form is hyperbolic** (computed).
   - The frame form (frame and connection as unknowns) is symmetric hyperbolic, with speeds
     $0, \pm\tfrac12, \pm1$.
   - The coordinate form with Christoffels from $\gamma$ is **not** strongly hyperbolic: there is a
     Jordan block whenever $E \ne 0$.
   - Note 17's "symmetric hyperbolic in this same chart" is now marked as true only for the frame
     form.
3. **Theorem 1 (proved).** For every conformally flat time-symmetric vacuum data set, the
   $\Delta\mathrm{Ric}$, $t^4$ and $t^6$ laws hold.
   - The proof is a reduction (naturality, Möbius/Kelvin normalisation, weighted homogeneity,
     polarisation) plus 72 exact rational evaluations.
   - Claim 98, previously a fitted law, is proved. So is claim 105.
   - Corollary: no term linear in curvature survives beyond $t^2$.
4. **Lemma 2 (lag, proved).**
   - Grade the system by curl degree, with ε multiplying every curl. The curl-degree-$2j$ sector of
     $\gamma$ is $O(t^{2j+4})$.
   - ε = 0 is exactly the silent ODE.
   - Each sector solves a linear transport ODE along the silent flow.
   - Conformal flatness enters only through curl Ric = 0.
5. **Lemma 3 (proved).** The curvature energy identity with flux bound, checked exactly on the jets.
6. **Experiment (exploratory).** Curl-degree sectors to $t^{12}$ at three generic points:
   - the lag lemma holds;
   - past $j = 1$ the sectors shrink roughly geometrically at fixed order;
   - near a hole, $j = 1$ is 3 % of the silent part;
   - at the waist, the first three sectors are comparable.

## Why it matters for the eventual metric

- The coupling now has a precise form: a symmetric hyperbolic curl system plus ODEs, in the
  evaluator's own chart. This is the formulation in which a lossless continuation (M2) can be
  attempted.
- Theorem 1 gives exact control of the first three nonzero orders for any conformally flat data.
- The lag lemma and transport form (M1) give a hierarchy in which every term is an ODE solution in
  $t$, not a Taylor truncation.

## Proved versus conjectured

| Status | Items |
| --- | --- |
| Proved | E–H system (derivation + check); Theorem 1; corollary on linear terms; Lemma 2; transport form; Lemma 3; flux bound |
| Inherited | symmetric hyperbolicity of the frame form and constraint propagation (Friedrich–Rendall) |
| Computed | Jordan block of the coordinate form |
| Exploratory | sector sizes; per-sector radii |
| Conjectural / open | convergence of the curl-degree expansion; commuted energies; any continuation beyond the CK certificate |

## B- and A-class status

- **B-class.** The certified construction is unchanged (same domains and errors). It gained a
  verified hyperbolic formulation, an organising hierarchy, and an exact covariant initial layer.
- **A-class.** No candidate. The sectors are ODE solutions, not closed forms.

## Before a certified merger

- an explicit $\delta_1$ (notes/15);
- MOTS formation inside a certified domain;
- global horizon information.

None was advanced here.

## Continue or abandon

- **Reason to continue M1.** Near the holes the curl-degree sectors fall quickly. If each sector
  stays analytic to the local collapse time, the hierarchy is a genuinely new representation.
- **Reason to abandon M1.** If the $j = 1$ sector already blows up before the silent collapse
  time, nonlocal singularities appear at first order in the coupling, and the hierarchy offers
  nothing beyond the Taylor series.

## Next single decisive task

Integrate the $j = 1$ sector as an ODE in $t$ at a near-hole point up to the silent collapse time
$\approx 2.1$, and check whether it stays bounded. The source is the silent solution together with
its spatial 2-jet, with curl² of it as the forcing.

## Reproduce

```bash
python -m research.calc.cf_universal
```

```bash
python -m research.calc.bianchi_symbol
```

```bash
python -m research.calc.bl_bianchi sectors
```

```bash
python -m pytest tests/test_bl_bianchi.py
```
