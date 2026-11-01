# Report: corrected and extended certified metric, session 13

*2026-10-31. Branch `claude/bl-metric-extension`, based on `aabde39`. Mathematics:
[`notes/17-certified-extension.md`](../notes/17-certified-extension.md).*

## What changed in the metric

The representation is unchanged: the geodesic-slicing development of equal-mass BL data,
$g = -dt^2 + \gamma_{ij}dx^idx^j$, as a convergent series with certified tails. What changed is
where it is certified, and how accurate it is there.

| | Certified time | Region | Status |
| --- | --- | --- | --- |
| session 11 (note 15) | $2.988\times10^{-4}$ claimed | 0.02-tube | **invalid**: misread constant; correct value $2.885\times10^{-4}$ |
| session 12 (note 16) | $5.557\times10^{-4}$ claimed | tube | **invalid**; correct value $5.115\times10^{-4}$ |
| session 13, `ext` | $8.683\times10^{-4}$ | whole tube $U$ | certified (×1.70 over the corrected baseline) |
| session 13, `waist` | $1.545\times10^{-3}$ | waist region $U_w$ ($z \le 0.15$) | certified (×3.02) |

In every row $d \in [d_c - 10^{-13}, 0.767578125]$. Units are $M_{\rm ADM}$; the arithmetic is
128-bit Arb balls with outward rounding.

**Errors at $t = T_0/2 = 2.78\times10^{-4}$** (generic point $(0.17, 0.11, 0.5)$, $d = 0.7662$):

| | $g$ | $\partial g$ | $\partial^2 g$ |
| --- | --- | --- | --- |
| before (corrected), order 5 | $9.3\times10^{-6}$ | $2.8\times10^{-3}$ | $17.5$ |
| after, order 7 | $2.8\times10^{-10}$ | $2.4\times10^{-7}$ | $2.6\times10^{-3}$ |

At a waist point after the change: $\partial^2g$ to $2.9\times10^{-7}$. The full table, with the
axis point and a box, is in notes/17 §4.

## How

1. **Audit.** Baouendi–Goulaouic's condition (12) was read on the page: $a < R/(16RC + 8M)$. A
   NaN-skipping hole was found and closed: `x ** 2` on balls with zero midpoint.
2. **Theorem CK′, proved.** BG's proof with a weight exponent $\beta$, and the Lipschitz bound split
   into a derivative-loss part and an output-level zero-order part.
3. **Coverings.** A symmetry-reduced imaginary quadrant with sub-boxes, run in parallel, with
   parameters searched.
4. **Evaluator.** It combines certificates and CK′ schedules, adds a Cauchy-in-$t$ bound for
   $\partial_t\partial_x\gamma$, and caches jets across times.

The energy-estimate continuation was inspected and not executed. The obstacles are weak
hyperbolicity of fixed-lapse ADM, and the $H^s$ residual and bootstrap requirements (notes/17 §3).

## B → A

**No category-A expression was found.** Two bounded experiments:

- **Coefficient structure.** For conformally flat time-symmetric vacuum data:
  - $\Delta\mathrm{Ric} = 3\mathrm{Ric}^2 - \lvert\mathrm{Ric}\rvert^2\gamma_0$, so $\gamma_2$ is algebraic in Ric;
  - $\gamma_3$ obeys a five-term covariant law (fitted, held-out, exact at 3 rational points);
  - pointwise-algebraic dependence on Ric is **excluded at $t^6$** (certified).
- **Weyl variables.**
  - $E(t) = E_0 - t^2(\tfrac32E_0^2 - \tfrac12\mathrm{tr}E_0^2)$;
  - the speciality $27J^2/I^3$ drifts at $t^2$;
  - the tidal eigenframe rotates at $t^4$ (certified). So every closed isotropic law $E' = f(E)$
    is excluded.

Excluded classes in total for BL:

- algebraically special (session 10);
- conformally flat geodesic slices and silent flow (session 12);
- pointwise algebraic in Ric, and isotropic tidal laws (session 13).

Nothing more general is excluded.

## Physical scope

- **A certified local development**, not a merger proof. No merger is certified at any explicit
  $d$: $\delta_1$ is still not explicit.
- **The single-hole control** shows that the series' convergence time is the time to the $r = 0$
  curvature singularity (corrected wording). For BL that time is not interpreted.

## Rerun versus inherited

- **Rerun:** all four certificates; the comparison; the demo; the structure experiments; the new
  and existing tests; the full fast suite. The extension certificates were also replayed from a
  clean process at the committed source.
- **Inherited:** Theorems C, F and F′; Proposition B; Corollary J; the embedding into the maximal
  development (cited).

## Next decisive calculation

One validated energy-estimate step for the Einstein–Bianchi system in Gauss gauge. It is symmetric
hyperbolic in the present chart, and needs a certified $H^4$ defect bound of the series and
explicit constants.
