# Route J report (session 11)

*2026-10-31. Branch `claude/route-j-ck-certificate`, based on research head `afe0bed`.
Mathematics: [`notes/15-route-j-convergence.md`](../notes/15-route-j-convergence.md).*

## Target and outcome

**Target.** A certified, convergent local representation of the near-critical BL development
that covers the formation of the common MOTS. It needs three explicit things:

- a neighbourhood $U \supset S_c$;
- an interval of $d$;
- a time $T_{\min}$.

It should also give an explicit $\delta$ for the birth time $t_*(d)$.

**Outcome.**

- *Achieved* (computer-assisted; candidate new):
  - $U$ is the closed tube of radius $0.02$ about $S_c$;
  - the interval is $d \in [d_c - 10^{-13}, 0.767578125]$;
  - $T_{\min} = 2.988\times10^{-4}\,M$ (Proposition T).

  On that region, the geodesic-slicing development equals its convergent Taylor series, with an
  explicit truncation bound.
- *Qualitative only*: Corollary J. For $d \in (d_c, d_c + \delta_1)$ the local common-MOTS
  formation lies inside this region. $\delta_1$ is not explicit.
- *Not achieved*: an explicit $\delta$ and $C_2$ with
  $\lvert t_*(d) - c(d - d_c)\rvert \le C_2(d - d_c)^2$. The blocking estimates are
  (notes/15 §7.3):
  - a certified lower bound on the second eigenvalue $\lambda_2(L_{S_c})$ of the Jacobi operator,
    over all Fourier modes;
  - certified bounds on $\partial_t^2\Phi$, $\partial_t\partial_d\Phi$, $\partial_t\partial_{z_0}\Phi$,
    $\partial_d^2\Phi$, $\partial_d\partial_{z_0}\Phi$ and $\inf\partial_{z_0}^2\Phi$ for the true
    shooting function on $\Sigma_t$. These need a validated MOTS solver for non-conformally-flat
    data with $K \ne 0$.

## Baseline audit (research head `afe0bed`)

| Category | Items |
| --- | --- |
| established external theorems | Choquet-Bruhat–Geroch (MGHD); Andersson–Metzger 2009 (MOTS existence and stability); AMMS 2009 (MOTT at a degenerate MOTS); Dain–Weinstein–Yamada 2011 (no common surface at large separation); Hintz 2024/2026 (EMR merger spacetimes); *new this session*: abstract CK (Nishida 1977, read via Baouendi–Goulaouic) and Frittelli 1997 (constraint propagation) |
| reproduced identities | curvature-code checks; KT reduction; AMMS identity; $d_c$ and critical area vs Hui–Lin; Schwarzschild free fall from the series; Hamiltonian constraint order by order; *new*: momentum constraint order by order, and the Ricci reduction $A + B = \mathrm{Ric}$ |
| candidate analytic results | Proposition B (Lyapunov–Schmidt fold; smooth dependence cited); Lemma A; Proposition P (Petrov type I); Theorems D/E (our reading of Hintz, self-reviewed); *new*: Corollary J |
| computer-assisted claims | Theorem C (fold, $d_c$, slopes); Theorem F ($d = 0.78$); Theorem F′ ($d \in (d_c, 0.767578]$); KT horizon certificates; individual-horizon barriers; *new*: Proposition T |
| numerical observations | $T \approx 2.7$ (ratio test); KT $\tau_c$; onset slopes to 6 digits; the shooting landscape; medium necks unstable in floating point; nonzero $t^3$ coefficient of $B_{ij}$ |
| unresolved gaps | explicit $\delta_0$; first common apparent horizon for $0 < t < t_*$; event horizons for BL; independent review of the computer-assisted code; $d \in (0.767578, 0.78)$; outside check of Hintz Lemmas 6–7 |

**Stale or contradictory statements reconciled.**

- notes/10, "Open for every other $d$": superseded by Theorem F′.
- status.md, the gap-3 item and the "Next concrete calculation (b)" plan: done.
- Ledger claim 35 (C → partly N), claim 72 (heuristic), claim 74 (its basis was not a proof).
- Route map E/J: "$T$ not computed", "missing explicit $T$".
- notes/14 §6, "$T$ analytic in $d$ ⇒ bounded below": not a proof, now replaced.
- notes/13: "certified time" corrected to "certified slope".

## Rerun versus stored

| Item | How it was checked this session |
| --- | --- |
| Proposition T | **recomputed** (2 m 35 s): profile boxes and every constant bit-identical to the stored JSON. Precision check: 128 vs 192 bits, coarse cells, $T_{\min}$ agrees to $7\times10^{-10}$ |
| fast test suite at `afe0bed` | **run**: 111 passed, 7 slow deselected, 13 min |
| new and changed tests | **run**: `tests/test_ck_bound.py` (6) and `tests/test_bl_series.py` (4), 2 m 19 s |
| Theorem C box and slopes | **stored evidence** (the slow certificate was not rerun); the profile run starts from its box |
| Theorem F / F′ | **stored evidence**; the fast tests check Theorem F's log coverage and F′'s components; the sweeps and the F′ `verify` replay were not rerun |
| Baouendi–Goulaouic, Frittelli | **read** (primary PDFs) |

## Tests and independent checks

1. The reduced Ricci tensor agrees with the closed form at real and complex points ($10^{-25}$).
2. The Lipschitz bounds hold against an independent NumPy implementation at random complex
   points; a mutation that shrinks them is caught.
3. Refusal cases: $q \ge 1$, and the puncture branch condition.
4. The precision is restored after a run.
5. The certified profile contains the floating-point critical curve.
6. Schwarzschild: $K_0$ brackets the exact $\lvert\mathrm{Ric}\rvert$, and $T_{\min}$ is below the
   exact crunch time.
7. Constraints (Hamiltonian and momentum) vanish order by order, for BL and Schwarzschild.

## Next decisive calculation

1. Validated enclosures of $\lambda_2(L_{S_c})$ for every Fourier mode $m$ (ODE eigenvalue
   problems on the profile).
2. A validated shooting solver for $(\gamma(t), K(t))$ that consumes the Proposition T bounds.

Together they would give an explicit $\delta$, and with it a theorem at a specific numerical
separation $d$.
