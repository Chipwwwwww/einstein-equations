# Status

*Updated 2026-10-31 (sessions 1–13). Start here.*

## Session 13: corrected and extended certificate, better derivatives, geometric structure

Branch `claude/bl-metric-extension`, from `aabde39` (PR #3's head). Note:
[`notes/17-certified-extension.md`](../notes/17-certified-extension.md).

- **A foundational error, repaired.**
  - The abstract Cauchy–Kovalevskaya condition of Baouendi–Goulaouic is $a < R/(16RC + 8M)$
    (read on the page). It is not the $\min(1/(16C), R/(8M))$ that sessions 11–12 reconstructed
    from garbled text.
  - So the certified times of notes 15–16 were too large. For the same $C$ and $K_0$ they are
    $2.885\times10^{-4}$, not $2.988\times10^{-4}$, and $5.115\times10^{-4}$, not $5.557\times10^{-4}$.
  - The baseline "$T_0 = 5.556\times10^{-4}$" was therefore not justified.
  - A second latent hole, also repaired: NaN-producing ball powers could have been skipped
    silently.
- **Theorem CK′ (proved, notes/17 §2).** BG's proof with a weight exponent $\beta$, and the Lipschitz
  bound split into derivative-loss and output-level zero-order parts.
- **Certified extension (computer-assisted, candidate):**
  - $T = 8.683\times10^{-4}$ on the whole tube $U$ (×1.70 over the corrected baseline);
  - $T = 1.545\times10^{-3}$ on the waist region $U_w$ (×3.02);
  - for every $d \in I$.

  The coverings use the conjugation and reflection symmetry and imaginary sub-boxes.
- **Evaluator.** It combines all certificates and CK′ schedules, has an extra $\partial_t\partial_x$
  bound, and caches jets. The before/after table is in notes/17 §4. Second derivatives at $T_0/2$:
  17.5 (corrected baseline) → 2.6e-03 at the generic point, order 7.
- **B → A.** No category-A expression.
  - For conformally flat time-symmetric vacuum data, $\gamma_2$ is algebraic in Ric. $\gamma_3$
    follows a five-term covariant law (fitted, held-out, exact at 3 points).
  - Pointwise-algebraic-in-Ric evolution is excluded at $t^6$.
  - The tidal tensor obeys $E = E_0 - t^2(\tfrac32E_0^2 - \tfrac12\mathrm{tr}E_0^2)$. Its speciality
    drifts at $t^2$ and its eigenframe rotates at $t^4$, which excludes every closed isotropic law
    $E' = f(E)$.
- **Crunch audit.** For one hole, the series' convergence time at a comoving point is the time to
  reach the $r = 0$ **curvature singularity**. Notes 14 and 16 had called it a gauge breakdown.
- **Not achieved.** An energy-estimate continuation (the feasibility analysis is in notes/17 §3).
  No merger at any explicit $d$.

Session-13 evidence labels:

| Item | Label | How checked |
| --- | --- | --- |
| BG constant | primary source re-read on the page | — |
| certificates classic, sharp, ext, waist | computer-assisted | **recomputed** this session (`python -m research.calc.ck_bound <mode>`) |
| CK′ | proof in notes/17 §2 | lemma integrals sampled in tests |
| evaluator comparison and demo | computer-assisted | **recomputed** |
| $t^4$, $t^6$, $E$ laws and exclusions | exact and ball arithmetic | **recomputed**; tests |
| Theorems C, F, F′; Proposition B; Corollary J | as before | inherited, not rerun |



## Session 12: an evaluable metric (category B), and a bounded search for category A

Branch `claude/bl-metric-evaluator`, from `7c4bff6` (PR #2's head). Note:
[`notes/16-bl-metric.md`](../notes/16-bl-metric.md). Report: [`report_bl_metric.md`](report_bl_metric.md).

- **The metric evaluator** (`research/calc/bl_metric.py`). `evaluate_metric(t, x, y, z, d, ...)`
  returns certified enclosures of the following, all in comoving Cartesian geodesic-slicing
  coordinates ($g_{00} = -1$, $g_{0i} = 0$, $M_{\rm ADM} = 1$):
  - every component of $g_{\mu\nu}$, and $g^{\mu\nu}$;
  - $\partial g$ and $\partial^2g$, including time derivatives.

  It returns an itemised error budget (coefficient rounding, no spatial approximation, time-tail
  bounds, parameters). Points outside $U \times (-T_{\min}, T_{\min})$ or $d \notin I$ are refused
  or marked uncertified.
  - Demonstration at $d = 0.7662$, order 5: metric errors $10^{-14}$, $7\times10^{-10}$ and
    $4\times10^{-6}$ at $t/T_{\min} = 0.1$, $0.25$, $0.5$, at points on the axis, on the plane and
    at a generic point, plus patches. About 5 s per evaluation.
  - Lorentzian on the whole certified domain.
- **Audit and repair of the session-11 certificate.** Its scalar bounds had used round-to-nearest
  floats: cells, the $d$-hull, $\sup\lvert\psi\rvert^{-4}$, the $G$-box radius, $C$ and $T_{\min}$.
  Everything is now outward rounded. The replay moves $T_{\min}$ by $-5.6\times10^{-12}$
  (relative), so $2.988\times10^{-4}$ stands.
- **Proposition T′ (computer-assisted, candidate).** $T_{\min} = 5.556\times10^{-4}\,M$, 86% more.
  It comes from exact norms of the principal part of the Ricci operator: $\lVert A(I, \cdot)\rVert = 2$,
  and $\tfrac12\lVert T\rVert = \tfrac12\sqrt{14 + 2\sqrt{33}}$ instead of 3.
- **Limits of this method.**
  - Localising to the waist of $S_c$ gives no gain.
  - $s_0$ is capped by the puncture terms' complex singular cones $\sum_i(w_i - p_{\pm,i})^2 = 0$.
  - For any weights, $C \ge 2\sqrt3\sup\lvert\psi\rvert^{-2}$, so $T \lesssim 10^{-3}$.
  - Recentering cannot help: the analytic width is used up at rate $16C$.
- **B → A experiments.** A category-A metric was **not found**.
  - *A:* the $t^4$ coefficient is the universal $\tfrac56\mathrm{Ric}\cdot\mathrm{Ric} - \tfrac16\lvert\mathrm{Ric}\rvert^2\gamma_0 - \tfrac1{12}\Delta\mathrm{Ric}$
    (proved, and checked exactly).
  - *C1:* the geodesic slices are not conformally flat; the Cotton tensor is $\ne 0$ at $t^2$.
  - *C2:* the flow is not silent; $B_{ij} \ne 0$ at $t^3$. Both are certified by enclosures that
    exclude 0, with Schwarzschild as the control.
  - *D:* Padé, exploratory.
- **Merger.** Not certified at any explicit $d$. At $d = 0.7662$ the linear law places the common
  MOTS at $t_* \approx 0.011\,T_{\min}$, inside the slab, but $\delta_1$ is still not explicit.

Session-12 evidence labels:

| Item | Label | How it was checked |
| --- | --- | --- |
| classic and sharp certificates | computer-assisted | **recomputed** this session (`ck_bound classic`, `ck_bound sharp`, about 75 s of cells each) |
| evaluator and its demo | computer-assisted | **recomputed** (`python -m research.calc.bl_metric`) |
| structure experiments | ball / exact arithmetic | **recomputed** (`python -m research.calc.bl_structure`) |
| Theorems C, F, F′; Proposition B; Corollary J | as before | **inherited**, not rerun |
| embedding of the CK region in the MGHD | standard | cited, not re-checked |


## Session 11: route J — a certified convergence time near the critical surface

Branch `claude/route-j-ck-certificate`, from the research head `afe0bed` of
`claude/black-hole-merger-solutions-rpa9hj`. Note: [`notes/15-route-j-convergence.md`](../notes/15-route-j-convergence.md).

- **Proposition T (computer-assisted; candidate new).**
  - The objects: $U$, the closed tube of radius $0.02$ about the critical surface $S_c$, and
    $I = [d_c - 10^{-13}, 0.767578125]$.
  - The claim: for every $d \in I$, the geodesic-slicing vacuum development of the BL data exists
    on $U \times (-T_{\min}, T_{\min})$, with $T_{\min} = 2.988\times10^{-4}\,M$. It is real
    analytic there and equal to its convergent Taylor series in $t$, with an explicit truncation
    bound.
  - Method: the abstract Cauchy–Kovalevskaya theorem (Nishida; Baouendi–Goulaouic), with
    re-derived constants. The Lipschitz and source bounds are certified by Arb balls on 2,815
    complex cells.
  - Record: `research/calc/data/ck_bound.json`. Replay: `python -m research.calc.ck_bound`
    (about 2.5 min). It was replayed in session 11 and is bit-identical.
- **Corollary J (qualitative).**
  - For $d \in (d_c, d_c + \delta_1)$, the local common-MOTS formation of Proposition B happens
    inside the certified region.
  - $\delta_1$ is **not explicit**. The heuristic from the linear law is
    $\delta_1 \approx T_{\min}/2.39 = 1.25\times10^{-4}$.
- **Not achieved:** an explicit $\delta$ with $\lvert t_* - c(d - d_c)\rvert \le C_2(d - d_c)^2$.
  The blocking estimates are in notes/15 §7.3:
  - a certified spectral gap $\lambda_2(L_{S_c})$, over all modes;
  - certified second derivatives of the true shooting function on $\Sigma_t$. These need a
    validated MOTS shooting solver for non-conformally-flat data with $K \ne 0$.
- **Scope.** Local branch only. It is not shown to be the first common apparent horizon of the
  slicing, and it says nothing about event horizons.
- **Audit of route E** (notes/15 §8):
  - The note-14 claim "T is analytic in d, hence bounded below" was not a proof; it is replaced
    by the uniform certificate.
  - $T \approx 2.7$ remains numerical.
  - Constraint propagation is now cited to Frittelli 1997 (read). Both constraints are also
    checked order by order.

Session-11 evidence labels:

| Item | Label | How it was checked in session 11 |
| --- | --- | --- |
| Proposition T | computer-assisted | **recomputed** twice (bit-identical); precision 128 vs 192 bits agrees to $7\times10^{-10}$ |
| tests `test_ck_bound.py`, momentum constraint | tests | **run** |
| Theorem C box, $c \in [2.3768, 2.3936]$ | computer-assisted | **read from stored evidence** (the slow test was not rerun); the profile run starts from its box |
| Theorem F′, Theorem F | computer-assisted | **read from stored evidence**. The 111 fast tests run at the session-11 baseline (13 min) exercise their components and the Theorem F log coverage; the F′ `verify` replay and both sweeps were not rerun |
| Proposition B, Corollary J | analytic | written out; $\delta_0$, $\delta_1$ not explicit |


## Session 10: baseline audit and the explicit-binary programme

**Where the research lives.** Branch `claude/black-hole-merger-solutions-rpa9hj`
(`main` holds only the original discussion notes 00–06). No pull request existed before
session 10.

**Audit of the main claims** (full list: [`claim_ledger.md`](claim_ledger.md)):

| Claim | Category | Replayable? |
| --- | --- | --- |
| Theorem C (fold, $d_c$, signs) | computer-assisted (Arb) | yes: `python -m research.calc.fold_proof` (~10 min), and a slow test |
| Theorem F ($d = 0.78$) | computer-assisted | partly: the slow test checks the **coverage** of the stored log `research/calc/data/gap3_sweep_d078.log`; replaying the computation is `python -m research.calc.gap3_sweep` (11.4 CPU-h, run once in session 7); the neck lemma is replayed by a slow test |
| Theorem F′ ($d \in (d_c, 0.767578]$) | computer-assisted | the sweep was run once (log `research/calc/data/gap3_interval_stage1.log`); `verify` replays the band certificate for every waist and checks coverage, but does **not** re-run the sweep's boxes (that is a resume of `python -m research.calc.gap3_interval 4` on an empty log, ~7 h); fold lemma and neck runs stored with replay commands |
| Proposition B / Theorem H (apparent-horizon merger, near-critical BL) | analytic + computer-assisted; $\delta_0$ not explicit | via Theorem C and F′ |
| Theorems D/E, Lemmas 6–7 (event horizons, extreme mass ratio) | **conditional**: our reading of Hintz 2024/2026, self-reviewed in session 9, not externally reviewed | analytic |
| Kastor–Traschen horizons | reproduced; the common-MOTS absence before $\tau_c$ is numerical | yes |

The historically flagged issues stand as follows.

- *Local vs global MOTS:* settled at $t = 0$ by Theorem F′. Open for $0 < t < t_*$ outside a
  neighbourhood of $S_c$.
- *$d = 0.78$ vs the neighbourhood of $d_c$:* settled by Theorem F′, which starts at $d_c$.
- *Apparent vs event horizons:* the BL results are apparent-horizon results only. The event-horizon
  results are EMR and conditional.
- *Log coverage vs replay:* see the table above.

**New in session 10** ([`notes/14-explicit-binaries.md`](../notes/14-explicit-binaries.md),
[`route_map.md`](route_map.md)):

- **No explicit dynamical binary metric was found.** None is known in the literature checked.
- **Proposition P (proved, exact).** The BL head-on spacetime is Petrov type I on an open dense
  set of the initial slice, for every $d$. So no algebraically special ansatz describes it:
  Kerr–Schild, Kerr–Schild–Kundt, Robinson–Trautman, Kundt.
- **Literature obstruction.** Integrable methods need two commuting Killing fields. AF radiative
  spacetimes with two Killing fields are boost-rotation symmetric, with struts, strings or negative
  masses.
- **Proposition S (route E).** The geodesic-slicing development of BL data is an explicit
  convergent series in $t$ (Cauchy–Kovalevskaya): category B, local. It is validated against exact
  Schwarzschild free fall. Near $S_c$ the radius is about $2.7M$ (numerical), against
  $t_* \le 0.0033$ on the interval of Theorem F′.
- **Most promising route (J):** a certified lower bound on that radius. *Done in session 11:*
  Proposition T, $T_{\min} = 2.988\times10^{-4}$ on a tube about $S_c$. Making $\delta_0$ explicit
  is still open (see above).


## Established so far

**Tools.** The SymPy curvature code passes independent checks: an
all-lower-index Riemann formula written from scratch, the algebraic
symmetries, the conformally flat Ricci formula, the 3+1 evolution identity
(with a mutation test), Reissner–Nordström with its Kretschmann scalar, and
Kerr Ricci-flatness as an exact rational identity.

**Benchmark (Kastor–Traschen; Einstein–Maxwell, Λ > 0).** For an arbitrary
conformal factor $U(\tau, \vec x)$ the field equations reduce to three scalar
expressions. The KT family is everything in the ansatz, and **Λ = 0 forces
staticity**. One centre is RNdS with $Q = m$; $M_{\rm ext} = 1/(4|H|)$; only
$Hm_i$ and $H|x_{ij}|$ matter. MOTSs on KT slices are constant-mean-curvature
surfaces. For $m = 0.05$, $d = 0.5$, $H = -1$, the common apparent horizon is
born by a fold at $\tau_c = -0.012884047$ (numerical). Two disjoint MOTSs at
$\tau = -0.05$ and a common one at $\tau = -0.009$ are **certified** by
interval arithmetic plus Andersson–Metzger. Event-horizon coalescence rests
on BHKT 1994 and was not re-derived.

**Hypotheses.** H1 (umbilic zero-shift slicing at Λ = 0) fails by Lemma A.
H2 (hyperboloidal KT) fails: the slice curvature forbids centres. H3 (event
horizons of Hintz's EMR mergers) is designed, with two unverified inputs.
**H4 survives**: for near-critical Brill–Lindquist data (vacuum, Λ = 0,
asymptotically flat) the common apparent horizon appears at

```math
t_* = (d - d_c)\,\frac{\langle\phi_0, \partial_dH\rangle}{\int N\phi_0|A|^2\,dA}
\quad\Rightarrow\quad t_* = 2.3852\,(d - d_c)\ (N = 1), \qquad 10.644\,(d - d_c)\ (N = \psi^{-2}),
```

with $d_c = 0.7661974282\,M$. The formula is verified against direct
solutions to 6 significant figures, and $d_c$ and the critical area match the
literature to all quoted digits. Proposition B in notes/09 turns this into an
apparent-horizon merger statement.

**Session 3: the fold hypothesis (F) is proved by computer (Theorem C,
notes/09).** A Krawczyk test in Arb ball arithmetic
(`research/calc/fold_proof.py`) gives three things:

- a unique fold point of the shooting function, with
  $d_c = 0.766197428240 \pm 1.1\times10^{-13}\,M$;
- the signs $\partial^2_{z_0}F \in [240.87, 242.70]$, $\partial_dF \in [1.2315, 1.2364]$
  and $\partial_tF \in [-0.51816, -0.51655]$;
- a positive Jacobi field and an embedded critical surface.

Consequently the onset slope for geodesic slicing is **rigorously** in
$[2.3768, 2.3936]$. The implementation is checked against exact ODE solutions,
SymPy jets, the independent floating-point code, and a second step schedule.

Two corrections to session 2. Its first "verified" run was invalid (the
python-flint series cap). Its diagnosis of an axis "barrier" was wrong: the
real cause was low-order steps at large step sizes, fixed by a piecewise
schedule.

**Session 4: event horizons of extreme-mass-ratio mergers (notes/11).**
Hypothesis H3 is executed. **Theorem D** treats Hintz's Λ > 0 merger
spacetimes, where a small Kerr hole of mass $\epsilon m$ plunges into a slowly
rotating Kerr–de Sitter hole. Take $\mathcal B$ to be the set of events that
cannot signal the region beyond the cosmological horizon. Then, for all small
$\epsilon$:

- both interiors lie in $\mathcal B$, with closed trapped surfaces;
- the event horizons sit within any prescribed margins of $r_H$ and of the
  small hole's own $\hat r_+$;
- the two holes lie in different components of $\mathcal B \cap \Sigma_t$
  before the plunge, and in the same one after it;
- $\mathcal B$ is connected.

Hintz listed the last point as a conjecture. The proof uses only his stated
estimates and elementary causal theory: escape routes are timelike curves
with margin, which survive $C^0$ perturbations, and traps are time functions.

**Session 5: Λ = 0.** Hintz's 2026 stability theorem for subextremal Kerr
(arXiv:2606.28253) works in fixed coordinates, and states the Λ = 0 merger
construction as its consequence (3). So Theorem D now holds for Λ = 0 too,
with $\mathcal B$ = events that cannot signal arbitrarily large distances.
**Theorem E** uses the textbook $\mathcal B = M \setminus J^-(\mathscr I^+)$.
It shows two black-hole components before the plunge and one after, which is
the event-horizon merger of notes/07. Its extra input, that far points reach
$\mathscr I^+$, is derived (Lemma 6) from the 2026 paper's gauge at null
infinity. There the outgoing null component of the metric decays
integrably. The derivation is our reading of that paper.

New SymPy checks (`tests/test_eh_merger.py`, sessions 4–5):

- Kerr–de Sitter solves $\mathrm{Ric} = \Lambda g$;
- $g^{rr} = \Delta_r/\rho^2$;
- the spheres inside the hole are trapped, for every spin and every
  Λ ≥ 0;
- outside Kerr, $\tilde u = t - r_* + a\arctan(r/a)$ is a time function;
- $\partial_t + \frac12\partial_r$ is uniformly timelike for $r \ge 10M$.

**Session 7: no common apparent horizon at $t = 0$, $d = 0.78$ (notes/12).**
**Theorem F** (computer-assisted): Brill–Lindquist data with equal masses at
$d = 0.78$ ($d_c = 0.7662$) contain no closed minimal surface enclosing both
punctures. The proof has three parts.

- **Reduction** (session 6). The outermost common surface would be a stable,
  symmetric sphere, generated by a shooting curve from $(0, z_0)$ with
  $z_0 \in [0.545, 1.05]$.
- **Thin-neck lemma.** A waist of radius $\lambda \le 0.00488$ makes the
  surface unstable. Near it the surface is a catenoid up to $O(\lambda)$, and
  the $\lambda$-variation Jacobi field changes sign on both sides of the
  waist. This is certified for all small $\lambda$ at once. An analytic
  **approach lemma** hands any curve heading into the axis to it.
- **Validated sweep.** All 2,452 $z_0$-intervals end in a certified
  exclusion event, using Lohner's method (`research/calc/lohner.py`):
  - plane crossing at a non-right angle;
  - leaving the ball;
  - entering a throat ball;
  - thin neck.

  It took 2.8 hours on 4 cores; the log is stored.

Falsification checks pass. At $d = 0.76 < d_c$, the intervals containing the
two real common surfaces are not excluded. A second integrator agrees on the
neck Jacobi field.

Prior work: non-existence at *large* separation is a theorem
(Dain–Weinstein–Yamada 2011, non-quantitative). The explicit near-critical
case appears to be **potentially new**.

## Unresolved gaps

1. **Proposition B: now written out in notes/10**, except for the standard
   ingredients, which are cited but not re-checked. These are smooth
   dependence of the development on parameters (harmonic gauge), and
   smoothness of the expansion map on $C^{2,\alpha}$ graphs. The equivariant
   Lyapunov–Schmidt step and the fold argument are written out.
2. **The certificate's trust base.** It rests on Arb and on about 600 lines of
   our own code. Tests and cross-checks are listed in notes/09, but no
   independent person has reviewed the code.
3. **Global absence of a common MOTS on $\Sigma_0$ for $d > d_c$.**
   - Proved at $d = 0.78$ (Theorem F, session 7) and for every $d \in (d_c, 0.767578]$
     (Theorem F′, sessions 9–10), both computer-assisted.
   - Open for $d \in (0.767578, 0.78)$ and $d > 0.78$, except at large $d$, where Dain–Weinstein–Yamada
     apply (non-quantitative).
   - Open for $0 < t < t_*$ outside a neighbourhood of $S_c$. So the common MOTS of Theorem H is not
     shown to be the first common apparent horizon.
   - The code's trust base: about 1,100 lines (F) plus about 1,500 (F′) plus Arb, tested and
     falsified, not independently reviewed.
4. **Event horizons in vacuum.** Status:
   - Theorem D (Λ ≥ 0, escape-defined black-hole region) is done.
   - Theorem E (Λ = 0, with $\mathscr I^+$) is proved given Lemma 6. Lemma 6
     is our reading of Hintz 2026's structure at $\mathscr I^+$: Lemma 3.19,
     Definition 5.3, and the support properties in the proof of Theorem
     13.1. An expert check of that reading is the most valuable review.
   - *Session 7:* localisation and connectedness hold for the textbook
     region too. Corollary E1 shows $\mathcal B_{\mathscr I} = \mathcal B$,
     via Lemma 7, which extends Lemma 6 to all late times. Lemma 7 reads Hintz
     2026's (13.2)–(13.3), (3.44) and the proof of Lemma 11.5 for
     $t_* \ge 1$, and Lemma 6 (hence the support of the Nash–Moser
     correction) for $t_* < 1$.
   - *Session 9:* Lemmas 6 and 7 were reviewed step by step against the
     paper's text (notes/11, review table). The frame, components, Kerr and
     boost terms, cut-offs, function-space normalisations and the metric
     for $t_* < 1$ are all confirmed. One error was corrected: Lemma 7's
     case split, and its claim not to need the support property. This was
     a self-review; an outside expert check is still advised.
   - Theorem D rests on our reading of the uniformity in Hintz III (1.2),
     (6.2), and of the smallness of the decay constants in Hintz 2026.
   - The horizons are located only up to margins.
   - Proposition B (Brill–Lindquist) still says nothing about event horizons.
5. **KT**: the absence of a common MOTS before $\tau_c$ is numerical; event
   horizons are taken from BHKT.
6. **Literature items marked *not re-checked*** in `literature.md` (Hawking–Ellis
   and Wald propositions, Penrose 1965, Wald–Iyer 1991, Bartnik 1986).
7. *(Closed in session 3.)* The slope for $N = \psi^{-2}$ is certified as well:
   $[10.606, 10.682]$.

## Next concrete calculation

*(Session 11.)* The decisive calculation toward an explicit $\delta$ in Theorem H (notes/15 §7.3):

1. **The spectral gap of the Jacobi operator of $S_c$, over every Fourier mode.** Certify
   $\lambda_2(L_{S_c}) \ge \underline\lambda > 0$. Each mode $e^{im\varphi}$ gives an ODE eigenvalue
   problem on the profile. The $m = 0$ ground state $\phi_0 > 0$ with eigenvalue 0 is Theorem C.
   This quantifies the Lyapunov–Schmidt step of Proposition B.
2. **A validated axisymmetric MOTS shooting solver for $(\gamma(t), K(t))$.** Its inputs are the
   Proposition T enclosures and Cauchy estimates (notes/15 §6). Its outputs are certified bounds on
   $\partial_t^2\Phi$, $\partial_t\partial_d\Phi$, $\partial_t\partial_{z_0}\Phi$, $\partial_d^2\Phi$,
   $\partial_d\partial_{z_0}\Phi$ and $\inf\partial^2_{z_0}\Phi$ on
   $[z_{0,c} \pm r] \times [d_c, d_c + \delta] \times [0, T_{\min}]$.
3. Optional: a sharper convergence time (better scale or constants). The full interval of
   Theorem F′ needs $T \gtrsim 0.0033$, eleven times $T_{\min}$.

*Earlier items.* (a) Locate the merger event in Hintz's EMR spacetimes as ε → 0: the first point
where the two horizon components touch. This needs the horizon at the scale where the margins
$c$, $\eta$ of Theorem D are removed, i.e. null geodesics rather than timelike curves. The
heuristic answer is Emparan–Martínez's exact extreme-mass-ratio horizon. (b) Gap 3 for a whole
interval of separations: *done* as Theorem F′ (sessions 9–10).

**Gap 3 at $d = 0.78$: done (session 7, notes/12).** The obstruction
recorded in session 6 was the necks of every size near $z_{\rm top}$. It is
handled by the thin-neck lemma and the approach lemma (notes/12 §§2–3),
not by isolating blocks.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
slow checks listed in `README.md`.
