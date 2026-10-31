# Status

*Updated 2026-10-31 (sessions 1–10). Start here.*

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
- **Most promising route (J):** a certified lower bound on that radius would give an exact
  convergent representation of a region containing a proved apparent-horizon merger, and could
  make $\delta_0$ explicit.


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
   - Proved at $d = 0.78$ (Theorem F, computer-assisted, session 7).
   - Every other $d$ is open. That includes the $d$ near $d_c$ where the onset
     theorem applies: its $\delta_0$ is unspecified.
   - The new code's trust base: about 1,100 lines plus Arb, tested and
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

*(Session 10: the next decisive calculation is route J. Run the route-E recursion with interval
(Arb) jets, and add a Cauchy–Kovalevskaya majorant tail. That gives a certified lower bound on
the convergence radius near $S_c$, and certified $\partial_t^2F$ toward an explicit $\delta_0$.
See `route_map.md` §J. The items below are from session 8.)*

Two candidates, in order of tractability.

**(a) Locate the merger event as ε → 0.** This is the first point where the
two horizon components touch. It would need the horizon at the scale where
the margins $c$, $\eta$ of Theorem D are removed: null geodesics, not
timelike curves. The heuristic answer is Emparan–Martínez's exact
extreme-mass-ratio horizon, so the natural first step is a rigorous version
of their rescaling limit.

**(b) Gap 3 for a whole interval of separations.** The next target is
"no common MOTS at $t = 0$ for every $d \in (d_c, 0.78]$", which would join
Theorem F to the onset theorem. Route:

1. Make $d$ a ball parameter throughout. The code is ready for it in
   `f_bl`, `start_box`, the barriers and the neck lemma.
2. Sweep $(z_0, d)$-boxes as in Theorem F, leaving out a small neighbourhood
   of the fold point $(z_{0,c}, d_c)$.
3. Cover that neighbourhood with the local uniqueness of Theorem C. Its
   Krawczyk box is tiny, so a larger certified neighbourhood is needed
   first: a Krawczyk test on a bigger box, or Taylor models.

Cost driver: the medium necks (radius $0.005$–$0.02$), where the sweep needs
$z_0$-intervals of width about $10^{-3}(z_0 - z_{\rm top})$. A Lohner
variant that keeps parameters out of the moving frame
($x = m + Ar + P(p - p_0)$) would let the thin-neck lemma reach radius
about $0.03$, where it holds in floating point. That would remove most of
this cost (`research/calc/gap3_neck_lohner.py`, idea log session 7).

**Gap 3 at $d = 0.78$: done (session 7, notes/12).** The obstruction
recorded in session 6 was the necks of every size near $z_{\rm top}$. It is
handled by the thin-neck lemma and the approach lemma (notes/12 §§2–3),
not by isolating blocks.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
slow checks listed in `README.md`.
