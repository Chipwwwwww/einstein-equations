# Status

*Updated 2026-10-29 (sessions 1–5). Start here.*

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

## Unresolved gaps

1. **Proposition B: now written out in notes/10**, except for the standard
   ingredients, which are cited but not re-checked. These are smooth
   dependence of the development on parameters (harmonic gauge), and
   smoothness of the expansion map on $C^{2,\alpha}$ graphs. The equivariant
   Lyapunov–Schmidt step and the fold argument are written out.
2. **The certificate's trust base.** It rests on Arb and on about 600 lines of
   our own code. Tests and cross-checks are listed in notes/09, but no
   independent person has reviewed the code.
3. **Global absence of a common MOTS on $\Sigma_0$ for $d > d_c$** is
   numerical only. Theorem C is local, near $S_c$.
4. **Event horizons in vacuum.** Status:
   - Theorem D (Λ ≥ 0, escape-defined black-hole region) is done.
   - Theorem E (Λ = 0, with $\mathscr I^+$) is proved given Lemma 6. Lemma 6
     is our reading of Hintz 2026's structure at $\mathscr I^+$: Lemma 3.19,
     Definition 5.3, and the support properties in the proof of Theorem
     13.1. An expert check of that reading is the most valuable review.
   - Localisation and connectedness of the textbook region remain open.
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

Two candidates, in order of tractability.

**(a) Locate the merger event as ε → 0.** This is the first point where the
two horizon components touch. It would need the horizon at the scale where
the margins $c$, $\eta$ of Theorem D are removed: null geodesics, not
timelike curves. The heuristic answer is Emparan–Martínez's exact
extreme-mass-ratio horizon, so the natural first step is a rigorous version
of their rescaling limit.

**(b) Gap 3, below.**

**Gap 3: prove that no common MOTS exists on $\Sigma_0$ for $d > d_c$
(partial, session 6).** Established at $d = 0.78$:

- **Barriers** (`research/calc/bl_global.py`, interval arithmetic):
  - spheres of radius $\le 0.155$ about each puncture have $H < 0$;
  - spheres $|x| \ge 1.05$ have $H > 0$.
- **Reduction.** A common minimal surface would have a symmetric outermost
  representative, generated by a shooting curve from $(0, z_0)$ with
  $z_0 \in [0.545, 1.05]$.
- **Landscape (floating point).**
  - For $z_0$ below the top of the individual horizon ($z_{\rm top} = 0.5896816$)
    the curves dive into a throat.
  - For $z_0 \in [0.605, 0.87]$ they meet the plane with $F \ge 0.017 > 0$.
  - Beyond that they leave the ball.

What blocks a proof is the family just above $z_{\rm top}$. Those curves hug
the individual horizon, pass a neck of radius $\approx 4(z_0 - z_{\rm top})$
near its bottom ($z = 0.2126$), and flare out. Validated integration cannot
follow necks that pinch to zero.

Proposed route:

1. Reparametrize with $ds = \rho\,d\sigma$ and
   $(u, v) = \rho(\cos\alpha, \sin\alpha)$. The system is then regular, and
   the axis becomes a line of hyperbolic fixed points (eigenvalues $\pm1$).
2. Handle the passages near $(0, z_{\rm bot})$ and $(0, z_{\rm top})$ with
   isolating blocks and cone conditions (Shilnikov / covering-relation
   style). The individual horizon's profile, read in both directions, is a
   heteroclinic cycle between these two saddles. Curves near it can wind
   several times, giving multi-layer necked spheres.
3. Do a validated sweep on the rest of the window.

This is multi-session work.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
slow checks listed in `README.md`.
