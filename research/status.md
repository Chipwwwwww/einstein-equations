# Status

*Updated 2026-10-29 (sessions 1–4). Start here.*

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
**Theorem D₀** gives the same merger statement for Λ = 0, conditional on two
properties of the Kerr-stability development. Klainerman–Szeftel prove these
in their own coordinates (§3.8, read). New SymPy checks
(`tests/test_eh_merger.py`):

- Kerr–de Sitter solves $\mathrm{Ric} = \Lambda g$;
- $g^{rr} = \Delta_r/\rho^2$;
- the spheres inside the hole are trapped, for every spin and every
  Λ ≥ 0.

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
4. **Event horizons in vacuum: Λ > 0 done (Theorem D), Λ = 0 conditional
   (Theorem D₀).** The two remaining inputs for Λ = 0 are (Ex), existence
   via Hintz's Remark 6.2 and KS Remark 3.4.3, and anchoring KS's $r$ to
   Kerr's $r$ on the initial layer. Theorem D also rests on our reading of
   the uniformity in Hintz's (1.2) and (6.2). It uses the cosmological region
   in place of $\mathscr I^+$, and locates the horizons only up to margins.
   Proposition B (Brill–Lindquist) still says nothing about event horizons.
5. **KT**: the absence of a common MOTS before $\tau_c$ is numerical; event
   horizons are taken from BHKT.
6. **Literature items marked *not re-checked*** in `literature.md` (Hawking–Ellis
   and Wald propositions, Penrose 1965, Wald–Iyer 1991, Bartnik 1986).
7. *(Closed in session 3.)* The slope for $N = \psi^{-2}$ is certified as well:
   $[10.606, 10.682]$.

## Next concrete calculation

Two candidates, in order of tractability.

**(a) Anchoring for Theorem D₀ (Λ = 0).** Read how KS23 and its cited
references construct the initial data layer from Cauchy data (KS Remark
3.4.3, refs. [36], [37], [8]). Check whether their $r$ on the layer is
$O(\epsilon_0)$-close to Kerr's $r$. If so, (K-esc) and (K-trap) follow from
KS §3.8, and Theorem D₀ becomes conditional only on (Ex).

**(b) Gap 3, below.**

**Gap 3: prove that no common MOTS exists on $\Sigma_0$ for $d > d_c$.**
On the time-symmetric slice a common MOTS is a minimal surface enclosing both
punctures. One route: build a continuous sweep-out of the region outside the
two individual horizons by surfaces with $H > 0$, ending in two disjoint
spheres around the individual horizons. By the comparison (maximum)
principle, an enclosing minimal surface would have to touch some leaf from
inside, which is a contradiction. The leaves can be explicit polar surfaces,
with $H > 0$ certified by interval arithmetic as in notes/08 §6. Near $d_c$
the margins vanish at $S_c$ itself, so this works at a definite $d > d_c$
(e.g. $d = 0.78$), not uniformly down to $d_c$. Combined with Theorem C,
that would make "two separate apparent horizons at $t = 0$, a common one at
$t_*$" a fully rigorous statement about the vacuum development, for those
separations.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
slow checks listed in `README.md`.
