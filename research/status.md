# Status

*Updated 2026-10-29 (sessions 1–3). Start here.*

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

## Unresolved gaps

1. **Proposition B's analytic steps are sketches.** These are the smooth
   dependence of the zero-shift slicing (and of the shooting function of the
   true slices) on $(d, t)$, and the equivariance step showing that every MOTS
   near $S_c$ is axisymmetric and reflection symmetric. Both are standard;
   neither is written out with function spaces.
2. **The certificate's trust base.** It rests on Arb and on about 600 lines of
   our own code. Tests and cross-checks are listed in notes/09, but no
   independent person has reviewed the code.
3. **Global absence of a common MOTS on $\Sigma_0$ for $d > d_c$** is
   numerical only. Theorem C is local, near $S_c$.
4. **Nothing about event horizons in vacuum.** H3's inputs (A2), (A3) are
   unchecked. Proposition B says nothing about $\mathscr I^+$, global
   existence or Kerr asymptotics.
5. **KT**: the absence of a common MOTS before $\tau_c$ is numerical; event
   horizons are taken from BHKT.
6. **Literature items marked *not re-checked*** in `literature.md` (Hawking–Ellis
   and Wald propositions, Penrose 1965, Wald–Iyer 1991, Bartnik 1986).
7. **Only geodesic slicing is certified.** The slope for $N = \psi^{-2}$
   (10.644) is numerical; certifying it needs only the corresponding
   $\dot K_{ij}$ in `t_coefficient`.

## Next concrete calculation

**Write out the analytic part of Proposition B** (gap 1) as a short proof
note. It needs three steps:

- (i) A gauge in which the vacuum development depends smoothly on $(d, t)$,
  e.g. harmonic coordinates for the development, then the zero-shift slicing
  obtained by solving $\partial_t x = N n$, which is smooth for small $t$.
- (ii) The expansion map $\Theta(u; d, t)$ on $C^{2,\alpha}$ normal graphs over
  $S_c$, smooth, with $D_u\Theta = L$ at the fold.
- (iii) The Lyapunov–Schmidt reduction, equivariant under rotations and
  $z \to -z$, using that the kernel (spanned by $\phi_0$, Corollary to
  Theorem C) is invariant.

Then add the $N = \psi^{-2}$ certificate (gap 7). After that, the most
interesting open direction is gap 3: a rigorous proof that no common MOTS
exists on $\Sigma_0$ for $d > d_c$. One route is a sweep-out by surfaces with
$H > 0$, verified by interval arithmetic as in notes/08 §6. That would make
"two black holes at $t = 0$, a common horizon at $t_*$" a fully rigorous
statement about the vacuum development.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
symbolic-parameter Kerr check.
