# Status

*Updated 2026-10-29 (session 1). Start here.*

## Established this session

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
apparent-horizon merger statement, conditional on the fold hypotheses (F).

## Unresolved gaps

1. **(F) is numerical.** The nondegenerate fold at $d_c$ (simple zero
   eigenvalue, $\phi_0 > 0$, $\alpha, \beta \ne 0$) has no proof yet.
2. **Proposition B's proof is a sketch.** The Lyapunov–Schmidt argument needs
   function spaces, smooth dependence of the zero-shift slicing on $(d, t)$,
   and remainder estimates written out.
3. **Global absence of a common MOTS on $\Sigma_0$ for $d > d_c$** is
   numerical only.
4. **Nothing about event horizons in vacuum.** H3's inputs (A2), (A3) are
   unchecked. Proposition B says nothing about $\mathscr I^+$, global
   existence or Kerr asymptotics.
5. **KT**: the absence of a common MOTS before $\tau_c$ is numerical; event
   horizons are taken from BHKT.
6. **Literature items marked *not re-checked*** in `literature.md` (Hawking–Ellis
   and Wald propositions, Penrose 1965, Wald–Iyer 1991, Bartnik 1986).

## Next concrete calculation

**A computer-assisted proof of (F).** Use validated interval integration
(an interval Taylor method, or Lohner-type enclosures) of the axisymmetric
minimal-surface ODE and its first and second variational equations in $z_0$
and $d$. The aim is to enclose the shooting map, its $z_0$-derivative and the
mixed derivatives near $(z_0, d) = (0.63266, 0.76620)$. Then prove, by the
interval Newton method on the system {mismatch = 0, ∂mismatch/∂z₀ = 0}, that
a unique degenerate solution exists in an explicit box, with
∂²mismatch/∂z₀² and ∂mismatch/∂d of fixed sign. That gives (F) rigorously and
makes Proposition B conditional only on writing out standard nonlinear
analysis. The existing `bl_fold.py` already has every formula; what is
missing is a validated ODE integrator (or the `vnode`/CAPD-style algorithm
written in mpmath intervals). The singular axis start needs a rigorous series
step.

After that: unequal masses (to test (*) off reflection symmetry), and
boosted Bowen–York data, where the fold happens in time on one slice (AMMS
directly) rather than in a family.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
symbolic-parameter Kerr check.
