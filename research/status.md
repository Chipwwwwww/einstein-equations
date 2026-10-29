# Status

*Updated 2026-10-29 (sessions 1 and 2). Start here.*

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

**Session 2.** A computer-assisted proof of (F) was attempted
(`research/calc/fold_proof.py`) and **not** achieved. A first "verified"
Krawczyk run was invalid: python-flint silently truncates power series at
length 10. Once that was fixed, the naive validated integrator overestimates
near the axis (see gap 1). The integrator pieces are tested on ODEs with known
solutions.

## Unresolved gaps

1. **(F) is numerical.** The nondegenerate fold at $d_c$ (simple zero
   eigenvalue, $\phi_0 > 0$, $\alpha, \beta \ne 0$) has no proof yet. The
   session-2 attempt is blocked by a precisely identified bottleneck. The
   zeroth-order axis enclosure has relative width $\propto s_1^2$, and ball
   radii then grow like $(s/s_1)^2$, so the final width is about
   $C s_{\rm end}^2$, independent of $s_1$.
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

**Finish the computer-assisted proof of (F)** in `research/calc/fold_proof.py`.
Already in place: the jets, the Krawczyk test, the segment structure and the
checks (Jacobi-field positivity, embeddedness), all exercised on the
problem. Missing is one ingredient: a **high-order rigorous expansion about
the axis**. Write $A(s) = \sum_{k \le n} a_k s^{2k} + s^{2n+2}\tilde A(s)$ (and
likewise $Q$, $Z$). Compute the $a_k$ by the series fixed point
$A_{\rm new} = \sum_k \Phi_k s^k/(k + 2)$ in `arb_series`. Enclose the
remainders by the same weighted-average argument with weight
$\sigma^{2n+3}$; that needs Taylor-model arithmetic for $\Phi$ along the
solution. This makes the axis width $\propto s_1^{2n+2}$ and breaks the
$C s_{\rm end}^2$ barrier. If the regular part then still wraps too much, add
Lohner (QR) steps. Expected payoff: (F) with rigorous enclosures of $d_c$
(about $10^{-14}$) and of the onset slope.

After that: unequal masses (to test (*) off reflection symmetry), and
boosted Bowen–York data, where the fold happens in time on one slice (AMMS
directly) rather than in a family.

## Reproduce

See [`README.md`](README.md). Expected outputs: `python -m research.calc.kt_horizons`
and `python -m research.calc.bl_fold` print the tables quoted in notes/08 §6
and notes/09 H4. `python -m pytest -m "not slow"` runs everything except the
symbolic-parameter Kerr check.
