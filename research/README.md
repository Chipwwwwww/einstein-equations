# Research project: exact and rigorous results on black-hole mergers

A bounded, reproducible attempt to find new routes to exact binary-black-hole
merger solutions or rigorous existence results, with the honest bookkeeping
needed to tell a result from a hope.

## Scope

**Primary target.** A 4-dimensional, asymptotically flat, vacuum spacetime
($R_{\mu\nu} = 0$, $\Lambda = 0$, $G = c = 1$) describing two initially distinct
black holes that merge. The precise meaning of every word in that sentence is
in [`notes/07-the-merger-problem.md`](../notes/07-the-merger-problem.md).
Late-time approach to Kerr is a separate claim.

**Benchmark.** The Kastor–Traschen spacetimes (Einstein–Maxwell, $\Lambda > 0$,
extremal charges): an exact merger with different physics.
[`notes/08-kastor-traschen.md`](../notes/08-kastor-traschen.md).

**Original work.** Hypotheses, tests and outcomes:
[`notes/09-merger-hypotheses.md`](../notes/09-merger-hypotheses.md). Proofs:
[`notes/10-onset-theorem-proof.md`](../notes/10-onset-theorem-proof.md)
(apparent horizons, Brill–Lindquist) and
[`notes/11-event-horizon-merger.md`](../notes/11-event-horizon-merger.md)
(event horizons, Hintz's extreme-mass-ratio mergers) and
[`notes/12-no-common-horizon.md`](../notes/12-no-common-horizon.md)
(no common apparent horizon at $t = 0$, Brill–Lindquist at $d = 0.78$).

## Files

| File | Contents |
| --- | --- |
| [`status.md`](status.md) | Where things stand: established findings, gaps, next calculation. Start here. |
| [`claim_ledger.md`](claim_ledger.md) | Every claim, labelled: established (with source), reproduced here, candidate new, conjecture/numerical, gap. |
| [`idea_log.md`](idea_log.md) | Hypotheses, the tests run on them, and how they failed or survived. |
| [`literature.md`](literature.md) | Primary sources, what was checked in each and what could not be. |
| [`calc/`](calc/) | Executable calculations (numerical and computer-assisted). |
| [`../tests/`](../tests/) | pytest checks of every symbolic claim. |

## Conventions

Those of [`notes/00-conventions.md`](../notes/00-conventions.md): signature
$(-,+,+,+)$, Riemann tensor of MTW/Carroll, $G = c = 1$. In addition:

- Extrinsic curvature $K_{ij} = \nabla_i n_j$ with $n$ the future unit normal
  ($K > 0$ expanding; Wald's sign, opposite to MTW). With lapse $N$ and zero
  shift, $\partial_t\gamma_{ij} = 2NK_{ij}$. The evolution equation used
  throughout is checked in `tests/test_validation.py::test_adm_evolution_equation_sign`.
- Null expansions of a 2-surface in a slice: $\theta_\pm = \mathrm{tr}_S K \pm H_S$,
  with $H_S$ the mean curvature for the outward normal (positive for round
  spheres in flat space).
- Electromagnetism in Gaussian units: $T_{\mu\nu} = \frac1{4\pi}(F_{\mu\alpha}F_\nu{}^\alpha - \frac14 g_{\mu\nu}F^2)$.

## Dependencies

Python ≥ 3.12 and `requirements.txt`: SymPy ≥ 1.12, pytest ≥ 8, NumPy, SciPy
(numerical calculations in `calc/`), python-flint (Arb ball arithmetic, for the
validated integrators in `calc/fold_proof.py` and `calc/lohner.py`). mpmath, used for the interval
certificates, ships with SymPy. Tested with Python 3.13, SymPy 1.14, NumPy 2.5, SciPy 1.18,
mpmath 1.3, python-flint 0.9.

## Reproduce

```bash
pip install -r requirements.txt
```

```bash
python -m pytest -m "not slow"
```

```bash
python -m pytest -m slow
```

The second command runs the slow checks: the symbolic Kerr check (about three
minutes), Kerr–de Sitter at rational parameters (about two and a half
minutes), the computer-assisted fold certificate (about five minutes) and the
certified onset slopes for two lapses (about seven minutes), and the
thin-neck certificate of notes/12 (15–25 minutes on 4 cores). The z0 sweep
of notes/12 is not a test, because it takes 2.8 hours on 4 cores. Run it as
`python -m research.calc.gap3_sweep 4`; it appends each box to `gap3_sweep.log`
in the current directory. The certified log is
`calc/data/gap3_sweep_d078.log`, and a fast test checks that its boxes cover
the window. The
calculations in `calc/` each run as `python -m research.calc.<name>` from the
repository root; each module's docstring says what it computes and what it
prints. The expected outputs are recorded in [`status.md`](status.md).

## What counts as proof here

- A SymPy identity reduced to `0` by `cancel` on rational expressions, or by
  `simplify` where noted, is treated as proof of that identity.
- Evaluating an expression at sample points, even in exact arithmetic, is a
  diagnostic, never a proof.
- A failure of `simplify` is never treated as evidence that an expression is
  nonzero.
- Interval arithmetic (mpmath `iv`) with outward rounding is used for the
  computer-assisted sign checks in `calc/`; the conclusions are rigorous up to
  the correctness of that library and of the theorems they are combined with,
  which are stated with their hypotheses.
