# Research project: exact and rigorous results on black-hole mergers

## 2026-11-01 global-exclusion audit checkpoint

Read the [current audit report](report_global_mots_audit_2026-11-01.md) and
[literature matrix](literature_matrix.md) before interpreting historical claims.
The new [note22](../notes/22-global-mots-audit.md) rejects the constant-residual
global shooting shortcut and supplies exact evolved-data replacements.
The critical all-mode gap from the independent branch is integrated as
[note21](../notes/21-horizon-spectral-gap.md) and independently replayed.
Arithmetic handoffs are repaired; corrected full M/M'/N verification remains
pending. The inherited local bracket is neither global nor newly revalidated.
No global equal-mass merger metric or event-horizon theorem is established.
Priorities and precise missing lemmas: [research_gaps.md](research_gaps.md).

The following sections retain the historical programme description, subject
to this audit and the current claim ledger.


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
(no common apparent horizon at $t = 0$, Brill–Lindquist, at $d = 0.78$ and for
$d \in (d_c, 0.767578]$). The head-on example is assembled in
[`notes/13-head-on-collision.md`](../notes/13-head-on-collision.md). The explicit-binary routes
are in [`notes/14-explicit-binaries.md`](../notes/14-explicit-binaries.md), and the certified
convergence of the series representation near the critical surface (route J) is in
[`notes/15-route-j-convergence.md`](../notes/15-route-j-convergence.md), and the evaluable metric built on it in
[`notes/16-bl-metric.md`](../notes/16-bl-metric.md).

## Files

| File | Contents |
| --- | --- |
| [`status.md`](status.md) | Where things stand: established findings, gaps, next calculation. Start here. |
| [`claim_ledger.md`](claim_ledger.md) | Every claim, labelled: established (with source), reproduced here, candidate new, conjecture/numerical, gap. |
| [`idea_log.md`](idea_log.md) | Hypotheses, the tests run on them, and how they failed or survived. |
| [`literature.md`](literature.md) | Primary sources, what was checked in each and what could not be. |
| [`route_map.md`](route_map.md) | Routes toward explicit dynamical binary metrics (A–J): question, decisive test, status (session 10). |
| [`../notes/19-bianchi-coupling.md`](../notes/19-bianchi-coupling.md) | Session 15: Einstein–Bianchi form in Gauss gauge (verified, hyperbolicity corrected), the $t^6$ law proved for all conformally flat vacuum data, the curl-degree hierarchy and lag lemma, the energy identity. |
| [`report_session15.md`](report_session15.md) | Session 15 report. |
| [`../notes/18-structure-discovery.md`](../notes/18-structure-discovery.md) | Session 14: bottleneck map, hypothesis portfolio, silent/nonlocal split (H1), conformally flat foliations for maximal-type lapses (H2). |
| [`report_session14.md`](report_session14.md) | Session 14 report: what was tested, what decided it, what failed and in what scope. |
| [`../notes/17-certified-extension.md`](../notes/17-certified-extension.md) | Session 13: corrected certificate constant, Theorem CK′ and the extension, derivative bounds, $t^6$ and Weyl structure. |
| [`report_session13.md`](report_session13.md) | Session 13 report: corrected and extended certificate, derivative errors before/after, B-to-A results. |
| [`report_bl_metric.md`](report_bl_metric.md) | Session 12 report: the metric evaluator, the sharp certificate, B-to-A experiments, what was rerun. |
| [`report_route_j.md`](report_route_j.md) | Session 11 report: Proposition T (certified convergence time near $S_c$), baseline audit, what was rerun, blocking estimate. |
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
validated integrators in `calc/fold_proof.py` and `calc/lohner.py`, and the convergence certificate
`calc/ck_bound.py`). mpmath, used for the interval
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

### Continuation (sessions 9–10)

Each certificate below was run once, and its record is in `calc/data/`. Replay commands:

```bash
python -m research.calc.gap3_neck interval 4
```

```bash
python -m research.calc.gap3_fold 4
```

```bash
python -m research.calc.gap3_fold fd 4
```

```bash
python -m research.calc.gap3_interval 4
```

```bash
python -m research.calc.gap3_interval verify gap3_interval.log 0.767578125 4
```

1. The neck runs: `interval` is one of nine; the names are in `gap3_neck._INT_RANGES`. Each run
   takes about 10–30 minutes.
2. Fold lemma part 1: about 2.5 CPU-hours.
3. Fold lemma part 2: about 7 CPU-hours.
4. The interval sweep: about 7 hours on 4 cores up to $d = 0.767578$. It resumes from
   `gap3_interval.log` in the current directory, lowest $d$ first, so stop it once
   `python -m research.calc.gap3_interval progress gap3_interval.log` reports the stage you want.
5. Verify a log: cover check, holes, and the band certificate, replayed. The stored log is
   `calc/data/gap3_interval_stage1.log`.

Route E (the BL time series) and the Petrov check:

```bash
python -m research.calc.bl_series
```

```bash
python -m research.calc.petrov_bl
```

The first validates against Schwarzschild free fall and prints the coefficients near $S_c$
(about 6 minutes at order $t^{12}$). The second is exact arithmetic.

### Route J (session 11)

```bash
python -m research.calc.ck_bound
```

This certifies Proposition T of notes/15: $T_{\min} = 2.988\times10^{-4}$ on the tube of radius
$0.02$ about $S_c$, for every $d \in [d_c - 10^{-13}, 0.767578125]$.

- **Runtime.** About 2.5 minutes on one core: 80 s for the Lohner profile run from Theorem C's box,
  60 s for the 2,815 cells.
- **Output.** It writes `calc/data/ck_bound.json` (constants, cells, profile boxes). A run on
  which any hypothesis check fails returns `"ok": false` and writes nothing.
- **Arguments.** Optional: `s0 R_gamma R_g R_K`.
- **Tests.** `python -m pytest tests/test_ck_bound.py` (about 40 s).

The next decisive calculation is in `route_map.md` §J and notes/15 §7.3.

### The metric evaluator (session 12)

```bash
python -m research.calc.ck_bound sharp
```

```bash
python -m research.calc.bl_metric
```

```bash
python -m research.calc.bl_structure
```

1. **The sharp certificate** of notes/16 §2: $T_{\min} = 5.556\times10^{-4}$. It writes
   `calc/data/ck_bound_sharp.json` in about 2.5 minutes. `ck_bound classic` replays note 15's
   certificate with outward rounding.
2. **The evaluator demonstration** ($d = 0.7662$; three points, three times, two patches). It
   writes `calc/data/bl_metric_demo.json` in about 1.5 minutes.
3. **The B-to-A experiments**. They write `calc/data/bl_structure.json` in about 10 minutes; the
   Padé part is most of it.

Use it from Python:

```python
from research.calc.bl_metric import evaluate_metric
r = evaluate_metric("1e-4", "0.17", "0.11", "0.5", "0.7662", tolerance=1e-10)
r["g"], r["dg"], r["ddg"], r["error_budget"], r["status"]
```

Tests: `python -m pytest tests/test_bl_metric.py tests/test_ck_bound.py` (about 1 minute).

### Session 13: corrected and extended certificates

```bash
python -m research.calc.ck_bound classic
```

```bash
python -m research.calc.ck_bound sharp
```

```bash
python -m research.calc.ck_bound ext
```

```bash
python -m research.calc.ck_bound waist
```

```bash
python -m research.calc.bl_metric compare
```

1. **The four certificates.** They write `calc/data/ck_bound{,_sharp,_ext,_waist}.json`, using 4
   processes. Each record carries:
   - the corrected Baouendi–Goulaouic time (`T_bg`);
   - the Theorem CK′ schedule of $(\beta, a, B)$ triples (notes/17 §2);
   - the source commit.

   Wall times and peak memory are in notes/17 §3.
2. **The before/after table** of notes/17 §4. It writes `calc/data/bl_metric_compare.json`.

The evaluator now loads every certificate and uses all that apply to a point.

Tests: `python -m pytest tests/test_bl_metric.py tests/test_ck_bound.py tests/test_bl_structure.py`.

### Session 14: structure experiments

```bash
python -m research.calc.bl_silent
```

```bash
python -m research.calc.bl_foliation
```

1. **H1** (notes/18 §6). It writes `calc/data/bl_silent.json`: the exact low-order checks, the
   $(\mathrm{curl}\,H)_3$ law at held-out points, and the coefficient comparison to $t^{16}$ at six
   generic points. It takes about 20 minutes on one core.

   `python -m research.calc.bl_silent long` adds the two runs to $t^{20}$
   (`calc/data/bl_silent_N10.json`). It takes about 20 minutes on two cores.
2. **H2** (notes/18 §7). It writes `calc/data/bl_foliation.json`: the two-point certificates at three
   separations, plus the gauge and Schwarzschild controls. It takes about one minute.

Tests: `python -m pytest tests/test_bl_silent.py tests/test_bl_foliation.py` (about 35 s).

### Session 15: Einstein–Bianchi form and the $t^6$ theorem

```bash
python -m research.calc.cf_universal
```

```bash
python -m research.calc.bianchi_symbol
```

```bash
python -m research.calc.bl_bianchi sectors
```

1. **`cf_universal`** (about 6 s). It runs the 72 exact evaluations of the proof of Theorem 1 and
   writes `calc/data/cf_universal.json`.
2. **`bianchi_symbol`**. It computes the principal symbols: the frame form is symmetric hyperbolic,
   and the coordinate form has a Jordan block.
3. **`bl_bianchi sectors`** (about 14 min on three cores). It writes the curl-degree sector norms to
   $t^{12}$, `calc/data/bl_bianchi_sectors.json`.

Tests: `python -m pytest tests/test_bl_bianchi.py` (about 75 s).

### Session 16: the explicit merger (notes/20)

```bash
python -m research.calc.bl_barrier
```

```bash
python -m research.calc.bl_foliate
```

1. **`bl_barrier`** (about 40 min on two cores). Theorem M at $d = 0.7662$, $t = 1.2\times10^{-5}$: the two
   barrier surfaces (validated secant refinement, sign change of $F$, $\theta_+$ on every box),
   their Lohner enclosures, nesting and $\Omega \subset U$. It writes `calc/data/bl_barrier.json`.
   `bl_barrier.main(t='6.5e-6', h_in=-4e-8, out_path=...)` gives the near-onset record
   `calc/data/bl_barrier_near.json`.
2. **`bl_foliate`** (about 77 min on four cores; needs `bl_barrier_near.json`). Theorem N: the foliation
   of a shell about $S_c$ by leaves with $\theta_+ > 0$, for every $t \in [0, 5\times10^{-6}]$, and
   $W_t \supset \Omega'$. It uses Lohner sets with the level-$h$ source, $t$ as a set parameter and
   sub-step enclosures; monotone $z_0$-pieces; and a mean-value bound on the fold piece. It writes
   `calc/data/bl_foliate.json`. `bl_barrier` also uses it for its thin enclosures.
3. **`bl_foliate.main_interval()`** (about 50 min on four cores). Theorem M′: a common MOTS at
   $t = 1.2\times10^{-5}$ for every $d \in [0.76619742824011, 0.7662]$. Here $d$ is a Lohner set parameter,
   the barriers come from fixed $z_0$-brackets by the IVT, and nesting follows from the extreme members.
   It writes `calc/data/bl_barrier_interval.json`. Run it with:

   ```bash
   python -c "from research.calc import bl_foliate as bf; bf.main_interval()"
   ```

Tests: `python -m pytest tests/test_bl_barrier.py` (about 10 s; includes the checks of the stored records).
The slow test reproduces `fold_proof.shoot` bit for bit (about 5 min).

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
