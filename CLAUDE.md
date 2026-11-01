# CLAUDE.md — einstein-equations

A discussion repo about Einstein's field equations: Markdown notes with math,
a small SymPy toolkit (`einstein/`) and tests that check the notes' claims.

## Layout

- `notes/NN-topic.md` — one topic per file, numbered in reading order.
- `einstein/tensors.py` — `Spacetime` class: Christoffel, Riemann, Ricci, Einstein, Kretschmann, divergence.
- `einstein/metrics.py` — the metrics from the notes, plus a `REGISTRY` for the CLI.
- `einstein/__main__.py` — `python -m einstein <metric>` prints the curvature.
- `tests/` — pytest; every exact solution in `notes/04` has a test. `-m slow` marks long checks.
- `research/` — black-hole merger project: `status.md` (start here), claim ledger, idea log, literature; `research/calc/` holds calculations run as `python -m research.calc.<name>`.
- `docs/decision_log.md` — decisions and current progress (keep it updated so another machine can pick up).

## Commands

```bash
pip install -r requirements.txt
pytest
python -m einstein schwarzschild
```

On Windows, if `python -m venv` or `pip install` fails with path errors, the
checkout is too deep for MAX_PATH: create the venv at a short path
(e.g. `%TEMP%\eeq-venv`) instead.

## Rules

- Conventions are fixed: signature (−,+,+,+), MTW/Carroll Riemann sign, code in
  G = c = 1. They are stated in `notes/00-conventions.md` and the
  `einstein/tensors.py` docstring; the two must stay identical.
- When quoting a source with different conventions, convert to ours and say so.
- A formula in the notes that SymPy can check gets a test in `tests/`.
- Display math goes in ```` ```math ```` fences (GitHub renders them reliably);
  inline math uses `$...$`. Avoid `$$`.
- Each note ends with a "Discussion prompts" section.
- Record decisions and progress in `docs/decision_log.md`.
