# Decision log

Newest first. Each entry: what was decided, why, and where things stand.

## 2026-10-29 (later) — Session 2: attempted computer-assisted proof of the fold

**Decisions**

- **python-flint (Arb) for validated numerics** (`requirements.txt`). Its
  ball arithmetic and power series are much faster than mpmath intervals.
  `ctx.cap` must be set explicitly: the default silently truncates series at
  length 10.
- **A computer-assisted claim needs known-solution tests first.** Any
  validated integrator gets tests against ODEs with exact solutions at several
  orders before it is used for a proof (`tests/test_fold_proof.py`).
- **`research/calc/fold_proof.py` is kept, marked INCOMPLETE.** It does not
  certify (F). Its docstring records the bottleneck and the design for whoever
  continues.

**Progress**

- Validated integrator, parameter jets, axis fixed-point enclosure and
  Krawczyk test implemented. A false certificate was caught (series cap). The
  near-axis overestimation was diagnosed precisely: final width about
  $C s_{\rm end}^2$, independent of $s_1$.

**Next**

- A high-order Taylor-model expansion about the axis (see `research/status.md`).

## 2026-10-29 — Research project: exact and rigorous results on black-hole mergers

**Decisions**

- **Layout.** Research notes continue the numbered series (`notes/07`–`09`).
  Project bookkeeping (status, claim ledger, idea log, literature) lives in
  `research/`, and executable calculations in `research/calc/` (run as
  `python -m research.calc.<name>`). Every symbolic claim has a test in
  `tests/`.
- **Proof standards.** Only exact identities (`cancel`/`simplify` to 0) count
  as symbolic proofs. Exact evaluation at sample points is labelled a
  diagnostic. Computer-assisted sign checks use mpmath interval arithmetic,
  with parameters entered as exact decimals.
- **Dependencies.** NumPy and SciPy added to `requirements.txt` for the
  numerical calculations; mpmath comes with SymPy.
- **Toolkit changes.** `Spacetime(simplify=...)` also accepts a function (for
  example `sp.cancel`, which decides identities between rational functions).
  Kerr uses $\chi = \cos\theta$ so that its components are rational. New:
  `einstein/fields.py` (Maxwell stress tensor), Reissner–Nordström, Kerr, and
  conformastatic/Kastor–Traschen metrics. A `slow` pytest marker covers the
  symbolic-parameter Kerr check (about 3 min); everything else runs in CI.
- **Conventions added** (stated in `research/README.md`): $K_{ij} = \nabla_in_j$
  (Wald's sign), $\theta_\pm = \mathrm{tr}_SK \pm H_S$, Gaussian units for
  Maxwell.
- **Event-horizon merger definition** tracks the components of
  $\mathcal B \cap \Sigma$ that contain each hole, rather than requiring the
  later slice to be connected.

**Progress**

- Problem definition (`notes/07`), Kastor–Traschen benchmark (`notes/08`),
  four hypotheses (`notes/09`). Main results: the reduced KT system and its
  classification (Λ = 0 ⇒ static); a certified KT apparent-horizon merger;
  Lemma A; the H2 obstruction; the near-critical onset law (*) for the
  Brill–Lindquist common horizon, verified two ways (Proposition B,
  conditional on the fold hypotheses (F)).

**Next**

- A computer-assisted proof of (F) (validated ODE enclosures plus an
  interval Newton step). See `research/status.md`.

## 2026-10-29 — Repository made public

**Decision:** switched visibility from private to public so anyone can read
the notes and take part in Discussions. Discussions stay enabled; Issues and
pull requests are open to all GitHub users.

## 2026-10-29 — Repository created

**Decisions**

- **Purpose.** A place to discuss Einstein's field equations: notes to argue
  with, code that checks the math, GitHub Discussions for the conversation.
- **Conventions: signature (−,+,+,+), MTW/Carroll Riemann sign.** The most
  common choice in current GR textbooks (Carroll, MTW, Wald), so readers can
  compare formulas without converting.
- **Code in geometrized units, G = c = 1.** Keeps SymPy expressions short; the
  notes restore G and c where it matters physically.
- **Notes as Markdown with GitHub math, not LaTeX/PDF.** Renders in the
  browser, diffs line by line, and reviewers can comment on a single equation
  in a pull request. Display math in ```` ```math ```` fences because GitHub
  parses them more reliably than `$$`.
- **Plain SymPy, not EinsteinPy/Cadabra/xAct.** One `pip install`, and the
  toolkit is short enough to audit against `notes/00-conventions.md`.
- **Every exact solution in `notes/04` gets a test.** Keeps the notes honest.
- **Kerr is in the notes but not yet in the code.** Simplifying the Kerr Ricci
  tensor in SymPy is slow; left as a first contribution (mark the test slow).
- **GitHub repo starts private, with Discussions enabled.** (Made public the
  same day; see the entry above.)

**Progress**

- Notes 00–06 written; toolkit and tests in place; CI runs pytest on push.

**Next**

- Add Kerr to `einstein/metrics.py` with a slow-marked vacuum test.
- Reissner–Nordström with the electromagnetic stress-energy tensor.
- A note on the Einstein–Hilbert action variation in full, including the
  Gibbons–Hawking–York boundary term.
- Seed Discussions with the prompts from `notes/06-open-questions.md`.

**Environment note**

On Windows, a checkout at a deep path can push pip's vendored files past
MAX_PATH and break `python -m venv`. Creating the venv at a short path
(e.g. `%TEMP%\eeq-venv`) works.
