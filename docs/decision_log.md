# Decision log

Newest first. Each entry: what was decided, why, and where things stand.

## 2026-10-29 (session 5) — Λ = 0 event-horizon merger via Hintz 2026

**Decisions**

- **Use Hintz's 2026 subextremal Kerr stability theorem (arXiv:2606.28253)
  instead of Klainerman–Szeftel for Λ = 0.** It is posed in fixed
  coordinates, with $g = g_b + h$ and $h$ decaying in every asymptotic
  region, so the causal arguments of Theorem D transfer directly. With KS we
  would have had to anchor their coordinates to Kerr's. The 2026 paper also
  states the Λ = 0 merger construction itself (its consequence (3)).
- **Two definitions of the black-hole region, kept separate.**
  - Theorem D uses "cannot signal arbitrarily large distances". It is robust
    under $C^0$ perturbations, and agrees with the textbook one in exact
    Kerr.
  - Theorem E uses $M \setminus J^-(\mathscr I^+)$ and is stated conditional
    on (N). Reaching $\mathscr I^+$ needs retarded-time control. We traced
    where it would come from (Hintz 2026, Definition 5.3), but did not check
    the two remaining steps.
- **Session 4's Theorem D₀ is marked superseded**, not deleted: the ledger
  keeps claim 44 with a pointer, and note 11 keeps a short account of the
  KS route.

**Progress**

- notes/11: Theorem D now covers Λ ≥ 0. Theorem E, Lemma 5 (an escape
  time function outside Kerr) and the far-field steering field were added.
- Two new exact SymPy checks in `tests/test_eh_merger.py`.
- Claim ledger 45–48; status, idea log, literature (Hintz 2026, Chen–Klainerman,
  Hintz Duke 2026, Shen) updated.

**Next**

- Prove (N) from Hintz 2026's structure at $\mathscr I^+$ (status.md, (a)).
- Gap 3 (Brill–Lindquist, no common MOTS at $t = 0$).

## 2026-10-29 (session 4) — Event horizons of extreme-mass-ratio mergers

**Decisions**

- **Pursue H3 rather than gap 3.** Gap 3 (no common MOTS at $t = 0$) needs a
  certified sweep-out and is a long computation. H3 could be settled by
  reading what Hintz and Klainerman–Szeftel actually state.
- **Measure escape against the cosmological region for Λ > 0.** Future null
  infinity of Kerr–de Sitter lies outside Hintz's domain, but the region
  beyond the cosmological horizon lies inside it. In exact Kerr–de Sitter the
  two definitions of the black-hole region agree (notes/11, Lemma 1).
- **Use timelike curves with margin instead of null geodesics.** This needs
  only $C^0$ control of the metric, which Hintz provides on both scales. It
  removes the session-1 inputs (A2) and (A3).
- **For Λ = 0, state hypotheses rather than claim a theorem.** Hintz's
  Remark 6.2 is a remark. KS prove the needed escape and trap properties
  only in their own coordinates. So Theorem D₀ is labelled conditional, with
  (Ex), (K-esc), (K-trap) spelled out.
- **Kerr–de Sitter added to `einstein/metrics.py`**, in Carter's form with
  $\chi = \cos\theta$, so identities reduce by `cancel`. The full
  $\mathrm{Ric} = \Lambda g$ check takes about 2.5 minutes at rational
  parameters and is marked `slow`. The Schwarzschild–de Sitter case runs in
  CI.

**Progress**

- notes/11: Theorem D (Λ > 0) proved from Hintz Theorem 1.1, Definition 4.1
  and eqs. (1.2), (2.19), (6.2). Conclusions: both interiors in
  $\mathcal B$, horizons localised, two components before the plunge and one
  after, $\mathcal B$ connected. Theorem D₀ (Λ = 0) conditional.
- `tests/test_eh_merger.py` covers Kerr–de Sitter vs Kerr and SdS, the
  formula for $g^{rr}$, the area density, and the positivity of
  $\partial_r E$ that makes the interior spheres trapped.
- Claim ledger 28 updated; claims 40–44 added. The idea log records H3 as
  surviving, and why the session-1 route was abandoned.

**Next**

- Anchoring KS's $r$ on the initial layer (would remove (K-esc), (K-trap)).
- Gap 3 via a certified sweep-out at a fixed $d > d_c$.
- An independent review of `research/calc/fold_proof.py`.

## 2026-10-29 (session 3) — The fold hypothesis (F), proved by computer

**Decisions**

- **Piecewise step schedule** for the validated integrator (`MU_SCHEDULE` in
  `research/calc/fold_proof.py`): low order with small $h/s$ near the axis,
  high order beyond. This replaced session 2's diagnosis of a structural
  barrier, which measurement showed to be wrong. Both the error and its
  correction are recorded in the module docstring, notes/09 and the idea log.
- **Acceptance criteria for a computer-assisted claim**, all met before
  claiming Theorem C:
  - known-solution tests;
  - jets checked against SymPy;
  - agreement with the independent floating-point code;
  - overlap of enclosures under a second step schedule;
  - literature agreement.
- **The certificate test is `slow`-marked** (about 5 minutes); everything else
  stays in CI.

**Progress**

- Theorem C (notes/09): unique fold point,
  $d_c = 0.766197428240 \pm 1.1\times10^{-13}\,M$; signs of $F_{zz}$, $F_d$,
  $F_t$; positive Jacobi field; embedded critical surface. Proposition B now
  depends only on standard analytic steps, written as sketches. The onset
  slope for geodesic slicing is rigorously in $[2.3768, 2.3936]$.

**Next**

- Write out the analytic steps of Proposition B. Then certify the
  $N = \psi^{-2}$ slope, and attempt a rigorous proof that no common MOTS
  exists at $t = 0$ (see `research/status.md`).

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
