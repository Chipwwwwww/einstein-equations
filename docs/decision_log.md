# Decision log

Newest first. Each entry: what was decided, why, and where things stand.

## 2026-10-31 (session 10) — Explicit binaries: route screening, Petrov obstruction, time series

**Decisions**

- **Starting point.** Research continues on `claude/black-hole-merger-solutions-rpa9hj`. This is
  the only research branch; `main` has the discussion notes only, and there were no PRs.
- **The interval sweep stops at $d = 0.767578$.** That is enough for Theorem H (note 13).
  Continuing to $0.78125$ costs about 18 more hours and is optional.
- **Routes screened before computing** (`research/route_map.md`, ten routes):
  - two were rejected by precise obstructions: algebraically special forms (proved here) and
    two-Killing integrable methods (literature);
  - one was executed (route E, the time series);
  - one was identified as most promising (J, near-critical BL);
  - the rest remain untested, each with its decisive test stated.
- **Two-vector Kerr–Schild was not attempted.** Its decisive test is a solvability analysis of
  an overdetermined linear system at $O(m_1m_2)$. It is recorded, not started.
- **Bugs found by validation, kept in the record.** The jet engine had two bugs:
  - a broadcasting error (constants added to every coefficient);
  - array aliasing ($\gamma_{\rho\rho}$ and $\gamma_{zz}$ shared storage).

  Both were found by the exact Schwarzschild comparison and the Hamiltonian-constraint check,
  before any BL number was used.
- **Monitoring pitfall.** mawk block-buffers piped input, so awk-based log watches missed stuck
  sweep boxes for an hour. Use `awk -W interactive`.

**Progress**

- Theorem F′ verified for $d \in (d_c, 0.767578125]$: notes/12 §7, the stored log, and the
  replayed band certificate. Theorem H (note 13) is no longer pending.
- New: note 14, `research/route_map.md`, `research/calc/petrov_bl.py`,
  `research/calc/bl_series.py`, tests `test_petrov_bl.py` and `test_bl_series.py`, ledger 64–75,
  a literature section on explicit binaries, and the status audit.

**Next**

- Route J: Arb jets for the route-E recursion plus a majorant tail. The goal is a certified
  convergence radius near $S_c$, then certified $\partial_t^2F$ and an explicit $\delta_0$.
- Optional: extend Theorem F′ to $0.78125$ (resume the sweep), joining it to Theorem F.

## 2026-10-30 (session 9) — Review of Lemmas 6–7; Theorem F toward d_c (in progress)

**Decisions**

- **Lemmas 6–7 re-checked step by step against Hintz 2026** (notes/11,
  review table). One correction: the case split at $t_* = 1$, and Lemma 7
  does use the support property, through Lemma 6. A self-review; an outside
  expert is still advised.
- **Gap 3 over an interval after all** (user request), reversing session 8's
  deprioritisation. The plan avoids the hundreds of CPU-hours estimated then:
  - One sweep over boxes in $(\delta, d)$ for $d \in [D_1, 0.78125]$,
    $D_1 = 0.766235$. It uses the relative height $z - d/2$, so $d$ enters
    only through the lower puncture, and Lohner sets with $d$ as a
    parameter.
  - A **waist event W** that hands medium necks ($\lambda \le 0.05$) to a
    $d$-independent **band certificate**: dyadic cells in $(\lambda, \zeta)$,
    each certified for all $d$ at once.
  - **A local fold lemma** for the last slab $d \in (d_c, D_1]$. At the fold
    point, $F = F_z = 0$ (Theorem C), so $F_{zz} > 0$ and $F_d > 0$ on a
    $z_0$-window give $F > 0$ there for every $d > d_c$. The sweep covers
    the rest of the slab.
- **How the fold-lemma enclosures were made tight** (`gap3_fold.py`). Each
  step was needed; each failure is recorded in the module docstring.
  - Block-diagonal Lohner frames, one block per jet slot. With a full QR
    frame the set blows up at σ ≈ 15.
  - $z_0$ as a dummy PSet parameter. Its first-order dependence then stays
    correlated across the jet blocks.
  - Near the plane, switch to $\tau$ with $z = z_{sw}(1-\tau)$. All curves
    cross at $\tau = 1$, and the jets there are the crossing derivatives.
  - Cross-check: a thin run reproduces Theorem C's $F_{zz}$ and $F_d$, which
    come from a different integrator.
- **Neck window widened to $\zeta \in [-0.234, -0.133]$.** The first sweep
  got stuck at $\delta \approx \pm0.0017$. A floating-point map of the necks
  shows why: the pinch is at $\delta = 0$ (the individual horizon itself),
  and the neck height moves off the old window as $|\delta|$ grows. The
  small-neck certificate is being extended; earlier exclusions stay valid
  (the window only grew), so the sweep resumes from its log.
- **`LAM_MAX` 0.04 → 0.05.** Band cells certify up to λ = 0.049, and the
  necks at $\delta \approx 0.0076$ have λ ≈ 0.039–0.041.
- **Pitfall recorded.** python-flint prints a ball as `[+/- |mid|+rad]` when
  the radius is not small relative to the midpoint. A positive enclosure
  [0.66, 0.80] printed as `[+/- 0.803]` and was misread for a while. Use
  `lower()` / `upper()`.

**Progress (when this entry was written)**

- Small-neck certificate over $d \in [0.765625, 0.78125]$, $\zeta \in [-0.2109, -0.1465]$:
  - λ ≤ 0.00488: 660 pieces, margin ≥ 0.511;
  - λ ∈ [0.00488, 0.0078]: 396 pieces, margin ≥ 0.297;
  - λ ∈ [0.0078, 0.0117]: 540 pieces, margin ≥ 0.0029.
  - Extensions to $\zeta \in [-0.234, -0.211]$ and $[-0.146, -0.133]$ are
    running.
- Fold lemma over $z_0 \in [0.62891, 0.63696]$, $d \in [d_c, D_1]$ (33 pieces):
  running; every piece so far passes with $F_{zz} \ge 170$ and $F_d \ge 0.65$.
- Main sweep (7,328 initial boxes) resumed; the slab sweep is queued after it.

**Next**

- Finish the runs. Then `python -m research.calc.gap3_interval verify gap3_interval.log`
  and `verify-slab gap3_slab.log` (cover checks plus the band certificate).
  Write up Theorem F′ in notes/12, with the ledger and status.

## 2026-10-30 (session 8) — Textbook black-hole region = escape-defined region (Λ = 0)

**Decisions**

- **Gap 3 for an interval of separations is deprioritised.** The cost
  estimate is hundreds of CPU-hours (idea log, session 8), and it would not
  locate the merger in time anyway.
- **Closed the cheapest open item in notes/11 instead.**
  - Localisation and connectedness for $\mathcal B_{\mathscr I}$ followed
    from extending Lemma 6 to all late times (Lemma 7). The 2026 paper
    itself records the needed decay of $h_{00}$ on the whole late region.
  - This gives Corollary E1, $\mathcal B_{\mathscr I} = \mathcal B$.
  - The session-6 plan, identifying $\partial\mathcal B_{\mathscr I}$ with
    a dynamical horizon, was not needed.
- **Kept the new Lohner tools, documented as tools.**
  - `lohner.PSet`/`step_param` keeps parameters out of the moving frame.
    It is tested.
  - `run_lohner(variational=True)` (event E4) did not help; it is kept as
    an option and is not used by the certificate.

**Progress**

- notes/11: Lemma 7, Corollary E1, the mass-mismatch identity (new exact
  test), intro and §7 table.
- `research/calc/lohner.py`: `Dl` (tangent numbers), `PSet`, `step_param`
  (new test).
- Ledger 62–63; status, literature (Hintz 2026 sections read), idea log,
  README.

**Next**

- An expert check of Lemmas 6–7 (our reading of Hintz 2026).
- Locating the merger event as ε → 0. This probably needs Hintz's
  intermediate-zone estimates to $o(\epsilon)$, because the big horizon's
  $O(\epsilon)$ teleological shift is $O(1)$ in rescaled units.

## 2026-10-30 (session 7) — Gap 3 closed at d = 0.78 (Theorem F, computer-assisted)

**Decisions**

- **Thin necks are excluded by instability, not followed.** Session 6
  planned isolating blocks at the axis saddles. Instead:
  - A certified thin-neck lemma: near a thin waist the surface is a
    catenoid up to $O(\lambda)$. Its $\lambda$-variation Jacobi field
    changes sign on both sides, so the surface is unstable. The outermost
    surface is stable, so it has no such neck.
  - An analytic approach lemma hands any curve heading into the axis to the
    neck lemma.

  This removes the heteroclinic-cycle and multi-layer analysis entirely.
- **Lohner's method for the sweep** (`research/calc/lohner.py`). Plain
  interval boxes blow up on curves that wind around the individual horizon.
- **Medium necks by brute force, with a stored partition.**
  - Necks of radius 0.005–0.02 are too wide for the certified lemma, and too
    narrow to integrate through cheaply.
  - Three fixes failed: mid-run splitting, a wider neck certificate with
    plain intervals, and a parametric Lohner neck certificate.
  - What worked: calibrated $z_0$-intervals in that window, stored in
    `research/calc/data/gap3_partition_d078.txt`, so that the certificate
    is one reproducible run.
- **Deliberate falsification before trusting it.**
  - At $d = 0.76 < d_c$, the code does not exclude the real common
    surfaces.
  - The neck certificate fails where it must ($\sigma = 1.0$).
  - Two integrators agree on the neck Jacobi field.
- **Novelty is stated narrowly.** Dain–Weinstein–Yamada (2011) already
  prove non-existence for large separation. Our result is the explicit
  near-critical case ($d/d_c = 1.018$), and only for one $d$.

**Progress**

- notes/12 (Theorem F); `lohner.py`, `gap3_neck.py`, `gap3_sweep.py`,
  `gap3_neck_lohner.py` (cross-check only); `tests/test_lohner.py`,
  `tests/test_gap3.py`.
- Thin-neck certificate: 640 pieces, margin ≥ 0.527, 25 min. CI reruns it
  as a slow test.
- Sweep: 2,452 intervals, all excluded, 2.8 h on 4 cores; the log is in
  `research/calc/data/`.
- Ledger 56–61; status, idea log, literature, notes/10 and the READMEs
  updated.

**Next**

- Gap 3 for all $d \in (d_c, 0.78]$ (status.md (b)). That joins Theorem F
  to the onset theorem.
- A parametric Lohner variant ($x = m + Ar + P(p - p_0)$) to certify the
  medium necks cheaply.
- Independent review of the new code.

## 2026-10-29 (session 6) — Gap 3 analysed; merger-timing corollary

**Decisions**

- **Gap 3 is not closed this session, by choice.** The sweep reduces to a
  window of axis heights, but just above the individual horizon the
  shooting curves pass necks of every size. A proof needs an isolating-block
  argument at the axis, which is multi-session work. We certified the two
  barrier lemmas (they are needed in any version) and recorded the
  obstruction, the scaling law and the failed shortcuts, so the next attempt
  starts from there.
- **Locating the merger event as ε → 0 was not attempted.** It needs light
  rays over long rescaled times. We recorded instead the cheap rigorous
  statement that on the slow scale the merger happens at the crossing
  (Corollary D1).

**Progress**

- `research/calc/bl_global.py` and `tests/test_bl_global.py`:
  - Lemmas G1 and G2 (interval arithmetic);
  - the floating-point landscape;
  - the desingularized $\sigma$-system, checked by SymPy.
- Corollary D1 in notes/11.
- Ledger 50–55; status and idea log updated.

**Next**

- Isolating block for the neck passage (gap 3).
- An expert check of note 11's Lemma 6.

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
  - Theorem E uses $M \setminus J^-(\mathscr I^+)$. Reaching
    $\mathscr I^+$ needs retarded-time control. Lemma 6 gets it from Hintz
    2026: the outgoing null component is integrable (Lemma 3.19,
    Definition 5.3, and the support properties of the proof of Theorem
    13.1). This is labelled as our reading.
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

- Expert check of Lemma 6. Locating the merger event (status.md, (a)).
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
