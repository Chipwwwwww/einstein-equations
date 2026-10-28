# Decision log

Newest first. Each entry: what was decided, why, and where things stand.

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
