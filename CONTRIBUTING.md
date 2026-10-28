# Contributing

This repo is for discussion first, so there are three ways in.

## 1. Ask or argue: Discussions

Use [Discussions](https://github.com/Chipwwwwww/einstein-equations/discussions) for anything that is a conversation rather
than a fix:

- **Q&A**: "Why is the Einstein tensor divergence-free?", "Where does the 8πG come from?"
- **Ideas**: proposals for new notes, worked examples, or code.
- **General**: open-ended debate. [`notes/06-open-questions.md`](notes/06-open-questions.md) has starters.

When you reference a formula, link to the line in the note
(click the line number on GitHub, then "Copy permalink").

## 2. Report a mistake: Issues

Open an [Issue](https://github.com/Chipwwwwww/einstein-equations/issues) for a sign error, a wrong factor, an unclear step,
or a failing test. Say which note and which equation, and what you expected.

## 3. Add something: pull requests

### Notes

- One topic per file in `notes/`, numbered in reading order (`07-…`).
- Follow [`notes/00-conventions.md`](notes/00-conventions.md). If a source uses
  other conventions, convert and say so.
- Display math in fenced blocks, which GitHub renders reliably:

  ````markdown
  ```math
  R_{\mu\nu} - \tfrac12 R\, g_{\mu\nu} = 8\pi G\, T_{\mu\nu}
  ```
  ````

  Inline math as `$G_{\mu\nu}$`. Avoid `$$`.
- State units: keep $G$ and $c$, or say "in $G = c = 1$".
- End with a **Discussion prompts** section: two to four questions the note raises.

### Code and tests

If a note claims something SymPy can check (a metric solves the equations, a
limit gives Newton, an identity holds), add a test in `tests/`. Run:

```bash
pip install -r requirements.txt
pytest
```

Tests run on every push and pull request via GitHub Actions.
