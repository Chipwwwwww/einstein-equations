# Einstein's Equations

A shared workspace for discussing the Einstein field equations of general
relativity: what they say, where they come from, what their solutions mean, and
what is still open. It has three parts: notes you can argue with, a
small SymPy toolkit that checks the math, and GitHub Discussions for the arguing.

```math
G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}\, T_{\mu\nu}
```

Curvature of spacetime on the left, energy and momentum of matter on the right.

## What's here

| Path | Contents |
| --- | --- |
| [`notes/00-conventions.md`](notes/00-conventions.md) | Sign conventions, index notation, units. Read this first when comparing with a textbook. |
| [`notes/01-the-field-equations.md`](notes/01-the-field-equations.md) | The equation term by term, counting degrees of freedom, trace-reversed form, the action, uniqueness. |
| [`notes/02-geometry-toolkit.md`](notes/02-geometry-toolkit.md) | Metric, Christoffel symbols, Riemann, Ricci, Bianchi identities. |
| [`notes/03-newtonian-limit.md`](notes/03-newtonian-limit.md) | How Newton's gravity falls out, and where it stops working. |
| [`notes/04-exact-solutions.md`](notes/04-exact-solutions.md) | Schwarzschild, de Sitter, FLRW (Friedmann equations), Kerr. |
| [`notes/05-gravitational-waves.md`](notes/05-gravitational-waves.md) | Linearized gravity, polarizations, the quadrupole formula. |
| [`notes/06-open-questions.md`](notes/06-open-questions.md) | Discussion starters: energy, Λ, uniqueness, singularities, quantum gravity. |
| [`einstein/`](einstein/) | SymPy code: Christoffel symbols, Riemann, Ricci, Einstein tensor, Kretschmann scalar. |
| [`tests/`](tests/) | Checks that the solutions in the notes really solve the equations. |

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pytest                           # verify every solution in notes/04
python -m einstein schwarzschild # print Christoffels, Ricci, Einstein tensor
```

`python -m einstein --help` lists the built-in metrics. To try your own metric:

```python
import sympy as sp
from einstein import Spacetime

t, r, th, ph = sp.symbols("t r theta phi")
A = sp.Function("A")(r)
st = Spacetime(sp.diag(-A, 1 / A, r**2, r**2 * sp.sin(th) ** 2), (t, r, th, ph))
st.einstein  # G_{mu nu} as a SymPy matrix
```

## How to take part

- **Questions and ideas** go in [Discussions](../../discussions): Q&A for "why does…",
  Ideas for new notes or topics, General for open-ended debate.
  [`notes/06-open-questions.md`](notes/06-open-questions.md) has starters.
- **Mistakes** in a note or the code go in [Issues](../../issues).
- **New notes or derivations** come in as pull requests. See
  [CONTRIBUTING.md](CONTRIBUTING.md).

## Conventions in one line

Signature $(-,+,+,+)$, Riemann and Einstein tensors as in Misner–Thorne–Wheeler
and Carroll, Greek indices $0..3$, Latin $1..3$, $x^0 = ct$. The code uses
geometrized units $G = c = 1$ unless stated. Details in
[`notes/00-conventions.md`](notes/00-conventions.md).

## References

- A. Einstein, "Die Feldgleichungen der Gravitation", *Sitzungsber. Preuss. Akad. Wiss.* (1915) 844–847.
- S. Carroll, *Spacetime and Geometry* (2004). Same conventions as this repo.
- C. Misner, K. Thorne, J. Wheeler, *Gravitation* (1973).
- R. Wald, *General Relativity* (1984).
- B. Schutz, *A First Course in General Relativity* (2nd ed., 2009).
