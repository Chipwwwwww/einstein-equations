"""Print the curvature of a named metric, e.g. `python -m einstein schwarzschild`."""

import argparse

import sympy as sp

from .metrics import REGISTRY
from .tensors import nonzero


def _section(title, rows):
    print(f"\n{title}")
    for label, value in rows:
        print(f"  {label} = {sp.sstr(value)}")
    if not rows:
        print("  (all zero)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m einstein",
        description="Print Christoffel symbols, Ricci and Einstein tensors of a metric (G = c = 1).",
    )
    parser.add_argument("metric", choices=sorted(REGISTRY))
    parser.add_argument("--kretschmann", action="store_true", help="also print R_abcd R^abcd")
    args = parser.parse_args(argv)

    st = REGISTRY[args.metric]()
    x = [sp.sstr(c) for c in st.coords]

    _section("Metric g_ab", [(f"g_{x[a]}{x[b]}", v) for (a, b), v in nonzero(st.metric) if a <= b])
    _section(
        "Christoffel symbols Gamma^a_bc (b <= c)",
        [(f"Gamma^{x[a]}_{x[b]}{x[c]}", v) for (a, b, c), v in nonzero(st.christoffel) if b <= c],
    )
    _section("Ricci tensor R_ab", [(f"R_{x[a]}{x[b]}", v) for (a, b), v in nonzero(st.ricci) if a <= b])
    _section("Ricci scalar", [("R", st.ricci_scalar)] if st.ricci_scalar != 0 else [])
    _section("Einstein tensor G_ab", [(f"G_{x[a]}{x[b]}", v) for (a, b), v in nonzero(st.einstein) if a <= b])
    if args.kretschmann:
        _section("Kretschmann scalar", [("K", st.kretschmann)])


if __name__ == "__main__":
    main()
