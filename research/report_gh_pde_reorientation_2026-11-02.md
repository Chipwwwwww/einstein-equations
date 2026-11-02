# PDE-first research checkpoint (2026-11-02)

**Parent inspected:** branch research/global-mots-audit-2026-11-01, commit fa5fbce2fc7fd6f071d3cff55af193789f62862b (PR #7). This inherited branch explicitly does not certify a global first-common-MOTS bracket or a global merger spacetime. This new branch leaves that work untouched.

**Primary mathematical question:** Can a smooth approximate generalized-harmonic metric with a certified continuum residual be upgraded, quantitatively and constructively, to a true Einstein vacuum metric across a strong-field slab? Directly relevant foundations: Hughes–Kato–Marsden local quasilinear Einstein well-posedness, Kato symmetric hyperbolic estimates, Lindblom et al. first-order GH and constraint propagation, and special exact gluing theorems of Hintz. Full finite-slab numerical-merger certificate remains open.

**Chosen missing hypothesis:** guaranteed continuum residual and slab-to-slab stability with explicit, computable constants. Grid point residual tests do not provide this, and generic high-Sobolev Newton inverses may lose derivatives.

**Completed proofs (not novelty claims):**

1. Explicit resonant Fourier family disproves bounded fixed-space wave inversion from $L^\infty_tH^{s-1}_x$ to $L^\infty_tH^{s+1}_x$.
2. A full 4D harmonic 2+2 metric subsector reduces Ricci-flatness to a scalar linear wave equation in $\log H$; a null-characteristic a posteriori theorem bounds metric, inverse and first derivatives by a recurrence with exactly evaluated constants.
3. A conditional symmetric-hyperbolic $L^2$ error inequality gives a standard computable Gronwall formula **provided** genuine-solution derivative bounds are separately available. It is **not** an existence theorem.

**Reproducible computer-assisted result:** For $H_{\rm app}=1+(x-t)/2+t^2/100$ on $0\le t\le1/4,\ |x|+t\le1$, exact rational enclosures give $H_{\rm app}\in[1/2,2401/1600]$, residual at most $3/25$, $\|\log H-\log H_{\rm app}\|_\infty\le3/800$, and metric-component discrepancy at most $7203/1275200$. The independent comparator has actual maximum error $1/1600$. The certificate holds on the **continuum domain**, not merely at mesh points. Five local tests passed (exact fractions, independent symbolic four-dimensional Ricci, harmonic gauge, slab splitting, Schwarzschild Kerr–Schild horizon sanity check).

**Precise scope:** The gauge-wave sector is flat vacuum in nontrivial harmonic coordinates. The Schwarzschild control is a separate exact spacetime, **not** a GH merger validation. The result does not handle all coupled GH metric fields, nontrivial strong curvature evolution, constraint-preserving excision, a binary merger, or a proof of a global exact metric. No extension of the certified BL evolution time is claimed.

**Next missing bridge:** explicit high-order commuted-energy and tame estimates, positive symmetrizer/metric tube, gauge/reduction subsidiary equations, incoming boundary fields and certified interpolation bounds from actual GH merger data. These need a proof, not a formal invocation of a semilinear or small-data theorem.

Proof and full hypothesis/limitation matrix: [notes/23-a-posteriori-gh-slab.md](../notes/23-a-posteriori-gh-slab.md).

Run from repository root:

    python -m research.calc.gh_slab_certificate
    python -m pytest tests/test_gh_slab_certificate.py -q


## Follow-on curvature-carrying proof (same session)

The initial flat-gauge test was extended to a **nonflat cylindrical vacuum**
Einstein–Rosen metric, including nonlinear constraint reconstruction. The full
separate proof and exact rational certificate are in
[notes/24-er-cylindrical-slab.md](../notes/24-er-cylindrical-slab.md).

On the rational annular causal diamond the radial-wave residual equals 1/50,
the quadratic gamma-time constraint residual is <=11/2500, and exact
characteristic/constraint bounds give psi error <=11/9800, gamma error
<=70079/24010000, and all relevant metric-component errors <=44/4889 (<0.01).
The independent analytic comparator has psi error <=1/1600 and gamma error
<=21/20000. SymPy verifies the independent full 4D Ricci tensor vanishes and
a nonzero Riemann component R^z_{tzt}=-209/1000 at a rational point. The
symmetry reduction is not a generic GH metric, and no 4D asymptotically flat
merger is certified.

    python -m research.calc.er_slab_certificate
    python -m pytest -q tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py

**Actual combined local result:** 8 passed (Python 3.13.5 / SymPy 1.14).
The exact rational continuum enclosures do not rely on mesh sampling.
Full repository regression/remote CI still pending.
