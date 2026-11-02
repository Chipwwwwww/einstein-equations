# 26 — Full 50-field GH principal structure, Schwarzschild geometry, and the obstruction map

**2026-11-02. Outcome:** verified symbolic 50-field generalized-harmonic (GH) principal symmetrizer, a valid exact Schwarzschild Kerr–Schild gauge source, and an exact-rational Lorentzian geometric tube on a curved shell. **NOT** a complete full-Einstein 50-field a posteriori evolution, constraint, or merger certificate. Read [the general H^3 theorem](25-full-hs-a-posteriori.md) first.

## 1. Exact GH principal structure

Primary source: Lindblom–Scheel–Kidder–Owen–Rinne, *A new generalized harmonic evolution system*, Class. Quantum Grav. 23 (2006) S447–S462, eqs. (27)–(39), https://arxiv.org/abs/gr-qc/0512093 (DOI 10.1088/0264-9381/23/16/S09). Full 50 dynamical fields are ten symmetric spacetime pairs $(ab)$, each with $(\psi_{ab},\Pi_{ab},\Phi_{1ab},\Phi_{2ab},\Phi_{3ab})$.

Fix their parameters **exactly** at $\gamma_1=-1$, $\gamma_2=0$; the algebraic $\gamma_0$ constraint damping terms are NOT evaluated here. Given positive lapse $N$, shift $\beta^i$, positive inverse spatial metric $\gamma^{ij}$, the principal matrices per 5-field pair are

\[
A^k=\begin{pmatrix}
0&0&0_{1\times3}\\
0&-\beta^k&N\gamma^{kj}\\
0_{3\times1}&N\delta_i{}^k&-\beta^k\delta_i{}^j
\end{pmatrix},\qquad
S_5=\operatorname{diag}(1,1,\gamma^{ij}),
\]
and $S_{50}=I_{10}\otimes S_5$, $A^k_{50}=I_{10}\otimes A^k$. For every spatial $k$,
\[
S_{50}A^k_{50}=(A^k_{50})^T S_{50},\qquad
\lambda_{\min}(S_{50})=\min(1,\lambda_{\min}(\gamma^{ij}))>0.
\]
These identities are independently checked on 50-by-50 SymPy matrices. A negative control rejects using the unweighted identity symmetrizer, which fails on a nonflat test point. There is no arbitrary field-number restriction in the abstract theorem in note25; however, the 50-field **algebraic sources** are not generated or bounded in this code.

The implementation additionally checks **nonzero $\gamma_2$**: with $\gamma_1=-1$, $\gamma_3=\gamma_1\gamma_2=-\gamma_2$, add $A^k_{\Pi\psi}=\gamma_2\beta^k$ and $A^k_{\Phi_i\psi}=-\gamma_2N\delta_i{}^k$. Use the genuine $S_{\psi\psi}=\Lambda^2$, $S_{\psi\Pi}=-\gamma_2$, $S_{\Pi\Pi}=1$ block with **$\Lambda^2>\gamma_2^2$**. A 50-field exact test at $\gamma_2=1/4$, $\Lambda^2=1$ passes, and an invalid $\gamma_2=3/2$ fails. The positive cross term condition is essential; this extends the *principal* check but does not implement the constraint-damping source terms. Nor should one substitute the repo's coordinate Einstein–Bianchi system of note19; that system was found to contain a Jordan block.

## 2. Valid generalized harmonic gauge source for nonflat Schwarzschild

In Cartesian horizon-penetrating Kerr–Schild coordinates, define $r=|\mathbf x|$, $h=2M/r$, $\ell_a=(1,x/r,y/r,z/r)$, $\eta=\operatorname{diag}(-1,1,1,1)$:
\[
g_{ab}=\eta_{ab}+h\ell_a\ell_b,\quad
g^{ab}=\eta^{ab}-h\ell^a\ell^b,\quad \det g=-1.
\]
The **fixed coordinate gauge source** (independent of unknown field derivatives) is
\[
H^a_{KS}=-\Gamma^a(g_{KS})
=\frac h r(1,-x/r,-y/r,-z/r),\qquad
H_a^{KS}=-\frac h r(1,x/r,y/r,z/r).
\]
The explicit identity $\Gamma^a=-\partial_b g^{ab}$ follows from $\det g=-1$; the code differentiates $g^{ab}$ to verify $H^a+\Gamma^a=0$ exactly. Also $\ell_aH^a=0$. This is a legal GH gauge source prescribed as a function of coordinates, away from $r=0$, and is a legitimate background gauge check. **It is not an evolution residual-to-vacuum certificate.**

At the Schwarzschild horizon $r=2M$ and $x=r,y=z=0$, the ADM values are $N=1/\sqrt2$, $\beta^x=1/2$, $\gamma^{ij}=\operatorname{diag}(1/2,1,1)$, so the actual 50-field GH symmetrizer eigenvalue lower bound is $1/2$. One 5-field $A^x$ block has characteristic eigenvalues $0,0,-1,-1/2,-1/2$. This checks the principal hyperbolic cone at a nonflat horizon-penetrating point. The source/constraint subsystems remain unverified.

## 3. A continuum, exact-rational geometric tube on a curved shell

Choose $M=1$, $3/2\le r\le5/2$, and a **hypothesized** pointwise metric error bound
$\|g-g_{KS}\|_{\mathrm{op}}\le 6R$ with $R=1/100$. A correctly justified global $H^3$ error bound could imply this using the Fourier embedding constant 6 from note25, but **the shell is not a torus**: we do not infer such a bound on the shell without a separate Sobolev extension/boundary theorem.

Over the entire continuum shell, $h\in[4/5,4/3]$. Hence the reference inverse norm is at most $B=1+2h_{\max}=11/3$ and the Neumann multiplier is
\[
z=B(6R)=\frac{11}{50}<1.
\]
The inverse metric is enclosed by
\[
\|g^{-1}-g_{KS}^{-1}\|_{\rm op}
\le \frac{B^2(6R)}{1-z}=\frac{121}{117}.
\]
Along the linear path $g_{KS}+\theta(g-g_{KS})$, $\theta\in[0,1]$, all matrices are nonsingular; their symmetric inertia is therefore unchanged. The entire pointwise tube is **Lorentzian**, not merely the central approximate metric. Spatial eigenvalues and the GH symmetrizer are bounded by
\[
\lambda_{\min}(\gamma)\ge\frac{47}{50},\quad
\lambda_{\max}(\gamma)\le\frac{359}{150},\quad
\boxed{\lambda_{\min}(S_{50})\ge\frac{150}{359}}.
\]
Also $448/585\le-g^{00}\le394/117$, ensuring a strictly positive ADM lapse. Exact-rational code and explicit failed-wide-tube tests reproduce these inequalities. They do not provide $\|u-v\|_{H^3}$ or any dynamic existence conclusion.

## 4. Exact gauge, reduction and Einstein constraints: unfinished

The general GH gauge constraint is $\mathcal C_a=H_a+\Gamma_a$, with first-order reduction constraints $\mathcal C_{iab}=\partial_i\psi_{ab}-\Phi_{iab}$ and $\mathcal C_{ijab}=\partial_i\Phi_{jab}-\partial_j\Phi_{iab}$. Hamiltonian and momentum constraints are additional required initial identities.

The contracted Bianchi identities yield a homogeneous **subsidiary hyperbolic system** for the gauge constraints if the **complete** reduced GH equations are used. Schematically, for undamped second-order GH,
\[
\Box_g\mathcal C_a+R_a{}^b\mathcal C_b=0
\]
(up to the chosen sign and convention, and the first/zero-order homogeneous damping terms when included). Exactly zero gauge constraints and their normal derivatives on an initial Cauchy surface, plus exact reduction constraints and compatible boundary fields, imply zero subsidiary solutions by uniqueness. Small gauge residuals and numerical damping are **not** sufficient to infer $\mathrm{Ric}(g)=0$.

**No complete 50-field nonlinear source, its $H^3$ composition bounds, Einstein initial-constraint check, or exact subsidiary-system computer-assisted proof was produced.** Thus we have **not** applied the note25 existence theorem to all GH Einstein equations. The Schwarzschild gauge-source identity is an independent exact benchmark only.

## 5. Missing hypotheses and precise next theorem

1. **Full 50-field GH coefficient/source evaluator** with finite-dimensional interval-certified inverse metric/lapse, coefficient derivatives, $H^3$ Moser constants over a verified tube, and a continuum residual $r_{GH}$.
2. **Admissible Cauchy domain:** the three-torus theorem of note25 does not directly apply to an annular black-hole exterior. Establish a boundaryless causal diamond by finite propagation *with its own local energy estimates*, or prove maximally dissipative, constraint-preserving boundary estimates up to $H^3$ with compatible incoming and excision characteristic fields. A principal symmetrizer is not a high-order boundary theorem.
3. **Constraint-exact initial GH data** including $\mathcal C_a=0$, normal derivative, and $\mathcal C_{iab}=\mathcal C_{ijab}=0$. Match the geodesic-slicing Brill–Lindquist metric to a GH chart with rigorous overlaps and inverse maps.
4. **Extended slabs and physical merger geometry:** certify nonlinear $H^3$ residuals, metric errors, boundary/excision geometry, curvature, gauge, and eventually the wave-extraction exterior. A MOTS birth-time bracket alone does not perform this task.
5. **Long-time estimates:** reduce globally pessimistic $K$ by localized/adapted energy, redshift and geometric multipliers *only after proving their hypotheses*; small-data Minkowski/Kerr estimates cannot be assumed through comparable-mass merger.

Reproduce the exact geometric and principal tests:

    python -m pytest -q tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py
    python -m research.calc.gh_ks_geometry_tube

**Classification of this note alone:** verified curved-geometry/principal-symbol result, **not** a verified full-Einstein a posteriori existence theorem or a nonflat 3+1 dynamical error certificate.
