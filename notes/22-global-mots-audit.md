# 22 · Global MOTS exclusion: corrected shooting and stability tools

*2026-11-01. Units $G=c=M_{\rm ADM}=1$. This note does not prove global
MOTS formation or extend the certified vacuum spacetime. It repairs a
proposed exclusion argument and supplies analytic tools for its replacement.*

The exact status and replay record are in
[the session report](../research/report_global_mots_audit_2026-11-01.md).
The local symmetry reduction already appears in note 10, Step 3; it is not
a new result of this session. Likewise, [note 21](21-horizon-spectral-gap.md)
is an integration of an independently developed branch, not a new discovery.

## 1 · What the audit changes

The barrier theorem used in note 20 has the correct sign and smoothness
hypotheses: inward-to-the-annulus normal on its inner boundary, outward on
its outer boundary, and strict negative/positive outgoing expansion.
Andersson–Metzger [AM09, Theorems 3.1 and 4.1] give a stable embedded MOTS
homologous to the outer boundary. No constraint equation is required for
this elliptic existence result. Calling the data **vacuum slice data** is
a separate obligation, discharged locally by the inherited ADM/CK argument.

For BL, $\Delta\psi=0$ off the punctures and $K=0$ imply the initial
constraints. The ADM evolution and contracted Bianchi identity give a
homogeneous constraint system. Analyticity and induction on time derivatives
then give identically zero constraints on the analytic development. A
finite-order residual test alone would not establish this fact. Embedding
in the maximal development uses local uniqueness on a domain of dependence;
it is not an assertion of a global Gaussian coordinate chart.

The two barrier spheres bound an annulus in the punctured coordinate space.
An embedded separating component of a homologous MOTS surrounds the inner
sphere and hence both punctures. This establishes a **common MOTS**. It
does not identify that surface with the outermost apparent horizon of the
entire slice. That identification additionally requires the exterior
trapped-region construction and control of all relevant ends.

The convex-ridge comparison in note 20 is sound for surfaces bounding the
specified interiors. At the outermost contact, normals agree. Mean
curvature of an inner touching surface is at least that of the leaf, and
$\operatorname{tr}_S K$ agrees because the tangent plane agrees. At a genuine
convex ridge a smooth inner surface would have to be tangent to two distinct
planes. This proof needs no axial symmetry of the candidate MOTS. It
excludes surfaces wholly in the shell, not surfaces crossing out of it.

### 1.1 Arithmetic and statement corrections

The old `bl_barrier.main` passed the decimal time to the Moore sign
calculation, but `float(t)` to the Lohner geometry calculation. At 128 bits
`arb('6.5e-6')` and `arb(float('6.5e-6'))` are disjoint. The corresponding
problem also occurs for `'1.2e-5'`. Thus the two calculations did not
literally enclose the same parameter, even though their difference is far
smaller than the displayed geometric margins. Continuity without a bound
does not repair a certificate handoff.

The code now passes the original parameter and encloses it by outward
binary64 endpoints. Pole-height brackets receive the same treatment.
Tail evaluation times round upward; angular bin assignment uses interval
division; radial gap subtraction and domain-bin endpoints are outward
enclosed. Disjoint supposed enclosures of one quantity now raise errors
instead of silently choosing one. Old JSON records are retained as history.
Only the explicitly recorded corrected replays upgrade their status.

Theorem N's printed pole-height endpoints had been rounded **outward** from
the tested window. The window is defined by the exact binary64 endpoints
in its JSON record; the shorter decimal endpoints are illustrations only.
The assertion that M' computes the entire interval $(d_c,d_c+\delta_0)$ is
also too strong: its tested interval starts at the fixed number
$d_\ell=0.76619742824011$. The small missing interval cannot be removed
merely by citing an unquantified local existence theorem.

## 2 · Why constant residual levels cannot give global exclusion

Write $\Theta=\Theta^{(1)}+R$. A true MOTS satisfies
$\Theta^{(1)}(S)(x)=-R(S)(x)$. The right side can vary along the meridian
and depends on its tangent and curvature. A sweep over **constant**
$h\in[-\varepsilon,\varepsilon]$ does not include all such surfaces, even
if $|R|\le\varepsilon$. Moreover the existing bound on $R$ is evaluated
on the chosen level surfaces, not on arbitrary candidate surfaces with
uncontrolled curvature. Note 20's existence and foliation uses do not
make this error; its proposed global shooting extension did.

**Falsification example.** For every $\varepsilon>0$ consider

```math
y''+y=h(s),\qquad y(0)=0,\quad y'(0)=\pi\varepsilon,
\qquad F=y'(2\pi).
```

Every constant $h$, of any size, gives $F=\pi\varepsilon>0$. The admissible
smooth nonconstant forcing $h=-\varepsilon\cos s$ gives
$y=\pi\varepsilon\sin s-(\varepsilon/2)s\sin s$ and $F=0$.
Thus a constant-parameter sweep can exclude every tested trajectory while
missing an allowed varying residual. This is a counterexample to the
inference, not evidence of an earlier BL horizon.

There are two sound replacements: a validated differential inclusion with
arbitrary admissible forcing, or the **true evolved-metric ODE** below.
A constant level remains entirely legitimate when it constructs actual
barrier leaves and a comparison theorem, rather than enumerating MOTSs.

## 3 · Exact shooting in an evolved, nonconformal metric

Use note 20's notation: $T=(\cos\alpha,0,\sin\alpha)$,
$n_i=(-\sin\alpha,0,\cos\alpha)$ (a covector), $n^i=\gamma^{ij}n_j$, $q_n=\sqrt{\gamma^{ij}n_i n_j}$,
$g_T=\gamma(T,T)$, and lower Christoffels $\Gamma_{i,jk}$.
Set $C_T=n^i\Gamma_{i,jk}T^jT^k$ and $C_y=n^i\Gamma_{i,yy}$.

**Proposition GS-ODE.** On $\rho>0$, in a positive definite rotationally
and meridionally reflection-invariant metric, a surface of revolution is
a MOTS with the chosen outward normal precisely when

```math
\begin{aligned}
\rho'&=\cos\alpha, & z'&=\sin\alpha,\\
\alpha'&=-C_T-\frac{g_T}{\gamma_{yy}}
                 \left(\frac{\sin\alpha}{\rho}+C_y\right)
       +q_n\left(K(T,T)+\frac{g_T}{\gamma_{yy}}K_{yy}\right).
\end{aligned}                                                    \tag{1}
```

*Proof.* Note 20's expansion is affine in $\alpha'$, with coefficient
$-1/(q_n g_T)<0$. Solving that exact identity gives (1); no division by
$K$, time, or a small expansion residual occurs. The other two equations
are the definition of flat arclength. Conversely substitution into the
expansion gives zero. ∎

At the upper regular axis pole, rotational smoothness gives
$\gamma_{xx}=\gamma_{yy}=a$, $\gamma_{zz}=b$, $\gamma_{xz}=0$,
$K_{xx}=K_{yy}$ and $\Gamma_{z,xx}=\Gamma_{z,yy}$.
Writing $\alpha(s)=c_0s+O(s^3)$ and $\rho=s+O(s^3)$ gives

```math
c_0=\frac{K_{xx}}{\sqrt b}
     -\frac{\partial_x\gamma_{xz}-\tfrac12\partial_z\gamma_{xx}}{b}.
                                                               \tag{2}
```

Indeed $\sin\alpha/\rho\to c_0$, so (1) contains $c_0$ on both sides
and gives (2). In conformal data it reduces to
$c_0=2\psi_z/\psi+\psi^2\operatorname{tr}_S K/2$, as required.
Equation (2) does **not** validate the omitted polar collar; a regular
axis integral equation and an inclusion test are still required.

`research/calc/mots_exact.py` implements these identities. It accepts
Arb data enclosures and checks against the existing expansion formula,
its conformal limit, anisotropic cross terms and a true CK-evaluator point.
The point test proves no surface existence or global-domain coverage.

## 4 · Removing the evolving stability drift in this symmetry class

Let $S$ be a smooth, connected, embedded MOTS sphere invariant under
rotations and reflections through meridian planes, in slice data sharing
those symmetries. Its spatial-deformation stability operator, with the
repository convention $K=+\tfrac12\dot\gamma$, has standard form

```math
L=-\Delta_q+2\omega^A D_A+
       (Q+\operatorname{div}_q\omega-|\omega|_q^2),
\qquad \omega_A=K(\nu,\partial_A).
                                                               \tag{3}
```

Here $Q=R_q/2-\mu+J(\nu)-|\chi^+|^2/2$, with
$J=\operatorname{div}K-d\operatorname{tr}K$ as in [AM09, §2];
$\chi^+=A+K|_{TS}$ and $\mu=J=0$ in vacuum. Equivalently $Q$ can be
defined by the zeroth-order term of (3). The result concerns this spatial
variation operator, not a different null deformation operator.

**Theorem GS-S (exact drift transformation).** Under these hypotheses
$\omega=d\chi$ for a smooth real, symmetry-invariant function on $S$, and

```math
e^{-\chi}Le^{\chi}=A_q:=-\Delta_q+Q.                    \tag{4}
```

Consequently $L$ is self-adjoint in $L^2(e^{-2\chi}dA_q)$, is similar to
$A_q$ on ordinary $L^2$, and has a real spectrum. This is a structural
result for the exact metric, not a small-time truncation.

*Proof.* Away from the poles an invariant one-form is
$a(\theta)d\theta+b(\theta)d\varphi$. Meridian reflection reverses
$d\varphi$ while preserving the other data, so $b=0$ and $d\omega=0$.
Smoothness extends this identity across the poles. Since
$H^1(S^2;\mathbb R)=0$, $\omega=d\chi$ globally. Normalize $\chi$ to
have mean zero; invariance follows by uniqueness of this primitive.
Expanding $\Delta(e^\chi f)$ in (3) cancels the drift and both
$\Delta\chi$ and $|D\chi|^2$, giving (4). Multiplication by $e^{-\chi}$
is unitary from the weighted space to $L^2(dA_q)$, proving the assertions. ∎

**Limitations and prior art.** General axisymmetric MOTS shooting and its
Jacobi/stability interpretation are already developed in [BHP21, §§III–IV,
equations (41), (81)–(83)]. Equation (1) specializes that geometry to this
repository's flat-arclength parameter and metric evaluator; it is not a new
shooting principle. This is a standard exact-drift/gauge
transformation, not a new spectral principle. Nonspinning head-on MOTS
spectra were studied in [PK20]. Axial symmetry by itself is insufficient:
$\sin^2\theta\,d\varphi$ is smooth and rotation-invariant but has nonzero
exterior derivative. Nor can (4) be applied to a nonsymmetric surface
without establishing the hypotheses. Local symmetry of all small graphs
is already proved in note 10. Global outermost symmetry requires the
trapped-region theorem and its exterior/barrier hypotheses.

### 4.1 A conditional all-mode estimate with checkable inputs

Choose any smooth positive $\beta$ on $S$ and a diffeomorphism from the unit
round sphere $(\mathbb S^2,\sigma)$ to $S$. Suppose, after pullback,

```math
\beta^2\frac{dA_q}{d\Omega}\,q^{-1}\ge a\,\sigma^{-1},\qquad
\beta^2\frac{dA_q}{d\Omega}\le b,\qquad
\left|\frac{A_q\beta}{\beta}\right|\le\varepsilon,
\quad a,b>0.                                             \tag{5}
```

**Corollary GS-G.** Index eigenvalues starting at zero. For every angular
mode, without a truncation,

```math
-\varepsilon\le\lambda_0(L)\le\varepsilon,\qquad
\lambda_k(L)\ge(a/b)\mu_k-\varepsilon,\qquad
\lambda_1(L)-\lambda_0(L)\ge2a/b-2\varepsilon.             \tag{6}
```

Here $\mu_k$ are the round-sphere Laplacian eigenvalues. On
$\beta^\perp\subset L^2(dA_q)$ the projected operator $P A_q P$ obeys
$P A_q P\ge(2a/b-\varepsilon)I$. This is an inverse bound for the
**projected** operator; $\beta^\perp$ need not be invariant under $A_q$.
For $L$, use the corresponding weighted projection onto
$(e^\chi\beta)^\perp$ in $L^2(e^{-2\chi}dA_q)$.

*Proof.* With $r=(A_q\beta)/\beta$, integration by parts gives

```math
\langle\beta v,A_q(\beta v)\rangle
 =\int\beta^2|Dv|_q^2dA_q+\int r\beta^2|v|^2dA_q.
```

The first two inequalities of (5), min–max and $r\ge-\varepsilon$ give
the ordered lower bounds. Test $A_q$ on $\beta$ for the upper bound on
$\lambda_0$. Weighted variance and the round-sphere Poincaré constant
2 give the projected estimate, exactly as in note 21. Finally (4)
preserves the eigenvalues. ∎

For example, verified inputs $a\ge24(1-\eta)$, $b\le768(1+\eta)$ give
$\lambda_1\ge(1-\eta)/(16(1+\eta))-\varepsilon$. The inherited
critical surface has $\eta=\varepsilon=0$. **No nonzero-time values of
$\eta,\varepsilon$ are certified here.** This isolates the necessary
pointwise geometric/residual inequalities without a separate sweep of
infinitely many nonsymmetric modes. It is not a Hölder inverse estimate.

### 4.2 Thin-neck instability at positive times

**Corollary GS-N.** Under GS-S, suppose a smooth one-parameter family of
local solutions of the **true** equation (1) varies a neck. If its physical
normal variation $f$ on a candidate closed MOTS has consecutive zeros on
the two boundary circles of an annulus, with no zeros between them, and
that annulus is a proper subdomain of the connected surface, the candidate
is unstable.

*Proof.* Differentiating $\theta_+=0$ along the family gives $Lf=0$;
tangential variation contributes zero because the expansion is identically
zero. The function $e^{-\chi}f$ solves $A_q(e^{-\chi}f)=0$, is zero
on the boundary and has one sign inside. It is the principal Dirichlet
eigenfunction on the annulus, with eigenvalue zero. Strict domain
monotonicity for the self-adjoint Schrödinger operator implies
$\lambda_0(A_q,S)<0$. By (4), $\lambda_0(L,S)<0$. ∎

This extends the *logic* of note 12's neck exclusion to the evolving
reflection-symmetric problem. Its old numerical neck bounds do not carry
over automatically. Variations of the first-order level equation satisfy
that equation's linearization, not $L$ of the true MOTS; substituting them
here would be circular.

## 5 · The quadratic birth-time correction

Use the **true** local shooting function $F(z,d,t)$ from note 10, or its
smooth reduced expansion function. At a simple fold,
$F=F_z=0$, $F_t F_{zz}\ne0$. Let $\delta=d-d_c$ and solve
$F(z(d),d,t_b(d))=F_z(z(d),d,t_b(d))=0$.

**Proposition GS-Q.** All derivatives below are at the critical point:

```math
c=t_b'(d_c)=-F_d/F_t,\qquad
z'(d_c)=-(F_{zd}+cF_{zt})/F_{zz},
```

```math
t_b''(d_c)=-\frac{1}{F_t}\left[
F_{dd}+2cF_{dt}+c^2F_{tt}
-\frac{(F_{zd}+cF_{zt})^2}{F_{zz}}\right].                \tag{7}
```

*Proof.* Differentiate both fold equations once. Differentiate $F=0$
again; the term $F_z z''$ vanishes at the fold. Substitute the formula
for $z'$ into the resulting quadratic expression. ∎

The final negative square is the motion of the minimizing pole height.
Omitting it amounts to freezing the critical surface. Moreover
$\gamma(t)=\gamma_0-t^2\operatorname{Ric}(\gamma_0)+O(t^4)$ and
$K(t)=-t\operatorname{Ric}(\gamma_0)+O(t^3)$ imply that the first-order
slice model computes $F_d,F_t$ correctly at zero, but **need not compute
$F_{tt}$ correctly**. Time reversal makes the metric even; it does not
make the outgoing expansion, or this second derivative, vanish. The
minimal ingredients for (7) are second derivatives of the true shooting
map, using the $t^2$ metric term and its spatial derivatives. This is a
testable structural target, rather than a request for indefinitely more
Taylor coefficients.

For an explicit **uniform** remainder one can use an ordinary validated
implicit-function argument. Put $G=(F,F_z)$, $y=(z-z_c,t)$, and
$B=G_y(0,0)^{-1}$. On $|\delta|\le D$, $\|y\|\le r$, suppose

```math
\sup\|I-BG_y\|\le q<1,\qquad
\sup_{|\delta|\le D}\|BG(\delta,0)\|\le(1-q)r.            \tag{8}
```

The map $y\mapsto y-BG(\delta,y)$ is a contraction of the ball into
itself, so it has a unique fixed point for every parameter. On this tube
let $\kappa=\|B\|/(1-q)$, $M_1\ge\sup\|G_\delta\|$, and
$M_2\ge\sup\|D^2G\|$, using the product norm
$\|(\delta,y)\|=|\delta|+\|y\|$. Then implicit differentiation yields

```math
\|y'\|\le\kappa M_1,\qquad
\|y''\|\le\kappa M_2(1+\kappa M_1)^2,
\quad
|t_b(d)-c\delta|\le\tfrac12\kappa M_2(1+\kappa M_1)^2\delta^2. \tag{9}
```

Choose a norm with $|t|\le\|y\|$. To identify this curve as the local
birth threshold, retain uniform $F_{zz}>0$, $F_t<0$ and control the
boundary of the local shooting interval. $G\in C^2$ requires $F\in C^3$.
No numerical bounds for (8)–(9) have been certified; neither M' nor the
$L^2$ gap by itself supplies them.

## 6 · The exact remaining global lemma

A sufficient, noncircular route at fixed $d=0.7662$ is to certify all of:

1. A smooth embedded geodesic-slicing slab on a common compact region
   $D=B_R\setminus(B_+\cup B_-)$, plus exterior/end barriers which force
   relevant common outermost MOTSs into $D$. Initial strict barriers alone
   do not give an explicit time interval of persistence.
2. The trapped-region construction in these data, including that a common
   candidate produces a connected stable outermost component surrounding
   both holes, and the topology/symmetry needed for shooting. Uniqueness
   gives symmetry of the outermost component, not of arbitrary MOTSs.
3. A covering of **all true** regular shooting starts: every admissible
   closing trajectory either lies wholly in Theorem N's shell, violates
   embeddedness/common enclosure, or has a true neck variation satisfying
   GS-N. Straddling trajectories must be covered too.
4. The certified local foliation on that same time slab.

If these hold, an assumed common MOTS has an enclosing stable outermost
common component by item 2. Items 1 and 3 classify it; item 4 or GS-N
then contradicts its existence or stability. This implication is rigorous;
items 1–3 have **not** been verified for $0<t\le5\times10^{-6}$.
This is the remaining global-exclusion lemma, not a claimed global theorem.

## 7 · Consequences for the metric programme

The silent-flow branch already establishes that pointwise conjugation
preserves the troublesome principal symbol and does not prove convergence
of the sum of curl sectors. Its results remain on the pinned branch
`3c13047`; they must not be counted as a new local-time improvement here.
The next genuinely different continuation tool is a frame/connection
symmetric-hyperbolic formulation, with geometric energy estimates and
propagation of its defining and vacuum constraints [FR00, equations
(4.48) and (4.80)–(4.94) in the linked arXiv version].

Even a successful finite-slab energy certificate must control the chart
(or replace Gaussian slicing before caustics), the exterior and asymptotic
ends, and reconstruction into one spacetime. A generic inequality
$E'\le C(E)E+\mathrm{defect}$ has no certified consequence without the
actual Sobolev, coefficient, boundary-flux and gauge constants. We have
not supplied those constants. Complete future null infinity, event-horizon
identification and entry into a final-state stability regime remain open
for these finite-mass-ratio data. The present work removes two specific
proof obstacles; it supplies no global metric and no event-horizon merger.

## References

- [BHP21] Booth–Hennigar–Pook-Kolb, *Ultimate fate of apparent horizons
  during a binary black hole merger I*,
  [arXiv:2104.11343v1](https://arxiv.org/abs/2104.11343v1), §§II–IV.
  Full PDF inspected; see the [literature matrix](../research/literature_matrix.md).

- [AM09] Andersson–Metzger, *The area of horizons and the trapped region*,
  [arXiv:0708.4252](https://arxiv.org/abs/0708.4252), §§2–4, 7.
- [AMMS09] Andersson–Mars–Metzger–Simon, *The time evolution of marginally
  trapped surfaces*, [CQG 26, 085018](https://doi.org/10.1088/0264-9381/26/8/085018).
  Existing local/global MOTT theory; no claim that it supplies our missing
  explicit global exclusion constants.
- [PK20] Pook-Kolb et al., *Horizons in a binary black hole merger II*,
  [arXiv:2006.03940](https://arxiv.org/abs/2006.03940), §II.B and §III.
  Numerical spectral analysis of nonspinning mergers, not interval-certified
  bounds. Its null-normal conventions/deformation operators must not be
  identified silently with (3).
- [FR00] Friedrich–Rendall, *The Cauchy Problem for the Einstein Equations*,
  [arXiv:gr-qc/0002074](https://arxiv.org/abs/gr-qc/0002074), frame
  reductions in §§4.2–4.3 of the linked version (section numbering differs
  from the citation inherited by note 19).

## Discussion prompts

- Can true-data shooting with GS-N exclude the complement of the local
  shell more economically than constructing a complete positive foliation?
- Can (5) be enclosed along a validated evolved MOTS family, and then
  combined with quantitative elliptic regularity for (8)?
- What finite exterior slab and frame-energy certificate would first
  remove a domain obstruction, rather than merely improve the local clock?
