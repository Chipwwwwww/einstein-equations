# 20 · An all-mode stability estimate for the critical common horizon

*Research branch `research/horizon-stability-2026-11-01`, based on `3c13047`.
The target is still a genuine, asymptotically flat, finite-mass-ratio vacuum
merger. This note concerns one initial minimal surface, not that spacetime.*

## 1 · The mathematical result and its status

The missing spectral estimate in [note 15 §7.3](15-route-j-convergence.md)
does not require solving infinitely many angular eigenvalue problems.
The positive Jacobi field already supplied by the fold converts it into a
weighted Dirichlet-form comparison with the **unit round sphere**.

The analytic comparison theorem below is proved directly. The explicit
Brill–Lindquist constants are verified by
`python -m research.calc.horizon_gap`; the output must report all assertions
true before the numerical corollary is used. The certificate files and current
run status are recorded in [the consolidated report](../research/status.md).

For the equal-mass critical common minimal sphere $S_c$, the resulting bounds,
in units of the total ADM mass $M$, are

```math
\lambda_{\rm gap}(L_{S_c})\ge\frac{1}{16M^2},\qquad
\left\|(L_{S_c}|_{\phi_0^\perp})^{-1}\right\|_{L^2\to L^2}\le16M^2.
```

The orthogonal complement is taken in physical surface area measure. The
estimate includes nonsymmetric perturbations, not just the axisymmetric
shooting sector. In the axisymmetric, equator-even sector it improves to
$3/(16M^2)$; for an azimuthal mode $m\ne0$ it is at least $m^2/(2M^2)$.

**Scope:** this removes the all-mode $L^2$ spectral-gap subproblem in note 15.
It does not supply the nonlinear Hölder/Sobolev inverse bounds, an explicit
formation time, a longer vacuum development, or an event horizon. The basic
ground-state transformation and spectral comparison are standard tools;
the contribution here is their explicit, validated application to this
repository's critical surface. Priority or publication novelty is not claimed.

## 2 · A ground-state comparison theorem

Let $(S,q)$ be a smooth closed connected surface, and

```math
L=-\Delta_q+V,\qquad L\phi=0,\qquad\phi>0.
```

Here $V$ is a smooth real potential. For the minimal-surface problem,
$V=-\operatorname{Ric}_{\gamma}(\nu,\nu)-|A|^2$; this is the repository's
Jacobi sign convention. Let $F:\mathbb S^2\to S$ be a diffeomorphism and
$\sigma$ the unit round metric. Suppose, after pullback, that for every
covector $\xi$,

```math
\phi^2\frac{dA_q}{d\Omega}\,q^{-1}(\xi,\xi)
   \ge a\,\sigma^{-1}(\xi,\xi),\qquad
\phi^2\frac{dA_q}{d\Omega}\le b,
\qquad a>0,\ b>0.                                      \tag{1}
```

**Theorem 1 (comparison).** Number eigenvalues with multiplicity starting at
zero: $\lambda_0(L)=0$. If $\mu_k$ are the eigenvalues of $-\Delta_\sigma$,
then

```math
\lambda_k(L)\ge\frac{a}{b}\mu_k\quad(k\ge0).
                                                               \tag{2}
```

In particular $\ker L=\operatorname{span}\{\phi\}$ and
$\lambda_{\rm gap}\ge2a/b$. If the metric, potential, and $\phi$ are
axisymmetric and equator-even, the first positive eigenvalue restricted to
axisymmetric equator-even functions is at least $6a/b$.

*Proof.* Put $u=\phi v$. Integration by parts and $L\phi=0$ give the exact
ground-state identity

```math
\int_S uLu\,dA_q
 =\int_S\bigl(|\nabla u|_q^2+Vu^2\bigr)dA_q
 =\int_S\phi^2|\nabla v|_q^2\,dA_q.                    \tag{3}
```

For complex functions replace squares by modulus squares. Positivity of
$\phi$ on the compact surface makes multiplication by $\phi$ an isomorphism
of both $H^1$ and $L^2$. The only zero-energy $v$ is constant, proving the
kernel assertion. By (1), for every nonzero $v$,

```math
\frac{\int\phi^2|\nabla v|_q^2dA_q}{\int\phi^2|v|^2dA_q}
 \ge\frac{a}{b}\,
       \frac{\int|\nabla v|_\sigma^2d\Omega}{\int|v|^2d\Omega}.
                                                               \tag{4}
```

Apply the min–max principle to all $(k+1)$-dimensional spaces of $v$.
Their images under multiplication by $\phi$ are exactly all such spaces
of $u$. This gives (2). On the unit sphere the eigenvalues are
$\ell(\ell+1)$, with multiplicity $2\ell+1$. In the axisymmetric
equator-even subspace only even $\ell$ occur, and the first nonconstant
one is $\ell=2$. Multiplication by $\phi$ preserves that subspace. ∎

**Orthogonality check.** One must not apply the ordinary zero-mean sphere
inequality directly to $v=u/\phi$: $u\perp\phi$ implies
$\int\phi^2v\,dA_q=0$, not $\int v\,d\Omega=0$. The min–max proof avoids
this error. Alternatively, use weighted variance:

```math
\int\phi^2v^2dA_q
 =\inf_c\int\phi^2(v-c)^2dA_q
 \le b\inf_c\int(v-c)^2d\Omega
 \le\frac b2\int|\nabla v|_\sigma^2d\Omega.
```

The first equality uses precisely the weighted mean-zero condition.

**Corollary 1 (residual correction).** If $f\perp\phi$ and an approximate
solution $w\perp\phi$ has residual $r=Lw-f$, then the unique exact solution
$u\perp\phi$ obeys

```math
\|u-w\|_{L^2(q)}\le\frac{b}{2a}\|r\|_{L^2(q)}.          \tag{5}
```

This is an $L^2$ assertion. It is not a pointwise or $C^{2,\alpha}$ estimate.

**Residual version (a posteriori comparison).** Exact knowledge of a zero
mode is not essential for the comparison tool. For any smooth $\varphi>0$
on a closed surface, set $r=(L\varphi)/\varphi$. If the weights formed with
$\varphi$ obey (1), then the same integration by parts gives

```math
Q[\varphi v]=\int\varphi^2|\nabla v|^2dA
                  +\int r\,|\varphi v|^2dA.
```

Consequently

```math
\lambda_k(L)\ge(a/b)\mu_k+\inf r,\qquad
\inf r\le\lambda_0(L)\le\sup r,\qquad
\lambda_1(L)-\lambda_0(L)\ge2a/b-\operatorname{osc}r.    \tag{5a}
```

The lower bounds follow by min–max; testing with $u=\varphi$ gives the
principal-eigenvalue upper bound. On $\varphi^\perp$, weighted variance
also proves $Q[u]\ge(2a/b+\inf r)\|u\|_2^2$. Thus the projected operator
$PLP$ there is coercive whenever $2a/b+\inf r>0$, even if $\varphi$ is
only an approximate eigenfunction. This is the standard ground-state/Barta
method with the comparison constants exposed, not a newly invented
spectral principle.

For a future application, **if** the same bounds $a=24,b=768$ are
reverified and $|r|\le\epsilon<1/32$ in $M=1$ units, at most one
eigenvalue lies below $1/16-\epsilon$, and the principal spectral gap is
at least $1/16-2\epsilon$. No such perturbed-surface residual has been
certified in this session. The statement assumes a self-adjoint operator;
it must not be applied automatically to a general evolving MOTS operator
with a first-derivative drift term.

## 3 · Geometric weights without singular polar divisions

Consider an axisymmetric star-shaped surface in a conformally flat spatial
metric $\gamma=\psi^4\delta$. Write its radial graph as $r=h(\theta)>0$,
and set

```math
c=\frac{h}{\sqrt{h^2+h_\theta^2}}\in(0,1].
```

Then the induced metric and its density are

```math
q=\psi^4h^2\bigl(c^{-2}d\theta^2+\sin^2\theta\,d\varphi^2\bigr),
\qquad\frac{dA_q}{d\Omega}=\frac{\psi^4h^2}{c}.          \tag{6}
```

The eigenvalues of $\phi^2(dA_q/d\Omega)q^{-1}$ relative to $\sigma^{-1}$
are $\phi^2c$ and $\phi^2/c$. Thus it suffices to certify the two scalar
inequalities

```math
A(\theta):=\phi^2c\ge a,\qquad
B(\theta):=\frac{\phi^2\psi^4h^2}{c}\le b.              \tag{7}
```

These can be evaluated directly on the existing meridian ODE. With flat
arclength $s$, tangent $(\rho_s,z_s)=(\cos\alpha,\sin\alpha)$, outward
flat unit normal $n=(-\sin\alpha,\cos\alpha)$, and pole-height parameter
$z_0$,

```math
h^2=\rho^2+z^2,\qquad
c=\frac{-\rho\sin\alpha+z\cos\alpha}{\sqrt{\rho^2+z^2}},
\qquad
\phi=\psi^2\left(-\sin\alpha\,\partial_{z_0}\rho
                    +\cos\alpha\,\partial_{z_0}z\right).        \tag{8}
```

The factor $\psi^2$ is essential: the Jacobi function is **physical normal
displacement**, whereas the shooting variation first gives Euclidean normal
displacement. Tangential changes introduced by switching from arclength to
the height parameter cancel in $n\cdot\partial_{z_0}X$.

The formula for $c$ contains no division by $\rho$ or $\sin\theta$ and
extends over the poles. Moreover,
$d\theta/ds=c/h>0$: positivity of this quantity, together with the
embedded reflected meridian, certifies the global radial-graph chart.
This is checked rather than assumed for $S_c$.

**Corollary 2 (azimuthal estimate).** If additionally
$q_{\varphi\varphi}=\psi^4\rho^2\le R_\varphi^2$ everywhere, then in every
nonzero Fourier sector $m$,

```math
\lambda\bigl(L|_m\bigr)\ge m^2/R_\varphi^2.             \tag{9}
```

Indeed the angular part of (3) already bounds the full denominator times
$m^2/R_\varphi^2$; the meridian part is nonnegative. Smoothness at the poles
is included by density in the form domain. No angular-mode cutoff is used.

## 4 · Application to critical Brill–Lindquist data

Use the same data and conventions as the earlier fold certificate:

```math
\gamma=\psi^4\delta,\quad K=0,\quad
\psi=1+\frac1{4r_+}+\frac1{4r_-},\quad
r_\pm=\sqrt{\rho^2+(z\mp d/2)^2},\quad M_{\rm ADM}=1.
```

The exact critical $(z_{0,c},d_c)$ is the unique root of $(F,F_{z_0})=0$
inside the original dyadic floating-point box of radius $10^{-12}$ around
$(0.6326595528058354,0.766197428240025)$. Those decimal centres label the
binary64 numbers specified by hexadecimal strings in the certificate;
they are not claimed to be the exact root.

At that root the meridian reflects smoothly through $z=0$. Differentiating
the minimal-surface family at fixed $d_c$ gives $L\phi=0$; the fold
condition $F_{z_0}=0$ supplies the reflected Jacobi boundary condition.
This is the existing Theorem C mechanism, now replayed along with the new
weights. Positive $\phi$ and $c$ are verified on every enclosure.

**Theorem 2 (computer-assisted numerical corollary).** The validated
geometric inequalities are

```math
A\ge24,\qquad B\le768,\qquad\psi^4\rho^2\le2.          \tag{10}
```

Therefore, with $\lambda_0=0$ and eigenvalues counted with multiplicity,

```math
\lambda_k(L_{S_c})\ge\frac{\ell(\ell+1)}{32M^2},
\qquad \ell^2\le k<(\ell+1)^2.                          \tag{11}
```

In particular the gap is $\ge1/(16M^2)$. In note 15's one-based notation
this is a bound on $\lambda_2$, not on the principal zero eigenvalue.
Equations (5) and (9) give the inverse and azimuthal estimates in §1.
The axisymmetric equator-even inverse norm is at most $16M^2/3$.

**Verification design.**

1. Recompute the centre shooting jet and its enclosure on the entire
   original parameter box. Recompute the Krawczyk inclusion and a strict
   contraction bound using an exactly represented fixed preconditioner.
2. Use the existing analytic pole collar, validated Taylor integrator and
   three regular coordinate segments. Evaluate (8) and (7) on the full
   **step enclosures**, not just the step endpoints.
3. Check the whole closed pole collar separately, including $s=0$.
   Reflect the results to the southern hemisphere by the exact symmetry.
4. Compare Arb balls directly with the exact integers in (10). No
   rounded decimal or floating-point eigenvalue enters those comparisons.
5. Store the fold box, Krawczyk image, all geometric enclosures, phase
   counts, arithmetic precision, library version and source hashes.

The trust base is the proof here, the inherited axis/ODE enclosure argument
in `fold_proof.py`, its implementation, and python-flint/Arb. This is a
computer-assisted result, not a proof checked by a formal proof assistant.
The inherited CK evolution and global no-common-horizon sweeps are not
inputs to Theorem 2 and are not replayed by this checker.

The diagnostic interval for $c$ can have an upper endpoint slightly above
1, and the interval for $\psi^4\rho^2$ can have a slightly negative lower
endpoint. These are enclosure slack at the extrema/poles. The exact
geometric identities give $c\le1$ and $\psi^4\rho^2\ge0$ independently;
the proof uses the certified positive lower bound on $c$ and upper bound
on $\psi^4\rho^2$, not those slack endpoints.

## 5 · Implications for the merger programme

The critical surface has one soft deformation, its positive zero mode.
Every orthogonal deformation now has a certified restoring coefficient;
nonsymmetric modes do not conceal a second zero eigenvalue. Qualitative
simplicity was already known from positivity. The new information is an
explicit bound and inverse estimate usable in an error budget.

This is about an elliptic surface operator at $t=0$. It is **not** a
stability estimate for time evolution of the spacetime or for gravitational
waves. The root separation $d_c$ is a threshold in the initial-data family,
not an event-horizon merger time.

The next bottleneck is now sharply located: control the nonlinear expansion
map on nearby graphs in the **true** evolved data. Its derivative is $L$ at
the critical point, but the graph remainder and the norms required by a
quantitative implicit-function argument remain uncertified. A possible
economy is to certify a strictly outer-trapped trial surface and an outer
untrapped barrier, instead of locating the exact fold time. Both require
explicit pointwise errors; (5) alone is insufficient.

Separately, a complete merger still needs controlled evolution over a
meaningful exterior interval and a global causal analysis connecting the
early and late horizon components to $\mathscr I^+$. No argument in this
note replaces those steps.

## 6 · Literature and novelty boundary

- Andersson–Mars–Metzger–Simon, [arXiv:0811.4721v2](https://arxiv.org/html/0811.4721v2),
  especially §5, supplies the established zero-principal-eigenvalue MOTS
  framework. Its qualitative horizon-birth mechanism is not new here.
- Pook-Kolb et al., [arXiv:2006.03940](https://arxiv.org/html/2006.03940v1),
  §§II.2–III, study horizon stability spectra numerically during head-on
  mergers. Our result is an interval-certified lower bound at the specified
  initial-data fold, not a new discovery that stability spectra matter.
- Ground-state transforms, min–max comparison, and bounded-density
  Poincaré comparison are standard spectral methods. The proof is included
  to expose every constant and the correct weighted orthogonality.
- Hintz, [arXiv:2408.06715](https://arxiv.org/html/2408.06715v1), §6 and
  Remark 6.2, constructs exact vacuum merger spacetimes in an extreme-mass-ratio
  existence sense. The paper explicitly leaves further event-horizon analysis
  open. This is distinct from an explicit finite-ratio merger metric.

A bounded literature search did not establish priority for the particular
explicit bound (10). It must not be advertised as previously unknown to all
researchers. It is a concrete addition to the present repository.

## Discussion prompts

- Can sharper weighted comparisons give useful nonlinear inverse constants
  without explicitly computing the low eigenfunctions?
- Can an outer-trapped witness be certified inside the existing spacetime
  tube with less work than certifying the exact local bifurcation time?
- Which exterior continuation estimates would turn that local horizon
  result into a statement about a complete merger?
