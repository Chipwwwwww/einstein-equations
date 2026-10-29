# 10 · When the common horizon is born: the proof

[09](09-merger-hypotheses.md) stated Proposition B. Take two black holes
momentarily at rest (Brill–Lindquist data), placed just too far apart for a
common apparent horizon. In the vacuum evolution, a common horizon appears
after a time proportional to the excess separation $d - d_c$. This note proves
it. The proof has two parts:

- a **computer-assisted** part (Theorem C, in 09), which locates the critical
  surface and fixes the signs of four derivatives;
- an **analytic** part, which turns those signs into a statement about the
  spacetime.

Each step says whether it is fully written out, standard and cited, or
checked by computer.

Conventions: [07](07-the-merger-problem.md) and
[`research/README.md`](../research/README.md). In particular
$\theta_\pm = \mathrm{tr}_S K \pm H_S$ with $K_{ij} = \nabla_i n_j$.

## Setting

- **Data.** $\Sigma = \mathbb R^3 \setminus \{p_+, p_-\}$ with $p_\pm = (0, 0, \pm d/2)$,
  $\gamma_d = \psi_d^4\delta$, $\psi_d = 1 + \frac{1}{4|x - p_+|} + \frac{1}{4|x - p_-|}$,
  and $K = 0$. These are time-symmetric vacuum data with ADM mass 1, smooth
  and real-analytic in $(x, d)$ away from the punctures.
- **Symmetry group.** $\mathcal G = O(2) \times \mathbb Z_2$: rotations
  about the $z$-axis, reflections through planes containing it, and
  $z \to -z$. It acts by isometries of every $\gamma_d$.
- **Critical surface.** $S_c$ and $\phi_0$ from Theorem C. $S_c$ is an
  embedded $\mathcal G$-invariant minimal sphere enclosing both punctures, and
  $\phi_0 > 0$ spans the kernel of its Jacobi operator
  $L = -\Delta_{S_c} - \mathrm{Ric}(\nu, \nu) - |A|^2$.
- **Slicing.** Fix a compact neighbourhood $V \subset \mathbb R^3$ of $S_c$
  that avoids the punctures for every $d$ near $d_c$. $\Sigma$ itself depends
  on $d$, but $V$ does not, so data on $V$ can be compared across $d$. The slices $\Sigma_t$ of the development have zero shift
  and lapse $N > 0$. We treat $N = 1$ (Gaussian normal coordinates) in
  detail; any smooth, $\mathcal G$-invariant lapse family works the same way.

## Step 1 · Smooth dependence (standard; cited)

The vacuum Einstein equations in harmonic gauge form a quasilinear
symmetric-hyperbolic system. Smooth data that depend smoothly on a parameter
give solutions that depend smoothly on it, on a common time interval, by the
standard theory: Choquet-Bruhat 1952; Hughes–Kato–Marsden 1977; Ringström,
*The Cauchy Problem in General Relativity* (2009) (*not re-checked this
session*). By finite speed of propagation, the development near $V$ for
$0 \le t \le t_0$ depends only on the data on a slightly larger compact set,
so the punctures play no role.

Gaussian normal coordinates $(t, x) \mapsto \exp_x(t\,n)$ are smooth in
$(t, x, d)$ for $|t| \le t_0$, $x \in V$, $|d - d_c| \le \delta_0$. In them the
induced data $(\gamma(t), K(t))$ on $\Sigma_t \cap V$ depend smoothly on
$(t, d)$. Moreover

```math
\gamma(t) = \gamma_d + O(t^2), \qquad K(t) = t\,\dot K + O(t^2), \qquad
\dot K_{ij} = D_iD_jN - N\,R_{ij}(\gamma_d),
```

by the evolution equations at a moment of time symmetry
([`tests/test_validation.py::test_adm_evolution_equation_sign`](../tests/test_validation.py)).
Everything is $\mathcal G$-invariant: the data are, Gaussian normal
coordinates are built from them, and uniqueness of the development carries
the symmetry along.

## Step 2 · The expansion map (standard)

For $u \in C^{2,\alpha}(S_c)$ small, let $S_u$ be the normal graph
$\{\exp^{\gamma(t)}_x(u(x)\nu(x))\}$ in $\Sigma_t$. Let $\Theta(u; d, t) \in C^{0,\alpha}(S_c)$
be its outward expansion $\theta_+$, pulled back to $S_c$. $\Theta$ is a
quasilinear second-order elliptic operator in $u$ whose coefficients depend
smoothly on $(x, d, t)$. So $\Theta : C^{2,\alpha} \times \mathbb R^2 \to C^{0,\alpha}$
is $C^\infty$ near $(0, d_c, 0)$, with $\Theta(0; d_c, 0) = 0$.

Since $K = 0$ at $t = 0$, $\Theta(u; d, 0)$ is the mean curvature of $S_u$, and

```math
D_u\Theta(0; d_c, 0) = L ,
```

the Jacobi operator. Equivariance: for $g \in \mathcal G$,
$\Theta(g \cdot u; d, t) = g \cdot \Theta(u; d, t)$, because $g$ is an
isometry of every slice preserving $S_c$ and its outward normal.

## Step 3 · Lyapunov–Schmidt reduction, with symmetry (written out)

$L$ is self-adjoint and elliptic on the closed surface $S_c$. Its kernel is
$\mathrm{span}\{\phi_0\}$: by Theorem C, $\phi_0 > 0$ lies in the kernel. A
positive eigenfunction belongs to the principal eigenvalue, which is simple.
Every other eigenvalue is therefore positive. Write

```math
C^{2,\alpha} = \mathrm{span}\{\phi_0\} \oplus W, \qquad
C^{0,\alpha} = \mathrm{span}\{\phi_0\} \oplus R, \qquad
W, R = \{\phi_0\}^{\perp_{L^2}} .
```

$L : W \to R$ is an isomorphism (Fredholm of index 0, no kernel in $W$).

Let $P$ be the $L^2$ projection onto $R$. By the implicit function theorem
there are neighbourhoods and a unique smooth $w(s, d, t) \in W$ with
$w(0, d_c, 0) = 0$ and $P\,\Theta(s\phi_0 + w; d, t) = 0$. Normal graphs
$u = s\phi_0 + w$ near $0$ are MOTSs exactly when the bifurcation function

```math
f(s, d, t) = \langle \phi_0, \Theta(s\phi_0 + w(s, d, t); d, t) \rangle
```

vanishes.

**Lemma (every MOTS near $S_c$ is $\mathcal G$-invariant).** For $g \in \mathcal G$:
$g \cdot \phi_0 = \phi_0$ (the principal eigenfunction is unique up to scale,
and $g$ maps it to a positive eigenfunction of the same eigenvalue);
$g$ preserves $W$ and $R$; and $\Theta$ is equivariant. So $g \cdot w(s, d, t)$
solves the same projected equation as $w(s, d, t)$, and by uniqueness the two
are equal. Hence every MOTS that is a small normal graph over $S_c$ is
$\mathcal G$-invariant: axisymmetric, and symmetric under $z \to -z$. ∎

## Step 4 · From MOTSs to the shooting function (written out)

Let $\tilde F(z_0, d, t)$ be the shooting function of the *true* slices $\Sigma_t$:

1. start the meridian curve at the axis point $(0, z_0)$ with horizontal
   tangent;
2. solve the ODE for $\theta_+ = 0$ in $(\gamma(t), K(t))$ until the curve
   meets $z = 0$;
3. record $\alpha + \pi/2$ there.

The regular solution from the axis exists and depends smoothly on
$(z_0, d, t)$. Near the axis it is the unique fixed point of the
averaging map $A \mapsto s^{-2}\int_0^s\sigma\Phi$ of notes/09, a contraction
for small $s$ whose data depend smoothly on the parameters. Away from the
axis the ODE is regular.

- A $\mathcal G$-invariant MOTS near $S_c$ meets the axis at two points.
  Its meridian is the regular solution through the upper one, $(0, z_0)$ with
  $z_0$ near $z_{0,c}$, and by the $z \to -z$ symmetry it crosses $z = 0$
  perpendicularly. So $\tilde F(z_0, d, t) = 0$.
- Conversely, take a zero of $\tilde F$ near $z_{0,c}$. Reflecting its meridian
  through $z = 0$ gives a $C^1$ closed curve, hence a smooth $\mathcal G$-invariant
  solution of $\theta_+ = 0$ by ODE uniqueness. By continuous dependence it is
  $C^2$-close to $S_c$, so it is a MOTS near $S_c$.

At $t = 0$, $\tilde F = F$, the function of Theorem C. Also
$\partial_t\tilde F = \partial_tF$ at $t = 0$, because the true ODE differs
from the one in Theorem C (data $(\gamma_d, t\dot K)$) by terms of order $t^2$.

## Step 5 · The fold (written out; signs from Theorem C)

At $(z_{0,c}, d_c, 0)$, Theorem C gives:

```math
\tilde F = \partial_{z_0}\tilde F = 0, \quad
\partial^2_{z_0}\tilde F \in [240.87, 242.70], \quad
\partial_d\tilde F \in [1.2315, 1.2364], \quad
\partial_t\tilde F \in [-0.51816, -0.51655] .
```

1. Since $\partial^2_{z_0}\tilde F > 0$, the implicit function theorem gives a
   smooth $\zeta(d, t)$ with $\partial_{z_0}\tilde F(\zeta, d, t) = 0$,
   $\zeta(d_c, 0) = z_{0,c}$. On a neighbourhood $I \ni z_{0,c}$, $\tilde F$ is
   strictly convex in $z_0$.
2. Let $m(d, t) = \tilde F(\zeta(d, t), d, t)$, the minimum over $I$. Then
   $m(d_c, 0) = 0$, $\partial_dm = \partial_d\tilde F > 0$ and
   $\partial_tm = \partial_t\tilde F < 0$ there (the $\partial_{z_0}$ term
   vanishes at $\zeta$).
3. On $I$, $\tilde F(\cdot, d, t)$ has two zeros if $m < 0$, one if $m = 0$,
   and none if $m > 0$. For $(d, t)$ near $(d_c, 0)$ the zeros stay in $I$,
   since $|z_0 - \zeta| = O(\sqrt{|m|})$.
4. Since $\partial_tm \ne 0$, $m(d, t) = 0$ defines a smooth curve
   $t = t_*(d)$ with $t_*(d_c) = 0$ and
   $t_*'(d_c) = \partial_d\tilde F / (-\partial_t\tilde F) \in [2.3768, 2.3936]$.
   Because $\partial_tm < 0$, $m(d, t) \le 0$ exactly when $t \ge t_*(d)$.

**Theorem (Proposition B, geodesic slicing).** There are $\delta_0, t_0 > 0$
and a neighbourhood $\mathcal U$ of $S_c$ such that, for every $d$ with
$d_c < d < d_c + \delta_0$ and every $0 \le t \le t_0$, the slice $\Sigma_t$
of the vacuum development contains a MOTS in $\mathcal U$ **if and only if**
$t \ge t_*(d)$. Here

```math
t_*(d) = c\,(d - d_c) + O\big((d - d_c)^2\big), \qquad c \in [2.3768, 2.3936],
\qquad d_c = 0.766197428240 \pm 1.1\times10^{-13}
```

in units $G = c = M_{\rm ADM} = 1$. Every such MOTS is axisymmetric,
reflection-symmetric, encloses both punctures, and is a small normal graph
over $S_c$.

*Proof.* Steps 3 and 4 identify the MOTSs in $\mathcal U$ with zeros of
$\tilde F$ near $z_{0,c}$; Step 5 counts them. ∎

## Step 6 · Before the common horizon (written out)

At $t = 0$ and $d > d_c$, $m(d, 0) > 0$, so $\Sigma_0$ has no MOTS in
$\mathcal U$. For every $d \in [0.76, 0.78]$, interval arithmetic certifies the
signs of $H_S$ on two concentric spheres about each puncture
([`tests/test_bl_fold.py`](../tests/test_bl_fold.py)). By Andersson–Metzger
(2009, Thm 3.1), each puncture is then enclosed by its own minimal surface,
and the two are disjoint. These persist for small $t$: the barrier
inequalities are strict, and $\theta_+$ of a fixed surface is continuous in
$t$. So along geodesic slicing the apparent horizon has two separate
components at $t = 0$ and a component enclosing both by $t = t_*(d)$.

## What is and is not proved

| Ingredient | Status |
| --- | --- |
| Theorem C (fold location, signs, positive Jacobi field, embedding) | computer-assisted (Arb), `research/calc/fold_proof.py` |
| Step 1 (smooth dependence) | standard; cited, not re-checked |
| Step 2 (smoothness of $\Theta$, $D_u\Theta = L$) | standard |
| Steps 3–6 | written out above |
| other lapses | same proof. For $N = \psi^{-2}$: $\partial_t F \in [-0.11612, -0.11574]$, so $c \in [10.606, 10.682]$ (`certify_slopes()`) |
| no common MOTS **anywhere** on $\Sigma_0$ (not just in $\mathcal U$) | **not proved**. At $d = 0.78$: barrier lemmas certified, and reduced to a shooting window; the remaining obstruction is a family of necked surfaces hugging the individual horizons (`research/calc/bl_global.py`, research/status.md gap 3) |
| anything about event horizons | **not addressed** |

So the vacuum development of near-critical Brill–Lindquist data contains an
apparent-horizon merger along geodesic slicing. The common MOTS forms at a
time linear in $d - d_c$, with a certified coefficient. As far as our bounded
literature search found (notes/09), this is **potentially new**.

## Discussion prompts

- The proof needs $\phi_0 > 0$ only to know the kernel is one-dimensional and
  symmetric. What could happen at a fold where the kernel were bigger, for
  example for unequal masses with an accidental degeneracy?
- Step 6 shows two separate apparent horizons before $t_*$. Is it plausible
  that the *event* horizon is already connected at $t = 0$ for these data?
- The coefficient $c$ depends on the lapse, but the weighted proper time
  $\int Nt_*\phi_0|A|^2/\int\phi_0|A|^2$ does not. Is there a slicing-independent
  way to phrase "when the common horizon is born"?
