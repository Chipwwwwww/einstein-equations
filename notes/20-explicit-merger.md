# 20 · An explicit apparent-horizon merger: barrier certificates in the certified development

Notes 13 and 15 left one gap between the two halves of the near-critical Brill–Lindquist (BL)
programme:

- **Theorem H** (notes/13) proves that a common apparent horizon forms in the geodesic slicing when
  $d$ is slightly above $d_c$, but only for $d \in (d_c, d_c + \delta_0)$ with $\delta_0$ not computed;
- **Proposition T** and its extension (notes/15, 17) give the spacetime metric on a tube $U$ about the
  critical surface $S_c$ for $\lvert t\rvert < 8.683\times10^{-4}$ as an explicit convergent series with
  certified error bounds.

Before this note, no separation was known for which the merger provably happens **inside** the
certified region (status: "no certified merger at any explicit $d$"). The route proposed in
notes/15 §7.3 needed a certified spectral gap of the Jacobi operator over all Fourier modes and a
validated MOTS solver with $C^{2,\alpha}$ bounds.

This note closes the gap by a different and much cheaper argument.

- **Existence needs only $C^0$ information on two surfaces.** By Andersson–Metzger, a MOTS exists
  between an outer-trapped and an outer-untrapped surface. We construct both as level sets of the
  *first-order* expansion, on which the true expansion differs from a known constant only by a
  remainder controlled pointwise by the CK certificate.
- **Non-existence in a shell needs only $C^0$ information on a foliation** (Lemma 3). By the
  comparison principle, no MOTS lies in a region foliated by surfaces with $\theta_+ > 0$. The
  leaves may have convex ridges, so no closing condition is needed. Its certificate is not
  completed (§4).

Neither argument uses the spectrum of the Jacobi operator.

*Status (session 16).*

| Item | Category |
| --- | --- |
| expansion of a surface of revolution in an axisymmetric metric (Lemma 1) | derived; checked against an independent SymPy level-set computation |
| remainder bound (Lemma 2) | proved (mean-value form) |
| barrier theorem | literature: Andersson–Metzger 2009, Thms 3.1 and 4.1, **read on the page** this session |
| foliation exclusion with convex ridges (Lemma 3) | proved here (elementary comparison argument) |
| **Theorem M**: common MOTS on $\Sigma_t$, $d = 0.7662$, $t = 1.2\times10^{-5}$ and $6.5\times10^{-6}$ | **computer-assisted** (new) |
| Theorem N (no MOTS in the shell for $t \le t_0$) | **not certified**: all hypotheses of Lemma 3 hold on test pieces except the ridge condition $F > 0$ over $z_0$-pieces (§4) |
| **Corollary**: explicit apparent-horizon merger at $d = 0.7662$ (two components at $t = 0$, a common MOTS by $t = 6.5\times10^{-6}$) | computer-assisted, combining F′ and M |
| first common apparent horizon of the slicing; event horizons | **open** (§6) |
| category-A metric | **no** |

Conventions are those of notes 00, 15 and 16: signature $(-,+,+,+)$, $G = c = M_{\rm ADM} = 1$,
geodesic slicing $g = -dt^2 + \gamma_{ij}(t, x)\,dx^idx^j$ in comoving Cartesian coordinates,
$K_{ij} = \tfrac12\partial_t\gamma_{ij}$ (Wald's sign), $\theta_\pm = \mathrm{tr}_SK \pm H$ with $H$ the
mean curvature for the outward normal. These are Andersson–Metzger's $\theta_\pm = P \pm H$.

## 1 · Baseline

- **Base.** `main` at `90f023b`, which contains PRs #1–#6 (sessions 1–15).
- **Rerun at the start:** the fast tests of the modules used here (`tests/test_bl_metric.py`,
  `tests/test_ck_bound.py`, `tests/test_fold_proof.py`, `tests/test_bl_fold.py`, fast parts: 43
  passed).
- **Inherited and not rerun:** the CK certificates (`research/calc/data/ck_bound*.json`, recomputed in
  session 13), Theorems C and F′.
- **Regression of the reused integrator.** `bl_barrier.shoot` (a copy of `fold_proof.shoot` with
  the source as a parameter) reproduces `fold_proof.shoot` at the fold point **bit for bit**, all six
  jet components (`tests/test_bl_barrier.py::test_integrator_copy_reproduces_fold_proof`).

## 2 · Three lemmas

### 2.1 The expansion of a surface of revolution

Let $\gamma$ be a Riemannian metric invariant under rotations about the $z$-axis and under
$y \to -y$, so that $\gamma_{xy} = \gamma_{yz} = 0$ on the half-plane $y = 0$, $x = \rho > 0$. Let
$S$ be the surface of revolution of a profile curve $(\rho(s), z(s))$ parametrised by flat
arclength, with tangent angle $\alpha$: $T = (\cos\alpha, 0, \sin\alpha)$, flat outward normal
$n = (-\sin\alpha, 0, \cos\alpha)$.

> **Lemma 1.** At the meridian point,
>
> ```math
> H = -\frac{1}{\lvert n\rvert_\gamma}\left[\frac{\alpha' + n\cdot\Gamma(T,T)}{\gamma_{TT}} + \frac{\sin\alpha/\rho + n\cdot\Gamma_{yy}}{\gamma_{yy}}\right],
> \qquad
> \mathrm{tr}_SK = \frac{K_{TT}}{\gamma_{TT}} + \frac{K_{yy}}{\gamma_{yy}},
> ```
>
> where $n\cdot\Gamma(V, W) = n^l\Gamma_{l,jk}V^jW^k$, $n^l = \gamma^{lk}n_k$,
> $\lvert n\rvert_\gamma^2 = n_ln^l$ and $\Gamma_{l,jk} = \tfrac12(\partial_j\gamma_{lk} + \partial_k\gamma_{lj} - \partial_l\gamma_{jk})$.
> The only components that enter are $\gamma_{xx}, \gamma_{yy}, \gamma_{zz}, \gamma_{xz}$, their $x$-
> and $z$-derivatives, $\partial_y\gamma_{xy}$, $\partial_y\gamma_{yz}$, and $K_{xx}, K_{yy}, K_{zz}, K_{xz}$.

*Proof.* Parametrise $X(s, \varphi) = (\rho\cos\varphi, \rho\sin\varphi, z)$ and evaluate at
$\varphi = 0$. The induced metric is diagonal, $h = \mathrm{diag}(\gamma_{TT}, \rho^2\gamma_{yy})$,
because $\gamma_{xy} = \gamma_{yz} = 0$. The second fundamental form is
$A_{ab} = -\nu_i(\partial_a\partial_bX^i + \Gamma^i_{jk}\partial_aX^j\partial_bX^k)$ with
$\nu_i = n_i/\lvert n\rvert_\gamma$. Use $\partial_s^2X = \alpha' n$ and $\partial_\varphi^2X = (-\rho, 0, 0)$,
and take the trace. ∎

*Checks.*

- For $\gamma = \psi^4\delta$ it reduces to $H = \psi^{-2}(-\alpha' - \sin\alpha/\rho + 4\,n\cdot\nabla\psi/\psi)$,
  the formula behind Theorems C and F.
- For a non-conformally-flat axisymmetric metric with $\gamma_{xz} \ne 0$ and a nonzero $K$, it
  agrees with $\mathrm{tr}_SK + \mathrm{div}_\gamma(\nabla f/\lvert\nabla f\rvert)$ computed by SymPy
  from a level-set function, to $10^{-15}$ (`tests/test_bl_barrier.py`).

### 2.2 Level surfaces of the first-order expansion, and the remainder

To first order in $t$ the slice data are $\gamma_0 = \psi^4\delta$ and $K = t\kappa_0$, with
$\kappa_0 = -\mathrm{Ric}(\gamma_0) = g_1$ (notes/15 §7.1). The **first-order expansion** is

```math
\theta_+^{(1)} = H[\gamma_0] + t\,\mathrm{tr}_S^{\gamma_0}\kappa_0 .
```

For a constant $h$, the profile of an axisymmetric surface with $\theta_+^{(1)} = h$ solves the
shooting ODE of Theorems C and F with one extra term:

```math
\alpha' = -\frac{\sin\alpha}{\rho} + 4\,\frac{n\cdot\nabla\psi}{\psi} + \psi^2\big(t\,p - h\big),\qquad p = \mathrm{tr}^{\gamma_0}_S\kappa_0 .
```

Start it on the axis at height $z_0$ with horizontal tangent, and let $F(z_0) = \alpha_{\rm end} + \pi/2$
where the curve meets $z = 0$.

- If $F$ changes sign on $[z_0 - \delta, z_0 + \delta]$, some $z_0^*$ in between closes the curve at a
  right angle. Its rotation and reflection is a smooth closed surface $S_h$ (the reflected curve
  solves the same ODE with the same Cauchy data on the plane). On $S_h$,
  $\theta_+^{(1)} = h$ **exactly**.
- A validated run with $z_0$ ranging over the whole ball encloses every such curve.

Write $p$ for the pointwise data $(\gamma_{ij}, \partial_l\gamma_{ij}, K_{ij})$ that enter Lemma 1, and
$\Theta(p)$ for the expansion as a function of them at a fixed point of the curve, with fixed
$\alpha$, $\alpha'$, $\sin\alpha/\rho$.

> **Lemma 2.** Let $p_0 = (g_0, \partial g_0, t g_1)$ and $p = p_0 + \delta p$ be the true data of
> $\Sigma_t$. Then on $S_h$ the true expansion is $\theta_+ = h + R$, with
>
> ```math
> R \in \Big\{\tfrac{d}{d\sigma}\Theta(p_0 + \sigma\,\delta p) : \sigma \in [0, 1]\Big\},
> ```
>
> and componentwise
>
> ```math
> \delta\gamma \in t^2g_1 + [-V, V],\qquad \delta\partial\gamma \in t^2\partial g_1 + [-D, D],\qquad
> \delta K \in 2t^3g_2 + [-E, E],
> ```
>
> where $V$ is the CK tail of $\gamma$ beyond $t^2$ (order $N = 1$), $D$ the tail of
> $\partial_x\gamma$ beyond $t^2$ ($N = 2$), and $E$ half the tail of $\partial_t\gamma$ beyond $t^3$ ($N = 2$).

*Proof.* $\Theta(p_0) = h$ on $S_h$ by construction, and the mean-value theorem applies to
$f(\sigma) = \Theta(p_0 + \sigma\delta p)$. The series $\gamma = \sum_k g_kt^{2k}$ and
$K = \tfrac12\partial_t\gamma$ converge on $U$, and the tails are bounded componentwise by the
Frobenius bounds of notes/16 §3 (`bl_metric.tail_bounds`). ∎

In the computation, $g_1$, $\partial g_1$ and $g_2$ are enclosed on each box by the exact jet
recursion (`bl_metric.series_jets`). $V$, $D$, $E$ are the minimum over the certificates whose domain
contains the box. One `arb_series` evaluation in $\sigma$, with $\sigma = [0, 1]$ and $\delta p$ the
balls above, encloses $f'$ on the whole box. The cancellation that makes $\theta_+$ small happens
only in $\Theta(p_0) = h$, which is exact. The enclosure of $R$ needs only crude bounds, because it
multiplies a derivative of size $O(\psi^{-6}) \approx 10^{-2}$ by $\delta p = O(t^2)$.

### 2.3 The two geometric principles

> **Barrier theorem** (Andersson–Metzger 2009, Theorem 3.1 with Theorem 4.1; read on the page).
> Let $(M, g, K)$ be a smooth compact initial data set with $\partial M = \partial^-M \sqcup \partial^+M$,
> both non-empty smooth closed surfaces. Suppose $\theta_+[\partial^-M] < 0$ for the normal pointing
> into $M$, and $\theta_+[\partial^+M] > 0$ for the normal pointing out of $M$. Then $M$ contains a
> smooth, embedded, stable MOTS homologous to $\partial^+M$.

No energy condition and no constraint equation is assumed.

> **Lemma 3 (foliation exclusion).** Let $\{S_\lambda\}_{\lambda\in[a,b]}$ be a $C^0$-continuous
> family of closed embedded surfaces bounding regions $I_\lambda$ with
> $\overline{I_\lambda} \subset I_\mu$ for $\lambda < \mu$. Each $S_\lambda$ is smooth except possibly
> along finitely many curves where two smooth sheets meet at a **convex ridge**: $I_\lambda$ is
> locally the intersection of the two sides, and the sheets' normals differ. Suppose
> $\theta_+(S_\lambda) > 0$ on the smooth part, for the normal pointing out of $I_\lambda$. Then the
> shell $W = \overline{I_b}\setminus I_a$ contains no closed embedded MOTS (expansion taken for the
> normal pointing out of the region the MOTS bounds).

*Proof.* Let $\Sigma \subset W$ be a MOTS, and let
$\lambda^* = \min\{\lambda : \Sigma \subset \overline{I_\lambda}\}$, which exists by continuity.

- **Case $\lambda^* = a$.** Then $\Sigma \subset \overline{I_a} \cap W = S_a$, so $\Sigma = S_a$. This is
  impossible, since $\theta_+(S_a) > 0$, or $S_a$ is not smooth.
- **Case $\lambda^* > a$.** If $\Sigma \cap S_{\lambda^*} = \emptyset$, then $\Sigma$ is a compact
  subset of $I_{\lambda^*}$, so it lies in $I_\lambda$ for some $\lambda < \lambda^*$. That
  contradicts the minimality of $\lambda^*$. So $\Sigma$ touches $S_{\lambda^*}$ from inside at some
  point $q$.
  - *$q$ on a ridge.* Then $\Sigma$ lies in both smooth sides, so its tangent plane at $q$ equals
    both sheets' tangent planes. That is impossible, since the normals differ.
  - *$q$ on the smooth part.* The outward normals agree at $q$. Write both surfaces as graphs over
    the common tangent plane, with $\Sigma$ below. The comparison principle for the mean-curvature
    operator gives $H_\Sigma(q) \ge H_{S_{\lambda^*}}(q)$. The term $\mathrm{tr}_SK$ depends only on
    the point and the tangent plane, so $0 = \theta_+(\Sigma)(q) \ge \theta_+(S_{\lambda^*})(q) > 0$,
    a contradiction. ∎

A leaf made from the level-$h$ shooting curve with $F(z_0) > 0$ meets the plane at the angle $F$ from
the vertical, leaning outward. Reflected, it has an outward equatorial ridge with interior angle
$\pi - 2F < \pi$, which is convex. So **the leaves need not close smoothly**: the closing condition
$F = 0$ is replaced by the open condition $F > 0$, to be certified over intervals of $z_0$. §4
explains why that last step is not yet done.

## 3 · Theorem M: a common MOTS at an explicit time

> **Theorem M (computer-assisted).** Let $(M, g)$ be the vacuum development of equal-mass BL data
> with separation $d = 0.7662$, so $d - d_c = 2.5718\times10^{-6}$ (units $M_{\rm ADM} = 1$, bare masses
> $\tfrac12$). Let $\Sigma_t$ be its geodesic slicing, and $t_1 = 1.2\times10^{-5}$.
>
> 1. $\Sigma_{t_1}$ contains two smooth embedded spheres $S_{\rm in}$ and $S_{\rm out}$, with $S_{\rm in}$
>    inside $S_{\rm out}$. Both are invariant under rotations about the axis and under $z \to -z$, are
>    radial graphs about the origin, and enclose both punctures. The closed shell $\Omega$ between them
>    lies in the certified domain $U$ of the CK certificate.
> 2. Everywhere on them,
>
>    ```math
>    \theta_+[S_{\rm in}] \in [-7.580, -7.420]\times10^{-7}, \qquad \theta_+[S_{\rm out}] \in [1.992, 2.008]\times10^{-6}.
>    ```
>
> 3. Hence $\Sigma_{t_1}$ contains a smooth, embedded, stable MOTS $\Sigma^* \subset \Omega$. One of its
>    components separates $S_{\rm in}$ from $S_{\rm out}$, so it encloses both punctures: a **common
>    MOTS**.
>
> The same holds at $t_1' = 6.5\times10^{-6}$, with
> $\theta_+[S_{\rm in}] \in [-4.162, -3.838]\times10^{-8}$ and $\theta_+[S_{\rm out}] \in [1.998, 2.002]\times10^{-6}$.
> That is $1.06$ times the linear-law birth time $t_*^{\rm lin} = c\,(d - d_c) \in [6.113, 6.156]\times10^{-6}$.

**Proof.**

| Step | How | Certified numbers ($t_1 = 1.2\times10^{-5}$ / $t_1' = 6.5\times10^{-6}$) |
| --- | --- | --- |
| the development on $U \times (-T, T)$ is the convergent series, with tail bounds | Proposition T and the CK′ certificates (notes/15, 17), inherited | $T = 8.683\times10^{-4}$ (tube), $1.545\times10^{-3}$ (waist); $t_1, t_1' \ll T$ |
| $S_{\rm in} = S_{h_{\rm in}}$, $S_{\rm out} = S_{h_{\rm out}}$ exist and are smooth closed surfaces | §2.2. A validated secant refinement of the centre $z_0$, then the sign change of $F$ on $[z_0 - 10^{-10}, z_0 + 10^{-10}]$ (Moore point runs). Levels: $h_{\rm in} = -7.5\times10^{-7}$ / $-4\times10^{-8}$; $h_{\rm out} = 2\times10^{-6}$ | $S_{\rm in}$: $F = \pm2.7\times10^{-12}$ / $\pm7\times10^{-13}$, radii $\le 5\times10^{-14}$ |
| $\theta_+ = h + R$ on every box | Lemma 2, on the Moore ball run's boxes ($6148$ per surface: axis segment plus 4 sub-boxes per step) | $\lvert R\rvert \le 7.95\times10^{-9}$ / $1.62\times10^{-9}$ |
| embedded, star-shaped | $\rho > 0$ and $\sin\alpha < 0$ after the axis segment, $A < 0$ on it; $x\cdot n \ge 0.238$ | — |
| nested, $\Omega \subset U$ | Lohner sets over the same brackets with sub-step enclosures (`bl_foliate.barrier_geometry`), 4000 angular bins | radial gap $\ge 2.24\times10^{-4}$ / $3.5\times10^{-5}$; $\Omega$ needs at most $0.0029$ of the $0.02$ inflation that defines $U$ |
| independent re-check of the barrier inequalities | Lemma 2 on the Lohner sub-step boxes ($\approx 21500$ per surface) | $\theta_+[S_{\rm in}] \in [-7.542, -7.458]\times10^{-7}$ / $[-4.053, -3.947]\times10^{-8}$; $\theta_+[S_{\rm out}] \in [1.996, 2.004]\times10^{-6}$ / $[1.9995, 2.0005]\times10^{-6}$ |
| MOTS | barrier theorem (§2.3) applied to $(\Omega, \gamma(t_1), K(t_1))$ | — |

*Why the MOTS encloses both holes.* $\Omega$ is a shell $S^2 \times [0, 1]$, so $H_2(\Omega) = \mathbb Z[S_{\rm out}]$.
A closed embedded surface that does not separate the two boundary spheres has intersection number 0
with a radial segment from $S_{\rm in}$ to $S_{\rm out}$. A surface homologous to $S_{\rm out}$ has
intersection number 1. So some component of $\Sigma^*$ separates the two spheres. $S_{\rm in}$ is a radial
graph with axis height $0.6326 > d/2 = 0.3831$, so both punctures lie inside it. ∎

*Why the data on $\Omega$ are the slice data of the spacetime.* By Proposition T(4) the CK solution
embeds isometrically into the maximal development, with $t$ the proper time along the normal
geodesics of $\Sigma_0$. Its domain of dependence contains every $(t, x)$ with
$\mathrm{dist}(x, \partial U) > \lvert t\rvert/4$. Here $t \le 1.2\times10^{-5}$, while $\Omega$ lies at distance $\ge 0.02 - 0.0029 > 0.017$ from $\partial U$.

**Size of the margins.** The barrier inequality holds with margin $\lvert h\rvert/\lvert R\rvert \approx 95$ (and
$25$ at $t_1'$). The certified evaluator gives the true $\theta_+ - h$ at sample points on $S_{\rm in}$ as
$\approx 10^{-11}$, so most of the bound on $R$ is overestimation. Two inputs share the bound: the crude box enclosures of
$t^2\partial g_1$ ($\lvert\partial g_1\rvert$ is enclosed to $\sim 10^3$ on a box, against $\sim 40$ at a point), and
the CK tail of $\partial\gamma$ ($6.3\times10^{-8}$ on the tube at $t_1$, $3\times10^{-9}$ on the waist region).
Both enter multiplied by $\partial\Theta/\partial(\partial\gamma) = O(\psi^{-6}) \approx 10^{-2}$.


## 4 · Toward Theorem N: excluding a MOTS from the shell before the merger (not completed)

The aim was the complementary statement:

> **Theorem N (target, not certified).** For $t_0 = 5\times10^{-6}$ and $h = 5\times10^{-8}$, the leaves
> $S_t(z_0)$, $z_0 \in [z_a, z_b] \supset [z_{\rm in}, z_{\rm out}]$, $t \in [0, t_0]$ (the level-$h$
> shooting curves, reflected), foliate a shell $W_t \supset \Omega'$ and satisfy the hypotheses of
> Lemma 3. Hence no MOTS lies in $\Omega'$ for $t \le t_0$.

Combined with Theorem M at $t_1'$, that would bracket the local birth time:
$t_b \in [0.81, 1.06]\,t_*^{\rm lin}$.

**What was built and checked.** `research/calc/bl_foliate.py` carries the leaves as Lohner sets.

- **Machinery.** It reuses `gap3_fold`'s $\sigma$-system, block frames, $z_0$ as a dummy parameter
  and $z$ as the clock before the plane. To that it adds the level-$h$ source, $t$ as a set
  parameter, and **sub-step enclosures in Lohner form**,
  $x(u) \in \varphi_u(m) + J_x(u)(Ar) + (J_x(u)P + J_p(u))q$. These are tested against an exact rotation
  flow with a parameter.
- **The conditions of Lemma 3**, checked on two test pieces of width $6\times10^{-5}$ and
  $1.35\times10^{-4}$ (the first next to the fold), for all $t \in [0, t_0]$:
  - $\theta_+ \in [4.90, 5.10]\times10^{-8}$ (so $\lvert R\rvert \le 1.03\times10^{-9}$);
  - the Jacobi field is positive;
  - $x\cdot n > 0$, $\rho > 0$, $\sin\alpha < 0$.

  For a single leaf $z_0 = 0.632632$ and all $t \in [0, t_0]$, the ridge condition holds as well:
  $u_{\rm end} \in [1.48, 8.01]\times10^{-7} > 0$.

**What fails, precisely.** The ridge condition $u_{\rm end} > 0$ ($F > 0$) cannot be certified on a
$z_0$-piece of width $6\times10^{-5}$ near the fold. The enclosure is $u_{\rm end} \in [-9.5, 10.5]\times10^{-6}$,
against the true range $[1.5, 2.3]\times10^{-7}$. A decomposition of the final set shows why:

- the $t$-interval is harmless;
- the $z_0$-spread inflates the Lohner remainder $Ar$ to $9.6\times10^{-6}$, about 350 times the
  true second-order term $\tfrac12\rho F_{zz}w^2 \approx 3\times10^{-8}$. The interval Jacobian over the
  set re-boxes the second-order dependence at every step, so the inflation grows like (width)² ×
  steps.

Pieces about ten times narrower near the fold would pass, which means about 100 runs of a quarter
of an hour each.

**The fix**, which is the next task, is a convexity argument in $z_0$, as `gap3_fold` used at
$t = 0$:

- carry the second jet ($z_0^2$, slot 4) on each piece, and get $\inf F_{zz} > 0$;
- get $F(z_c, t)$ and $F_z(z_c, t)$ from tight point-$z_0$ runs with the $t$-interval;
- conclude $F \ge F(z_c) - F_z(z_c)^2/(2\inf F_{zz})$ on the piece.

With the measured values ($F(z_c) \approx 5.9\times10^{-7}$, $F_z^2/(2F_{zz}) \le 2.2\times10^{-7}$) this
closes on pieces of width $6\times10^{-5}$.

## 5 · The explicit merger

> **Corollary (an explicit apparent-horizon merger).** In the geodesic-slicing development of
> equal-mass Brill–Lindquist data with $d = 0.7662$:
>
> - $\Sigma_0$ contains no closed embedded MOTS enclosing both punctures (Theorem F′);
> - $\Sigma_{t}$ contains a smooth, embedded, stable MOTS enclosing both punctures, for
>   $t = 6.5\times10^{-6}M$ and for $t = 1.2\times10^{-5}M$ (Theorem M).
>
> This happens where the spacetime metric is the certified convergent series of route E: an explicit
> metric representation with certified errors. So the apparent horizon of the geodesic slicing
> changes from two components to one enclosing both holes within $6.5\times10^{-6}M$, that is
> $1.06$ times the linear-law onset time.

What this adds to Theorem H:

- at this separation, $\delta_0$ is no longer needed for the existence half;
- the upper bound on the birth time is explicit and within 6 % of the linear law;
- the lower bound $t_b > t_0$ is not certified (§4);
- the exclusion of other common MOTSs at $0 < t < t_1'$ is open everywhere on the slice, as before.

## 6 · What is and is not proved, and what blocks the rest

**What the result is.**

- The spacetime is the maximal vacuum development of explicit, asymptotically flat, time-symmetric
  two-black-hole data (equal-mass BL, $d = 0.7662$).
- On the region where the merger is proved, its metric is the convergent series of route E, with
  certified errors and a certified evaluator (category B, notes/16).
- The merger statement is a theorem about that metric (category C). In the geodesic slicing, the
  apparent horizon has two components at $t = 0$. By $t_1 = 1.2\times10^{-5}M$ it has a component
  enclosing both holes.
- Before this session, both halves existed but did not meet at any explicit separation.

**What it is not.**

1. *Not a closed-form metric.* The category-A obstructions of notes 14–19 are unchanged: Petrov type I,
   hence no algebraically special form; no conformally flat geodesic or maximal-type leaves; not
   silent; and integrable methods need two commuting Killing fields.
2. *Apparent, not event, horizons.* The event horizon is teleological and needs the global
   development: complete $\mathscr I^+$ and the final state. For near-critical data it is expected
   to be connected already on $\Sigma_0$ (numerical relativity; not proved here).
3. *Gauge.* "Two components at $t = 0$, one by $t_1$" refers to the geodesic slicing. The common
   MOTS on $\Sigma_{t_1}$ is a spacetime object. Which slice first carries one is gauge dependent.
4. *Not the first common apparent horizon.* For $0 < t < t_1'$, no common MOTS is excluded
   anywhere: not near $S_c$, since Theorem N is not completed (§4), and not elsewhere, where
   Theorem F′ covers only $t = 0$. So $6.5\times10^{-6}$ is an upper bound for the birth time of the
   common apparent horizon, not a bracket.
5. *Not trapped.* The inner barrier is outer trapped but not trapped. On $S_c$,
   $\mathrm{tr}_SK = tp + O(t^3)$ with $p = \mathrm{Ric}(\nu,\nu)$ positive on a band (floating point:
   $p \in [-0.58, 0.11]$). So $\theta_- = 2tp - \theta_+ > 0$ there, and no surface $C^1$-close to
   $S_c$ is trapped at small $t$.

**What blocks the next steps**, in order of value.

1. **A lower bound for the birth time.**
   - *Near $S_c$:* complete Theorem N with the convexity argument of §4. This needs second
     $z_0$-jets on Lohner pieces plus point-$z_0$ runs; roughly 15 runs of 20–30 minutes.
   - *Globally* (so that the birth time is bracketed for the whole slice), two ingredients are
     needed:
     - a CK certificate on the region where Theorem F′ works, $B(0, 1.05)$ minus the throat balls,
       instead of the tube $U$ of radius $0.02$. Only $t \lesssim 10^{-5}$ is needed, so crude
       constants suffice.
     - a perturbed Theorem F′ sweep. The outermost common MOTS of $\Sigma_t$, if any, is unique,
       hence axisymmetric and reflection symmetric, so the reduction of notes/12 §1 still applies.
       Run the same exclusion events for the level-$h$ curves with $\lvert h\rvert \le \sup\lvert R\rvert$,
       and $t$ as a parameter. Every exclusion event of F′ is an open condition, so it survives the
       perturbation away from $S_c$. Near $S_c$, Theorem N takes over.

     The cost is of the order of the F′ sweep (hours on four cores).
2. **An interval of separations.** Theorem M holds at one $d$. Wide $z_0$- or $d$-brackets blow up
   in the Moore integrator (radius $5\times10^{-5}$ fails; $10^{-8}$ fails in the last segment), so
   an interval needs Lohner sets (`gap3_fold`'s machinery, extended in `bl_foliate`), with the same
   convexity treatment of the $z_0$-spread as Theorem N. Alternatively
   a smoothing lemma for shallow creases would remove the closing condition (idea log, session 16).
3. **Event horizons.** These need global existence and the final state, i.e. nonlinear stability of
   Schwarzschild for these large data. Out of reach. The only rigorous event-horizon mergers known
   to us are the extreme-mass-ratio ones (notes/11, conditional on our reading of Hintz).
4. **A closed form.** No candidate. A proof that none exists would first need a precise class of
   "closed forms" (for example, differentially algebraic over elementary functions) and a
   differential-Galois-type obstruction for the reduced 2+1 system. Nothing of this kind is
   available.


## References

- [AM09] L. Andersson and J. Metzger, *The area of horizons and the trapped region*, Commun. Math.
  Phys. 290 (2009) 941–972, arXiv:0708.4252. §2 conventions ($\theta_\pm = P \pm H$), Theorem 3.1,
  Theorem 4.1, Theorem 5.1 and Lemma 2.14 (Kriele–Hayward smoothing) **read in the PDF** this session.
- [KH97] M. Kriele and S. A. Hayward, *Outer trapped surfaces and their apparent horizon*, J. Math.
  Phys. 38 (1997) 1593. Cited through [AM09] Lemma 2.14 only.
- Notes 09 (Theorem C), 12 (Theorem F′), 15–17 (Proposition T, CK′, the extended certificates), 16
  (the evaluator and its tail bounds).

## Discussion prompts

- Existence used only $C^0$ control of $\theta_+$ on two surfaces, and exclusion only $C^0$ control
  on a foliation. Which other horizon properties (area, position, stability of the new MOTS) can be
  certified from $C^0$ information alone, without the spectrum of the Jacobi operator?
- The newborn common MOTS is marginally *outer* trapped, but its inward expansion is positive on a
  band near the waist, where $\mathrm{Ric}(\nu,\nu) > 0$. Is that generic for head-on mergers? What
  does it say about the signature of the marginally trapped tube at its birth (AMMS)?
- Both certificates use $t \lesssim 10^{-5}M$, while the CK certificate reaches $8.7\times10^{-4}M$.
  How far in $d - d_c$ can the barrier method go before the CK tails or the crude $t^2\partial g_1$
  enclosures dominate? What would a CK certificate on a larger domain buy (§6, item 1)?

