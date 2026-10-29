# 12 · No common horizon at $t = 0$: a computer-assisted proof at one separation

Note [10](10-onset-theorem-proof.md) proved that, for time-symmetric
Brill–Lindquist data slightly beyond the critical separation $d_c$, a common
apparent horizon appears near the critical surface $S_c$ at a time linear in
$d - d_c$. That theorem is local. It says that no common horizon exists
*near* $S_c$ at $t = 0$, but not that none exists *anywhere*. If one existed
elsewhere, the "birth" of the common horizon at $t_*$ would not be the first
moment the two holes share a horizon. This was gap 3 of
[`research/status.md`](../research/status.md).

This note closes gap 3 for one separation, $d = 0.78$ (in units
$G = c = M_{\rm ADM} = 1$; $d_c = 0.7662$).

> **Theorem F (computer-assisted).** Take the time-symmetric Brill–Lindquist
> data $\gamma = \psi^4\delta$, $\psi = 1 + \frac{1}{4r_+} + \frac{1}{4r_-}$
> (two punctures of bare mass $\frac12$ at $z = \pm d/2$), with $d = 0.78$.
> Then no closed embedded minimal surface encloses both punctures. The
> apparent horizon of $\Sigma_0$ consists of the two individual horizons.

On a time-symmetric slice, MOTSs are minimal surfaces ($\theta_\pm = \pm H$),
so "common MOTS" and "minimal surface enclosing both punctures" mean the same
thing here.

*Status.* The proof combines:

- analytic lemmas, written out below;
- interval and ball-arithmetic certificates (Arb, through python-flint) for
  inequalities and ODE enclosures.

Its trust base is Arb plus about 700 lines of our own code. That code is
tested against flows with exact solutions, but nobody outside this project
has reviewed it. The result is for $d = 0.78$ only. The same programme
should work for any fixed $d > d_c$, but runs longer as $d \to d_c$, where a
common surface is about to appear.

## 1 · Reduction to a family of curves

Session 6 established three facts ([`research/calc/bl_global.py`](../research/calc/bl_global.py),
claims 51–53).

1. **Barriers.** Coordinate spheres of radius $s \le 0.155$ about each
   puncture have $H < 0$ (normal pointing away from the puncture). Spheres
   $|x| \ge 1.05$ have $H > 0$.
2. **Outermost surface.** If a common minimal surface existed, the outermost
   minimal surface would have a component $S^*$ enclosing both punctures
   (Andersson–Metzger). The outermost surface is unique, so $S^*$ is
   invariant under rotations about the axis and under $z \to -z$. It is also
   **stable**: its Jacobi operator has no negative eigenvalue. Otherwise,
   pushing it outward along the first eigenfunction would give a trapped
   surface outside it.
3. **Profile.** $S^*$ is a sphere of revolution. Its profile curve runs from
   an axis point $(0, z_0)$ to the plane $z = 0$. It meets that plane once,
   at a right angle. A second crossing would meet its own mirror image, so
   the surface would not be embedded. The curve stays inside $B(0, 1.05)$ and
   outside the balls $B(p_\pm, 0.155)$, so $z_0 \in [0.545, 1.05]$.

The profile solves the ODE for $H = 0$. With $\rho$ the distance from the
axis, $s$ flat arclength, and $\alpha$ the angle of the tangent, we change
the time to $ds = \rho\,d\sigma$ and use $(u, v) = \rho(\cos\alpha, \sin\alpha)$:

```math
\rho' = u,\qquad z' = v,\qquad u' = \rho - vP,\qquad v' = uP,\qquad
P = \frac{4\,(u\,\psi_z - v\,\psi_\rho)}{\psi}.
```

This system is regular up to the axis (`tests/test_bl_global.py`), and
$u^2 + v^2 = \rho^2$ is conserved. The axis points are fixed points with
eigenvalues $\pm 1$: saddles. A curve therefore never reaches the axis in
finite $\sigma$-time, only in the limit.

**The difficulty.** For $z_0$ just above the top of the individual horizon,
$z_{\rm top} = 0.5896816$, the shooting curve runs down along that horizon.
It passes close to the saddle at its bottom, $(0, 0.2126)$, through a neck
of radius about $4(z_0 - z_{\rm top})$, and flares out. As
$z_0 \to z_{\rm top}$ the necks pinch to zero. No finite computation follows
all of them, so the necks are handled by a lemma instead.

## 2 · Thin necks are unstable

Near a waist of radius $\lambda$ at height $z_w$ (a point where the curve is
vertical, $u = 0$), rescale:

```math
\rho = \lambda R,\quad z = z_w + \lambda Z,\quad u = \lambda U,\quad v = \lambda V .
```

The system keeps its form, with $P = 4\lambda(U\psi_z - V\psi_\rho)/\psi$.
At $\lambda = 0$, $P = 0$ and the waist is a **catenoid**:
$R = \cosh\sigma$, $Z = \pm\sigma$. This is the usual picture of a thin neck
in a minimal surface. Gravity (the conformal factor) is invisible at the
scale of the neck.

Changing $\lambda$ at fixed $z_w$ moves the curve through minimal surfaces.
The normal component of that motion is therefore a Jacobi field. Up to a
positive factor it is

```math
\phi = -V\,(R + \lambda R_\lambda) + U\,(Z + \lambda Z_\lambda),
\qquad \phi\big|_{\lambda = 0} = \mp\,(\cosh\sigma - \sigma\sinh\sigma).
```

On the catenoid this is the dilation field. It vanishes where
$\sigma\tanh\sigma = 1$, at $\sigma = \pm 1.1997$
([`tests/test_gap3.py`](../tests/test_gap3.py) checks that it solves the
catenoid's Jacobi equation). That is the classical fact that a catenoid
piece longer than this is unstable.

**Thin-neck lemma (certified).** Let $\lambda \in [0, 0.00488]$ and
$z_w \in [0.1797, 0.2422]$, for a waist moving up or down. Then $\phi(0)$ and
$\phi(\pm 1.5)$ have opposite signs
([`research/calc/gap3_neck.py`](../research/calc/gap3_neck.py)).

The certificate splits the parameter rectangle into pieces and integrates
the rescaled system on each piece. The integration uses Taylor steps in ball
arithmetic, with $\lambda$ and $z_w$ as balls and $\partial_\lambda$ carried
as a jet. (Certified pieces and smallest margin: see §5.)

**Consequence.** If a profile curve has such a waist, $\phi$ has a zero on
each side of it. Between the two nearest zeros lies an annulus on which
$\phi$ is a Dirichlet eigenfunction of the Jacobi operator with eigenvalue
$0$. The annulus is a proper subdomain of the closed surface. By strict
domain monotonicity, the first eigenvalue of the closed surface is
therefore negative: it is **unstable**, and cannot be $S^*$.

## 3 · Getting into the neck: the approach lemma

The lemma above needs the waist radius and height. These come from a state
on the way in. Let $K$ bound $4(|\psi_z| + |\psi_\rho|)/\psi$ on the
rectangle $[0, 0.02] \times [z_{\rm lo}, z_{\rm hi}]$, the certified $z_w$
range. Interval arithmetic gives $K = 16.53$. Since $|u|, |v| \le \rho$, we
get $|P| \le K\rho$ there.

Suppose a curve is at a state with $u < 0$ (moving toward the axis) and
$\rho \le \rho_e$. While $u < 0$, $\rho$ is a parameter along the curve, and
$dv/d\rho = P$, so $|dv/d\rho| \le K\rho$. Two cases are possible.

- **$u < 0$ forever.** Then $\rho \to 0$, because otherwise
  $u' \ge \rho(1 - K\rho) > 0$ would bring $u$ to $0$. The curve closes up on
  the axis at a point with $z > 0$. It is a sphere around $p_+$ alone, not
  the profile of $S^*$.
- **$u$ reaches $0$ at a waist.** Its radius is

  ```math
  \rho_w = |v_w| \le |v_e| + \tfrac12 K\rho_e^2 .
  ```

  Its height moves by at most

  ```math
  |z_w - z_e| \le \frac{2\rho_e/e + K\rho_e^2/2}{\sqrt{1 - K\rho_e - K^2\rho_e^2/4}} .
  ```

  This uses $\rho^2 - v^2 \ge c\,(\rho^2 - \rho_w^2)$, together with
  $\rho_w\,\mathrm{arccosh}(\rho_e/\rho_w) \le \rho_w\log(2\rho_e/\rho_w) \le 2\rho_e/e$.
  The inequalities are checked in `tests/test_gap3.py`.

If both bounds land inside the thin-neck rectangle, the curve is excluded
either way. This is event **N** below.

## 4 · The sweep

Each $z_0$-interval is integrated as a set: an enclosure of all its curves,
from the axis outward, in the $\sigma$-system. Plain interval boxes fail
here. A box carried around the individual horizon and past a saddle is
re-boxed at every step, and it blows up (the *wrapping effect*). We use
**Lohner's method** ([`research/calc/lohner.py`](../research/calc/lohner.py)).
The enclosure is a parallelepiped $m + Ar$ whose frame $A$ follows the flow.
Each step is a Taylor polynomial with a Lagrange remainder, the step map's
Jacobian comes from dual numbers, and the frame is re-orthogonalised by QR.
`tests/test_lohner.py` checks it against rotation, a hyperbolic saddle and a
Riccati equation, all with exact solutions.

A $z_0$-interval is **excluded** when the enclosure of all its curves shows
one of these events:

| Event | Condition (whole enclosure) | Why the curve is not $S^*$ |
| --- | --- | --- |
| E1 | $z < 0$; every step that allowed $z = 0$ had $u \ne 0$ | first plane crossing is not at a right angle |
| E2 | outside $B(0, 1.05)$ | $S^*$ lies inside (barrier) |
| E3 | inside a throat ball $B(p_\pm, 0.155)$ | $S^*$ avoids them (barrier) |
| N | $u < 0$, $\rho \le 0.02$, approach-lemma bounds inside the neck rectangle | closes on the axis, or has an unstable thin neck |

Intervals start at width $2^{-8}$ on $[139/256, 269/256] \supset [0.545, 1.05]$
and are bisected until they are excluded. (Results: §5.)

## 5 · Results

**In progress.** The certificate runs are under way. Until this section
records them, Theorem F is a target, not a result.

## 6 · What is and is not proved

| Ingredient | Status |
| --- | --- |
| barriers, $K$, thin-neck lemma, sweep | computer-assisted (Arb / mpmath intervals) |
| approach lemma, catenoid Jacobi field, reduction | analytic, written out above and in session 6 |
| stability of the outermost minimal surface | literature (Andersson–Mars–Simon, Andersson–Metzger); not re-derived |
| other separations $d > d_c$ | **not proved**. The same code applies at any fixed $d$; no uniform argument |
| independent review of the code | **not done** |

**Prior work.**

- Non-existence at **large** separation is a theorem. Dain, Weinstein and
  Yamada (*Class. Quantum Grav.* 28 (2011) 085015, arXiv:1012.4190,
  Proposition 1, after Chruściel–Mazzeo 2003, Thm 3.2) show that for
  sufficiently small mass-to-separation ratio, no connected minimal surface
  encloses both punctures. Their proof plays the Penrose inequality against
  curvature estimates for stable minimal surfaces. It gives no explicit
  threshold.
- Near the critical separation, non-existence had only numerical support
  (Brill–Lindquist 1963, Gibbons 1972, and later horizon finders).

Theorem F is explicit, at $d/d_c = 1.018$, within 2% of the critical
separation, where compactness arguments say nothing. As far as our bounded
search found, it is **potentially new**, as a computer-assisted proof.

## Discussion prompts

- The thin-neck lemma is a gravitational version of the instability of long
  catenoid pieces. How thin must a neck be before the conformal factor stops
  mattering, and could this be made quantitative without a computer?
- The proof treats one separation. What would a proof for all $d > d_c$ need:
  a continuity argument in $d$, or a monotonicity of the shooting landscape?
- Near $z_{\rm top}$ the excluded surfaces look like a sphere hugging the
  individual horizon, connected through a thin neck to an outer shell. Do
  such surfaces exist as (unstable) minimal surfaces for some $d$, and what
  is their area compared with $16\pi M^2$?
