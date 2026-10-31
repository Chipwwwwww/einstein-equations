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
$G = c = M_{\rm ADM} = 1$; $d_c = 0.7662$). §7 (session 9) extends it to every $d$ in
$(d_c, 0.767578]$, an interval starting at the critical separation (Theorem F′).

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

Its trust base is Arb plus about 1,100 lines of our own code (with
docstrings):

- `lohner.py`;
- `gap3_sweep.py`;
- `gap3_neck.py`;
- the jet, axis and Taylor-step parts of `fold_proof.py`;
- the barrier part of `bl_global.py`.

That code is tested against flows with exact solutions and deliberately
falsified (§4), but nobody outside this project has reviewed it. Theorem F
is for $d = 0.78$ only. Theorem F′ (§7) covers $(d_c, 0.767578]$, with a
local lemma at the fold, where a common surface is about to appear.

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
as a jet. All 640 pieces ($\lambda$-width $2^{-11}$, $z_w$-width $2^{-9}$,
two directions) passed on the first try. The certified values satisfy
$|\phi| \ge 0.527$ at $\sigma = 0$ and $\sigma = \pm1.5$, with the required
signs; the catenoid values are $1$ and $0.84$. That took 25 minutes on
4 cores.

**Consequence.** If a profile curve has such a waist, $\phi$ has a zero on
each side of it. Between the two nearest zeros lies an annulus on which
$\phi$ is a Dirichlet eigenfunction of the Jacobi operator with eigenvalue
$0$. The annulus is a proper subdomain of the closed surface. By strict
domain monotonicity, the first eigenvalue of the closed surface is
therefore negative: it is **unstable**, and cannot be $S^*$. The annulus
does lie on the surface. The curve cannot reach the axis within
$|\sigma| \le 1.5$, because axis points are fixed points of the system.

Wider necks look unstable too. In floating point the sign pattern holds up
to radius $0.04$, using $\sigma = \pm 2$. But certifying that with
plain interval steps needs impractically small parameter pieces. So the
sweep integrates through the medium necks instead (§4).

## 3 · Getting into the neck: the approach lemma

The lemma above needs the waist radius and height. These come from a state
on the way in. Let $K$ bound $4(|\psi_z| + |\psi_\rho|)/\psi$ on the
rectangle $0 \le \rho \le 0.02$, $0.1797 \le z \le 0.2422$ (the certified $z_w$
range). Interval arithmetic gives $K = 16.53$. Since $|u|, |v| \le \rho$, we
get $|P| \le K\rho$ there.

Suppose a curve is at a state with $u < 0$ (moving toward the axis) and
$\rho \le \rho_e$. While $u < 0$, $\rho$ is a parameter along the curve, and
$dv/d\rho = P$, so $|dv/d\rho| \le K\rho$. Two cases are possible.

- **$u < 0$ forever.** Then $\rho \to 0$, because otherwise
  $u' \ge \rho(1 - K\rho) > 0$ would bring $u$ to $0$. The curve ends on the
  axis at a point with $z > 0$, without reaching the plane, so it is not
  the profile of $S^*$. (That profile meets the axis only at
  $(0, \pm z_0)$.)
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

The range $[139/256, 269/256] \supset [0.545, 1.05]$ is covered by a stored
partition into $z_0$-intervals,
[`research/calc/data/gap3_partition_d078.txt`](../research/calc/data/gap3_partition_d078.txt).
Each interval is run, and any interval that is not excluded is bisected.

The partition comes from a first sweep: uniform intervals of width $2^{-8}$,
bisected on failure. That sweep left two gaps, $[0.59077, 0.59254]$ and
$[0.59322, 0.59364]$. These are the curves whose necks are too wide for the
thin-neck lemma (radius $0.005$ to $0.02$) but still narrow enough to make
the outgoing family spread fast. Lohner's first-order enclosure then
explodes just before the curve leaves the ball, unless the $z_0$-interval is
narrow. In the gaps, the stored intervals therefore have width about
$1.2\times10^{-3}\,(z_0 - z_{\rm top})$, between $2^{-20}$ and $2^{-17}$, a
rule calibrated on a few runs.

**Falsification checks.**

- At $d = 0.76 < d_c$, two common minimal surfaces exist, at axis heights
  $0.62181$ and $0.63778$. The same code leaves the width-$2^{-14}$
  intervals containing them unexcluded, and excludes an interval between
  them.
- The thin-neck certificate fails, as it must, when asked for a sign change
  at $\sigma = 1.0$, before the catenoid's Jacobi zero.
- The neck Jacobi field computed by a second integrator agrees with the
  first at test points: Lohner steps on the 10-dimensional variational
  system ([`research/calc/gap3_neck_lohner.py`](../research/calc/gap3_neck_lohner.py)).

All three are in `tests/test_gap3.py`.

## 5 · Results

Both certificates ran on 2026-10-29/30, on 4 cores.

- **Thin-neck lemma**: 640 pieces, all certified on the first try, margin
  $\ge 0.527$. That took 25 minutes; the CI run of the slow test repeats it.
- **Sweep**: all 2,452 $z_0$-intervals of the stored partition were excluded
  on the first try, with no bisection. That took 2.8 hours (11.4 CPU-hours).
  The log is
  [`research/calc/data/gap3_sweep_d078.log`](../research/calc/data/gap3_sweep_d078.log).
  The excluded intervals cover $[139/256, 269/256]$ without gaps.

The events fall into five consecutive ranges of $z_0$:

| $z_0$ | Event | Intervals | What the curve does |
| --- | --- | --- | --- |
| $[0.5430, 0.5884]$ | E3 | 182 | dives into the throat around $p_+$ |
| $[0.5884, 0.5908]$ | N | 38 | thin neck at the bottom of the individual horizon (on both sides of $z_{\rm top} = 0.58968$) |
| $[0.5908, 0.5973]$ | E2 | 1,579 | wider neck, then flares out of the ball |
| $[0.5973, 0.8730]$ | E1 | 593 | meets the plane at a non-right angle |
| $[0.8730, 1.0508]$ | E2 | 60 | leaves the ball |

This matches the floating-point landscape of session 6. Together with §§1–3,
it proves Theorem F. ∎

## 6 · What is and is not proved

| Ingredient | Status |
| --- | --- |
| barriers, $K$, thin-neck lemma, sweep | computer-assisted (Arb / mpmath intervals) |
| approach lemma, catenoid Jacobi field, reduction | analytic, written out above and in session 6 |
| stability of the outermost minimal surface | literature (Andersson–Mars–Simon, Andersson–Metzger); not re-derived |
| other separations $d > d_c$ | *(session 9)* proved for every $d \in (d_c, 0.767578]$: Theorem F′, §7. Not yet for $(0.767578, 0.78)$ |
| combination with the onset theorem of note 10 | *(session 9)* done: Theorem F′ starts at $d_c$, so it combines with the onset theorem whatever $\delta_0$ is (Theorem H, note 13) |
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

## 7 · An interval of separations reaching $d_c$: Theorem F′

> **Theorem F′ (computer-assisted).** For every separation $d$ with
> $d_c < d \le 0.767578125$, no closed embedded minimal surface of the
> Brill–Lindquist data encloses both punctures. The apparent horizon of
> $\Sigma_0$ consists of the two individual horizons. Here
> $d_c = 0.766197428240 \pm 1.1\times10^{-13}$ is the fold separation of
> Theorem C ([note 09](09-merger-hypotheses.md)), in units
> $G = c = M_{\rm ADM} = 1$.

The interval starts at the critical separation itself. At $d = d_c$ a
degenerate common minimal surface exists; just below $d_c$ there are two
(the fold). So $d_c$ is exactly where the common horizon of the initial data
appears. Combined with the onset theorem of [note 10](10-onset-theorem-proof.md),
this gives the head-on collision of [note 13](13-head-on-collision.md).

The code reaches up to $0.78125$. Each further $2^{-9}$ in $d$ costs
about 2.5 hours on 4 cores, so the theorem is stated for the part that has
been run. Joining it to Theorem F at $0.78$ takes about 18 more hours of
computing and is optional.

The proof follows §§1–5 with $d$ an interval parameter everywhere. Three
things are new:

- a waist event, with a band certificate that holds for all $d$ at once;
- a calibrated partition;
- a local lemma at the fold, where the sweep cannot work.

### 7.1 · What carries over, with $d$ an interval

- **Barriers** (`tests/test_bl_global.py`), for all
  $d \in [0.765625, 0.78125]$ at once:
  - throat spheres of radius $s \le 0.15$ have $H < 0$ (bound $-0.037$);
  - spheres $|x| \ge 1.05$ have $H > 0$ (bound $0.060$).

  The throat radius had to shrink from $0.155$, which fails at $d = 0.766$.
- **Relative height.** The state is $(\rho, \zeta, u, v)$ with
  $\zeta = z - d/2$, the height above the upper puncture. Then $d$ enters
  the equations only through the lower puncture, at $\zeta = -d$. With $z$
  itself as the state, a box in $d$ decorrelates $z$ from $d$ inside $\psi$,
  and the sets blow up at $\sigma \approx 16$.
- **Boxes in $(\delta, d)$.** The start is
  $z_0 = z_{\rm ref} + \delta + 0.543\,(d - 0.78)$. The slope follows the
  top of the individual horizon, so a box stays in the same place relative
  to the delicate region as $d$ changes. Each box is one Lohner set, with $d$
  kept out of the moving frame (`lohner.PSet`).
- **Thin necks (event N).** The thin-neck lemma is certified in the
  variables $(\lambda, \zeta_w)$ for all $d \in [0.765625, 0.78125]$ at once,
  on $\lambda \le 0.0117$ and $\zeta_w \in [-0.234, -0.133]$, both directions
  ([`research/calc/gap3_neck.py`](../research/calc/gap3_neck.py),
  `certify_region_rel`). Nine runs, 2,664 pieces, worst margin $0.0029$
  ([`research/calc/data/gap3_neck_interval.txt`](../research/calc/data/gap3_neck_interval.txt)).

  The $\zeta$-window had to be wider than at $d = 0.78$. As $z_0$ moves away
  from the top of the individual horizon, the neck height moves too, and the
  approach lemma's drift (up to $0.015$) carries the window past its old
  ends.

### 7.2 · Medium necks: the waist event W and the band certificate

At $d = 0.78$, the medium necks ($\lambda \approx 0.005$–$0.04$) needed
$z_0$-intervals as narrow as $2^{-20}$ (§4). Over an interval of $d$ that
would mean millions of boxes. Instead, the sweep stops at the waist.

- **Event W.** Take the window from the last step after which $u < 0$ on the
  whole set to the first step after which $u > 0$ on it. Every curve has its
  waist ($u = 0$) in one of the window's step boxes. The waist box is the
  union of their $\rho$-ranges times the union of their $\zeta$-ranges, and
  the direction is $\mathrm{sign}\,v$. W is recorded when $v$ has one sign
  and $\lambda < 0.05$.
- **Band certificate.** The sign pattern of the neck Jacobi field (§2) is
  certified on dyadic cells in $(\lambda, \zeta_w)$, each for all
  $d \in [0.765625, 0.78125]$ at once
  ([`research/calc/gap3_neck_lohner.py`](../research/calc/gap3_neck_lohner.py),
  `certify_rel`). The method is Lohner, with $(\lambda, \zeta_w, d)$ as
  parameters. Each direction stops at the first $\sigma \ge 1.25$ where
  $s\phi > 0$ is certified, up to $\sigma = 3$. Cells start at $2^{-10}$,
  and failures are quartered. A waist box counts once every cell it meets is
  certified.

The instability argument of §2 applies unchanged. A neck whose Jacobi field
changes sign on both sides makes the surface unstable, so it cannot be
$S^*$.

### 7.3 · The partition

A first sweep with coarse boxes stalled almost everywhere within $0.05$ of
the top of the individual horizon. Lohner's first-order enclosures explode
where the curves turn.

Test boxes at both ends of the $d$-range then fixed a table of bands
(`gap3_interval.BANDS`). The $\delta$-width runs from $2^{-10}$ far from the
top down to $2^{-15}$ next to it, and the $d$-width is $2^{-10}$ near the
top. This is close to what the single-$d$ sweep at $0.78$ needed. Boxes that
still fail are bisected, and the sweep runs lowest $d$ first, so the
certified range $(d_c, d_k]$ grows steadily.

### 7.4 · Near the fold: a local lemma

As $d \downarrow d_c$, the curves near the fold point $z_c = 0.63266$ meet
the plane almost at a right angle, since $F = \alpha + \pi/2 \to 0$ there.
Event E1 needs $u \ne 0$ at the crossing. Up to $d = 0.7676$ it fails even
for boxes of $\delta$-width $2^{-15}$, because $F \approx F_d\,(d - d_c)$ is
below $2\times10^{-3}$. So the sweep leaves out the *fold window*
$z_0 \in [0.62891, 0.63696]$ for $d \le 0.767578$. A lemma covers it instead.

**Local fold lemma.** Let $Z$ be an interval of $z_0$ containing the fold
point $z^*$, and $d_1 > d_c$. Suppose that for all $z_0 \in Z$ and
$d \in [d_c, d_1]$ the curve reaches the plane with $v < 0$, and that

- $\partial^2_{z_0}F > 0$ for $d$ at $d_c$;
- $\partial_dF > 0$ on all of $[d_c, d_1]$.

Then $F > 0$ on $Z \times (d_c, d_1]$. *Proof.* Theorem C gives
$F = \partial_{z_0}F = 0$ at $(z^*, d_c)$. By Taylor's theorem,

```math
F(z_0, d_c) = \tfrac12\,\partial^2_{z_0}F(\xi, d_c)\,(z_0 - z^*)^2 \ge 0 .
```

Then

```math
F(z_0, d) = F(z_0, d_c) + \int_{d_c}^{d}\partial_dF(z_0, s)\,ds > 0 .
```

∎

**Enclosures.** Jets in $(z_0, d)$ are carried through Lohner's method
([`research/calc/gap3_fold.py`](../research/calc/gap3_fold.py)). Three
choices were each needed to make them tight; the failures are recorded in
the module docstring and the idea log.

1. **A block-diagonal frame**, one block per jet order. The base curve does
   not see its jets, but a full QR frame lets the jets' large widths leak
   into it.
2. **$z_0$ as a parameter of the set**, next to $d$. Its first-order effect
   is then carried exactly and correlated across all jet orders.
3. **$z$ as the clock near the plane.** The last stretch uses
   $\tau \in [0, 1]$ with $z = z_{\rm sw}(1 - \tau)$, as in Theorem C.
   Every curve crosses at $\tau = 1$, where the jets are the derivatives of
   the crossing state:

   ```math
   F_z = \alpha_x\cdot x_z,\qquad
   F_{zz} = \alpha_{xx}(x_z, x_z) + \alpha_x\cdot x_{zz},\qquad
   F_d = \alpha_x\cdot x_d .
   ```

As a check, a thin run at the fold point reproduces Theorem C's
$\partial^2_{z_0}F$ and $\partial_dF$. Those were computed by a different
integrator in a different parametrisation (`tests/test_gap3_fold.py`).

**Results.** Pieces of width $2^{-12}$ in $z_0$; wider ones blow up.

| Part | Region | Pieces | Bounds |
| --- | --- | --- | --- |
| 1 | $d \in [d_c - 10^{-13},\ 0.7662354]$ | 33, none split | $\partial^2_{z_0}F \ge 125.1$, $\partial_dF \ge 0.591$ |
| 2 | $d \in [0.7662354,\ 0.7677002]$ | 198, $2^{-12}$ in $d$, none split | $\partial_dF \ge 0.177$ |

The records are
[`gap3_fold_window.txt`](../research/calc/data/gap3_fold_window.txt) and
[`gap3_fold_fd.txt`](../research/calc/data/gap3_fold_fd.txt). The
holes the sweep leaves out map into the fold window, which is checked by
interval arithmetic (`gap3_interval.holes_inside`).

### 7.5 · Results

Run on 2026-10-30/31, on 4 cores. The stored log is
[`research/calc/data/gap3_interval_stage1.log`](../research/calc/data/gap3_interval_stage1.log).

- **Coverage.** The partition of
  $[\delta_{\min}, \delta_{\max}] \times [d_c - 10^{-13},\ 0.767578]$, minus the fold
  holes, has 2,644 boxes. All are excluded, 33 of them after bisection: 2,677 excluded boxes in
  total. `check_cover` confirms coverage, the holes map into the fold window, and the
  $z_0$-window $[d/2 + s_*, R_*]$ is covered.
- **Events** (merged across $d$), in the same order as at $d = 0.78$ (§5). The waist event W
  replaces the medium-neck E2 range there:

| $\delta$ (offset from the top of the individual horizon) | Event | Boxes |
| --- | --- | --- |
| $[-0.0625, -0.0027]$ | E3 | 575 |
| $[-0.0027, +0.0022]$ | N | 160 |
| $[+0.0022, +0.0095]$ | W | 291 |
| $[+0.0094, +0.2861]$, minus the fold window | E1 | 1,464 |
| $[+0.2861, +0.4688]$ | E2 | 187 |

- **Band certificate** for all 382 waist boxes recorded with $d \le 0.767578$ (291 of them in
  the final cover): 609 cells of $2^{-10}$, of which 573 certified at once, with the other 36
  quartered and all 144 quarters certified. That makes 717 cells, none failing, and every waist
  box is covered
  ([`gap3_interval_stage1_verify.txt`](../research/calc/data/gap3_interval_stage1_verify.txt)).

**Replay status.**

- The sweep was run once.
- `verify` replays the band certificate and checks coverage. It does not re-run the sweep's
  boxes; that is a resume of `python -m research.calc.gap3_interval 4` on an empty log, about
  7 hours.
- The neck and fold-lemma certificates were each run once, and their records are stored.

Together with the barriers, the reduction of §1 and the lemmas of §§2–3, these prove
Theorem F′. ∎

**Falsification.** At $d = 0.76595 < d_c$, two common minimal surfaces
exist, at $z_0 = 0.630956$ and $0.634135$. The interval code leaves the
boxes containing them unexcluded, and excludes a box between them by E1
(`tests/test_gap3_interval.py`).

**Trust base.** Everything in §6, plus `gap3_interval.py`, `gap3_fold.py`
and `certify_rel`. As for Theorem F, nobody outside this project has
reviewed the code.

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
