# 11 · Event horizons of extreme-mass-ratio mergers

Hintz (2024, *Gluing small black holes along timelike geodesics III*,
arXiv:2408.06715) constructs vacuum spacetimes in which a light Kerr black
hole, of mass $\epsilon m$, falls into a slowly rotating unit-mass Kerr or
Kerr–de Sitter black hole. Afterwards the spacetime settles down to a single
black hole. He leaves "the construction of the event horizon, the structure of
the black hole region—conjecturally connected—" to future work (§1, p. 3–4).
Hypothesis H3 of [09](09-merger-hypotheses.md) asked whether the event
horizons in these spacetimes merge in the sense of [07](07-the-merger-problem.md).

This note answers it for Λ > 0 and reduces Λ = 0 to two properties of the
Kerr-stability spacetime.

- **Theorem D (Λ > 0).** In Hintz's merger spacetimes (his Theorem 6.1) take
  $\mathcal B$ to be the events that cannot send a signal into the expanding
  region beyond the cosmological horizon. Then, for all small $\epsilon$:
  - $\mathcal B$ contains the interiors of both holes, each with closed
    trapped surfaces;
  - both horizons are located to within any prescribed margin;
  - before the plunge the two interiors lie in **different** components of
    $\mathcal B\cap\Sigma_t$, and afterwards in the **same** one;
  - $\mathcal B$ is **connected**.

  The only inputs are statements of Hintz's Theorem 6.1 and his Definition 4.1
  (both read), plus elementary causal theory and the background facts checked
  in [`tests/test_eh_merger.py`](../tests/test_eh_merger.py).
- **Theorem D₀ (Λ = 0).** The same merger statement for Hintz's Remark 6.2.
  It is *conditional* on two properties (K-esc), (K-trap) of the Kerr-stability
  development of Klainerman–Szeftel. They prove these properties in their own
  coordinates. What is missing is the comparison of those coordinates with
  Kerr's on the initial slice.

**The idea in one paragraph.** Two facts do all the work.

1. Inside a black hole the area radius $r$ is a clock: it decreases along
   every future-directed causal curve.
2. Whether a curve is timelike changes only if the metric changes by a finite
   amount. If a curve is timelike with some margin, it stays timelike under
   any metric perturbation that is small in $C^0$, with no derivatives needed.

Hintz's metrics are $C^0$-close to the big black hole away from the small one,
and $C^0$-close to the small Kerr metric (in rescaled units) near it. So escape
routes can be drawn by hand in the exact backgrounds and survive. Inside the
holes, nothing escapes because $r$ (or the small hole's $\hat r$) runs down.
No light ray is followed precisely. That is why H3's old inputs drop out:
(A2), $C^1$ control of null geodesics, and (A3), strong asymptotic
predictability.

Conventions: [00](00-conventions.md), [07](07-the-merger-problem.md). Hintz and Klainerman–Szeftel also use signature $(-,+,+,+)$, so nothing needs converting.

## 1 · The spacetimes (Hintz, Theorem 6.1)

**Background.** Slowly rotating Kerr–de Sitter $(M_\bullet, a_\bullet, \Lambda > 0)$,
with Cauchy, event and cosmological horizon radii $r_{\rm CH} < r_H < r_C$. It
lives on $M = \mathbb R_{t_*} \times X$ with $X = \{r_- < r < r_+\}$, where
$r_{\rm CH} < r_- < r_H$ and $r_+ > r_C$. The time function is
$t_* = t_{\rm BL} - T(r)$, with $\varphi_* = \varphi - \Phi(r)$, and its level
sets cross both horizons ([HV18, (3.9)], as cited by Hintz; *not re-checked*).

**The plunge.** $C$ is a future-inextendible timelike geodesic through
$t_* = 0$ with $r(C(0)) \in (r_H, r_C)$. It crosses $r = r_H$ in finite time.
By Hintz's proof of 6.1 it then reaches $r = r_-$, the inner boundary of the
domain, before some time $t_* = T$. On $C$ we use proper time $s$, and
**Fermi coordinates** $(s, x)$ on a tube $\{|x| < \rho_0\}$ around $C$. Let
$F(s) = t_*(C(s))$; this is increasing. For Λ = 0 (§6) the same $t_*$ is
taken equal to $t_{\rm BL}$ far out, so its level sets go to spacelike
infinity.

Two facts about $C$ are used below.

- **$C$ never reaches the cosmological horizon.** Beyond $r_C$, $r$ increases
  along every future causal curve (Lemma 2), and $C$ must still reach
  $r_H < r_C$. So $\sup r\circ C < r_C$; set
  $c_0 := r_C - \sup r\circ C > 0$.
- **After the crossing $C$ only goes deeper.** Once $r(C) < r_H$, $r\circ C$
  decreases.

**The small hole.** Subextremal Kerr parameters $b = (m, \mathbf a)$, with
$|\mathbf a| < m$ and $\hat r_\pm = m \pm \sqrt{m^2 - |\mathbf a|^2}$, so
$\hat r_- < m < \hat r_+$. In rescaled coordinates $\hat s = s/\epsilon$,
$\hat x = x/\epsilon$, Hintz's model $\hat g_b$ is Kerr on $\{\hat r = |\hat x| \ge m\}$.
Here $\hat r$ is the Boyer–Lindquist radius and $\hat s = \hat t_{\rm BL} - T(\hat r)$,
with angle $\phi_* = \phi - \Phi(\hat r)$; $d\hat s$ is timelike everywhere
(Hintz (2.19); **read**). The spacetime is
$M_\epsilon = (M \cap \{t_* \ge 0\}) \setminus \{|x| \le \epsilon m\}$ with metric
$g_\epsilon$. The boundary of the excised ball lies between the small hole's
two horizons.

**What we use from Hintz.** Each item below is our reading of the statement
cited.

- **(G1) Far field.** For $t_* \le T + 1$, $g_\epsilon \to g$ in $C^\infty$ on
  every compact set disjoint from $C$ (Theorem 1.1; Definition 4.1(2)).
- **(G2) Near field.** In Fermi coordinates,
  $(g_\epsilon)_{\mu\nu}(s, x) = (\hat g_b)_{\mu\nu}(x/\epsilon) + e_{\epsilon,\mu\nu}(s, x)$,
  where $e_\epsilon = (g - \eta) + h_\epsilon$ (eq. (1.2); by scaling,
  $\hat g_{\epsilon b}(x) = \hat g_b(x/\epsilon)$ componentwise). Here
  $|g - \eta| \le C|x|^2$, and $h_\epsilon \to 0$ together with its
  $\partial_s$ and $(\epsilon^2 + |x|^2)^{1/2}\partial_x$ derivatives. We read
  the convergence as uniform on the tube for $s$ in compact sets; Definition
  4.1(3) gives this. Since $\partial_{\hat x} = \epsilon\partial_x$ and
  $\partial_{\hat s} = \epsilon\partial_s$, two things follow:
  - **(G2a)** on $\{|x| \le \rho\}$, $\limsup_{\epsilon\to0}\sup|e_\epsilon| \le C\rho^2$;
  - **(G2b)** for each $\hat R$, $e_\epsilon \to 0$ in $C^1$, in the variables
    $(\hat s, \hat x)$, on $\{|\hat x| \le \hat R\}$.
- **(G3) Late times.** For $t_* \ge T + 1$, $g_\epsilon = g_{b_{\bullet,\epsilon}} + \tilde g_\epsilon$
  with $|\partial^i_{t_*}\partial^\beta_x \tilde g_\epsilon| \le C\epsilon e^{-\alpha t_*}$
  and $b_{\bullet,\epsilon} = b_\bullet + O(\epsilon)$ (eq. (6.2)). We read
  this as uniform over $X$.

## 2 · A black-hole region without null infinity

Future null infinity of Kerr–de Sitter lies beyond $r_+$, outside Hintz's
spacetime. What the spacetime does contain is the expanding region beyond the
cosmological horizon, so we measure escape against that region.

**Definition.** For $0 < c < r_+ - r_C$ let
$\mathcal C_c = \{p \in M_\epsilon : r(p) \ge r_C + c\}$ and

```math
\mathcal B = M_\epsilon \setminus J^-(\mathcal C_c),
```

with $J^-$ computed in $M_\epsilon$: the events from which no causal signal
reaches the cosmological region.

**Lemma 1.**

- (i) $\mathcal B$ is a future set: if $p \in \mathcal B$ then
  $J^+(p) \subset \mathcal B$.
- (ii) For exact Kerr–de Sitter, $\mathcal B = \{r \le r_H\}$, which is the
  usual black-hole region, because $J^-(\mathscr I^+) = J^-(\mathcal C_c)$
  there.
- (iii) For small $\epsilon$, $\mathcal B$ does not depend on
  $c \in [c_1, c_2] \subset (0, r_+ - r_C)$.

*Proof.* (i) If $q \in J^+(p)$ reaches $\mathcal C_c$, so does $p$.

(ii) Inside, $r$ decreases (Lemma 2(a)). On $r = r_H$, the horizon is an
achronal null hypersurface. Outside, $W_{\rm out}$ carries every point to
$\mathcal C_c$, and beyond $r_C$ the radius increases up to $\mathscr I^+$.

(iii) $C$ stays at $r \le r_C - c_0$, so the region $r \ge r_C + c_1$ stays
away from it. There $W_{\rm out}$ of Lemma 2 is $g_\epsilon$-timelike, by
(G1), (G3) and Lemma 3, and carries every point to $r \ge r_C + c_2$. ∎

(Extending $g_\epsilon$ beyond $r_+$ to $\mathscr I^+$ would need a stability
theorem for the expanding region, e.g. Fournodavlos–Schlue; *not used*.)

## 3 · Background facts

**Lemma 2 (Kerr–de Sitter, any $\Lambda \ge 0$).**

- **(a) $r$ is a clock wherever $\Delta_r < 0$.** There
  $g^{rr} = \Delta_r/\rho^2 < 0$, so $dr$ is timelike. In the black-hole
  interior $r_{\rm CH} < r < r_H$, $r$ decreases along future causal curves;
  beyond $r_C$ it increases. *(test: `g^{rr}`, signs of $\Delta_r$)*
- **(b) Steering fields.** There are $\partial_{t_*}$-invariant, future timelike
  vector fields:
  - $W_{\rm in}$ on $\{r_- \le r \le r_C - c'\}$ with $dr(W_{\rm in}) \le -1$;
  - $W_{\rm out}$ on $\{r_H + c' \le r \le r_+\}$ with $dr(W_{\rm out}) \ge 1$.

  *Proof.* At each point a future timelike $W$ with $dr(W) < 0$ exists:
  - if $\nabla r$ is future causal ($r \le r_H$), any $W$ works;
  - if $\nabla r = \beta n + v$ is spacelike, with $|v| > |\beta|$, take
    $W = n - \kappa\, v/|v|$ with $\kappa < 1$ close to 1.

  Such $W$ form a convex cone, so a $\partial_{t_*}$-invariant partition of
  unity over the compact $r$-range glues them; then rescale. $W_{\rm out}$
  is built the same way. ∎
- **(c) Trapped spheres.** If $r_{\rm CH} < r_0 < r_H$, every sphere
  $\{t_* = \text{const}, r = r_0\}$ is a closed trapped surface. With
  $t_* = t_{\rm BL} - T(r)$, the sphere is a Boyer–Lindquist sphere rotated
  in $\varphi$. Its area density squared is $E/(\Xi^2\Delta_\chi)$, with
  $\chi = \cos\theta$ and

  ```math
  E = \Delta_\chi (r^2 + a^2)^2 - \Delta_r a^2 (1 - \chi^2),
  ```

  ```math
  \partial_r E = 4r^3 + 2a^2 r (1 + \chi^2) + 2Ma^2(1 - \chi^2)
     + \tfrac{4}{3}\Lambda r a^2\chi^2(r^2 + a^2)
     + \Lambda a^2(1 - \chi^2)\big(\tfrac43 r^3 + \tfrac23 a^2 r\big) > 0 .
  ```

  This holds for all $r > 0$, every spin, and every $\Lambda \ge 0$. A normal
  vector $V$ changes the area at rate $V^r\,\partial_r\log\sqrt{E}$. Both
  future null normals have $V^r = dr(V) < 0$ by (a), so $\theta_\pm < 0$.
  *(test: the identity, exactly)*
- **(d) The same for the small Kerr hole ($\Lambda = 0$)** with $\hat r$, on
  $\{m \le \hat r\}$:
  - $\hat r$ decreases along future causal curves in $\{m \le \hat r < \hat r_+\}$;
  - $W^{\rm small}_{\rm in}$ (with $d\hat r \le -1$) exists on
    $\{m \le \hat r \le \hat R\}$, and $W^{\rm small}_{\rm out}$ (with
    $d\hat r \ge 1$) on $\{\hat r_+ + \eta \le \hat r \le \hat R\}$;
  - the spheres $\{\hat s = \text{const}, \hat r = \hat r_1\}$ with
    $m < \hat r_1 < \hat r_+$ are trapped. This uses Hintz's
    $\hat s = \hat t_{\rm BL} - T(\hat r)$ and $\phi_* = \phi - \Phi(\hat r)$.

Two things make (a)–(d) uniform. The backgrounds are stationary, so every
constant is uniform in $t_*$ (respectively $\hat s$). And $\hat g_b \to \eta$
as $|\hat x| \to \infty$, so on all of $\{|\hat x| \ge m\}$ the margin
$\hat g_b^{-1}(d\hat s, d\hat s) \le -\kappa_b < 0$ is uniform.

**Lemma 3 ($C^0$ stability).** Fix a compact set $K$ and a Riemannian
reference metric $h$. Suppose $g'$ is close to $g$ on $K$.

- If $g(W, W) \le -\kappa|W|_h^2$ on $K$ and $|g' - g|_{C^0(K)}$ is small
  enough (depending on $\kappa$), then $W$ is $g'$-timelike, with the same
  time orientation.
- If $d\tau$ is $g$-timelike with a margin, the same holds for $g'$.
- A compact surface with $\theta_\pm < -\kappa$ for $g$ stays trapped for $g'$
  if $|g' - g|_{C^1}$ is small. Expansions involve one derivative of the
  metric, and a constant rescaling $g \to \epsilon^{-2}g$ only rescales
  them. ∎

## 4 · Adapted slices

**Lemma 4.** Fix $\rho_0$ small and $K_0$ large, and set
$\delta_1 = \rho_0 e^{-K_0}$. There is a function $\tau$ on $M_\epsilon$ with:

- $\tau = F(s)$ on $\{|x| \le \delta_1\}$, so the slices are Fermi slices,
  i.e. $\hat s = $ const near the small hole;
- $\tau = t_*$ on $\{|x| \ge \rho_0\}$;
- $d\tau$ is $g_\epsilon$-timelike everywhere, for all $\epsilon < \epsilon_0$.

Write $\Sigma_t = \{\tau = t\}$.

*Proof.* Put $\tau = t_* + \chi\,(F(s) - t_*)$, where the cut-off $\chi$
depends on $\log|x|$, equals 1 for $|x| \le \delta_1$ and 0 for
$|x| \ge \rho_0$, and has $|d\chi/d\log|x|| \le 2/K_0$. Then

```math
d\tau = (1 - \chi)\,dt_* + \chi F'(s)\,ds + (F(s) - t_*)\,\frac{d\chi}{d\log|x|}\,\frac{d|x|}{|x|}.
```

The first two terms form a convex combination of two future timelike
covectors, so they are timelike with a margin $\mu$. Since $F(s) = t_*$ on
$C$, we have $|F(s) - t_*| \le A|x|$, so the last term has size at most
$2A/K_0 < \mu$ for $K_0$ large. Hence $d\tau$ is $g$-timelike.

- On $\{|x| \ge \delta_1\}$, which stays away from $C$, use (G1), (G3) and
  Lemma 3.
- On $\{|x| \le \delta_1\}$, $d\tau \propto ds \propto d\hat s$. By (G2a),
  $g_\epsilon = \hat g_b(\hat x) + O(\delta_1^2) + o(1)$. Also $\hat g_b$ is
  bounded with bounded inverse and has $\hat g_b^{-1}(d\hat s, d\hat s) \le -\kappa_b$
  on $\{|\hat x| \ge m\}$. ∎

A cut-off with $|d\chi/d|x|| \sim 1/\rho_0$ need not work. Fermi slices and
$t_*$-slices are tilted by the boost of $C$, so $F(s) - t_*$ is of first order
in $|x|$, and the product is of order 1, not small. The logarithmic cut-off
makes it small.

## 5 · Theorem D and its proof

**Constants.** $L$ is a Lipschitz constant on the tube: it bounds
$|r(s, x) - r(C(s))| / |x|$ and $|dr(C(s))/ds|$. It also bounds
$|r - r(C(s(t)))|/\rho$ on $\Sigma_t \cap \{|x| \le \rho\}$, $\rho \le \rho_0$.
Here $s(t) := F^{-1}(t)$ is the Fermi time of the slice $\Sigma_t$ near $C$.
Choose them in this order:

1. $c \in (0, \min(c_0, r_+ - r_C)/2)$ with $r_H + 2c < r(C(0))$ and
   $r_- + 4c < r_H - 4c$, and $\eta \in (0, (\hat r_+ - m)/2)$;
2. $\rho_0 \le c/(4L)$ and $\delta_1$ (Lemma 4);
3. the tube radius $\delta' = \min(\delta_1, c/(8L))$;
4. $\hat R$ with $|\hat g_b - \eta| \le \kappa$ for $|\hat x| \ge \hat R$;
5. finally $\epsilon_0$.

Let $N := \{|x| \le \delta'\}$, the closed tube in which $g_\epsilon$ is *not*
close to $g$. Let $P := \{m < \hat r \le \hat r_+ - \eta\}$, the small hole's core.

**Theorem D.** For $0 < \epsilon < \epsilon_0$, with $\mathcal B$ as in §2 and
$\Sigma_t$ from Lemma 4:

1. **The small hole is a black hole.** $P \subset \mathcal B$. Whenever $C$
   is inside the domain, $\Sigma_t \cap P$ contains a closed trapped surface
   $T_1(t) = \{s = s(t), \hat r = \hat r_1\}$, with
   $m < \hat r_1 < \hat r_+ - \eta$.
2. **The large hole is a black hole.** $\{r \le r_H - c\} \setminus N \subset \mathcal B$.
   For every $t$, $\Sigma_t$ contains a closed trapped sphere $T_2(t)$ in this
   set.
3. **The horizons are localised.**
   - $\mathcal B \cap \{r \ge r_H + c\} \setminus N = \emptyset$.
   - At Fermi times when $r(C(s)) \ge r_H + 2c$,
     $\mathcal B \cap N \subset \{\hat r < \hat r_+ + \eta\}$.

   So, up to the margins $c$ and $\eta$ (which are arbitrary), the event
   horizons sit at the Kerr–de Sitter horizon $r_H$ and, in the small hole's
   own units, at its Kerr horizon $\hat r_+$.
4. **Merger.**
   - If $r(C(s(t))) \ge r_H + 2c$, then $T_1(t)$ and $T_2(t)$ lie in different
     connected components of $\mathcal B \cap \Sigma_t$.
   - If $r_- + 4c \le r(C(s(t))) \le r_H - 3c$, they lie in the same
     component.
5. **$\mathcal B$ is path-connected.**

By 1, 2 and 4, $g_\epsilon$ describes an event-horizon merger in the sense of
[07](07-the-merger-problem.md): two black holes, initially distinct, end in
one component. The definition there is phrased with $\mathscr I^+$ and Cauchy
slices; here $\mathscr I^+$ is replaced by $\mathcal C_c$ (§2), and the slices
are the level sets of $\tau$.

### Proof

**Step 1: the small hole's core lies in $\mathcal B$.** By (G2b) and Lemma 3,
$\hat r$ is a $g_\epsilon$-time function on $\{m \le \hat r \le \hat r_+ - \eta/2\}$
and decreases to the future. A future causal curve from $P$ therefore stays
in $P$ until it leaves $M_\epsilon$ through the excision or the inner
boundary. All the while it lies within $\epsilon\hat r_+$ of $C$. There
$r \le r_C - c_0 + L\delta' < r_C + c$, so it never reaches $\mathcal C_c$.

The surfaces $T_1(t)$ are trapped for $\hat g_b$ (Lemma 2(d)). By (G2b) they
are trapped for $\epsilon^{-2}g_\epsilon$, hence for $g_\epsilon$.

**Step 2: a gain lemma for the large hole.** Let $\gamma$ be a future causal
curve from $p$, where either

- (i) $p \notin N$ and $r(p) \le r_H - c$, or
- (ii) $p \in N$, at a Fermi time $s_p$ with $r(C(s_p)) \le r_H - c - L\delta'$.

*Claim:* $r \le r_H - c + 2L\delta' \le r_H - 3c/4$ along $\gamma$.

Outside $N$, $r$ decreases. This holds as long as $r \le r_H - c/2$, by
Lemma 2(a), Lemma 3, (G1) and (G3). On $N$, the Fermi time is
$s = F^{-1}(\tau)$ (Lemma 4, since $\delta' \le \delta_1$), and $\tau$
increases along $\gamma$. So each visit to $N$ happens at a later Fermi time
than the one before.

Suppose $\gamma$ first enters $N$ at Fermi time $s_e$. Then
$r(C(s_e)) \le r(\gamma) + L\delta' < r_H$. So $C$ is already inside, and
$r\circ C$ decreases from then on. During that visit, and any later one,

```math
r(\gamma) \le r(C(s)) + L\delta' \le r(C(s_e)) + L\delta' \le r(p) + 2L\delta' .
```

Between visits $r$ decreases. So $\gamma$ never reaches $r_H - c/2$, let
alone $\mathcal C_c$. The same argument gives case (ii) directly.

For $T_2(t)$, take $\Sigma_t \cap \{r = r_2\}$, with $r_2 \in \{r_H - 2c, r_H - 3c\}$
chosen at distance $\ge c/2$ from $r(C)$. Then it avoids
$\{|x| \le \rho_0\}$, so it is the sphere $\{t_* = t, r = r_2\}$. That
sphere is trapped by Lemma 2(c) and (G1)/(G3).

**Step 3: escape.**

- **(a) Compact part.** Consider the compact set
  $Q = \{0 \le t_* \le T + 1,\ r_H + c \le r \le r_C + c\} \setminus \{|x| < \delta'\}$.
  Each $p \in Q$ has a $g$-timelike curve to $\mathcal C_{2c}$ that avoids
  $C$: an integral curve of $W_{\rm out}$, nudged off $C$ if needed, which is
  possible because a curve can always be moved off another curve in four
  dimensions.

  Timelike connections to an interior point form an open condition. So
  finitely many such curves, together with short connecting geodesics, serve
  all of $Q$, and their union is a compact set disjoint from $C$, with a
  common timelike margin. By (G1), (G3) and Lemma 3 these curves are
  $g_\epsilon$-timelike, so $Q \subset J^-(\mathcal C_c)$.
- **(b) Late part.** On $\{t_* \ge T + 1,\ r \ge r_H + c\}$, $W_{\rm out}$ is
  $g_\epsilon$-timelike with $dr \ge 1/2$, by (G3).
- **(c) Near the small hole while it is outside.** Take $q \in N$ with
  $\hat r(q) \ge \hat r_+ + \eta$ at a time when $r(C) \ge r_H + 2c$. From
  $q$:
  - if $\hat r(q) < \hat R$, follow $W^{\rm small}_{\rm out}$ out to
    $\hat r = \hat R$ (by (G2b));
  - then move radially away from $C$ at Fermi speed $1/2$, up to $|x| = 2\delta'$.
    This path is timelike, because on $\{\hat R \le |\hat x|,\ |x| \le 2\delta'\}$
    the metric is within $\kappa + 4C\delta'^2 + o(1)$ of Minkowski (G2a).
  - It arrives outside $N$. The trip takes Fermi time at most
    $4\delta' + O(\epsilon\hat R)$, so the arrival point has
    $r \ge r_H + 2c - 7L\delta' \ge r_H + c$. It is therefore a point of $Q$,
    or of the late part, and escapes by (a) or (b).

This proves 3.

**Step 4: separation before the plunge.** Let $r(C(s(t))) \ge r_H + 2c$. The
coordinate sphere $S = \{s = s(t), |x| = \delta'/2\}$ lies in $\Sigma_t$
(since $\delta' \le \delta_1$). It has $\hat r = \delta'/(2\epsilon) \ge \hat r_+ + \eta$,
so $S \subset J^-(\mathcal C_c)$ by Step 3(c).

$S$ splits $\Sigma_t$ into an inside, which contains $T_1(t)$, and an
outside, which contains $T_2(t)$ (since $T_2$ avoids $|x| \le \rho_0$). A
connected subset of $\mathcal B \cap \Sigma_t$ cannot cross $S$.

**Step 5: joined after the plunge.** Let
$r_- + 4c \le r(C(s(t))) \le r_H - 3c$ and
$Y = \Sigma_t \cap \{r \le r_H - c\} \cap M_\epsilon$.

- **$Y \subset \mathcal B$**, by Step 2: case (i) off $N$, case (ii) on $N$.
- **$Y$ is connected.** The part of $\Sigma_t$ with $|x| \le \rho_0$ has
  $|r - r(C)| \le L\rho_0 \le c/4$. So it lies inside $\{r_- < r \le r_H - c\}$.
  It is a graph over the Fermi disc, minus the excision, which is connected.
  Outside it, $\Sigma_t$ is the $t_*$-slice, where $\{r_- < r \le r_H - c\}$
  is a spherical shell. The two pieces meet along $|x| = \rho_0$. The shell
  minus a small interior ball, glued to the disc minus a smaller ball, is
  connected.
- **$T_1(t), T_2(t) \subset Y$.**

**Step 6: $\mathcal B$ is path-connected.** Let

```math
G = \big(\{r \le r_H - c\} \setminus N\big) \;\cup\; \big(N \cap \{s : r(C(s)) \le r_H - 3c\}\big) \;\cap\; M_\epsilon .
```

By Step 2, $G \subset \mathcal B$. $G$ is the connected four-manifold
$\{r_- < r \le r_H - c\}$ minus closed tubular neighbourhoods of two arcs:
the excised worldtube, and a piece of $\bar N$. An arc has codimension 3, so
removing such neighbourhoods leaves $G$ path-connected.

Now take any $q \in \mathcal B$. Every curve used below is future causal, so
by Lemma 1(i) it stays in $\mathcal B$.

- **$q \notin N$.** By Step 3, $r(q) < r_H + c$. Follow $W_{\rm in}$, which
  is $g_\epsilon$-timelike off $N$, with $dr \le -1/2$. Either the curve
  reaches $r = r_H - c$ off $N$, which puts it in $G$, or it first meets
  $N$ (on $|x| = \delta'$), which is one of the next two cases. It cannot
  reach $r_-$ first.
- **$q \in N$ at a time with $r(C) \le r_H - 3c$.** Then $q \in G$.
- **$q \in N$ otherwise.** Here the tube lies inside the domain. Move
  radially *toward* $C$ at Fermi speed $1/2$ to $\hat r = \hat R$, then
  follow $W^{\rm small}_{\rm in}$ into $P$.
  - $P$ over Fermi times up to the first time with $r(C) \le r_H - 3c$ is an
    interval times a shell, so it is connected, and its end lies in $G$.
  - Points of $P$ after that time lie in $G$ already. ∎

### What made this possible

- No null geodesic is integrated and no curvature decay is used. Every
  escape route is a *timelike curve with margin*, and every trap is a *time
  function*. Both survive $C^0$ perturbations.
- Derivatives of $g_\epsilon$ enter only through the trapped surfaces.
- Strong asymptotic predictability ((A3) of note 09) is never invoked:
  membership in $\mathcal B$ is shown directly, not inferred from trapped
  surfaces.
- The price is the margins $c$, $\eta$. The theorem says nothing about the
  horizon *at* the scale where its exact shape matters, such as its
  creases, the moment of merger, or the toroidal phases seen by
  Emparan–Martínez–Zilhão.

## 6 · Λ = 0: Theorem D₀, conditional

Hintz's Remark 6.2 replaces Hintz–Vasy by the Kerr stability theorem of
Klainerman–Szeftel (KS23), with GKS22 and She23.

- Glue up to a time $T_1 \ge T + 1$.
- Take the induced data on $\Sigma_{T_1} = \{t_* = T_1\}$. These are
  asymptotically flat, and compactly supported perturbations of Kerr data by
  finite speed of propagation.
- Evolve them with the stability theorem.

A remark is not a theorem, so we state what we need as hypotheses.

- **(Ex)** For small $\epsilon$, these data have a future development
  $\mathcal M^+_\epsilon$ with complete $\mathscr I^+$, and
  $M_\epsilon := (\text{Hintz region} \cap \{t_* \le T_1\}) \cup \mathcal M^+_\epsilon$
  is a spacetime.
- **(K-esc)** There is $R_0$, independent of $\epsilon$ and $T_1$, such that
  every point of $\Sigma_{T_1}$ with $r \ge R_0$ lies in $J^-(\mathscr I^+)$.
- **(K-trap)** There is $c_K \in (0, (r_H - r_-)/16)$ such that, for small
  $\epsilon$, no point of $\Sigma_{T_1}$ with $r \le r_H - c_K$ lies in
  $J^-(\mathscr I^+)$.

**Theorem D₀.** Assume (Ex), (K-esc), (K-trap). Let $c = 2c_K$ and
$\delta' \le c_K/(4L)$. Let $\Theta$ be the largest $t_*$-time a $W_{\rm out}$
curve needs to go from $\{r_H + c \le r \le R_0 + 1\}$ to $r = R_0 + 1$, and
take $T_1 \ge T + 1 + \Theta$. Then, for small $\epsilon$, with
$\mathcal B = M_\epsilon \setminus J^-(\mathscr I^+)$, the following
conclusions of Theorem D hold:

- 1;
- 2 for $t_* \le T_1$;
- 3 for $t_* \le T + 1$;
- 4.

How each step changes:

- Conclusion 1 needs only (Ex). The curves from $P$ stay in the tube, which
  leaves the domain before $t_* = T$.
- Step 3 ends at $\Sigma_{T_1}$ and uses (K-esc). For $r \ge R_0 > 2M_\bullet$,
  the stationary Killing field of Kerr is timelike. It keeps $r$ fixed until
  $t_* = T_1$. Far out the metric is exactly Kerr, by finite speed of
  propagation.
- Step 2 ends at $\Sigma_{T_1}$, at $r \le r_H - 2c_K + 2L\delta' < r_H - c_K$,
  and uses (K-trap).

Steps 4 and 5 are unchanged. We do not claim connectedness (Step 6) for
Λ = 0.

**How close KS come to (K-esc) and (K-trap)** (KS23, arXiv:2104.11857,
**read** §3.4, §3.8). Their admissible future complete spacetime is
$\mathcal M = {}^{(ext)}\mathcal M \cup {}^{(int)}\mathcal M$, split by
$\{r = r_0\}$ (Definition 3.4.5).

- **Escape side.** Every point of ${}^{(ext)}\mathcal M$ lies on an outgoing
  null geodesic with $e_4(r) = 1$ that reaches $\mathscr I^+$. KS state this
  for $\mathcal T = \{r = r_0\}$ (§3.8, "Existence of a future event horizon").
- **Trap side.** On ${}^{(int)}\mathcal M \cap \{r \le r_+(1 - \delta_H/2)\}$,
  $e_3(r) \le -1/2$ and $e_4(r) \lesssim -\delta_H$. So $r$ decreases along
  future null geodesics, "so that, in particular, γ(s) cannot reach
  $\mathscr I^+$" (same section).

  KS state this for null geodesics. For a causal vector
  $\lambda e_4 + \mu e_3 + f^b e_b$ with $\lambda, \mu \ge 0$, $r$ still
  decreases, because $e_b(r) = 0$. So it extends to all causal curves.

Here $r$, $e_3$, $e_4$ are KS's own, defined by their PG structures. The
missing link is **anchoring**: on the initial layer, KS's $r$ must be
uniformly close to the Kerr $r$ of our $\Sigma_{T_1}$, so that points with
$r \ge R_0$ and points with $r \le r_H - c_K$ lie in, or flow into, those two
regions. Two further caveats:

- KS's theorem starts from an admissible *initial data layer*. That Cauchy
  data generate one is cited from other work (KS Remark 3.4.3; refs. [36],
  [37], [8] there). We did not read it, nor check that Hintz's data at
  $T_1$ qualify.
- Their $\delta_H$ is a fixed small constant. So (K-trap) is expected with
  $c_K \approx r_H\delta_H/2$, which is enough here.

## 7 · What is and is not proved

| Statement | Status |
| --- | --- |
| Lemma 2(a): $g^{rr} = \Delta_r/\rho^2$ | **checked** (`tests/test_eh_merger.py`, exact) |
| Lemma 2(c): area identity, $\partial_r E > 0$ for all $r > 0$, all spins, Λ ≥ 0 | **checked** (exact identity); positivity read off term by term |
| Kerr–de Sitter solves $\mathrm{Ric} = \Lambda g$; reduces to Kerr and SdS | **checked** (rational parameters exactly, slow test; SdS symbolic) |
| Lemmas 1, 2(b), 2(d), 3, 4 | written out here (elementary) |
| (G1)–(G3) | Hintz Theorem 1.1, Definition 4.1, (1.2), (6.2): **read**. Uniformity on the tube and over $X$ is *our reading* |
| $t_* = t_{\rm BL} - T(r)$ form of Hintz–Vasy's time function | cited (HV18 (3.9)), *not re-checked* |
| **Theorem D** (Λ > 0): both holes in $\mathcal B$, localisation, merger, $\mathcal B$ connected | **proved here, given (G1)–(G3)**; potentially new (Hintz lists it as open) |
| $\mathcal B$ defined with $\mathscr I^+$ instead of $\mathcal C_c$, for Λ > 0 | not proved: needs the spacetime beyond $r_+$ |
| **Theorem D₀** (Λ = 0) | **conditional** on (Ex), (K-esc), (K-trap) |
| (K-esc), (K-trap) | proved by KS in their coordinates; anchoring to Kerr's $r$ *not checked* |
| Exact location, smoothness and merger time of the horizon | not addressed (only to within $c$, $\eta$) |

**Novelty.** Hintz states the structure of the black-hole region ("conjecturally
connected") as future work. A bounded search (2026-10-29) found no later
treatment. Theorem D is therefore **potentially new**. The argument is
elementary given Hintz's estimates, and experts may regard it as routine; its
content is that $C^0$ control suffices. For Λ = 0 we offer a reduction, not a
theorem.

## 8 · Physical picture

Before the plunge, the small hole is a tiny Kerr black hole moving through
nearly empty space near the big one.

- **Inside its horizon,** the small hole's own radius $\hat r$ only runs
  down, so nothing gets out.
- **A little outside its horizon,** a rocket can climb out of the small
  hole's well. There space looks almost flat, so it can coast away, and then
  outrun the big hole into the cosmological region. Its worldline is
  *timelike with room to spare*, which is why small metric errors cannot
  block it.

So a shell of escaping events surrounds the small hole and separates it from
the big hole's interior: two black holes.

After the small hole has fallen well inside the big hole, every event near
it has big-hole $r$ below $r_H$. The big hole's $r$ is a clock that runs
down, and the small hole only ever moves deeper. Its neighbourhood can push
$r$ up by at most the tube's width. So everything near the small hole is
trapped together with the big hole's interior: one black hole.

The merger of event horizons is forced by the big hole being a trap at late
times. This is the teleological nature of event horizons, made into an
estimate.

## Discussion prompts

- Lemma 3 says only $C^0$ closeness of the metric matters for "can a signal
  get out". Why, then, do the full stability theorems need so much more
  (decay of derivatives, gauge choices)? Which of their conclusions are
  causal, and which are not?
- Step 6 shows that $\mathcal B$ is connected. Does this rule out a
  transient toroidal horizon cross-section in some slicing? (Compare
  [07](07-the-merger-problem.md) and Emparan–Martínez–Zilhão.)
- The margins $c$ and $\eta$ are arbitrary but fixed before $\epsilon \to 0$.
  What would it take to locate the merger *event*, the first point where the
  two horizon components touch, as $\epsilon \to 0$?
- For Λ = 0 everything hinges on anchoring KS's $r$ to Kerr's $r$ on the
  initial slice. Is that a statement about the gauge, or about physics?
