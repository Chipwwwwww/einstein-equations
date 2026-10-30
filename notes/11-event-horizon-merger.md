# 11 · Event horizons of extreme-mass-ratio mergers

Hintz (2024, *Gluing small black holes along timelike geodesics III*,
arXiv:2408.06715) constructs vacuum spacetimes in which a light Kerr black
hole, of mass $\epsilon m$, falls into a unit-mass Kerr–de Sitter black hole.
Afterwards the spacetime settles down to a single black hole.

- For Λ > 0 this is his Theorem 6.1.
- For Λ = 0 he sketches it in his Remark 6.2. Hintz's 2026 stability theorem
  for the whole subextremal Kerr family (*Nonlinear stability of subextremal
  Kerr black holes*, arXiv:2606.28253, consequence (3) of §1) now supplies
  the missing late-time step.

He leaves "the construction of the event horizon, the structure of the black
hole region—conjecturally connected—" to future work (III, §1). Hypothesis H3
of [09](09-merger-hypotheses.md) asked whether the event horizons in these
spacetimes merge in the sense of [07](07-the-merger-problem.md).

- **Theorem D (Λ ≥ 0).** Take $\mathcal B$ to be the events that cannot send
  a signal to the *escape region*. That region is the expanding region
  beyond the cosmological horizon if Λ > 0, and arbitrarily large distances
  if Λ = 0. Then, for all small $\epsilon$:
  - $\mathcal B$ contains the interiors of both holes, each with closed
    trapped surfaces;
  - both horizons are located to within any prescribed margin;
  - before the plunge the two interiors lie in **different** components of
    $\mathcal B\cap\Sigma_t$, and afterwards in the **same** one;
  - $\mathcal B$ is **connected**.

  The inputs are the statements of Hintz's Theorem 6.1 and Definition 4.1
  (gluing), and, for Λ = 0, of Theorem 1.1/13.1 of the 2026 paper (late
  times). We read all of them; where we rely on our reading of their
  uniformity, §7 says so. The rest is elementary causal theory and background
  facts checked in [`tests/test_eh_merger.py`](../tests/test_eh_merger.py).
- **Theorem E (Λ = 0, null infinity).** With the textbook black-hole region
  $M \setminus J^-(\mathscr I^+)$, the merger holds in the exact sense of
  note 07: two components before the plunge, one after. This needs one more
  input: points far out on the initial slice of the stability development
  can signal to $\mathscr I^+$. Lemma 6 derives it from the 2026 paper's
  description of null infinity. That derivation is our reading of the
  paper, not an independent verification.
- **Corollary E1 (session 7).** The textbook region and the escape-defined
  region coincide, $\mathcal B_{\mathscr I} = \mathcal B$, so all of
  Theorem D holds for the textbook black-hole region: localisation and
  connectedness included. It rests on Lemma 7, the all-late-times version
  of Lemma 6, read from the same paper.

*History.* Session 4 stated the Λ = 0 case conditionally (as "Theorem D₀"),
through Klainerman–Szeftel's stability theorem. That needed their coordinates
to be anchored to Kerr's. Hintz's 2026 theorem works in fixed coordinates, so
anchoring is no longer needed (§6).

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

## 1 · The spacetimes

**Background.**

- **Λ > 0:** slowly rotating Kerr–de Sitter $(M_\bullet, a_\bullet, \Lambda)$,
  with Cauchy, event and cosmological horizon radii $r_{\rm CH} < r_H < r_C$.
  It lives on $M = \mathbb R_{t_*} \times X$ with $X = \{r_- < r < r_+\}$,
  where $r_{\rm CH} < r_- < r_H$ and $r_+ > r_C$ (Hintz III, Thm 6.1).
- **Λ = 0:** subextremal Kerr $(M_\bullet, a_\bullet)$, with $X = \{r > r_-\}$
  and $r_- = M_\bullet$, which lies between $r_{\rm CH}$ and $r_H$. This is
  the inner radius of the domain $\Omega = \{r \ge m_0\}$ in the 2026 paper.

In both cases the time function is $t_* = t_{\rm BL} - T(r)$, with
$\varphi_* = \varphi - \Phi(r)$. Its level sets cross the horizons ([HV18,
(3.9)] as cited by Hintz; for Λ = 0 the function $\tilde t$ of the 2026
paper). For Λ = 0, $T = 0$ for large $r$, so the level sets go to spacelike
infinity.

**The plunge.** $C$ is a future-inextendible timelike geodesic through
$t_* = 0$ with $r(C(0)) > r_H$ (and $< r_C$ if Λ > 0). It crosses $r = r_H$ in finite time.
By Hintz's proof of 6.1 (the same argument works for Λ = 0) it then reaches
$r = r_-$, the inner boundary of the
domain, before some time $t_* = T$. On $C$ we use proper time $s$, and
**Fermi coordinates** $(s, x)$ on a tube $\{|x| < \rho_0\}$ around $C$. Let
$F(s) = t_*(C(s))$; this is increasing.

Two facts about $C$ are used below.

- **$C$ stays below the escape region.** For Λ > 0, $r$ increases beyond
  $r_C$ along every future causal curve (Lemma 2), and $C$ must still reach
  $r_H < r_C$. So $\sup r\circ C < r_C$; set
  $c_0 := r_C - \sup r\circ C > 0$. For Λ = 0, $\sup r\circ C < \infty$,
  since $C$ is inside the domain only for a compact time.
- **After the crossing $C$ only goes deeper.** Once $r(C) < r_H$, $r\circ C$
  decreases.

**The small hole.** Subextremal Kerr parameters $b = (m, \mathbf a)$, with
$|\mathbf a| < m$ and $\hat r_\pm = m \pm \sqrt{m^2 - |\mathbf a|^2}$, so
$\hat r_- < m < \hat r_+$. In rescaled coordinates $\hat s = s/\epsilon$,
$\hat x = x/\epsilon$, Hintz's model $\hat g_b$ is Kerr on $\{\hat r = |\hat x| \ge m\}$.
Here $\hat r$ is the Boyer–Lindquist radius and $\hat s = \hat t_{\rm BL} - T(\hat r)$,
with angle $\phi_* = \phi - \Phi(\hat r)$; $d\hat s$ is timelike everywhere
(Hintz (2.19); **read**). The boundary of the excised ball $\{|x| \le \epsilon m\}$
lies between the small hole's two horizons.

**The spacetime $(M_\epsilon, g_\epsilon)$.**

- **Λ > 0:** Hintz's Theorem 6.1 gives $g_\epsilon$ on
  $M_\epsilon = (M \cap \{t_* \ge 0\}) \setminus \{|x| \le \epsilon m\}$.
  Set $T_{\rm late} := T + 1$.
- **Λ = 0:**
  1. Fix $T_{\rm late} \ge T + 1$ (chosen in §6). Glue with Hintz's
     Theorem 5.12 on $\{0 \le t_* \le T_{\rm late}\}$.
  2. The data on $\Sigma_{T_{\rm late}} = \{t_* = T_{\rm late}\}$ equal the Kerr data
     outside a compact set, by finite speed of propagation (III, Theorem
     5.1(3)). They are $O(\epsilon)$-close to them in $C^k$ (III,
     Definition 4.1(2)), because the small hole has left the domain by then.
  3. Evolve these data with Theorem 1.1/13.1 of the 2026 paper, on
     $\{t_* \ge T_{\rm late}\} \cap \{r \ge r_-\}$.
  4. Both pieces are developments of the same data near
     $\Sigma_{T_{\rm late}}$, so they glue (local uniqueness). This is consequence
     (3) in §1 of the 2026 paper.

**What we use from Hintz.** Each item below is our reading of the statement
cited.

- **(G1) Far field.** For $t_* \le T_{\rm late}$, $g_\epsilon \to g$ in
  $C^\infty$ on every compact set disjoint from $C$ (Theorem 1.1;
  Definition 4.1(2)). For Λ = 0, also $g_\epsilon = g$ outside a compact set
  in this range (Theorem 5.1(3)).
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
- **(G3) Late times.** For $t_* \ge T_{\rm late}$, $g_\epsilon$ is uniformly
  $C^0$-close to the background, with $C^1$-closeness on compact sets. The
  sources differ by Λ.
  - **Λ > 0:** $g_\epsilon = g_{b_{\bullet,\epsilon}} + \tilde g_\epsilon$
    with $|\partial^i_{t_*}\partial^\beta_x \tilde g_\epsilon| \le C\epsilon e^{-\alpha t_*}$
    and $b_{\bullet,\epsilon} = b_\bullet + O(\epsilon)$ (III, eq. (6.2)). We
    read this as uniform over $X$.
  - **Λ = 0:** the 2026 paper (Theorems 1.1, 13.1) gives
    $g_\epsilon = g_{b} + h$ in the fixed coordinates $(\tilde t, x)$, with
    $|b - b_\bullet|$ small. Here $h$ decays in every asymptotic region:
    $\tilde t^{-2-\epsilon_K}$ on compact sets, $r^{-1}$ near null and spatial
    infinity, and $t_*^{-1}\log t_*$ in between. We read the constants as
    small when the data are small. The theorem is a small-data Nash–Moser
    construction, and its exterior part states this smallness explicitly
    (eq. (5.56) there).

## 2 · Black-hole regions defined by escape

For Λ > 0, future null infinity lies beyond $r_+$, outside Hintz's
spacetime. What the spacetime does contain is the expanding region beyond the
cosmological horizon. For Λ = 0 the natural target is "arbitrarily far away".

**Definition (escape region).**

- **Λ > 0:** $\mathcal E = \{r \ge r_C + c_E\}$, for some
  $0 < c_E < r_+ - r_C$.
- **Λ = 0:** $\mathcal E = \{r \ge R\}$, for some
  $R > \max(r_H, \sup r\circ C) + 2$.

Then

```math
\mathcal B = M_\epsilon \setminus J^-(\mathcal E),
```

with $J^-$ computed in $M_\epsilon$: the events from which no causal signal
reaches the escape region.

**Lemma 1.**

- (i) $\mathcal B$ is a future set: if $p \in \mathcal B$ then
  $J^+(p) \subset \mathcal B$.
- (ii) For exact Kerr–de Sitter, $\mathcal B = \{r \le r_H\}$, which is the
  usual black-hole region, because $J^-(\mathscr I^+) = J^-(\mathcal E)$
  there.
- (iii) For small $\epsilon$, $\mathcal B$ does not depend on the
  parameter: $c_E \in [c_1, c_2] \subset (0, r_+ - r_C)$, or $R$ in a
  compact range.
- (iv) For Λ = 0,
  $\mathcal B \subseteq \mathcal B_{\mathscr I} := M_\epsilon \setminus J^-(\mathscr I^+)$.
  A causal curve that reaches $\mathscr I^+$ passes through every
  $\{r \ge R\}$.

*Proof.* (i) If $q \in J^+(p)$ reaches $\mathcal E$, so does $p$.

(ii) Inside, $r$ decreases (Lemma 2(a)). On $r = r_H$, the horizon is an
achronal null hypersurface. Outside, $W_{\rm out}$ carries every point to
$\mathcal E$. For Λ > 0, beyond $r_C$ the radius increases up to
$\mathscr I^+$. For Λ = 0, the outgoing principal null geodesics from
$r > r_H$ reach $\mathscr I^+$.

(iii) $C$ stays below the escape region, so the region $r \ge r_C + c_1$
(respectively $r \ge R_1$) stays away from it. There $W_{\rm out}$ of
Lemma 2 is $g_\epsilon$-timelike, by (G1), (G3) and Lemma 3, and carries
every point further out. ∎

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
  - $W_{\rm in}$ on $\{r_- \le r \le r_C - c'\}$ with $dr(W_{\rm in}) \le -1$
    (for Λ = 0, on $\{r_- \le r \le R'\}$ for any $R'$);
  - $W_{\rm out}$ on $\{r_H + c' \le r \le r_+\}$ with $dr(W_{\rm out}) \ge 1$.

  *Proof.* At each point a future timelike $W$ with $dr(W) < 0$ exists:
  - if $\nabla r$ is future causal ($r \le r_H$), any $W$ works;
  - if $\nabla r = \beta n + v$ is spacelike, with $|v| > |\beta|$, take
    $W = n - \kappa\, v/|v|$ with $\kappa < 1$ close to 1.

  Such $W$ form a convex cone, so a $\partial_{t_*}$-invariant partition of
  unity over the compact $r$-range glues them; then rescale. $W_{\rm out}$
  is built the same way.

  For Λ = 0 the range of $W_{\rm out}$ is not compact. For $r \ge 10M$ we
  take $W_{\rm out} = \partial_t + \tfrac12\partial_r$, which has bounded
  components and $g(W, W) \le -\tfrac14$ for every $|a| \le M$
  *(test: polynomial with negative coefficients)*; then glue as before. ∎
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
  Lemma 3. For Λ = 0 this includes $r \to \infty$, where $d t_* = dt$ has a
  uniform margin.
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

1. $c > 0$ (for Λ > 0 also $c < c_0/2$) with $r_H + 2c < r(C(0))$ and
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
   horizons sit at the background horizon $r_H$ and, in the small hole's
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
slices; here $\mathscr I^+$ is replaced by the escape region $\mathcal E$
(§2; Theorem E in §6 restores $\mathscr I^+$ for Λ = 0), and the slices
are the level sets of $\tau$.

**Corollary D1 (when the horizons merge).** Let $s_{\rm cross}$ be the proper
time at which $C$ crosses the background horizon $r = r_H$. For every
$\eta > 0$ there is $\epsilon_0 > 0$ with the following property. For
$\epsilon < \epsilon_0$, on the adapted slices through $C(s)$:

- the two holes lie in different components of $\mathcal B$ for
  $s \le s_{\rm cross} - \eta$;
- they lie in the same component for $s$ in
  $[s_{\rm cross} + \eta, s_{\rm cross} + \eta_1]$, where $\eta_1 > \eta$ is
  any value at which $C$ is still inside the domain.

On the slow time scale, then, the event horizons merge at the moment the
small hole's worldline crosses the background horizon. The details of the
merger (the Emparan–Martínez horizon) live on the fast scale, in rescaled
time $\hat s = s/\epsilon$. The same holds for $\mathcal B_{\mathscr I}$ in
Theorem E.

*Proof.* $C$ crosses $r = r_H$ transversally and to the future. There
$\nabla r$ is future null and $\dot C$ is future timelike, so
$\tfrac{d}{ds}r(C(s)) = g(\nabla r, \dot C) < 0$. Choose $c$ so small that
$r(C(s)) \ge r_H + 2c$ for $s \le s_{\rm cross} - \eta$ and
$r(C(s)) \le r_H - 3c$ for $s \ge s_{\rm cross} + \eta$. Then apply
conclusion 4 with that $c$. ∎

### Proof

**Step 1: the small hole's core lies in $\mathcal B$.** By (G2b) and Lemma 3,
$\hat r$ is a $g_\epsilon$-time function on $\{m \le \hat r \le \hat r_+ - \eta/2\}$
and decreases to the future. A future causal curve from $P$ therefore stays
in $P$ until it leaves $M_\epsilon$ through the excision or the inner
boundary. All the while it lies within $\epsilon\hat r_+$ of $C$. There
$r \le \sup r\circ C + L\delta'$, which is below the escape region, so it
never reaches $\mathcal E$.

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
alone $\mathcal E$. The same argument gives case (ii) directly.

For $T_2(t)$, take $\Sigma_t \cap \{r = r_2\}$, with $r_2 \in \{r_H - 2c, r_H - 3c\}$
chosen at distance $\ge c/2$ from $r(C)$. Then it avoids
$\{|x| \le \rho_0\}$, so it is the sphere $\{t_* = t, r = r_2\}$. That
sphere is trapped by Lemma 2(c) and (G1)/(G3).

**Step 3: escape.**

- **(a) Compact part.** Consider the compact set
  $Q = \{0 \le t_* \le T_{\rm late},\ r_H + c \le r \le r_{\mathcal E}\} \setminus \{|x| < \delta'\}$,
  where $r_{\mathcal E}$ is the inner radius of $\mathcal E$. Each $p \in Q$
  has a $g$-timelike curve into the interior of $\mathcal E$ that avoids
  $C$: an integral curve of $W_{\rm out}$, nudged off $C$ if needed, which is
  possible because a curve can always be moved off another curve in four
  dimensions.

  Timelike connections to an interior point form an open condition. So
  finitely many such curves, together with short connecting geodesics, serve
  all of $Q$, and their union is a compact set disjoint from $C$, with a
  common timelike margin. By (G1), (G3) and Lemma 3 these curves are
  $g_\epsilon$-timelike, so $Q \subset J^-(\mathcal E)$.
- **(b) Late part.** On $\{t_* \ge T_{\rm late},\ r \ge r_H + c\}$, $W_{\rm out}$ is
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
so $S \subset J^-(\mathcal E)$ by Step 3(c).

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

## 6 · Λ = 0 with null infinity: Theorem E

For Λ = 0, Theorem D measures escape to arbitrarily large distances. The
textbook black-hole region is
$\mathcal B_{\mathscr I} = M_\epsilon \setminus J^-(\mathscr I^+)$. Here
$\mathscr I^+$ is the null infinity of the late-time development, the face
$\mathscr I^+$ of the compactified spacetime in the 2026 paper. By Lemma
1(iv), $\mathcal B \subseteq \mathcal B_{\mathscr I}$.

Reaching $\mathscr I^+$ is more than reaching large $r$: the signal's
retarded time must also stay bounded. In Schwarzschild, for comparison, it
drifts like $2M\log r$ unless one uses the right retarded time. What is
needed is:

- **(N)** There is $R_N$, depending only on the background, such that every
  point of the late development near $\Sigma_{T_{\rm late}}$ with $r \ge R_N$
  lies in $J^-(\mathscr I^+)$.

We derive (N) from the 2026 paper's description of null infinity (Lemma 6
below).

**Theorem E.** Let Λ = 0. Choose $T_{\rm late}$ large enough
that the curves of Step 3(a), started in
$\{0 \le t_* \le T + 1,\ r_H + c \le r \le R_N + 1\} \setminus \{|x| < \delta'\}$,
reach $r \ge R_N + 1$ before $t_* = T_{\rm late}$. Then, for small $\epsilon$,
conclusions 1, 2 and 4 of Theorem D hold with $\mathcal B_{\mathscr I}$ in
place of $\mathcal B$. That is:

- the small hole's core $P$ and the large hole's interior lie in
  $\mathcal B_{\mathscr I}$, with closed trapped surfaces on every slice;
- before the plunge they lie in different components of
  $\mathcal B_{\mathscr I} \cap \Sigma_t$;
- after it they lie in the same component.

So these asymptotically flat vacuum spacetimes describe the merger of two
initially distinct black holes in the sense of
[07](07-the-merger-problem.md), with $\mathscr I^+$.

*Proof.* Conclusions 1, 2 and the "same component" half of 4 follow from
Theorem D, because $\mathcal B \subseteq \mathcal B_{\mathscr I}$: a connected
subset of $\mathcal B$ is a connected subset of $\mathcal B_{\mathscr I}$.

For the "different components" half we need the separating sphere $S$ of
Step 4 to lie in $J^-(\mathscr I^+)$. By Steps 3(c) and 3(a), each point of
$S$ has a $g_\epsilon$-timelike curve that reaches $r \ge R_N + 1$ before
$t_* = T_{\rm late}$. From there follow $\partial_t$ up to
$\Sigma_{T_{\rm late}}$. This keeps $r$ fixed and stays timelike, because for
$r \ge 10M$ we have $g(\partial_t, \partial_t) \le -4/5$ with $C^0$
closeness. The endpoint has $r \ge R_N$, so it lies in $J^-(\mathscr I^+)$
by (N), that is, by Lemma 6. ∎

*Session 7: localisation (3) and connectedness (5) for
$\mathcal B_{\mathscr I}$ as well.* Lemma 7 below extends (N) from the
neighbourhood of $\Sigma_{T_{\rm late}}$ to all late times. Corollary E1
then gives $\mathcal B_{\mathscr I} = \mathcal B$, so all of Theorem D
holds for the textbook black-hole region. (The session-6 version of this
note left (3) and (5) open. It suggested identifying
$\partial\mathcal B_{\mathscr I}$ with the dynamical horizon of
Chen–Klainerman and Hintz (Duke 2026); that is not needed.)

**Lemma 5 (escape where the data are exactly Kerr).** Work in Kerr, outside
the horizon ($\Delta > 0$). Let $dr_*/dr = (r^2 + a^2)/\Delta$ and
$\tilde u = t - r_* + a\arctan(r/a)$.

- (a) $\rho^2 g^{-1}(d\tilde u, d\tilde u) = -a^2r^2/(r^2 + a^2) - (\text{a nonnegative term})$.
  So $\tilde u$ is nondecreasing along future causal curves in $r > r_+$
  *(test, exact)*.
- (b) Along an outgoing principal null geodesic, $t - r_*$ is constant, so
  $\tilde u$ changes by less than $\pi|a|/2$.
- (c) Let $\Sigma$ be a $t_*$-slice and $r_1$ large enough compared with
  $r_0$. Then the outgoing principal null geodesic from a point of $\Sigma$
  with $r \ge r_1$ never meets $J^+(\Sigma \cap \{r < r_0\})$.

*Proof of (c).* Suppose a causal curve runs from $\Sigma \cap \{r < r_0\}$ to
a point $q$ of the geodesic. Take its last crossing of $\{r = r_0\}$; after
it the curve stays where $\tilde u$ is a time function. That gives
$\tilde u(q) \ge t_*(\Sigma) - r_*(r_0) + O(1)$. By (b),
$\tilde u(q) \le t_*(\Sigma) - r_*(r_1) + \pi|a|$. These are incompatible
for $r_1$ large. ∎

So the geodesic stays in the domain of dependence of
$\Sigma \cap \{r \ge r_0\}$. If the data on $\Sigma$ are exactly Kerr there,
the spacetime is exactly Kerr along it (local uniqueness), and it reaches
$\mathscr I^+$. Thus (N) holds outside the support of the perturbation.

**Why Lemma 5 is not enough.** Our signals leave from next to the small
hole, inside the future of the region the gluing has disturbed. By (a)–(b),
$\tilde u$ grows along them, so they never enter the exactly-Kerr region of
Lemma 5. They cross $\Sigma_{T_{\rm late}}$ where the data have been
perturbed, at size $O(\epsilon)$. There (N) must come from the 2026 paper's
description of null infinity.

**Lemma 6 ((N) holds).** This is our reading of the 2026 paper. Near
$\mathscr I^+$ it uses:

- the null coordinates of the mass-$m_0$ Schwarzschild model, $x^0$ and
  $x^1 = t_*$, where $t_*$ is the hyperboloidal time;
- near $\mathscr I^+$, $t_*$ equals the model's retarded time $u$;
- the null frame $\partial_0$ (outgoing, $\partial_0 r = \frac12(1 - 2m_0/r)$)
  and $\partial_1$ (ingoing), from its (3.30)–(3.33).

Then (N) holds.

*Proof.*

1. **Where the final metric is explicit.** In the proof of its Theorem 13.1
   (Steps 2–3), the final metric on $\{t_* < 1\}$ is the boosted exterior
   solution $\phi^*g_{b_0} + h'$. Two facts give this:
   - the Nash–Moser correction $\Phi(U)$ is supported in $\{t_* \ge 1\}$.
     Its Definition 12.1 places the unknowns on $\Omega^* = \mathrm{cl}\{t_* \ge 1\}$,
     and Theorem 12.3 maps into spaces of distributions supported there
     (the superscript $\bullet$);
   - the cut-off $\chi_K$ vanishes for $t_* \le 10$ (preamble of Theorem
     13.1). Here $\phi$ is the small asymptotic boost, and
   $h'$ lies in the space of its Definition 5.3, with small norm (its
   (5.56) and Step 2).
2. **The outgoing null component $g(\partial_0, \partial_0)$**, term by term:
   - *Kerr background.* The paper states
     $(g_{m,a} - g_m)_{\bar\mu\bar\nu} \in r^{-2}C^\infty$ (after its
     Definition 3.16), and $g_m(\partial_0, \partial_0) = 0$. So
     $g_{b_0}(\partial_0, \partial_0) = O(r^{-2})$. *(test: in these
     coordinates it is $-a^2\sin^2\theta/(4r^2) + O(r^{-3})$, and exactly 0
     for Schwarzschild.)*
   - *The boost.* By its Lemma 3.19,
     $\phi^*g_{m,a} - g_{m,a} = r^{-1}h_{S,11}(dx^1)^2 + O(\rho_0\rho_I^2)$.
     Since $dx^1(\partial_0) = 0$, the boost adds only $O(\rho_I^2)$ here.
   - *The perturbation.* $h'(\partial_0, \partial_0) = h'_{00}$ is part of
     $\pi^C h'$. Its index set at $\mathscr I^+$ has real parts above
     $1 + 2\gamma^\Upsilon$, and its remainder weight is $\ell_I > 1$
     (Definition 5.3, (5.8)). The spaces are weighted $b$-Sobolev spaces with
     respect to a $b$-density (its (5.2)), and $H_b^k \subset L^\infty$ for
     $k \ge 3$. So $|h'_{00}| \le C\epsilon'\rho_I^{1+\delta}$ (times a power
     of $\rho_0$), with $\delta > 0$ and $\epsilon' \to 0$ as $\epsilon \to 0$.

   Near $\mathscr I^+$ we have $\rho_I \sim |t_*|/r$. So at fixed $t_*$,
   $|g(\partial_0, \partial_0)| \le C(t_*)\,r^{-1-\delta}$. Also, for
   $r \ge R_N$, $g(\partial_0, \partial_1)$ is within $1/8$ of $-\frac12$ and
   $|g(\partial_1, \partial_1)| \le 1/8$.
3. **The escape curve.** Take $V = \partial_0 + f\,\partial_1$ with
   $f = 4|g(\partial_0, \partial_0)| + r^{-2}$. Then
   $g(V, V) \le g_{00} - \tfrac34 f + \tfrac18 f^2 < 0$, and $V$ is future
   directed. Along its integral curves:
   - $dr/d\lambda \ge 1/8$ and $dt_*/d\lambda = f$;
   - so $t_*$ increases by at most $\int 8f\,dr < \infty$. Near spatial
     infinity the bound on $h'_{00}$ carries a factor $\rho_0^\alpha$, so the
     integral is finite for each starting radius. For fixed $R_N$ and small
     $\epsilon$ it is below 1;
   - the curve reaches $r = \infty$ at finite $t_*$ and fixed angle, that is,
     at a point of $\mathscr I^+$.

   A starting point $q$ at $r \ge R_N$ near $\Sigma_{T_{\rm late}}$ has
   $t_*(q) \approx -r \ll 0$, so the whole curve stays in $\{t_* < 1\}$,
   where step 1 applies. ∎

**Lemma 7 ((N) at all late times).** There is $R_N$ such that every point
of the late development with $r \ge R_N$ lies in $J^-(\mathscr I^+)$. Here
$r$ and the hyperboloidal time $t_*$ (its (11.3)) are the coordinates of
the 2026 paper. This is our reading of that paper, like Lemma 6, and it
contains Lemma 6. *(Corrected in session 9, see the review below: the
session-7 proof split the cases at $t_* = 2$ instead of $t_* = 1$, and
wrongly said the support property of the Nash–Moser correction was not
needed.)*

*Proof.* Use the escape curve of Lemma 6, step 3: $V = \partial_0 + f\,\partial_1$
with $f = 4|g(\partial_0, \partial_0)| + r^{-2}$, started at a point $q$
with $r(q) \ge R_N$ and $t_*(q) = T$, at any time. Three facts are needed.

- **Frame, everywhere on $\{r \ge R_N\}$.** $g(\partial_0, \partial_1)$
  is within $1/8$ of $-\frac12$, $|g(\partial_1, \partial_1)| \le 1/8$,
  and $f \le 1/2$. These follow because:
  - $g = g_{b_0, b, -S} + h$ by its (13.3);
  - $h$ is small in $L^\infty$ (small weighted $b$-Sobolev norm);
  - the reference metric $(1 - \chi_K)(\phi_S)^*g_{b_0} + \chi_K g_b$ is
    within $O(|b - b_0| + |S|)$ of $g_{b_0}$;
  - the Kerr components tend to $-\frac12$ and $0$.

  So $dr/d\lambda \ge 1/8$ and $dt_*/d\lambda = f$ (since
  $dx^1 = dt_*$, its (3.29)).
- **Bounds for $t_* < 1$.** By Steps 2–3 of the proof of its Theorem
  13.1, $g = (\phi_S)^*g_{b_0} + h_{-S}$ there:
  - $\chi_- = 1$ on $(-\infty, \frac74]$;
  - the Nash–Moser part $\Phi(U)$ is supported in $\{t_* \ge 1\}$
    (Definition 12.1, Theorem 12.3);
  - $h_{-S}$ is an exterior solution, with the component structure of
    Definition 5.3.

  This is Lemma 6, step 2: $|g(\partial_0, \partial_0)| \le C\rho_I^{1+\delta}$,
  with a power of $\rho_0$ near $I^0$.
- **Bounds for $t_* \ge 1$**, on $\Omega^*$.
  - $h$ lies in the space of its (13.2). Restricted to $\Omega^*$ this
    is the space $X^k$ of its Lemma 11.5. The proof of that lemma
    records $h_{00} \in C_b^{k-3,(1+2\gamma^\Upsilon, 1-\epsilon, 2+\epsilon_K)}$,
    because $h_{00}$ is a component of $\pi^C h$ ((11.8), Definition 5.2).
  - Near $\mathscr I^+$ we have $\chi_K = 0$ ($\mathrm{supp}\,\chi_K \subset \{r < t/2\}$,
    Definition 3.17). So the reference metric is $(\phi_S)^*g_{b_0}$,
    whose outgoing component is $O(r^{-2}) + O(\rho_0\rho_I^2\rho_+)$
    by its (3.44); $\chi_{\rm IVP} = 1$ on $\Omega$.
  - So at finite $t_*$, $|g(\partial_0, \partial_0)| \le C_T\, r^{-1-2\gamma^\Upsilon}$
    for $r \ge C'T$. Here $\rho_I = (t_* + 2)/r$, its Lemma 3.15.
  - On the compact range $R_N \le r \le C'T$ it is bounded. There it
    contains the mass-mismatch term $(m - m_0)/r$ from $\chi_K g_b$, and
    $h_{00} = O(\rho_+^{1-\epsilon})$.

So $\int_{R_N}^\infty f\,dr < \infty$ along the curve. It grows at most like
$\log T + T^{\epsilon}$, but only finiteness is used. $t_*$ increases by
at most $8\int f\,dr$, which closes the continuity argument. The curve
reaches $r = \infty$ at finite $t_*$, at a point of $\mathscr I^+$. ∎

*Coordinates.* The 2026 paper's coordinates are boosted relative to the
background coordinates by $\phi_S$. Near $\Sigma_{T_{\rm late}}$ the boost
shifts time by $O(|S|r)$ and space by $O(|S|t)$, so the two radii differ by
$O(1)$ there. At late times they differ by $O(|S|t)$. In Theorem D and
Corollary E1, $r$ on the late development means the paper's $r$. Escape
curves increase it at a rate $\ge 1/8$, so the choice is immaterial for
"reaching arbitrarily large $r$".

Why the reference mass matters: in the outgoing null frame of mass $m_0$,
Schwarzschild of mass $m$ has

```math
g(\partial_0, \partial_0) = \frac{(m - m_0)\,\bigl(1 - (m + m_0)/r\bigr)}{r - 2m} \sim \frac{m - m_0}{r},
```

which is not integrable: retarded time would drift like
$(m - m_0)\log r$ (`tests/test_eh_merger.py`). The final black hole has a
different mass from the ADM mass $m_0$, but in the 2026 paper's construction
the switch to it ($\chi_K$) stays away from $\mathscr I^+$. Along the curve,
$t_*$ can grow like $\log T$ for late starting points. That is harmless: only
finiteness is needed.

**Corollary E1 ($\mathcal B_{\mathscr I} = \mathcal B$).** Let Λ = 0. Take
the escape radius $R \ge R_N + 1$ (allowed by Lemma 1(iii)). Then
$\mathcal B_{\mathscr I} = \mathcal B$. Hence all five conclusions of
Theorem D hold for the textbook black-hole region
$M_\epsilon \setminus J^-(\mathscr I^+)$:

- both holes lie in it, with closed trapped surfaces;
- localisation within the margins;
- two components before the plunge and one after;
- connectedness.

*Proof.* $\mathcal B \subseteq \mathcal B_{\mathscr I}$ by Lemma 1(iv). For
the converse, let $p \notin \mathcal B$. A causal curve runs from $p$ to a
point $q$ with $r(q) \ge R$. There are two cases.

- If $t_*(q) \ge T_{\rm late}$, then $q \in J^-(\mathscr I^+)$ by Lemma 7.
- If $t_*(q) < T_{\rm late}$, follow $\partial_t$ from $q$ up to
  $\Sigma_{T_{\rm late}}$, as in the proof of Theorem E. This keeps $r$
  fixed and stays timelike for $r \ge 10M$. Then apply Lemma 7.

Either way $p \in J^-(\mathscr I^+)$, so $p \notin \mathcal B_{\mathscr I}$. ∎

**What Lemmas 6 and 7 rest on.** We read these statements but did not verify them
independently:

- the support properties in step 1, from Definition 12.1, the target
  spaces of Theorem 12.3, and the proof of Theorem 13.1 (session 6: read);
- that Definition 5.3 describes $h'$ uniformly up to $\mathscr I^+$ near
  its corner with spatial infinity. Definition 5.3 is stated on the exterior
  domain, which contains that corner.

For Lemma 7, add: (13.2)–(13.3), Definitions 3.17 and 5.2, (3.29)–(3.36),
(3.44), Lemma 3.15, and the proof of Lemma 11.5. Lemma 7 uses the support
property of the Nash–Moser correction through Lemma 6, for $t_* < 1$. The
session-7 text cited Remark 13.4(1), Lemma 11.38 and Definitions 11.47 and
11.55. Those confirm that the correction terms vanish near $\mathscr I^+$,
but (13.2) already contains this.

**Independent review of Lemmas 6 and 7 (session 9).** Each step was
re-checked against the text of the 2026 paper, not against our notes.

| Step | Check | Outcome |
| --- | --- | --- |
| null frame | (3.29)–(3.33): $x^1 = t_*$, $\partial_0 = \frac12(1 - \frac{2m}{r})\partial_r$ at fixed $x^1$, $g_m = -(1 - \frac{2m}{r})dx^0dx^1 + r^2\slashed g$ | confirmed: $\partial_0$ is outgoing null for $g_m$; $g_m(\partial_0, \partial_1) = -\frac12(1 - \frac{2m}{r})$; $dt_*(\partial_0) = 0$ |
| "$h_{00}$" | (3.35) splits $S^2T^*$ by $(dx^0)^2$, $dx^0dx^1$, …, so $h_{00} = h(\partial_0, \partial_0)$; Definition 5.2: $\pi^C h = (h_{00}, h_{0\bar b}, \frac12\slashed{\mathrm{tr}}\, h)$ | confirmed |
| Kerr part | $(g_{m,a})_{\bar\mu\bar\nu} - (g_m)_{\bar\mu\bar\nu} \in r^{-2}C^\infty(M)$, after Definition 3.16 | confirmed (also checked exactly in `tests/test_eh_merger.py`) |
| boost | (3.44) holds on $\mathrm{supp}\,\chi_{\rm IVP}(1 - \chi_K)$ | confirmed applicable: $\chi_{\rm IVP} = 1$ for $t_{\rm IVP} > -r/10$, and $\mathrm{supp}\,\chi_K \subset \{t > 50m_0,\ r < t/2\}$ (Definition 3.17) |
| weights in §5 | (5.2): weighted $b$-Sobolev spaces relative to a $b$-density, so weights are pointwise rates | confirmed, so Lemma 6's $\rho_I^{1+\delta}$ is right |
| weights in §§11–13 | different normalisation ($3 + \epsilon_I$ for the same objects, footnote 103) | handled by citing the paper's own conversion (proof of Lemma 11.5) instead of converting ourselves |
| metric for $t_* < 1$ | proof of Theorem 13.1, Steps 2–3: $\chi_- = 1$ on $(-\infty, \frac74]$, $\Phi(U)$ supported in $t_* \ge 1$ | confirmed; **correction**: Lemma 7 needs this too |
| split of Lemma 7's proof | session 7 used Remark 13.4(1) for $t_* \ge 2$ and Lemma 6 for $t_* \le 2$, but Lemma 6 covers only $t_* < 1$ | **corrected**: use (13.2) and Lemma 11.5 on $\Omega^* = \{t_* \ge 1\}$ |
| region near $\iota^+$ | the escape curve crosses $R_N \le r \le C'T$, where $g_{00}$ is only $O(1/r)$ (mass mismatch) | added to the proof: the range is compact for fixed $T$; $f \le \frac12$ uniformly |
| coordinates | the paper's coordinates are boosted by $\phi_S$ | added: displacement $O(\lvert S\rvert t)$ at late times; immaterial for escape |

What this review does *not* do:

- It does not re-prove the paper's statements (13.2), (3.44), Lemma 11.5.
- It does not check the identification of our $\Sigma_{T_{\rm late}}$ data
  with the paper's initial-value setting (Remark 13.3, §5.7.1).
- It was done by the same reasoning agent that wrote the lemmas, against the
  primary text. An outside expert check is still the most valuable review.

Lemmas 6 and 7 are the least robust link in Theorem E, and the place where an
expert check of the 2026 paper would help most.

**The Klainerman–Szeftel route (session 4, superseded).** Before the 2026
paper, the Λ = 0 case went through KS23 (arXiv:2104.11857; §3.4 and §3.8
read). Their spacetime has the needed structure:

- every point of ${}^{(ext)}\mathcal M$ lies on an outgoing null geodesic
  that reaches $\mathscr I^+$;
- on ${}^{(int)}\mathcal M \cap \{r \le r_+(1 - \delta_H/2)\}$, $r$ decreases
  along causal curves.

But both hold in KS's own coordinates, which would have to be anchored to
Kerr's on the initial layer, and we could not check that. Hintz's fixed
coordinates remove this step for Theorem D and reduce Theorem E to (N).

## 7 · What is and is not proved

| Statement | Status |
| --- | --- |
| Lemma 2(a): $g^{rr} = \Delta_r/\rho^2$ | **checked** (`tests/test_eh_merger.py`, exact) |
| Lemma 2(c): area identity, $\partial_r E > 0$ for all $r > 0$, all spins, Λ ≥ 0 | **checked** (exact identity); positivity read off term by term |
| Kerr–de Sitter solves $\mathrm{Ric} = \Lambda g$; reduces to Kerr and SdS | **checked** (rational parameters exactly, slow test; SdS symbolic) |
| Lemma 2(b) for Λ = 0: $\partial_t + \frac12\partial_r$ has $g(W, W) \le -\frac14$ for $r \ge 10M$ | **checked** (polynomial with negative coefficients) |
| Lemma 5(a): $\tilde u$ is a time function outside a Kerr hole | **checked** (exact identity) |
| Lemma 6, Kerr part: $g(\partial_0, \partial_0) = -a^2\sin^2\theta/(4r^2) + O(r^{-3})$, and 0 for Schwarzschild | **checked** (exact, limit) |
| Lemmas 1, 2(b), 2(d), 3, 4, 5(b,c) | written out here (elementary) |
| (G1)–(G3), Λ > 0 | Hintz III Theorem 1.1, Definition 4.1, (1.2), (6.2): **read**. Uniformity on the tube and over $X$ is *our reading* |
| (G3), Λ = 0 | Hintz 2026 Theorems 1.1 and 13.1, consequence (3), Remark 13.3, (5.56): **read**. Smallness of the decay constants for small data is *our reading* |
| $t_* = t_{\rm BL} - T(r)$ form of Hintz–Vasy's time function | cited (HV18 (3.9)), *not re-checked* |
| **Theorem D** (Λ > 0 and Λ = 0): both holes in $\mathcal B$, localisation, merger, $\mathcal B$ connected | **proved here, given (G1)–(G3)**; potentially new (Hintz lists the event horizon and connectedness as open) |
| Corollary D1: on the slow time scale, the merger happens when $C$ crosses $r_H$ | **proved here** (from Theorem D) |
| **Theorem E** (Λ = 0, $\mathscr I^+$): both holes in $\mathcal B_{\mathscr I}$, two components then one | **proved here, given (G1)–(G3) and Lemma 6**; potentially new. Strengthened by Corollary E1 below |
| Lemma 6 = (N) | derived from Hintz 2026 Lemma 3.19, Def. 5.3, (5.2), (5.56), proof of Thm 13.1 (**read**). Support properties and uniformity up to the corner of $\mathscr I^+$ are *our reading* |
| Lemma 7 ((N) at all late times) | derived from Hintz 2026 (13.2)–(13.3), Defs. 3.17, 5.2, (3.29)–(3.36), (3.44), Lemma 3.15, proof of Lemma 11.5, and Lemma 6 for $t_* < 1$ (**read**). Reviewed against the source in session 9, which corrected the proof's case split; outside expert check still advised |
| Corollary E1: $\mathcal B_{\mathscr I} = \mathcal B$, so localisation and connectedness hold for $\mathcal B_{\mathscr I}$ | **proved here, given (G1)–(G3) and Lemma 7** (session 7) |
| Mass mismatch: $g_m(\partial_0^{m_0}, \partial_0^{m_0}) = (m - m_0)(1 - (m + m_0)/r)/(r - 2m)$ | **checked** (exact) |
| $\mathcal B$ defined with $\mathscr I^+$ for Λ > 0 | not proved: needs the spacetime beyond $r_+$ |
| Exact location, smoothness and merger time of the horizon | not addressed (only to within $c$, $\eta$) |

**Novelty.** Hintz states the structure of the black-hole region ("conjecturally
connected") as future work. A bounded search (2026-10-29) found no later
treatment. Theorem D is therefore **potentially new**. The argument is
elementary given Hintz's estimates, and experts may regard it as routine; its
content is that $C^0$ control suffices. Theorem E (Λ = 0, with
$\mathscr I^+$) is potentially new in the same sense. It leans further on
our reading of the 2026 paper (Lemma 6).

## 8 · Physical picture

Before the plunge, the small hole is a tiny Kerr black hole moving through
nearly empty space near the big one.

- **Inside its horizon,** the small hole's own radius $\hat r$ only runs
  down, so nothing gets out.
- **A little outside its horizon,** a rocket can climb out of the small
  hole's well. There space looks almost flat, so it can coast away, and then
  outrun the big hole's pull to arbitrarily large distances (into the
  cosmological region, if Λ > 0). Its worldline is
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
- Theorem E needs (N): light from far out on a slice must reach null
  infinity even after crossing gravitational waves of size $\epsilon$. Why is
  reaching large $r$ not enough? (Think of the $2M\log r$ drift of retarded
  time in Schwarzschild.) Lemma 6 shows the gauge of the 2026 paper keeps
  the outgoing component $g(\partial_0, \partial_0)$ integrable. Why is that
  exactly the component that matters?
