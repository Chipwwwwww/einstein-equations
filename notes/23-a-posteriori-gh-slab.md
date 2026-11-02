# A posteriori GH Einstein evolution: mathematical reorientation and a finite-slab proof

**2026-11-02. Status:** proved a special-case finite-slab residual theorem and a derivative-loss counterexample; verified exact-rational certificate. **Not** a generic nonlinear Einstein validation theorem, first-common-horizon theorem, new globally evaluable merger, or claimed mathematical novelty.

## 1. Precise problem and literature positioning

The desired construction starts with constrained asymptotically flat Cauchy data, a smooth approximate metric on a strong-field spacetime slab, a **continuum** residual, and interval-enclosed hyperbolicity and geometry. It must produce a genuine Lorentzian vacuum metric with quantified metric and derivative errors and preserve the gauge and Einstein constraints. A pointwise numerical residual alone does not guarantee existence.

Relevant results, with their restrictions:

- **Hughes, Kato, Marsden (1977)**, *Well-posed quasi-linear second-order hyperbolic systems with applications to nonlinear elastodynamics and general relativity*, DOI 10.1007/BF00251584, https://authors.library.caltech.edu/records/yqzjw-rqj49 . In the Einstein harmonic-reduction setting the short-time Cauchy solution uses data in $H^{s+1}(\mathbb R^3)\times H^s(\mathbb R^3)$ for $s>3/2$, an appropriately regular Lorentzian background and a sufficiently small metric-component neighborhood preserving the hyperbolicity cone. The result is **local well-posedness**, not a published numerical a posteriori bound with evaluated constants for a merger slab. Its physical data must also satisfy the constraints.
- **Kato (1975)**, *The Cauchy Problem for Quasi-linear Symmetric Hyperbolic Systems*, Arch. Rational Mech. Anal. 58, 181–205, DOI 10.1007/BF00280740. An appropriate first-order symmetric hyperbolic reduction supplies local existence and energy estimates under regularity and positivity assumptions. It does not by itself certify discrete numerical data.
- **Lindblom, Scheel, Kidder, Owen, Rinne (2006)**, *A new generalized harmonic evolution system*, https://arxiv.org/abs/gr-qc/0512093 , equations (25)–(52). Their first-order variables $\{\psi_{ab},\Pi_{ab},\Phi_{iab}\}$ form a symmetric-hyperbolic system with a positive symmetrizer under the stated parameter and Lorentzian-geometry conditions, together with reduction, gauge and constraint subsystems and boundary fields. It is the closest direct formulation for a computer-assisted Einstein certificate. The subsidiary GH constraint propagation equation is derived from the contracted Bianchi identity; it is **not** valid to infer Einstein vacuum from a small reduced-equation residual without controlling these constraints.
- **Hintz (2024)**, *Gluing small black holes along timelike geodesics III*, https://arxiv.org/abs/2408.06715 . Assumes a smooth globally hyperbolic vacuum background and an inextendible timelike geodesic, constructs true solutions by perturbing an $O(\epsilon^\infty)$-accurate small Kerr-hole gluing on precompact subsets with specialized uniform linearized-gauge-fixed estimates. Includes an extreme-mass-ratio merger application, **not** generic comparable-mass merger data or an arbitrary residual-based numerical validator.
- **Christodoulou–Klainerman (1993)**, *The Global Nonlinear Stability of Minkowski Space*, and **Lindblad–Rodnianski (2010)**, *The Global Stability of Minkowski Space-Time in Harmonic Gauge*, https://annals.math.princeton.edu/2010/171-3/p01 : global control for suitably small asymptotically flat perturbations, in the authors' norms, **not** generic strong-field merger data.
- **Christodoulou (2009)**, *The Formation of Black Holes in General Relativity*, https://arxiv.org/abs/0805.3880 : geometric trapped-surface formation for prepared short-pulse data, not arbitrary binary-metric validation.
- **Klainerman–Szeftel Kerr-stability program** and **Szeftel (2026)**, *Kerr stability in the full subextremal range*, https://arxiv.org/abs/2609.40068 : extends nonlinear stability to $|a|<m$ for sufficiently small perturbations of Kerr. The full 166-page proof and all norm-level hypotheses were **not** exhaustively checked in this session, and should not be treated as a black-box theorem for strong-field comparable-mass merger slabs.
- Rigorous a posteriori methods for hyperbolic PDE (e.g. reconstruction/relative entropy for conservation laws and validated fixed points for nonlinear wave problems) require their exact regularity, inverse bounds and boundary hypotheses. They are **not** automatically available for ten coupled quasilinear Einstein tensor components.

There is no mathematically unambiguous single “strongest theorem” across all these settings: the function spaces, geometry, gauge, and boundary conditions differ. HKM gives a sharp directly relevant *local Cauchy* result; Hintz gives an unusually strong *specific gluing/existence* result; GH gives a practical *first-order symmetric hyperbolic formulation*. None directly supplies a general finite-slab numerical-merger a posteriori existence certificate.

**Selected missing hypothesis:** a uniform **continuum residual enclosure** and an exactly computable slab-to-slab stability mechanism. Ordinary merger code output does not provide guaranteed interval enclosures of second metric derivatives, symmetrizers, all commutators, boundary fluxes and the gauge/reduction-constraint propagation on every point of a slab. We attack this restricted verification problem first rather than pretending a small-data theorem covers the merger.

## 2. Proved no-go for a naive inverse estimate (standard resonant example)

**Proposition.** For the zero-Cauchy-data wave equation on $[0,2\pi]\times\mathbb T$, there is no bounded inverse from forcing $L^\infty_t H^{s-1}_x$ to solution $L^\infty_t H^{s+1}_x$.

**Proof.** Put

$$
f_n=n^{-(s-1)}\sin(nt)\sin(nx),\qquad
v_n=n^{-(s-1)}\left(\frac{\sin nt}{2n^2}
-\frac{t\cos nt}{2n}\right)\sin(nx).
$$

Then $(\partial_t^2-\partial_x^2)v_n=f_n$ and both initial data vanish. At $t=2\pi$,

$$
\frac{\|v_n(2\pi)\|_{H^{s+1}}}
{\|f_n\|_{L^\infty_tH^{s-1}}}
=\pi(n+n^{-1})\to\infty.
$$

This invalidates a naive fixed-Sobolev-space Newton–Kantorovich inverse, even on flat spacetime. It does **not** invalidate well-balanced energy-space inverses, validated tame estimates or Nash–Moser methods. QED.

## 3. Exact 4D Einstein subsector and rigorous a posteriori theorem

On the shrinking causal diamond

$$
D_T=\{(t,x):0\le t\le T,\ |x|+t\le X\},\qquad T\le X,
$$

consider the **four-dimensional** metric

$$
g_H=H(t,x)(-dt^2+dx^2)+dy^2+dz^2,\qquad H>0.
$$

Direct computation of the full Ricci tensor and density gives

$$
\partial_\mu(\sqrt{-g_H}\,g_H^{\mu\nu})=0,\quad
R_{tt}=-\tfrac12\Box\log H,\quad
R_{xx}=+\tfrac12\Box\log H,
$$

with all other Ricci components zero and $\Box=\partial_t^2-\partial_x^2$.
Thus this **special gauge-wave class** is exactly harmonic, and the
**full four-dimensional vacuum equation restricted to it** becomes
$\Box\phi=0$, $\phi=\log H$. Exact solutions in this class are flat
coordinate-distorted vacuum metrics. This is a rigorous manufactured
benchmark, **not a nontrivial-curvature merger model**.

**Theorem (proved, restricted class).** Suppose $H_{\rm app}>0$ is $C^2$ on
$D_T$, with $m\le H_{\rm app}\le M$, and let
$r=\Box\log H_{\rm app}$. On each subslab $[t_j,t_{j+1}]$, assume the
**continuum** bound $|r|\le\delta_j$. Let $\phi$ be the exact wave-equation
Cauchy solution with prescribed smooth initial data, set
$w=\phi-\log H_{\rm app}$, and suppose at $t_0=0$ we know

$$
\|w(0)\|_\infty\le E_0,\qquad
\max_\pm\|(\partial_t\pm\partial_x)w(0)\|_\infty\le P_0.
$$

For $h_j=t_{j+1}-t_j$, define the **computable** recurrence

$$
P_{j+1}=P_j+h_j\delta_j,\qquad
E_{j+1}=E_j+h_jP_j+\tfrac12h_j^2\delta_j.
$$

On the full nested subdiamond at $t=t_{j+1}$, the true solution satisfies
$\|w\|_\infty\le E_{j+1}$ and
$\max_\pm\|(\partial_t\pm\partial_x)w\|_\infty\le P_{j+1}$.
The reconstructed $g=e^\phi(-dt^2+dx^2)+dy^2+dz^2$ is smooth,
Lorentzian, harmonic, and satisfies **all** components of Ric$(g)=0$.
If $E_{j+1}<1$, then for the nonzero differing components,

$$
|g-g_{\rm app}|\le {M E\over1-E},\qquad
|g^{-1}-g_{\rm app}^{-1}|\le {E\over m(1-E)},
$$

and with $K_a=\|\partial_a\log H_{\rm app}\|_\infty$,
$a\in\{t,x\}$,

$$
|\partial_a(g-g_{\rm app})|
\le {M(EK_a+P)\over1-E}.
$$

**Proof.** The one-dimensional d'Alembert Cauchy problem gives an exact
$\phi$. With $p_\pm=(\partial_t\pm\partial_x)w$,
$(\partial_t\mp\partial_x)p_\pm=-r$. The past null
characteristic of any point in $D_T$ remains in $D_T$, yielding the
recurrence for $P$ by direct integration. Since
$w_t=(p_++p_-)/2$, integrating in time yields the recurrence for $E$.
Vacuum and gauge follow from the **full** Ricci identities above.
The bound $e^E-1\le E/(1-E)$ for $0\le E<1$ and the product rule give
all displayed metric, inverse and derivative estimates. QED.

This theorem is intentionally *restricted*. Its linearization in
$\log H$ is **not** valid for generic ten-field Einstein equations;
no nonlinear full-Einstein existence theorem is claimed.

## 4. Exact-rational reproduction and controls

Set $X=1$, $T=1/4$,

$$
H_{\rm app}=1+\tfrac12(x-t)+\tfrac1{100}t^2.
$$

By analytic inequalities holding on **every point** of $D_T$,

$$
\tfrac12\le H_{\rm app}\le\tfrac{2401}{1600},\qquad
r=\frac{1+\tfrac12(x+t)-t^2/100}{50H_{\rm app}^2},
\qquad |r|\le\tfrac3{25}.
$$

The comparator $H_{\rm true}=1+\tfrac12(x-t)$ has **exactly
the same initial position and velocity**. Two rational slabs of width
$1/8$ yield

$$
P_1=\tfrac3{200},\quad E_1=\tfrac3{3200},\quad
P_2=\tfrac3{100},\quad E_2=\tfrac3{800},
$$

and the certified metric-component bound

$$
|g-g_{\rm app}|\le\frac{7203}{1275200}\approx0.0056485.
$$

The independently known actual bound is $T^2/100=1/1600=0.000625$;
it lies within the certificate. No floating-point sampling or
unverified interpolation is used to assert the continuum residual
bound. Symbolic tests independently compute the **full four-dimensional
Ricci tensor**, harmonic condition and comparison with the exact
solution. A separate Schwarzschild Kerr–Schild Ricci test confirms
regular geometry at $r=2M$, but is **not** a residual-to-error estimate
for the Schwarzschild or merger backgrounds.

To reproduce:

    python -m research.calc.gh_slab_certificate
    python -m pytest tests/test_gh_slab_certificate.py -q

Local baseline at submission: 5 passed. The arithmetic certificate
is exact-rational, and the sign and Ricci formulas are checked by SymPy.

## 5. Additional conditional general stability lemma and missing bridge

**Lemma (proved, conditional L2 stability, standard energy method).**
On a spatial torus, let smooth fields $u$ and $v$ solve a first-order
symmetric system

$$
A^0(u)u_t+\sum_i A^i(u)u_i=F(u),\quad
r=A^0(v)v_t+\sum_i A^i(v)v_i-F(v),
$$

with symmetric matrices, $\lambda I\le A^0(u)\le\Lambda I$,
and on the solution/approximation tube operator Lipschitz constants
$L_\mu$ for $A^\mu$ and $L_F$ for $F$. Assume certified bounds

$$
D\ge\left\|\partial_tA^0(u)+\sum_i\partial_i A^i(u)\right\|_{L^\infty},
\quad
B\ge\sum_{\mu=0}^3L_\mu\|\partial_\mu v\|_{L^\infty}.
$$

Then for $e=u-v$,
$Y=(\int e^\top A^0(u)e\,dx)^{1/2}$,

$$
Y(t)\le e^{\kappa(t-t_0)/2}Y(t_0)
+\frac1{\sqrt\lambda}\int_{t_0}^t
e^{\kappa(t-s)/2}\|r(s)\|_{L^2}\,ds,\qquad
\kappa=\frac{D+2L_F+2B}{\lambda}.
$$

**Proof.** Subtract the equations, contract with $e$, integrate the
symmetric spatial principal part by parts (periodic boundary), use
the operator bounds and $E\ge\lambda\|e\|_2^2$, then apply the scalar
integrating factor to $\sqrt{E}$. QED.

**Important limitation:** $D$ depends on derivatives of the unknown
true $u$. This is a conditional stability estimate, **not** an
a posteriori existence theorem. Closing this gap for GH Einstein
requires an explicit $H^s$ high-norm bootstrap with commutator
control, a computable symmetrizer/hyperbolicity tube, compatible
incoming boundary or excision flux estimates, and a proof that
the GH gauge and reduction constraints propagate. An $L^2$
estimate alone does not control pointwise hyperbolicity.

## 6. Next theorem to attempt

Derive a **two-level tame, quantitatively computable finite-slab
existence theorem for first-order GH Einstein**, tracking
$H^s$/lower-energy residuals and all symmetrizer/commutator constants,
whose closure criterion is checked on validated data from a numerical
binary merger. Work first on a finite causal slab without incoming
boundary flux, then extend to moving inner boundaries/excision and
asymptotic exterior. Small-data Minkowski and near-Kerr theorems cannot
be inserted into the merger step without proving their hypotheses.

No certified extension of the existing Brill–Lindquist merger metric has
been obtained in this session, and the present work does **not**
constitute a previously unknown merger metric.
