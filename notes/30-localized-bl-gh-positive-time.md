# 30 — A quantitative local full-GH Einstein development of equal-mass BL data

**Research date:** 2026-11-02. **Base:** PR #10 at `f4e56812cba4c057fc89af775641b244e318ec1f`.  
**Status:** a constructive *extremely short-time* local vacuum certificate using published symmetric-hyperbolic existence/uniqueness and GH subsidiary theorems, analytic exact-rational majorants and an explicit GH two-jet evaluator. **No merger, no extension of prior spacetime domain, and no claim to a new PDE principle.** The numerical lifespan is astronomically pessimistic and has no physical usefulness. Independent peer review and remote CI are pending.

## 1. Convention and precise physical region

Signature $(-,+,+,+)$; $K_{ij}=-(2N)^{-1}(\partial_t\gamma_{ij}-\mathcal L_\beta\gamma_{ij})$; $\Pi_{ab}=-n^c\partial_cg_{ab}$; $\Phi_{iab}=\partial_i g_{ab}$. Use Lindblom–Scheel–Kidder–Owen–Rinne (2006), first-order parameters $\gamma_1=-1$, $\gamma_2=\gamma_0=0$ and $H_a=0$. The genuine evolution is $\partial_tu+A^k(u)\partial_ku=F(u)$ for 50 fields $u=(g_{ab},\Pi_{ab},\Phi_{iab})$, implemented in `gh_full_source.py` with all terms of their equations (35)–(39), including the actual gravitational quadratic contractions. In contrast to the geodesic-slicing evaluator from PR #7, the GH initial lapse is $N=p^{-2}$ and the initial shift vanishes.

Let $d=3831/5000=0.7662$, $c_\pm=(0,0,\pm3831/10000)$, $r_\pm=|x-c_\pm|$,
\[
 p=1+\frac1{4r_+}+\frac1{4r_-},\quad
 g_0=\operatorname{diag}(-p^{-4},p^4,p^4,p^4),\quad
 \Pi_0=0,\quad \Phi_{i0}=\partial_i g_0.
\]
Write
\[
 \Omega_{\rm phys}=\{|x|<11/10,\ r_+>1/5,\ r_->1/5\},\quad
 \Omega_*=\{|x|\le9/10,\ r_+\ge3/10,\ r_-\ge3/10\}.
\]
The Euclidean separation from $\Omega_*$ to the complement of $\Omega_{\rm phys}$ is **at least $1/10$** by the 1-Lipschitz property of each radial distance. $x=0$ lies in $\Omega_*$ and $p(0)=1+5000/3831>2.3$, $g_{ii}(0)>28$; this is near the binary, not a weak-field far-exterior substitution. This region is not certified to contain a complete common horizon.

## 2. A global *smooth nonvacuum* auxiliary torus extension

Work on $\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3$ with canonical $[-\pi,\pi]^3$. Put $\theta(s)=e^{-1/s}$ for $s>0$, $\theta(s)=0$ otherwise, and $\chi(s)=\theta(s)/(\theta(s)+\theta(1-s))$ for $0<s<1$, extended by $0$ for $s\le0$ and $1$ for $s\ge1$. Then $\chi\in C^\infty(\mathbb R)$, $0\le\chi\le1$. With $a=3/25$, $b=9/50$, $R_0=11/10$, $R_1=3/2$, define
\[
 I_\pm(x)=\chi\left(\frac{r_\pm^2-a^2}{b^2-a^2}\right),\quad
 O(x)=\chi\left(\frac{R_1^2-|x|^2}{R_1^2-R_0^2}\right),
\quad p_E(x)=1+\frac{O(x)}4\left(\frac{I_+(x)}{r_+(x)}+\frac{I_-(x)}{r_-(x)}\right).
\]
Set $I_\pm/r_\pm=0$ where $r_\pm\le a$. This is $C^\infty$ at the puncture centers. It equals $p$ **exactly** on $\Omega_{\rm phys}$ and equals $1$ near the torus seams because $R_1<\pi$; hence it is a legitimate positive smooth periodic function. The auxiliary 50 fields are $u_{E,0}=(\operatorname{diag}(-p_E^{-4},p_E^4\delta),0,\partial_ig_E)$.

**Critical integrity condition:** $\Delta p_E$ need not vanish outside $\Omega_{\rm phys}$. The auxiliary torus solution is merely a solution of the *gauge-reduced first-order GH PDE*; no global vacuum claim is made. We recover vacuum only inside the quantitative protected domain of dependence.

## 3. Analytic continuum derivative and Sobolev bounds

No grid sampling is used. If $r=s^{-1}$, the $n$th derivative of $\theta(s)$ for $s>0$ is $e^{-r}P_n(r)$, where $P_0=1$ and $P_{n+1}=r^2(P_n-P_n')$. Since $e^{-r}r^m\le m^m$ and $\theta(s)+\theta(1-s)\ge e^{-2}>1/8$ on $0\le s\le1$, quotient-rule recursion yields **integer** upper bounds $S_n\ge\sup|\chi^{(n)}|$. Mixed derivatives of $\chi(Q(x))$ for quadratic $Q$ are enclosed using explicit multivariate Faà di Bruno coefficients, computed as coefficients of $(3D_1z+\tfrac32D_2z^2)^j$ with exact rational arithmetic. For $r\ge a$, the elementary recursion for $r^{-1}$ gives $|\partial^\alpha r^{-1}|\le6^{|\alpha|}|\alpha|!a^{-1-|\alpha|}$. Leibniz gives a global bound on $|\partial^\alpha p_E|$, $|\alpha|\le8$, with analytic positivity $1\le p_E\le31/6$. Explicit power/reciprocal Faà di Bruno bounds then control derivatives of $p_E^{\pm4}$ and $p_E^{-2}$ through order 8.

In the **normalized Fourier norm** $\|f\|_{H^s(\mathbb T^3)}^2=\sum_{k\in\mathbb Z^3}(1+|k|^2)^s|\hat f(k)|^2$, for integer $s\le7$,
\[
 \|f\|_{H^s}\le2^s\max_{|\alpha|\le s}\|\partial^\alpha f\|_{L^\infty},
\]
from expanding $(1+|k|^2)^s$ into positive derivative terms. The independently executable rational derivative recursion establishes the certified power-of-ten envelopes
\[
 \max_{|\alpha|\le8}\|\partial^\alpha\chi(\cdots)\|_\infty<10^{36},\quad
 \max_{|\alpha|\le8}\|\partial^\alpha p_E\|_\infty<10^{95},\quad
 \max_{|\alpha|\le8}\|\partial^\alpha p_E^{\pm4}\|_\infty<10^{789}.
\]
Thus **each reference field, inverse metric entry and $p_E^{-2}$ has $H^7$ norm less than $10^{795}$**; the latter bound also absorbs the finite 50-field dimension. Numbers are deliberately overlarge but are algebraically derived, not fitted. See `bl_local_cutoff.py` for the exact recurrence and `test_bl_local_certification.py` for its cross-checks.

## 4. Full 50-field nonlinear source calculus about a strong nonconstant metric

The previous PR #10 near-Kasner small-Minkowski-deviation graph cannot directly enclose $g_E$ because its spatial metric is $p_E^4\delta$ and is $>28\delta$ at the midpoint. We instead center at the **actual variable BL reference** $g_E$.

Let $U=795$, $P=192$ (a valid conservative multiplication constant in each of $H^3,H^4,H^5$), and choose a candidate tube $\|u-u_E(t)\|_{H^3}\le R=10^{-8050}$. Once $\|u-u_{E,0}\|_{H^3}\le2R$, since $\|g_{E,0}^{-1}\|_{H^3}$ is bounded, the matrix Neumann series for
\[
 g^{-1}=(I+g_{E,0}^{-1}(g-g_{E,0}))^{-1}g_{E,0}^{-1}
\]
converges in the **entire** $H^3$ ball. There is a uniform bound on the inverse and its Fréchet derivative. To bound lapse, put $-g^{00}=p_E^4(1+w)$ with $w=p_E^{-4}(-g^{00}-p_E^4)$; the $H^3$ binomial series $(1+w)^{-1/2}$ and its derivative converge uniformly for $P\|w\|_{H^3}<1/2$. This yields pointwise Lorentzian signature, positive lapse, positive spatial metric, and positive symmetrizer by a straight-line metric-inertia argument. Nothing depends on evaluating finitely many perturbation states.

**Exact power-of-ten bound graph:** a scalar node `E(n,l)` means $\|f\|_{H^s}\le10^n$ and $\|Df\|_{H^s\to H^s}\le10^l$ throughout the tube. An exact sum takes $\max(n_1,n_2)+1$, and a product takes $n_1+n_2+3$; its derivative takes $\max(n_1+l_2,n_2+l_1)+4$. These rules dominate the Banach product and chain rules for $P=192$, including repeated tensors/matrix index contractions, and are applied to **all 50 original GH source expressions** by calling `full_source()` with reference-centered certified metric-inverse/lapse/shift entries. No source term is omitted. All graph computations are integer operations and work under Python `-O`.

The resulting conservative bounds for every GH state in the full tube are:

| Quantity | Certified bound |
|---|---:|
| $\|F(u)\|_{H^3}$ | $10^{8871}$ |
| full source Fréchet $H^3$ Lipschitz constant | $10^{11287}$ |
| $\|A^i(u)\|_{H^3}$, each $i$ | $10^{4013}$ |
| $H^3$ Lipschitz constant of $A^i$, each $i$ | $10^{6425}$ |
| $\|DS(g)\|$ pointwise | $10^{5620}$ |
| $\|D(SA^i)(g)\|$ pointwise, each $i$ | $10^{8029}$ |
| $\|g_t\|_{L^\infty}$ on tube | $10^{3209}$ |

The proof is especially conservative. $S(g)$ and $SA^i(g)$ depend only on the **ten metric components**, so their divergence is bounded by $D=J_S\sup|g_t|+\sum_i J_{B_i}\sup|\partial_i g|$ without derivatives of unknown $\Pi_t,\Phi_t$. Their contribution plus the high-order commutators yields a valid global $H^3$ rate
\[
 K<10^{11302}, \qquad \lambda_{\min}(S)>10^{-4},\quad \lambda_{\max}(S)<2.
\]
The full source itself is genuinely 50-field, with no symmetry reduction. The **generic shift and off-diagonal metric** regression checks the independent second-order reduced Einstein identity using nontrivial off-shell jets.

## 5. Evaluable positive-time GH jet from real BL data

Define the exact global auxiliary first time jet (all 50 components) by
\[
 u_1(x)=F(u_{E,0})-A^i(u_{E,0})\partial_i u_{E,0}.
\]
This evaluates from the explicit $p_E$, its first and second derivatives and the **full** GH operator, without importing geodesic-slicing time coefficients. It has vanishing $g$ and $\Phi$ components, with $\Pi$ component generally nonzero. Set $u_2^g=-p_E^{-2}u_1^\Pi$, and $u_2^\Pi=u_2^\Phi=0$. Then the globally smooth auxiliary approximate evolution is
\[
 v(t,x)=u_{E,0}(x)+t u_1(x)+\tfrac12t^2u_2(x).
\]
It obeys the exact full-GH initial data, and its 50-field GH residual vanishes at $t=0$. Derivative majorants prove $\|u_1\|_{H^4}<10^{8878}$, $\|u_2\|_{H^4}<10^{9684}$.

**Explicit physical metric evaluator** ($\Delta p=0$ only within $\Omega_{\rm phys}$): with $p_i=\partial_i p$ and $D=\sum p_i^2$, the complete nonlinear GH initial equations yield
\[
 \dot\Pi_{00}(0)=\frac{4D}{p^{12}},\quad \dot\Pi_{0i}(0)=0,\quad
 \dot\Pi_{ij}(0)=\frac{4(D\delta_{ij}-4p_i p_j)}{p^4}.
\]
Thus on $\Omega_*$, the ten-component metric approximation is **explicit**:
\[
\boxed{g_{00}^{app}(t,x)=-p^{-4}-\frac{2t^2D}{p^{14}},\quad
 g_{0i}^{app}=0,\quad
 g_{ij}^{app}(t,x)=p^4\delta_{ij}+\frac{t^2}{p^6}
 (8p_i p_j-2D\delta_{ij}).}
\]
This is not the exact metric; it is an *evaluation formula with a rigorous error*. $\Pi^{app}=t\dot\Pi(0)$ and $\Phi^{app}_i=\partial_i g_0$ complete the 50 approximate fields. The symbolic evaluator `bl_gh_taylor_evaluator.py` checks exact rational coordinate/time inputs against the certified physical domain, and independently cross-checks selected jet values against all the original GH source contractions. At the reflection-symmetric midpoint $\nabla p=0$, the $t^2$ metric coefficient vanishes in this **particular harmonic lapse choice**; it is nonzero at generic neighbouring points. The nontriviality is independently tested.

## 6. Continuum residual and strict high-order bootstrap

Let $G(u)=F(u)-A^i(u)\partial_i u$. Since $v_t(0)=u_1=G(u_{E,0})$,
\[
 r_{GH}(t)=v_t-G(v)
 =tu_2-[F(v)-F(u_0)]+\sum_i[A^i(v)\partial_i v-A^i(u_0)\partial_i u_0].
\]
This identity is exact on the **entire** torus; the closed polynomial $v$ contains the spatially varying cutoff data outside the physical region too. Using the *full-ball* source/coefficients from §4 and the explicit $H^4$ bounds of $u_1,u_2$ gives
\[
 \|r_{GH}(t)\|_{H^3}\le t\,10^{20985},\quad
 \rho:=\sup_{0\le t\le T}\|r_{GH}(t)\|_{H^3}\le10^{-8060},
 \qquad T=10^{-29045}.
\]
This is a **continuum analytic majorant**, not a mesh-residual test. The true reduced-GH initial error is *exactly zero*. From the symmetric first-order $H^3$ estimate, the metric-sector-adapted divergence and the strictly positive symmetrizer,
\[
 \|u(t)-v(t)\|_{H^3}
 \le\frac{\Lambda\rho T}{\lambda(1-KT)}
 <10^{-37098}<R=10^{-8050},
\]
where $KT<10^{-17743}<1/2$. Kato-type local smooth symmetrizable-hyperbolic existence provides a nonzero starting interval; the strict tube bound excludes a first exit, and the $H^3$ continuation criterion excludes termination before $T$. Hence a unique smooth **reduced full 50-field GH** solution exists on the whole auxiliary torus for $0\le t\le T$.

No unknown exact-solution derivative appears in $K$; on the bootstrap interval the metric component $g_t=F_g(u)$ and $\nabla g$ are enclosed by the known full-tube graph. The estimate does not rely on spherical or axial symmetry.

## 7. A genuinely quantitative finite-speed local vacuum theorem

**Theorem 30.1 (local computable GH vacuum recovery; conditional only on the standard precisely stated symmetric-hyperbolic and GH subsidiary existence/uniqueness theorems).** Let smooth physical vacuum Cauchy data and initial gauge-compatible 50-field GH data be given on an open $\Omega$. Suppose a smooth periodic auxiliary extension agrees there; a full-neighbourhood $H^3$ residual/positive-symmetrizer bootstrap produces an exact reduced GH solution $u$ on $[0,T]\times\mathbb T^3$; and every evolution and **subsidiary constraint** characteristic speed in Euclidean spatial directions has magnitude $\le C$. For a closed target $\Omega_*$ with Euclidean distance $b=\mathrm{dist}(\Omega_*,\mathbb T^3\setminus\Omega)>CT$, the auxiliary solution restricted to $[0,T]\times\Omega_*$ solves $\mathrm{Ric}(g)=0$, is unique there as a local vacuum development up to the standard geometric isometry, and inherits the computed $H^3$ metric error bound and pointwise first-derivative consequences. No auxiliary vacuum constraint is needed outside $\Omega$.

**Proof.** The symmetric first-order GH system has a strictly positive symmetrizer $S$ and symmetric flux matrices $SA^i$. For two solutions with equal data on $\Omega$, their difference $e$ solves a linearized-in-difference symmetric hyperbolic equation. On the shrinking interior $\Omega_t=\{x:\mathrm{dist}(x,\mathbb T^3\setminus\Omega)>Ct\}$, integration of $e^T S e$ gives an outgoing-boundary quadratic form $(V_n S-SA^n)$ with boundary normal velocity $V_n=-C$. Since $\lambda_{\max}(A^n)\le C$, $-CS-SA^n\preceq0$; the usual lower-order Gronwall estimate with zero data proves $e=0$ in the corresponding lens domain. Lipschitz level sets are handled by smooth approximating distance functions/coarea limits, or the standard domain-of-dependence theorem of Friedrich–Rendall, §3.1, whose boundary flux argument is precisely this one. The same argument holds for the *homogeneous* GH subsidiary first-order system because its characteristic speeds are $0,-\beta\cdot n,-\beta\cdot n\pm N\sqrt{\gamma^{ij}n_in_j}$ for $\gamma_1=-1$, $\gamma_2=0$ (Lindblom et al. §4.1). On $\Omega$, the reduction constraints vanish initially by construction; $C_a=H_a+\Gamma_a=0$ because $\sqrt{-g}g^{ab}=\mathrm{diag}(-p^8,1,1,1)$ and initial time derivatives of $g$ vanish. The normal derivative $\nabla_n C_a$ vanishes because Hamiltonian/momentum constraints are exactly zero. Independently, the explicit full-GH first jet above gives $\dot C_i=N^{-1}\dot\Pi_{0i}=0$ and $\dot C_0=N^{-1}\dot\Pi_{00}+(N/2)g^{bc}\dot\Pi_{bc}=0$ for arbitrary $p>0$, $\Delta p=0$; the latter is a **direct exact algebraic check**, not a reliance on a heuristic damping argument. Thus all subsidiary data vanish on $\Omega$ and finite propagation keeps them exactly zero on the target causal diamond. There $C_a=0$ and the entire first-order reduction constraints vanish, so the reduced equations imply $R_{ab}=0$. Local geometric Cauchy uniqueness identifies this solution with the original BL development up to diffeomorphism, not as numerically identical geodesic-slicing components. QED.

**Reference-adapted characteristic bound for THIS concrete application.** At $t=0$, $\beta=0$, $N=p_E^{-2}\le1$, $\gamma^{ij}=p_E^{-4}\delta^{ij}$, hence the physical coordinate wave speed is $p_E^{-4}\le1$. On the *entire certified true-solution tube*, the spacetime inverse perturbation satisfies $\|g^{-1}-g_0^{-1}\|_{op}<10^{-8040}$: the reference inverse has op norm $<1000$, metric perturbation op norm $<50R$, and the matrix Neumann estimate gives $\Delta g^{-1}\le10^8R=10^{-8042}$. Consequently $N<2$, $|\beta|_2<4\times10^{-8040}$, and $\lambda_{\max}(\gamma^{ij})<2$, so **every GH and subsidiary characteristic speed is strictly below $C=4$**. This geometric estimate is vastly sharper than the generic full-expression H3 graph bound. Since $4T=4\times10^{-29045}<b=1/10$, all causal pasts from $[0,T]\times\Omega_*$ stay in the exact physical initial region.

**Concrete conclusion.** For every $(t,x)\in[0,10^{-29045}]\times\Omega_*$, the true BL vacuum development has a unique GH metric with the above explicit evaluator $g^{app}$ and
\[
 \|u-v\|_{L_t^\infty H_x^3(\mathbb T^3)}<10^{-37098},\quad
 \max_{a,b}|g_{ab}(t,x)-g^{app}_{ab}(t,x)|<6\cdot10^{-37098},\quad
 \max_{a,b,i}|\partial_i(g_{ab}-g^{app}_{ab})|<12\cdot10^{-37098}.
\]
These are GH coordinate-component errors. H3 in three dimensions does **not** control pointwise second derivatives and curvature errors; no curvature error certificate is claimed. The auxiliary torus metric is NOT necessarily Ricci-flat outside the protected region.

## 8. Relation to previous work, novelty and failure tests

- The earlier certified geodesic-slicing BL metric is valid on its own specified tube for about $8.683\times10^{-4}M$ in that clock. Our $10^{-29045}M$ harmonic-time interval is vastly shorter; **no region has been proved geometrically beyond the earlier spacetime**, and there is no certified gauge-coordinate overlap map. Local geometric uniqueness proves only local equivalence, not a new merger or extension beyond a known common horizon.
- The theorem exploits **standard** finite propagation/domain-of-dependence theory and GH homogeneous constraint propagation, applied with explicit constant and exact auxiliary data. The resulting whole-ball BL certificate is new to this repository as a computable mathematical checkpoint; no genuinely new PDE method or publishable novelty is claimed.
- There is no numerical relativity mesh residual. The approximate fields are explicit smooth functions with analytic continuum residual majorants. A finite set of 37+ positive and negative tests cannot substitute for this proof; independent mathematical peer review remains necessary.
- Falsification: evaluator refuses times larger than certified $T$, punctures, and coordinates outside the protected target. The auxiliary extension is never labeled globally vacuum. An independent generic off-shell metric jet with **nonzero shift and off-diagonal spatial metric** is checked against the second-order reduced Einstein wave tensor. The initial normal GH constraint is checked symbolically for arbitrary $p$ and its three gradients. All certificate conditions fail closed under Python `-O`.
- The absurdly tiny $T$ is the main practical obstruction, caused by global cutoff $H^7$ constants $\sim10^{795}$ and source-Jacobian bounds $\sim10^{11287}$; further decimal precision will not fix it. The next mathematical task is a **localized, background-relative high-order energy and functional-calculus estimate** avoiding global $C^\infty$ cutoff derivatives, then compactly supported spatial jets and validated characteristic-domain decomposition. Only that could make a strong-field merger-scale slab conceivable.

## 9. Exact reproduction

From the research repository root:

```bash
python -m research.calc.bl_local_cutoff
python -m research.calc.bl_torus_majorant
python -O -m research.calc.bl_torus_majorant
python -m research.calc.bl_gh_taylor_evaluator
python -m pytest -q tests/test_bl_local_certification.py tests/test_gh_generic_shift_jet.py
python -m pytest -q tests/test_gh_full_source_majorant.py tests/test_gh_reference_neumann.py tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py tests/test_full_hs_bootstrap.py tests/test_gh_50_principal.py tests/test_gh_ks_source.py tests/test_gh_ks_geometry_tube.py tests/test_bl_gh_initial_data.py tests/test_bl_local_certification.py tests/test_gh_generic_shift_jet.py
```

**Primary literature:** Kato (1975), *The Cauchy Problem for Quasi-linear Symmetric Hyperbolic Systems*, Arch. Ration. Mech. Anal. 58; Hughes–Kato–Marsden (1977), *Well-posed quasi-linear second-order hyperbolic systems*, Arch. Ration. Mech. Anal. 63; Friedrich–Rendall (2000), *The Cauchy Problem for the Einstein Equations*, §3.1, https://arxiv.org/abs/gr-qc/0002074; Lindblom–Scheel–Kidder–Owen–Rinne (2006), eqs. (35)–(39), (40)–(58), https://arxiv.org/abs/gr-qc/0512093; Rinne (2006), https://arxiv.org/abs/gr-qc/0606053, and Kreiss–Reula–Sarbach–Winicour (2009), https://doi.org/10.1007/s00220-009-0788-2, for physically imposed boundary conditions **not used** in this no-artificial-incoming-boundary construction. Literature results are invoked only under their exact smooth periodic/symmetric-hyperbolic/constraint assumptions. Full literature novelty search is not exhaustive.