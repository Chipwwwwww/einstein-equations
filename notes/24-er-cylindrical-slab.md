# Nonflat Einstein–Rosen vacuum: computable nonlinear slab validation

**2026-11-02. Status:** A proved symmetry-restricted residual-to-exact-vacuum theorem with nonlinear constraint propagation and an exact-rational proof certificate. This is a **nonflat**, curved Einstein test, but **not** a generic Einstein–GH continuation theorem, asymptotically flat black-hole binary or new mathematical discovery.

## 1. A genuinely curved vacuum Einstein sector

In coordinates $(t,r,z,\varphi)$ with $r>0$ and periodic $\varphi$, set

$$
g=e^{2(\gamma-\psi)}(-dt^2+dr^2)+e^{2\psi}dz^2+r^2e^{-2\psi}d\varphi^2.
$$

Direct Ricci calculation reduces **all four-dimensional vacuum equations** to the standard Einstein–Rosen system

$$
\psi_{tt}-\psi_{rr}-r^{-1}\psi_r=0,\qquad
\gamma_r=r(\psi_t^2+\psi_r^2),\qquad
\gamma_t=2r\psi_t\psi_r.
$$

References: Einstein and Rosen, *On Gravitational Waves*, J. Franklin Inst. 223 (1937), 43–54, DOI 10.1016/S0016-0032(37)90583-0; see also https://www.mdpi.com/2073-8994/13/8/1425 for the displayed metric/reduction. This is a standard integrable polarized cylindrical-wave reduction: the nonlinear quadratic $\gamma$ constraints matter even though $\psi$ satisfies a linear radial wave equation. No Cartesian GH condition is imposed on these coordinates. Cylinder symmetry means that 4D asymptotic flatness in the standard isolated-binary sense is **not** available.

**Domain:** $D_T=\{0\le t\le T,\ r_-+t\le r\le r_+-t\}$ with $0<r_-<r_+$ and $2T\le r_+-r_-$. Its two moving boundaries are null for the reduced radial wave operator; past characteristic lines remain within the domain. The axis $r=0$ is excluded to avoid a singular radial coefficient.

## 2. Constructive a posteriori theorem with nonlinear coupling

Suppose smooth approximate fields $(\psi_a,\gamma_a)$ are given on $D_T$. Define their **continuum** reduced-system defects

$$
\rho=\psi_{a,tt}-\psi_{a,rr}-r^{-1}\psi_{a,r},\quad
q_t=\gamma_{a,t}-2r\psi_{a,t}\psi_{a,r},\quad
q_r=\gamma_{a,r}-r(\psi_{a,t}^2+\psi_{a,r}^2).
$$

Assume on each slab $D_j$, with $h_j=t_{j+1}-t_j<r_-$,

$$
|\rho|\le\delta_j,\quad |q_t|\le q_j,\quad
|\psi_{a,t}|,|\psi_{a,r}|\le B_j.
$$

The prescribed smooth *exact* initial $\psi$ and $\psi_t$ are within
$E_0$ in function value and $P_0$ in both null derivatives of $\psi_a$
on the initial annulus. Fix smooth initial $\gamma(0,r)$ within $G_0$
of $\gamma_a(0,r)$, satisfying the exact radial Einstein constraint

$$
\partial_r\gamma(0,r)
=r\big(\psi_t(0,r)^2+\psi_r(0,r)^2\big).
$$

**Theorem.** There exists a genuine smooth Ricci-flat four-dimensional Einstein–Rosen metric with those initial data on $D_T$. Put
$w=\psi-\psi_a$, $h=\gamma-\gamma_a$, and on each slab define

$$
Q_j=\frac{P_j+h_j\delta_j}{1-h_j/r_-},\qquad
P_{j+1}=Q_j,\qquad E_{j+1}=E_j+h_jQ_j,
$$

$$
G_{j+1}
=G_j+h_j\left[q_j+2r_+(2B_j+Q_j)Q_j\right].
$$

Then the exact solution differs from its approximation by

$$
\|w(t_{j+1})\|_\infty\le E_{j+1},\quad
\max_\pm\|(\partial_t\pm\partial_r)w(t_{j+1})\|_\infty\le P_{j+1},
\quad
\|h(t_{j+1})\|_\infty\le G_{j+1}.
$$

The same bounds, with the corresponding partial-slab recurrences, hold on the interior of each slab. All recurrences are **algebraic rational operations** on verifiable continuum bounds. If the approximate metric coefficients admit upper envelopes $A\ge e^{2(\gamma_a-\psi_a)}$, $Z\ge e^{2\psi_a}$, then with $\xi=2(E+G)<1$, $\eta=2E<1$,

$$
|g_{tt}-g_{a,tt}|,\ |g_{rr}-g_{a,rr}|
 \le A\frac{\xi}{1-\xi},
\qquad
|g_{zz}-g_{a,zz}|\le Z\frac{\eta}{1-\eta},
$$

$$
|g_{\varphi\varphi}-g_{a,\varphi\varphi}|
\le r_+^2 e^{-2\inf\psi_a}\frac{\eta}{1-\eta}.
$$

These are coordinate-component error bounds on $D_T$, not coordinate-invariant spacetime distances.

**Proof.** The linear radial-wave Cauchy problem is smooth on this
causal diamond because $r\ge r_->0$. For
$p_\pm=w_t\pm w_r$, subtraction gives

$$
(\partial_t\mp\partial_r)p_\pm=-\rho+r^{-1}w_r.
$$

Since $|w_r|\le \max(|p_+|,|p_-|)$, the supremum $Q_j^*$ on a whole
slab satisfies
$Q_j^*\le P_j+h_j(\delta_j+Q_j^*/r_-)$.
Algebra yields $Q_j^*\le Q_j$. Integrating
$w_t=(p_++p_-)/2$ at fixed $r$ gives the formula for $E$.

Construct the true $\gamma$ by integrating
$\gamma_t=2r\psi_t\psi_r$ at fixed $r$ from the constrained initial data.
The identity

$$
\partial_t\left[\gamma_r-r(\psi_t^2+\psi_r^2)\right]
=2r\psi_t\left(\psi_{rr}+r^{-1}\psi_r-\psi_{tt}\right)=0
$$

proves that the exact radial constraint propagates. Thus the full 4D
Ricci tensor vanishes. For approximate and true fields,

$$
|\psi_t\psi_r-\psi_{a,t}\psi_{a,r}|
\le (2B_j+Q_j)Q_j,
$$

and subtraction of the time constraints gives
$|\partial_t(\gamma-\gamma_a)|
\le q_j+2r_+(2B_j+Q_j)Q_j$. Integrating at fixed $r$
derives $G$. Finally, differences of each exponential metric coefficient
are bounded by $e^\zeta-1\le \zeta/(1-\zeta)$ for $0\le\zeta<1$,
giving the displayed metric estimates. QED.

**Important qualification:** all these steps exploit an exact
cylindrical symmetry reduction, a **scalar** wave principal part and
integrable quadratic constraints. They cannot be imported unchanged
into a first-order 50-field GH Einstein evolution. In particular,
the gamma-constraint propagation shown here is *not* the general
GH subsidiary system.

## 3. Certified nonflat rational test

Let $a=1/10$, $\epsilon=1/100$, $T=1/4$, $r_-=1$, $r_+=2$, split
into two slabs of width $1/8$. Choose

$$
\psi_a=(a+\epsilon)t^2+\tfrac12ar^2,\qquad
\gamma_a=2(a+\epsilon)^2t^2r^2+\tfrac14a^2r^4.
$$

The **independent exact vacuum solution** matching the initial data is

$$
\psi=a(t^2+r^2/2),\qquad
\gamma=a^2(2t^2r^2+r^4/4).
$$

The analytic residual identities over the entire (uncountable) domain are

$$
\rho=2\epsilon=\tfrac1{50},\quad
q_r=0,\quad
q_t=4\epsilon(a+\epsilon)tr^2,\quad
|q_t|\le\tfrac{11}{2500}.
$$

Uniform derivative majorant $B=1/5$ and the two-slab rational
recurrences yield

$$
P_1=\tfrac1{350},\ E_1=\tfrac1{2800},\
P_2=\tfrac3{490},\ E_2=\tfrac{11}{9800},\
G_2=\tfrac{70079}{24010000}.
$$

Bounding all exponentials with rational inequalities gives the
**componentwise metric error certificate**

$$
\boxed{
\|g-g_a\|_{\text{entrywise }\infty}
\le \frac{44}{4889}<\frac1{100}.
}
$$

The independent exact errors satisfy
$\|\psi-\psi_a\|_\infty\le1/1600$ and
$\|\gamma-\gamma_a\|_\infty\le21/20000$, safely inside
the certified bounds. All quoted numbers are **exact fractions**;
no floating samples are substituted for a continuum enclosure.

The metric is **not flat**: independent four-dimensional symbolic
Riemann calculation gives at $t=0,r=1,a=1/10$

$$
R^z{}_{tzt}=-209/1000\ne0.
$$

Full symbolic Ricci components vanish for the exact comparison
metric. The proof does not depend on knowing that comparator:
its explicit formula serves as an independent regression oracle.

## 4. Exact reproduction and limits

Run:

    python -m research.calc.er_slab_certificate
    python -m pytest -q tests/test_gh_slab_certificate.py tests/test_er_slab_certificate.py

Baseline run on Python 3.13 / SymPy 1.14: **9 passed**.

These tests certify a **nonflat nonlinear constraint propagation
subsector** and improve the previous flat gauge-wave benchmark.
They do **not** establish a new previously unknown vacuum
solution, verify arbitrary nonzero curvature backgrounds, or
increase the Brill–Lindquist metric's certified evolution time.

**Missing full-Einstein theorem:** a strong-field
$H^s$ error-tube bootstrap for first-order GH with an evaluated
positive symmetrizer, commuted-energy inequalities, reduction
and gauge constraints, and certified excision/incoming boundary
fluxes. This demands genuine nonlinear PDE analysis rather
than the exact symmetry reduction here.

Independent curvature regression also evaluates R^z_{tzt} directly from the four-dimensional Christoffel symbols, separately from the closed-form Riemann expression.
