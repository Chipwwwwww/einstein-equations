# 25 — Explicit H^3 a posteriori existence/continuation for coupled symmetric-hyperbolic PDEs

**2026-11-02; branch research/full-gh-hs-validation.**
**Status: theorem proved (standard symmetric-hyperbolic + explicit Fourier estimates, not claimed mathematically novel), exact-rational 2-field 3D nonlinear certificate completed.** A 50-field generalized-harmonic Einstein principal-symbol test is separate. **No full 50-field Einstein source/residual certificate, no vacuum Einstein slab, no Brill–Lindquist extension, and no merger solution are certified.**

## 1. Exact norm, elementary constants and derivative-loss repair

On the flat three-torus $(\mathbb R/2\pi\mathbb Z)^3$, use normalized Fourier coefficients $\widehat f(k)=(2\pi)^{-3}\int f e^{-ik\cdot x}\,dx$ and
\[
\|f\|_3^2=\sum_{k\in\mathbb Z^3}(1+|k|^2)^3|\widehat f(k)|^2.
\]
For vector/matrix fields use the Euclidean/Frobenius fiber norm; inequalities between matrices use operator norm. We write $\Lambda^3=(1-\Delta)^{3/2}$, with $H^3=\{f:\|f\|_3<\infty\}$.

**Lemma 25.1 (explicit constants).**
\[
\|f\|_{L^\infty}\le6\|f\|_3,\quad
\|\nabla f\|_{L^\infty}\le12\|f\|_3,\quad
\|fg\|_3\le48\|f\|_3\|g\|_3,
\]
and for $3\times3$ spatial transport directions, any matrix coefficient $A$ and vector field $e$,
\[
\|[\Lambda^3,A]\partial_i e\|_{L^2_{\rm norm}}
   \le144\|A-C\|_3\|e\|_3
\quad (i=1,2,3),
\tag{25.1}
\]
where $C=C(t)$ is any matrix independent of the spatial variables. The constants are conservative and **universal**, independent of the number of field components.

**Proof.** The sup Fourier embedding constants are bounded by
\[
C_0^2=\sum_k(1+|k|^2)^{-3}<36,\quad
C_1^2=\sum_k |k|^2(1+|k|^2)^{-3}<144.
\]
For the cubes $|k|_\infty=n\ge1$, the number of integer lattice points is $(2n+1)^3-(2n-1)^3=24n^2+2$ and $n^2\le|k|^2\le3n^2$. Hence
\[
C_0^2\le1+24\zeta(4)+2\zeta(6)<1+24(4/3)+2(6/5)<36,
\]
\[
C_1^2\le72\zeta(2)+6\zeta(4)
\le72(7/4)+6(4/3)=134<144.
\]
(The elementary bounds are from integral comparison, with $\zeta(2)\le7/4$ after separating the first two terms.) Weighted convolution and
\[
(1+|p+q|^2)^{3/2}\le4[(1+|p|^2)^{3/2}+(1+|q|^2)^{3/2}]
\]
give the product bound $2^3 C_0=48$ by Young's inequality.

For the commutator, put $w_s(k)=(1+|k|^2)^{s/2}$. The mean value theorem yields
\[
|w_s(p+q)-w_s(p)|
\le s\,2^{s-2}|q|\,[w_{s-1}(p)+w_{s-1}(q)]\quad(s\ge2).
\]
Multiplying by $|p_i|\,|\widehat A(q)|\,|\widehat e(p)|$ splits the convolution into
(i) $\||q|\widehat A(q)\|_{\ell^1}\||p_i|w_{s-1}(p)\widehat e(p)\|_{\ell^2}$ and
(ii) $\||q|w_{s-1}(q)\widehat A(q)\|_{\ell^2}\||p_i|\widehat e(p)\|_{\ell^1}$.
Each is at most $C_1\|A\|_s\|e\|_s$. At $s=3$, the total constant is $2s\,2^{s-2}C_1\le144$. A spatially constant $C(t)$ commutes with $\Lambda^3\partial_i$. All matrix/vector estimates use only submultiplicativity of the Frobenius norm. **No fourth spatial derivative of the exact error is needed.** QED.

**Contrast with PR #8.** A naïve inverse wave estimate $L^\infty H^{s-1}\to L^\infty H^{s+1}$ is unbounded by the resonant Fourier example there. The commutator here is of transport type in a *symmetric first-order* formulation. No inverse of the wave operator in the invalid fixed-scale spaces is used.

## 2. Fully coupled system and a computable mathematical certificate

Let $u(t,x)\in\mathbb R^N$ solve the **general** nonlinear symmetric-hyperbolic system
\[
\partial_t u+\sum_{i=1}^3A_i(u)\partial_i u=F(t,x,u),\qquad u(0)=u_0.
\tag{25.2}
\]
Here $A_i:\Omega\to\mathbb R^{N\times N}$ and a **positive symmetric symmetrizer** $S:\Omega\to\mathrm{Sym}_N^+$ are smooth on an open state set $\Omega\subset\mathbb R^N$, and $B_i(w)=S(w)A_i(w)$ are symmetric for all admissible $w$. The source $F(t,x,w)$ is sufficiently smooth and periodic in $x$. This includes arbitrarily coupled $N$-component systems. $S$ and $A_i$ are assumed to depend on the state, not explicitly on $(t,x)$; the source may depend explicitly on $(t,x)$ as in GH gauge-source terms. Smooth extensions with $x$-dependent principal coefficients require additional known partial-derivative bounds.

Let $v$ be a given smooth approximate solution in $C^1([0,T];H^4)$ and
\[
r=v_t+\sum_iA_i(v)\partial_i v-F(t,x,v),\quad
\rho\ge\sup_{0\le t\le T}\|r(t)\|_3,\quad
\epsilon_0\ge\|u_0-v(0)\|_3.
\]
Pick a proposed radius $R>0$. The **entire** H^3 tube is
$\mathcal U_t=\{w\in H^3:\|w-v(t)\|_3\le R,\ \operatorname{range}(w)\subset\Omega\}$.
It must be verified *before* invoking this theorem that the geometric
range and pointwise positivity conditions hold uniformly in the tube, not just for the approximate field.

The following constants must be guaranteed **for every** $t\in[0,T]$ and every $w\in\mathcal U_t$ (these are input proof obligations, not unknown true-solution bounds):
\[
\begin{array}{ll}
\lambda I\preceq S(w(x))\preceq\Lambda I,\quad \lambda>0;
&\|A_i(w)\|_3\le M_{A_i},\quad\|F(t,\cdot,w)\|_3\le M_F,\\
\|A_i(w)-C_i(t)\|_3\le M_{c_i}
& C_i(t)\ {\rm spatially\ constant},\\
\|D S(z)\|_{\mathrm{op}}\le J_S,\quad\|D(SA_i)(z)\|_{\mathrm{op}}\le J_{B_i}
&z=w(x),\\
\|A_i(w)-A_i(v)\|_3\le L_{A_i}\|w-v\|_3,
&\|F(t,\cdot,w)-F(t,\cdot,v)\|_3\le L_F\|w-v\|_3 .
\end{array}
\tag{25.3}
\]
The finite-dimensional derivative norms in the third row are operator norms from Euclidean vectors to matrices in operator norm. Define
\[
V=\sup_t\|v(t)\|_3,\quad V_i=\sup_t\|\partial_i v(t)\|_3,
\]
and the **explicit** constants
\[
\begin{aligned}
U_t^*&=6M_F+72(V+R)\sum_iM_{A_i},\\
D&=J_SU_t^*+12(V+R)\sum_iJ_{B_i},\\
G&=L_F+48\sum_iL_{A_i}V_i+144\sum_iM_{c_i},\\
K&=\frac D{2\lambda}+\frac{\Lambda}{\lambda}G.
\end{aligned}\tag{25.4}
\]
The meaning of $U_t^*$ is a **tube-computed** bound on $\|u_t\|_\infty$ *if* an exact solution enters the tube; it is never read from an unknown solution.

**Theorem 25.2 (closed nonlinear a posteriori existence and uniqueness).**
Assume $u_0$ is smooth and belongs to the admissible state range, that the validated tube covers $v(0)$ and the exact initial state, and the classical symmetric-hyperbolic local existence/continuation theorem holds for (25.2) (which it does for smooth coefficients on an open hyperbolicity cone with the above smoothness and compact torus). If
\[
\boxed{
\mathcal B_T=
\frac{e^{KT}}{\sqrt\lambda}
\left(\sqrt\Lambda\,\epsilon_0+
\frac{\Lambda}{\sqrt\lambda}\,\rho T\right)
<R,
}\tag{25.5}
\]
then the **actual nonlinear solution** $u$ of (25.2) exists uniquely on all of $[0,T]$, stays within the tube, and
\[
\|u(t)-v(t)\|_3\le\mathcal B_T,\qquad 0\le t\le T.
\]
In particular, the coefficient hyperbolicity and all pointwise geometric tube inequalities persist. For $KT<1$, a fully rational computable sufficient condition substitutes $(1-KT)^{-1}$ for $e^{KT}$. A sequence of finite slabs may instead propagate the certified energy $Y_j$ by
\[
Y_{j+1}\le\frac{Y_j+h_j(\Lambda/\sqrt\lambda)\rho_j}{1-K_jh_j},
\quad K_jh_j<1,
\tag{25.6}
\]
provided the energy norms at the joins are compatible; when the symmetrizer changes, one must re-enclose $Y_{j+1}$ using the relevant change-of-energy factor. All the estimates in (25.5)–(25.6) are valid only for source/initial data satisfying the stated conditions.

**Proof.** Standard Kato theory supplies a local exact solution; this is the sole nonconstructive starting point, and does not assume existence to time $T$. Let $e=u-v$. Subtracting the equations *before* commuting yields
\[
e_t+\sum_iA_i(u)\partial_i e
=-r+[F(t,x,u)-F(t,x,v)]-
\sum_i[A_i(u)-A_i(v)]\,\partial_i v.
\tag{25.7}
\]
Assume temporarily $\|e\|_3<R$ until first exit. Put $z=\Lambda^3e$ and $\mathcal E=\int z^\top S(u)z\,dx/(2\pi)^3$. The exactly symmetric matrices $B_i=S A_i$ allow integration by parts with **no lost derivative**:
\[
\begin{split}
\mathcal E'&=
\int z^\top\left(\partial_t S(u)+\sum_i\partial_iB_i(u)\right)z\\
&\quad+2\langle z,S(u)\Lambda^3[
-r+F(u)-F(v)-\sum_i(A_i(u)-A_i(v))v_i]\rangle\\
&\quad-2\sum_i\langle z,S(u)[\Lambda^3,A_i(u)]\partial_i e\rangle.
\end{split}\tag{25.8}
\]
By Lemma 25.1, $\|\partial_i u\|_\infty\le12(V+R)$ and $\|u_t\|_\infty\le6M_F+72(V+R)\sum_iM_{A_i}=U_t^*$, using the **exact equation only on the bootstrap interval**. Chain rule gives
$\|\partial_tS+\sum_i\partial_iB_i\|_{\mathrm{op},\infty}\le D$.
The two Lipschitz bounds and product estimate control the source difference by $(L_F+48\sum_iL_{A_i}V_i)\|e\|_3$. Lemma 25.1, with each spatial constant $C_i(t)$ subtracted, bounds the commutator by $144\sum_iM_{c_i}\|e\|_3$. Because
$\lambda\|e\|_3^2\le\mathcal E\le\Lambda\|e\|_3^2$,
$Y=\sqrt{\mathcal E}$ satisfies (by a standard $\sqrt{\mathcal E+\delta}$ regularization at zeros)
\[
Y'\le \left(\frac D{2\lambda}+\frac{\Lambda G}{\lambda}\right)Y+
\frac{\Lambda}{\sqrt\lambda}\rho=KY+\frac{\Lambda}{\sqrt\lambda}\rho.
\tag{25.9}
\]
With $Y(0)\le\sqrt\Lambda\,\epsilon_0$, integrating the scalar inequality gives (25.5) for as long as the bootstrap interval exists. The strict bound $\mathcal B_T<R$ contradicts any first exit. In addition, the uniform H^3 bound $\|u\|_3\le V+R$ and strictly interior compact hyperbolicity range exclude breakdown of the local Kato solution before $T$ by its continuation criterion. Thus existence to $T$ follows, not just conditional stability. Uniqueness follows from the standard symmetric-hyperbolic energy argument (or local uniqueness and continuation). The use of a smooth approximation and smooth exact initial data justifies the commutations classically; extension to $H^3$ data follows by Friedrichs regularization and continuous dependence. QED.

**Caveat:** The theorem requires certifying (25.3) **over a full H^3 ball**, not just the pointwise values of $v$. General GH source rational/composition bounds on that ball have not been built or executed for actual merger data. To verify the Lipschitz bounds constructively, a sufficient sub-certificate is
$\|DA_i(v+\theta e)\|_3\le J^3_{A_i}$ for all $\theta\in[0,1]$ and admissible $e$; then $L_{A_i}=48J^3_{A_i}$ by the product lemma and fundamental theorem of calculus (similarly $F$). Polynomial maps admit exact Banach-algebra bounds; inverse metric and lapse square roots require separately verified inverse/positive-domain functional calculus. No sampling of a finite list of candidate states can certify an infinite-dimensional tube.

## 3. Exact-rational, genuinely coupled 3D test

Use $N=2$, $S=I$, $F=0$ on the full 3-torus, with
\[
A_1(u)=\begin{pmatrix}u_1&1/20\\1/20&u_2\end{pmatrix},\qquad
A_2(u)=\begin{pmatrix}u_2&0\\0&-u_1\end{pmatrix},\quad A_3=0.
\]
The matrices are symmetric and do **not** commute in general, so the PDE has nontrivial field coupling and genuinely 3D function spaces. Let exact initial state $u_0=0$ and
\[
v(t,x)=a t(\cos x_1,\sin x_1),\quad
a=10^{-4},\quad T=10^{-2},\quad R=10^{-3}.
\]
The independent exact comparator is $u\equiv0$, but its knowledge is **not** used in the a posteriori theorem.

Every continuum bound is analytic and exact rational:
\[
V=V_1\le3aT=3/10^6,\quad V_2=V_3=0;
\quad M_{A_1}=M_{A_2}=1/10;
\]
\[
\rho\le3a+48(1/10)(3aT)=393/1250000.
\]
The *exact affine* coefficient Lipschitz constants satisfy $L_{A_1}=L_{A_2}=1$; set $J_S=0$, $J_{B_1}=J_{B_2}=1$, $\lambda=\Lambda=1$, and choose spatial constants $C_i=A_i(0)$, giving $M_{c_i}\le V+R$.

Hence $D=3009/125000$, $G=18063/62500$, $K=75261/250000$, and the **closed** rational $H^3$ error bound is
\[
\boxed{\|u-v\|_{L^\infty_tH^3_x}
 \le\frac{393}{124623695}\approx 3.1535\times10^{-6}<R=10^{-3}.}
\tag{25.10}
\]
The independent exact error is $\sqrt8\,aT<3\times10^{-6}$, inside the certificate. The test is not a special scalar equation: the principal matrices include off-diagonal coupling and fail to commute. It is nevertheless a **manufactured** hyperbolic example, not Einstein or a black-hole merger.

A rational multi-slab proof splits the time interval without using unknown solutions; each subslab closes its tube. A zero-but-false residual, missing initial-data compatibility, excessively wide metric/coefficient tube, and failure of the nonlinear inequality are all rejected in negative tests.

Reproduce:
    python -m research.calc.full_hs_bootstrap
    python -m pytest -q tests/test_full_hs_bootstrap.py

## 4. Relation to existing mathematics and what remains

This is an explicit, conservative a posteriori **repackaging** of Kato-type local existence/continuation and standard commutator energy bounds. The Fourier constants are elementary and usable in exact arithmetic; neither the PDE theory nor the commutator lemma is claimed original. Relevant primary sources: Kato (1975, *The Cauchy Problem for Quasi-linear Symmetric Hyperbolic Systems*, Arch. Rational Mech. Anal. 58, 181–205); Hughes–Kato–Marsden (1977, *Well-posed quasi-linear second-order hyperbolic systems*, Arch. Rational Mech. Anal. 63, 273–294); Kato–Ponce (1988, *Commutator estimates and the Euler and Navier-Stokes equations*, CPAM 41, 891–907).

Full Einstein GH suitability is analyzed separately in note26. **It is not legitimate to promote the 2-field coefficient certificate to an evaluated 50-field Einstein source estimate.** A full-GH application must compute validated $H^3$ norms of the 50-field source, its Frechet derivatives and symmetrizer in a genuinely curved 3+1 region, deal with the GH gauge and first-order reduction constraints, supply a boundaryless causal domain or a proved constraint-preserving boundary estimate, and compare charts to the existing BL data. These jobs have not been completed.
