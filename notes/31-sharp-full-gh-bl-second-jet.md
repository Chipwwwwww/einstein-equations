# Note 31 — Sharper full-50-field BL GH continuation and dynamic resolution

**2026-11-02.** Parent PR #11 head: 9ff9dfe1ef1b3ea1c6b3f4ba632c84ac8687051f.  
**Scope:** Certified local full-GH vacuum time slab; no merger or global exterior. Methods are applications of known Kato symmetric-hyperbolic, Sobolev, Bernstein and characteristic finite-propagation estimates. No mathematical-priority novelty is claimed.

## 1. Independent audit: the PR #11 second jet was incomplete

Use signature (-,+,+,+), undamped gamma1=-1, gamma2=gamma0=0, fixed harmonic source H_a=0. The full GH fields are ten g_ab, ten Pi_ab and thirty Phi_iab; Pi=-n^c partial_c g, Phi_i=partial_i g. Write G(u)=F(u)-A^i(u)partial_i u for the entire nonlinear 50-field first-order system (Lindblom et al. 2006, equations 35–39).

At the initial conformally flat time-symmetric slice, set p>0, p_i=partial_i p, p_ij=partial_ij p, D=sum_i p_i² and L=Delta p, with initial lapse N=p^-2, shift=0, Pi=0 and Phi_i=partial_i g. Full initial jets u1=G(u0), u2=DG(u0)[u1] obey

\[
g_1=0,\quad \Phi_{1,i}=0,\quad
(\Pi_1)_{00}=4D/p^{12}-4L/p^{11},\quad
(\Pi_1)_{0i}=0,\quad
(\Pi_1)_{ij}=4(D\delta_{ij}-4p_ip_j)/p^4-4L\delta_{ij}/p^3,
\]

\[
\boxed{g_2=-p^{-2}\Pi_1,\qquad
(\Pi_2)_{0i}=\frac{16p_i(pL+5D)}{p^{13}},
\qquad(\Pi_2)_{00}=(\Pi_2)_{ij}=0,\qquad
\Phi_{2,i}=\partial_i g_2.}
\tag{31.1}
\]

The complete Phi2 formulas require p_ij and partial_i L; these are implemented exactly in the source. For physical BL, L=0 and (Pi2)_0i=80D p_i/p^13, generally nonzero. For the artificial extension L need not vanish. An independent full four-dimensional Christoffel-variation calculation with all three p_i nonzero and L nonzero agrees in all ten symmetric Pi2 components. The prior PR #11 set Pi2=Phi2=0, so its polynomial was **not** a complete second-order 50-field Taylor polynomial. This need not disprove its conservative O(t) residual bound: any first-jet-consistent polynomial can have O(t) residual.

## 2. Quantitative finite-C9 Brill–Lindquist auxiliary extension

Equal-mass d=3831/5000 in total-mass M=1 units. Let r_+/- denote distances to (0,0,+/-3831/10000). Preserve
\[
\Omega_{\rm phys}=\{|x|<1.1,\ r_\pm>0.2\},
\quad\Omega_*=\{|x|\le0.9,\ r_\pm\ge0.3\},
\]
with Euclidean initial-data buffer at least b=0.1.

Replace PR #11's flat exponential smooth cutoff by the finite-order polynomial
\[
S_9(z)=\frac{19!}{(9!)^2}\int_0^z t^9(1-t)^9dt,\qquad0\le z\le1,
\]
extended as 0 to the left and 1 to the right. The extended scalar is C9, **not C-infinity**, because the first nine derivatives vanish at each glue interface while the tenth does not. The piecewise positive auxiliary conformal factor is
\[
p_E=1+\frac14\chi_{\rm out}(|x|^2)
\sum_{\sigma=\pm}\chi_{\rm in}(|x-x_\sigma|^2)/r_\sigma.
\]
Each puncture summand is DEFINED as identically zero for r_sigma <=0.12, so no undefined 0/0 is introduced. Inner transition occurs at 0.12 to 0.18, and outer transition at 1.1 to 1.5. Therefore p_E equals the actual BL factor exactly throughout Omega_phys, equals 1 near every torus seam, is strictly positive, and gives C9 global auxiliary GH initial fields on T3. They are generally NOT Einstein vacuum outside Omega_phys.

**Lemma 31.2 (Bernstein derivative bound).** If a polynomial Q(z)=sum_j a_j z^j of degree d has Bernstein coefficients b_i=sum_{j<=i}a_j binom(i,j)/binom(d,j), then Q(z)=sum_i b_i binom(d,i)z^i(1-z)^(d-i), so |Q(z)|<=max_i |b_i| for all 0<=z<=1. Apply to each derivative of S9, then use the exact multivariate Faà di Bruno coefficients for the quadratic radial arguments, the elementary bound |partial^alpha(1/r)|<=6^n n! r^(-n-1) on r>=0.12, and convolution of exponential generating functions for p_E^4, p_E^-4 and p_E^-2. Since Phi_i=partial_i g, one more spatial derivative is included. Exact rational bounds, not grid samples, give
\[
\boxed{\|u_0\|_{H^s(T^3;R^{50})}<10^{U_s},
\quad(U_3,U_4,U_5,U_6,U_7)=(22,27,32,37,42).}
\tag{31.2}
\]
This proof uses actual derivatives at each order, not PR #11's single largest derivative bound raised to the eighth power. C9 suffices to produce initial H7 data; smooth Kato approximants and high-regularity persistence can be used for the H3 energy proof. Implemented in the exact-rational Bernstein and composition recurrences.

## 3. Full H3 source, Jacobian, Hessian, symmetrizer and transport bounds

Use the same full 50 GH system as PR #10 but center Neumann inversion at the variable conformal reference g0=diag(-p_E^-4,p_E^4,p_E^4,p_E^4), not Minkowski. Its inverse and initial lapse p_E^-2 have known Hs bounds from (31.2), and the full error tube is set to R=10^-320. The matrix-inverse Neumann series ratio, normalized binomial inverse-square-root series for lapse, Lorentzian inertia, and spatial-metric positivity all close on the whole infinite-dimensional H3 ball. For the known initial H4/H5 fields, the same graph is evaluated at the reference state, not by assuming unknown H4/H5 errors.

Apply the H3 Banach product bound 48 and the H4/H5 product bound 192<10^3. A symbolic logarithmic majorant node stores the value, first and second Frechet derivative norms. Product rules use conservative exponents n_ab <= n_a+n_b+3, l_ab <= max(n_a+l_b,l_a+n_b)+5, h_ab <= max(n_a+h_b,h_a+n_b,l_a+l_b)+7, with extra integer margin for sums. Explicit Neumann/binomial Hessian bounds close inversion and lapse. The full 50 GH source graph gives the following **upper** exponents on the entire H3 tube:
\[
\begin{array}{c|rrr}
 & \text{H3 value}&\text{Frechet 1}&\text{Frechet 2}\\ \hline
F & 368&547&724\\
A^i &148&318&487
\end{array}.
\tag{31.3}
\]
The first and second Frechet bounds apply to the algebraic maps from H3 fields to H3 outputs; the spatial derivative in G=F-A^i partial_i u is accounted for separately with known H4 jets.


**Independent-source falsification and repair.** An independently assembled graph through the actual production 50-field source found that a hand-duplicated contraction graph underestimated value, Jacobian and Hessian exponents by 5, 6 and 7, respectively, before the previous +3 final margin. The certifier was changed to call the actual full source with rigorous geometry majorant objects, rather than duplicating the contractions. All exponents below are from that corrected graph. This correction matters to proof integrity; the superseded 10^(-2832) time is withdrawn.

In the GH energy, S(g) and S(g)A^i(g) depend **only on g**, not directly on Pi/Phi. Therefore
\[
\partial_t S=DS(g)[g_t],\
\partial_i(SA^i)=D(SA^i)(g)[\partial_i g].
\]
To enclose the *pointwise* finite-dimensional matrix derivatives, use the stronger pointwise reference bound U3+2=24 rather than silently identifying H3 derivative bounds with pointwise DS. The resulting pointwise graph exponents are JS<=304 and JB<=401; this gives metric-sector energy divergence D<10^430. With source/transport G<10^551, symmetrizer eigen bounds lambda>=10^-4, Lambda<=10^4, one obtains
\[
\boxed{K<10^{560}}.
\tag{31.4}
\]
No unknown exact-solution derivative appears outside the strictly verified H3 tube. These are intentionally loose exact-arithmetic upper bounds, not fitted data.

## 4. Complete second jet yields continuum O(t²) residual

The full 50-field known coefficients satisfy
\[
\|u_1\|_{H^{3,4,5}}\le(10^{372},10^{427},10^{482}),
\quad
\|u_2\|_{H^{3,4}}\le(10^{923},10^{1058}).
\tag{31.5}
\]
Define v2(t)=u0+t u1+(t²/2)u2 with u1=G(u0), u2=DG(u0)[u1]. Put delta=v2-u0. Banach-space Taylor's integral formula, in the correct H4 to H3 mapping, gives
\[
r_2(t)=v_{2,t}(t)-G(v_2(t))
=-\tfrac12t^2 DG(u_0)[u_2]
-\int_0^1(1-\theta)D^2G(u_0+\theta\delta)
[\delta,\delta]\,d\theta.
\tag{31.6}
\]
The complete second jet eliminates ALL O(t) residual contributions, including Pi and Phi. The H3 source Hessian, principal derivative/commutator and known H4 jets yield on the ENTIRE continuum torus and slab
\[
\boxed{\|r_2(t)\|_{H^3}\le10^{2856}t^2}.
\tag{31.7}
\]
The auxiliary jet includes nonzero Delta p_E and grad Delta p_E outside Omega_phys. It would be invalid to use the physical harmonic p formula globally.

## 5. Closed nonlinear existence and a dynamically resolved metric

Choose
\[
T=10^{-2876}M,\qquad R=10^{-320},\qquad K<10^{560}.
\]
The exact majorant code verifies that v2 remains in the reference tube, the GH coefficient/symmetrizer bounds are valid throughout it, KT<1/2, and all eight closure/falsification inequalities are strict. Kato's H3 local existence continues by the symmetrized energy first-exit argument. The residual at T is <=10^-2896. Writing the error energy as Y, one has
\[
Y'\le KY+\frac{\Lambda}{\sqrt\lambda}\|r_2\|_{H^3},
\quad Y(0)=0.
\]
The integrating-factor estimate and lambda>=10^-4, Lambda<=10^4 yield
\[
\boxed{\|u-v_2\|_{L^\infty([0,T];H^3)}
 <10^{-5762}<R.}
\tag{31.8}
\]
For metric components the pointwise error is <6*10^-5762, first spatial derivative error <12*10^-5762. No pointwise curvature error has been established.

At the concrete protected point x_*=(0.1,0,0), set r=sqrt(0.1²+0.3831²). Since 0.39<r<0.4 exactly,
\[
(g_2)_{xx}=12p_x^2/p^6>1/20,
\qquad g^{app}_{xx}(T,x_*)-g_{xx}(0,x_*)>T²/40.
\]
The metric change divided by its certified error bound therefore exceeds 10^10/240>4*10^7. A **coordinate** dynamical change is resolved, but no gauge-invariant curvature or radiation signal is certified.

## 6. Protected vacuum and the unchanged geometric obstruction

The exact BL initial constraints, GH gauge and first-order reduction constraints vanish on Omega_phys. The normal gauge subsidiary datum vanishes by the reduced full GH equations and the exact Hamiltonian/momentum constraints; independently checked for the conformal initial values. The original homogeneous GH subsidiary equations propagate zero constraints locally. Direct reference-metric and tiny-tube bounds give lapse <2, shift speed <1 and spatial inverse eigen upper <2, implying that every GH and subsidiary coordinate characteristic speed is <4. Since 4T<0.1, all backward causal characteristics from [0,T]xOmega_* meet t=0 within Omega_phys. A localized symmetric-hyperbolic energy flux with a contracting domain has sign -CS-SA^n<=0 for C above the characteristic speeds; hence the artificial nonvacuum extension cannot contaminate the protected region. There Ric(g)=0 exactly. The vacuum metric is identified with the original BL development by local geometric uniqueness up to diffeomorphism, not by a proved explicit gauge map.

The new GH time is still vastly smaller than the old geodesic-slicing certified interval (~8.683e-4 M), and no additional geometric spacetime coverage, common MOTS-merger passage, radiative exterior, curvature error or global metric is claimed. The dominant next obstacle is not correctness of the second jet but the remaining 10^547/10^560 full-source energy overestimates and the absence of sharp strong-field localized energy/multislab continuation.

## 7. Reproduction and relation to mathematical literature

Run from repository root:

    python -m research.calc.bl_sharp_cutoff
    python -m research.calc.bl_sharp_certificate
    python -O -m research.calc.bl_sharp_certificate
    python -m research.calc.bl_full_second_jet
    python -m pytest -q tests/test_bl_sharp_local.py

The revised calculator must be judged against explicit rational inequalities, not only unit-test success. Independent symbolic checks of all ten second Pi jets and Phi compatibility, falsification of accelerated time, Bernsteins' convex-hull proof and exact cutoff matching are included. The earlier local prototype passed seven tests, and its independent Christoffel extension adds another. Remote CI and external mathematical peer review remain to be checked.

Primary sources: Kato (1975), Hughes-Kato-Marsden (1977), Kato-Ponce (1988), Friedrich-Rendall (2000) on domains of dependence, and Lindblom-Scheel-Kidder-Owen-Rinne (2006), especially full GH equations (35)–(39) and homogeneous subsidiary system (40)–(52). The PDE and polynomial tools are standard; the verified numerical constants and complete BL GH second jet are project-specific contributions, **not** original PDE literature discoveries.
