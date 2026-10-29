# Literature baseline

Bounded on purpose: enough to know what is established and what a result here
would add. "Checked" says what was actually read this session. Statements
marked *not re-checked* come from memory of standard texts and should be
verified before being relied on in a proof.

## Initial-value problem

| Source | What it gives | Checked |
| --- | --- | --- |
| Y. Fourès-Bruhat (Choquet-Bruhat), Acta Math. 88 (1952) 141–225 | local existence and uniqueness for vacuum Einstein in harmonic gauge | bibliographic only |
| Y. Choquet-Bruhat, R. Geroch, *Global aspects of the Cauchy problem in general relativity*, Comm. Math. Phys. 14 (1969) 329–335 | existence and uniqueness (up to isometry) of the maximal globally hyperbolic development | citation verified (journal page); theorem statement from secondary sources (Sbierski 2015) |
| R. Bartnik, *The mass of an asymptotically flat manifold*, CPAM 39 (1986) 661–693 | ADM mass well defined for decay $q > 1/2$ | *not re-checked* |
| R. Penrose, PRL 14 (1965) 57 | trapped surface + NEC + non-compact Cauchy surface ⇒ null incompleteness | *not re-checked* |

## Horizons

| Source | What it gives | Checked |
| --- | --- | --- |
| S. Hawking, G. Ellis, *The Large Scale Structure of Space-Time* (1973), §9.2 | trapped surfaces lie in the black hole region (Prop. 9.2.1); black holes do not bifurcate | *not re-checked* (proposition numbers from memory) |
| R. Wald, *General Relativity* (1984), §12.2 | same, Prop. 12.2.2 | *not re-checked* |
| R. Wald, V. Iyer, PRD 44 (1991) R3719 | Schwarzschild has Cauchy slices with no trapped surfaces | *not re-checked* |
| L. Andersson, J. Metzger, *The area of horizons and the trapped region*, CMP 290 (2009) 941–972, arXiv:0708.4252 | boundary of the trapped region is a smooth stable MOTS; existence of a stable MOTS between an outer-trapped inner boundary and an outer-untrapped outer boundary | abstract and secondary summary; **precise hypotheses of the barrier theorem to be re-read before use** |
| L. Andersson, M. Mars, J. Metzger, W. Simon, *The time evolution of marginally trapped surfaces*, CQG 26 (2009) 085018, arXiv:0811.4721 | MOTS persist while strictly stable; at zero principal eigenvalue, with $W = \lvert\chi^+\rvert^2 + \mathrm{Ric}(\ell^+,\ell^+) \not\equiv 0$, $D_t\vartheta = -\int\alpha W\psi < 0$ and the MOTT is tangent to the slice and lies to its future (Prop. 5.1, Thm 5.2); coalescence Thm 4.1 | read (ar5iv), statements and formula extracted |
| D. Pook-Kolb et al., *The interior of a binary black hole merger*, PRL 123 (2019), arXiv:1903.05626; follow-ups arXiv:2104.11344 | numerical: common MOTS bifurcates into inner and outer branches in head-on BL evolutions; self-intersecting MOTS | abstracts |

## Exact merger benchmark

| Source | What it gives | Checked |
| --- | --- | --- |
| D. Kastor, J. Traschen, *Cosmological multi-black-hole solutions*, PRD 47 (1993) 5370, arXiv:hep-th/9212035 | metric $-\Omega^{-2}dt^2 + a^2\Omega^2 d\vec x^2$, $A_t = \Omega^{-1}$, $H = \pm\sqrt{\Lambda/3}$; $H < 0$ "contracting": black holes merge; $\lvert Q\rvert = M$ | read (ar5iv): metric, potential, branch convention |
| D. Brill, G. Horowitz, D. Kastor, J. Traschen, *Testing cosmic censorship with black hole collisions*, PRD 49 (1994) 840, arXiv:gr-qc/9307014 | $U = H\tau + \sum M_i/r_i$ form; $M_{\rm ext} = 1/(4\lvert H\rvert)$; coalescence for $M < M_{\rm ext}$ via analytic null-geodesic bounds (equal masses); naked singularity beyond a Cauchy horizon for $M > M_{\rm ext}$; extension non-unique, $C^2$ / $C^1$ smoothness | read (ar5iv), statements extracted; proofs **not re-derived** |
| K. Nakao, T. Shiromizu, S. Hayward, *Horizons of the Kastor–Traschen multi-black-hole cosmos*, PRD 52 (1995) 796 | marginal surfaces for two equal masses | **could not access**; citation from secondary sources (McNutt–Coley 2018, arXiv:1811.02931) |

## Extreme mass ratio

| Source | What it gives | Checked |
| --- | --- | --- |
| R. Emparan, M. Martínez, *Exact event horizon of a black hole merger*, CQG 33 (2016) 155003, arXiv:1603.00712 | in the EMR limit the merger event horizon is the congruence of null geodesics in Schwarzschild that approach a null plane at infinity; elliptic functions | abstract |
| R. Emparan, M. Martínez, M. Zilhão, *Black hole fusion in the extreme mass ratio limit*, PRD 97 (2018) 044004, arXiv:1708.08868 | spinning small hole, arbitrary orientation; transient toroidal horizons; critical exponents | abstract |
| P. Hintz, *Gluing small black holes along timelike geodesics I–III*, arXiv:2306.07409, 2408.06712, 2408.06715 | given a vacuum (±Λ) globally hyperbolic spacetime and a timelike geodesic, a family $g_\epsilon$ of true solutions with a mass-$\epsilon m$ subextremal Kerr hole along it, on compact sets; **Thm 6.1** (KdS) and **Rem. 6.2** (Kerr, via KS23/GKS22/She23 stability): EMR merger followed by relaxation to Kerr(–dS). Leaves open the event horizon and the structure of the black hole region ("conjecturally connected"). | III: §1 and §6 read in full (PDF) |

## Stability, formation, gluing (context; read only as needed)

| Source | What it gives | Checked |
| --- | --- | --- |
| D. Christodoulou, S. Klainerman, *The global nonlinear stability of the Minkowski space* (1993) | small data: complete, no black hole | *not re-checked* |
| D. Christodoulou, *The formation of black holes in general relativity* (EMS 2009) | trapped surfaces form from regular vacuum data (short pulse) | *not re-checked* |
| P. Hintz, A. Vasy, Acta Math. 220 (2018) | nonlinear stability of slowly rotating Kerr–de Sitter | cited by Hintz III |
| S. Klainerman, J. Szeftel, *Kerr stability for small angular momentum*, PAMQ 19 (2023); E. Giorgi, S. Klainerman, J. Szeftel, PAMQ 20 (2024) 2865–3849, arXiv:2205.14808; D. Shen (She23) | nonlinear stability of slowly rotating Kerr | GKS citation verified; theorem statements not read |
| M. Dafermos, G. Holzegel, I. Rodnianski, M. Taylor, arXiv:2104.08222 (2021) | Schwarzschild stability, codimension 3 | *not re-checked* |
| G. Neugebauer, J. Hennig, arXiv:0905.4179; J. Hennig, G. Neugebauer, arXiv:1103.5248; review arXiv:1105.5830 | no stationary, axisymmetric vacuum two-black-hole configurations with regular axis (double-Kerr–NUT) | abstracts |
| F. Pretorius, PRL 95 (2005) 121101; M. Campanelli et al., PRL 96 (2006) 111101; J. Baker et al., PRL 96 (2006) 111102 | numerical binary black-hole evolutions through merger | *not re-checked* |

## Not covered

Gluing constructions of initial data (Corvino, Chruściel–Delay,
Carlotto–Schoen, Mao–Oh–Tao, Kehle–Unger), Robinson–Trautman convergence to
Schwarzschild (Chruściel), and higher-dimensional exact coalescence (e.g. on
Eguchi–Hanson base spaces) were not reviewed this session. They become
relevant if the hypotheses in notes/09 are pursued further.
