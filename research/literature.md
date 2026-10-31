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
| S. Dain, G. Weinstein, S. Yamada, *A counterexample to a Penrose inequality conjectured by Gibbons*, CQG 28 (2011) 085015, arXiv:1012.4190 | Prop. 1: equal-mass BL data with sufficiently small mass-to-separation ratio contain no closed connected minimal surface enclosing both punctures (Penrose inequality against curvature estimates for stable minimal surfaces; also from Chruściel–Mazzeo, CQG 20 (2003) 729, Thm 3.2). Not quantitative. Says non-existence near the threshold rests on numerics (Brill–Lindquist 1963, Gibbons 1972) | read (PDF text, session 7): statement and proof sketch. Calibrates notes/12: our Theorem F is the explicit, near-critical case |

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
| P. Hintz, *Gluing small black holes along timelike geodesics I–III*, arXiv:2306.07409, 2408.06712, 2408.06715 | given a vacuum (±Λ) globally hyperbolic spacetime and a timelike geodesic, a family $g_\epsilon$ of true solutions with a mass-$\epsilon m$ subextremal Kerr hole along it, on compact sets; **Thm 6.1** (KdS) and **Rem. 6.2** (Kerr, via KS23/GKS22/She23 stability): EMR merger followed by relaxation to Kerr(–dS). Leaves open the event horizon and the structure of the black hole region ("conjecturally connected"). | III: §1 and §6 read in full (PDF); session 4 also read Def. 4.1, eq. (2.19) (coordinates of $\hat g_b$), proof of Thm 6.1, Rem. 6.2 |

## Stability, formation, gluing (context; read only as needed)

| Source | What it gives | Checked |
| --- | --- | --- |
| P. Hintz, *Nonlinear stability of subextremal Kerr black holes*, arXiv:2606.28253 (v2 Aug 2026) | global nonlinear stability of the full subextremal Kerr family, for data close to Kerr (polyhomogeneous decay; e.g. equal to Kerr outside a compact set), in generalized wave map gauge on fixed coordinates: $g = g_b + h$ with decay rates in every asymptotic region. Consequence (3): with Hintz III, Λ = 0 EMR merger spacetimes | session 5: Thm 1.1, Thm 13.1, Rems 1.3, 13.2–13.4, Lemma 5.16, Def. 5.3, (3.31)–(3.37), (5.56), Thm 5.23, Rem. 5.24 **read**; session 7: Def. 11.3 and (11.8), proof of Lemma 11.5 ($h_{00} \in C_b^{k-3,(1+2\gamma^\Upsilon,\dots)}$ on $\Omega^*$), Defs. 11.37, 11.47, 11.55 and Lemma 11.38 (correction terms supported in $\mathrm{supp}\,d\chi_K$, disjoint from $\mathscr I^+$), (11.228)–(11.229), §1.5 ($\chi_K = 0$ near $\mathscr I^+$) **read** |
| X. Chen, S. Klainerman, *Regularity of the future event horizon in perturbations of Kerr*, arXiv:2409.05700 | event horizons of the KS stability spacetimes are smooth null hypersurfaces | abstract |
| P. Hintz, *Horizons of some asymptotically stationary spacetimes*, Duke Math. J. 175 (2026), arXiv:2411.12568 | a unique null hypersurface asymptotic to the Kerr horizon in dynamical spacetimes that approach Kerr | abstract |
| D. Shen, *Kerr stability in external regions*, arXiv:2303.12758 | stability of Kerr in external regions $r \ge R_0$ from Cauchy data (after Caciotta–Nicolò) | abstract |
| D. Christodoulou, S. Klainerman, *The global nonlinear stability of the Minkowski space* (1993) | small data: complete, no black hole | *not re-checked* |
| D. Christodoulou, *The formation of black holes in general relativity* (EMS 2009) | trapped surfaces form from regular vacuum data (short pulse) | *not re-checked* |
| P. Hintz, A. Vasy, Acta Math. 220 (2018) | nonlinear stability of slowly rotating Kerr–de Sitter | cited by Hintz III |
| S. Klainerman, J. Szeftel, *Kerr stability for small angular momentum*, PAMQ 19 (2023); E. Giorgi, S. Klainerman, J. Szeftel, PAMQ 20 (2024) 2865–3849, arXiv:2205.14808; D. Shen (She23) | nonlinear stability of slowly rotating Kerr | KS23 (arXiv:2104.11857): Thm 1.2.1, Defs 3.4.2–3.4.5, Main Theorem v2, Rem. 3.4.3 and §3.8 ("existence of a future event horizon", horizon location) **read** in session 4; GKS citation verified; She23 not read |
| M. Dafermos, G. Holzegel, I. Rodnianski, M. Taylor, arXiv:2104.08222 (2021) | Schwarzschild stability, codimension 3 | *not re-checked* |
| G. Neugebauer, J. Hennig, arXiv:0905.4179; J. Hennig, G. Neugebauer, arXiv:1103.5248; review arXiv:1105.5830 | no stationary, axisymmetric vacuum two-black-hole configurations with regular axis (double-Kerr–NUT) | abstracts |
| F. Pretorius, PRL 95 (2005) 121101; M. Campanelli et al., PRL 96 (2006) 111101; J. Baker et al., PRL 96 (2006) 111102 | numerical binary black-hole evolutions through merger | *not re-checked* |

## Explicit binaries: obstructions and known exact classes (session 10)

Inspection levels:

- **read**: the relevant section or theorem was read;
- **abstract**: the abstract or catalogue record only;
- **cited**: the citation was located, the source was not read.

| Source | What it gives | Checked |
| --- | --- | --- |
| V. Pravda, A. Pravdová, *Boost-rotation symmetric spacetimes – review*, Czech. J. Phys. 50 (2000) 333–376, gr-qc/0003067 | axisymmetric, locally AF (a piece of $\mathscr I^+$), electrovacuum, with a second Killing field forming a 2D Lie algebra with the axial one: boost (radiative) or translation (not radiative). Sources are uniformly accelerated: "The acceleration can be caused by nodal singularities or by mutual gravitational interaction. However, in the case with no nodal singularities, negative masses occur." | read (statement and assumptions; refs [6] vacuum, [7] electrovacuum) |
| J. Bičák, B. Schmidt, *Isometries compatible with gravitational radiation*, J. Math. Phys. 25 (1984) 600–606 | the vacuum theorem behind the above | cited; *not read* |
| J. Bičák, B. Schmidt, *Asymptotically flat radiative space-times with boost-rotation symmetry: the general structure*, Phys. Rev. D 40 (1989) 1827 | construction of all AF boost-rotation symmetric spacetimes from flat-space wave solutions | abstract (OSTI 5551002) |
| J. Bičák, *Selected solutions of Einstein's field equations*, Lect. Notes Phys. 540 (2000) 1–126, gr-qc/0004016, §11 | boost-rotation symmetric solutions: nodal singularities ("strings") or self-acceleration with negative mass; global smooth sections of $\mathscr I$ except where the sources start and end | read (§11, part) |
| W. Bonnor, N. Swaminarayan, Z. Phys. 177 (1964) 240 | uniformly accelerated (positive/negative mass) particles | cited (via the reviews) |
| R. Bach, H. Weyl (1922); W. Israel, K. Khan, Nuovo Cim. 33 (1964) 331 | static collinear vacuum black holes held apart by struts | secondary sources (search results) |
| M. Gürses, F. Gürsey, *Lorentz covariant treatment of the Kerr–Schild geometry*, J. Math. Phys. 16 (1975) 2385 | Einstein equations linear (mixed form) for Kerr–Schild metrics | abstract and citing papers (e.g. Gergely, gr-qc/0203101) |
| Kerr–Schild 1965; Debney–Kerr–Schild 1969; H. Stephani et al., *Exact Solutions* (2nd ed., 2003), ch. 32 | vacuum KS: geodesic shear-free $k$, algebraically special; classification | cited; *book not re-checked* |
| superposed Kerr–Schild black holes: arXiv:2405.03289 (*Some consequences of the Kerr–Schild gauge*), gr-qc/0305071 (*Physics and initial data for multiple black hole spacetimes*), gr-qc/0205060 (*Kerr–Schild type initial data for black holes with angular momenta*) | sums of KS black holes are not vacuum; used as approximate data in numerical relativity | search summaries only; authors not re-checked |
| E. Wessels, *Exact analytic characteristic initial data for axisymmetric, non-rotating, vacuum spacetimes, with an application to the binary black hole problem*, CQG 15 (1998) 2509, gr-qc/9811074 | exact data on a null hypersurface, not a spacetime | abstract |

A broad search (extended web search, plus the reviews above) found no claimed exact, regular, AF,
dynamical two-black-hole vacuum solution. A 2021 Phys. Lett. B paper (search result, not read)
states that no exact analytic solution of the merger phase is available.

## Cauchy–Kovalevskaya with constants, constraint propagation (session 11)

| Source | What it gives | Checked |
| --- | --- | --- |
| T. Nishida, *A note on a theorem of Nirenberg*, J. Differential Geom. 12 (1977) 629–633 (cited by Baouendi–Goulaouic as "to appear") | abstract nonlinear CK theorem in a scale of Banach spaces, hypotheses (H1) Lipschitz with loss $1/(s - s')$ and (H2) $\lVert F(0)\rVert_s \le M/(s_0 - s)$ | cited; the statement is read via Baouendi–Goulaouic |
| M. S. Baouendi, C. Goulaouic, *Le théorème de Nishida pour le problème de Cauchy abstrait par une méthode de point fixe*, Journées Équations aux dérivées partielles (1977), 1–8, numdam JEDP_1977____1_0 | short proof of Nishida's theorem: weighted norm, Lemmas 1–3, inequalities (10)–(11), smallness condition (12) $a < \min(1/(16C), R/(8M))$ for $s_0 = 1$; holomorphic solutions (Remark) | **read** (PDF text; (12) garbled in extraction, re-derived in notes/15 §3) |
| L. Nirenberg, *An abstract form of the nonlinear Cauchy–Kowalewski theorem*, J. Differential Geom. 6 (1972) 561–576 | the original abstract theorem | cited |
| S. Frittelli, *Note on the propagation of the constraints in standard 3+1 general relativity*, Phys. Rev. D 55 (1997) 5992 | eqs. (12)–(13): the Hamiltonian and momentum constraints evolve by a linear homogeneous first-order system when the ADM evolution equations hold | **read** |

## Not covered

Gluing constructions of initial data (Corvino, Chruściel–Delay,
Carlotto–Schoen, Mao–Oh–Tao, Kehle–Unger), Robinson–Trautman convergence to
Schwarzschild (Chruściel), and higher-dimensional exact coalescence (e.g. on
Eguchi–Hanson base spaces) were not reviewed this session. They become
relevant if the hypotheses in notes/09 are pursued further.
