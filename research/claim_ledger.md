# Claim ledger

Every claim made in notes/07–09, with its status. Categories:

- **E** established in the literature (source and hypotheses given; "checked" says what was read);
- **R** reproduced here independently (how);
- **N** candidate new result (with what it rests on);
- **C** conjecture or numerical observation;
- **G** unresolved gap.

## Foundations and definitions (notes/07)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 1 | Smooth vacuum data satisfying the constraints have a unique MGHD | E | Choquet-Bruhat–Geroch 1969 (citation verified; statement from secondary source) |
| 2 | MGHD of data with a closed trapped surface is future null incomplete | E | Penrose 1965 (NEC, non-compact Cauchy surface; *not re-checked*) |
| 3 | Closed trapped surfaces lie in $\mathcal B$ under strong asymptotic predictability + NEC | E | Hawking–Ellis Prop. 9.2.1, Wald Prop. 12.2.2 (*not re-checked*) |
| 4 | Boundary of the trapped region is a smooth stable MOTS; barrier existence theorem (no energy condition) | E | Andersson–Metzger 2009, Thms 1.3, 3.1, 5.1 (**read**: statements and conventions $\theta_\pm = P \pm H$) |
| 5 | Outermost MOTS can jump; at zero principal eigenvalue the MOTT is tangent to the slice and lies to its future, with $D_t\vartheta = -\int\alpha W\psi$ | E | AMMS 2009 Prop. 5.1, Thm 5.2 (**read**, ar5iv) |
| 6 | EMR merger spacetimes exist (small subextremal Kerr crossing the horizon of slowly rotating Kerr(–dS)), settling to Kerr(–dS); event horizon left open | E | Hintz 2024 III, Thm 6.1 (Λ > 0), Rem. 6.2 (Λ = 0, via KS23/GKS22/She23) (**read**, §1 and §6) |
| 7 | No stationary axisymmetric vacuum two-black-hole equilibria (double-Kerr–NUT) | E | Neugebauer–Hennig 2009, Hennig–Neugebauer 2011 (abstracts) |

## Curvature code (tests/test_validation.py)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 8 | `Spacetime.riemann` agrees with an independently coded all-lower-index formula | R | exact (`cancel`), generic 3-metric with 3 free functions |
| 9 | Riemann has pair, antisymmetry and first-Bianchi symmetries | R | exact |
| 10 | Conformally flat Ricci formula (Wald App. D) holds with our signs | R | exact |
| 11 | 3+1 evolution identity with $K_{ij} = \frac1{2N}\partial_t\gamma_{ij}$ | R | exact, 4 free functions; mutation test fails as it should |
| 12 | Reissner–Nordström solves Einstein–Maxwell; Kretschmann $8(6M^2r^2 - 12MQ^2r + 7Q^4)/r^8$ | R | `simplify` |
| 13 | Kerr is Ricci flat | R | exact rational identity in $\chi = \cos\theta$ (symbolic $M, a$: slow test; $M = 1, a = 3/5$: fast) |

## Kastor–Traschen benchmark (notes/08)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 14 | Reduced field equations for arbitrary $U(\tau, \vec x)$ (three scalar expressions) | R | exact identity (test) |
| 15 | Within the ansatz, solutions are exactly $U = H\tau + V$, $\Delta V = 0$, $H^2 = \Lambda/3$; Λ = 0 ⇒ static | R | identity + elementary argument (test of each step) |
| 16 | Two-centre KT solves Einstein–Maxwell–Λ | R | via 15 (proof) and exact evaluation of the unreduced residual at three rational points (diagnostic) |
| 17 | One centre = RNdS with $Q = m$ | R | scalar identity $\lvert\nabla R\rvert^2 = (1 - m/R)^2 - H^2R^2$ + generalised Birkhoff (Birkhoff step *not re-checked*) |
| 18 | $M_{\rm ext} = 1/(4\lvert H\rvert)$ | R | discriminant (matches BHKT) |
| 19 | Homothety; only $Hm_i$, $H\lvert x_i - x_j\rvert$ physical | R | exact |
| 20 | KT slices umbilic, $K = 3H$; MOTSs are CMC surfaces $H_S = 2\lvert H\rvert$ | R | exact (probably known; no novelty claimed) |
| 21 | Equal masses $m = 0.05$, $d = 0.5$, $H = -1$: common MOTS born at $\tau_c = -0.012884047$ by a fold | C | shooting; tolerance-stable to $2\times10^{-10}$; cross-checked by polar formula ($\lvert\theta_+\rvert \le 2\times10^{-8}$) |
| 22 | Same spacetime: two disjoint MOTSs at $\tau = -0.05$; a common MOTS at $\tau = -0.009$ | N (computer-assisted) | interval arithmetic + Andersson–Metzger Thm 3.1. New only as a certificate; the fact is expected |
| 23 | Event horizons coalesce in KT for $M < M_{\rm ext}$ | E | BHKT 1994 (analytic null-geodesic bounds, equal masses); **not re-derived** |
| 24 | Extensions non-unique; $C^2$ at dS horizon, $C^1$ at inner horizons; naked singularity for $M > M_{\rm ext}$ | E | BHKT 1994; **not re-derived** |

## Hypotheses (notes/09)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 25 | Lemma A: zero-shift umbilic foliations of AF Λ = 0 spacetimes with zero momentum density have $K \equiv 0$, frozen slices | R/N | elementary proof + identity test; probably folklore |
| 26 | Same + uniqueness of normal KIDs ⇒ static | G | uniqueness not proved |
| 27 | Hyperboloidal KT ansatz admits no centres, any Λ | N (small) | exact identity + elementary argument; negative result within one ansatz |
| 28 | EH merger in Hintz's EMR spacetimes | N / G | **Proved for Λ > 0 and Λ = 0** with escape-defined $\mathcal B$ (Theorem D, claims 40, 45). **With $\mathscr I^+$ (Λ = 0): proved given our reading of Hintz 2026 at $\mathscr I^+$** (Theorem E, claims 46, 49). The session-1 inputs (A2), (A3) are no longer used |
| 29 | AMMS identity $\langle\phi_0, p\rangle = -\int N\phi_0\lvert A\rvert^2$ at a degenerate time-symmetric minimal surface | E + R | AMMS 2009; re-derived; numerically to $3\times10^{-7}$, two lapses |
| 30 | BL $d_c = 0.7661974282\,M$ (1.53239486 bare-mass units); $A_{\rm crit} = 49.101988\,M^2$ | E + R | Hui & Lin 2025 (1.5323948; 196.40795): agreement to all quoted digits |
| 31 | Individual BL horizons exist and are disjoint for every $d \in [0.76, 0.78]$ | N (computer-assisted) | interval arithmetic over the $d$-interval + Andersson–Metzger |
| 32 | Onset law (*): $t_* = (d - d_c)\langle\phi_0, \partial_dH\rangle / \int N\phi_0\lvert A\rvert^2$ | N | Lyapunov–Schmidt (sketch); verified numerically: slopes 2.385183 ($N = 1$), 10.644372 ($N = \psi^{-2}$) vs direct Richardson limits 2.385183, 10.644370; for $N = 1$ also a **rigorous** enclosure $[2.3768, 2.3936]$ (Theorem C) |
| 33 | Proposition B (apparent-horizon merger in near-critical vacuum BL developments) | N | full proof in notes/10. (F) computer-assisted (Theorem C); equivariant Lyapunov–Schmidt and fold argument written out; smooth dependence and smoothness of the expansion map standard and cited (not re-checked). Slopes certified: $N = 1$: $[2.3768, 2.3936]$; $N = \psi^{-2}$: $[10.606, 10.682]$ |
| 34 | (F): nondegenerate fold at $d_c$ | N (computer-assisted) | Theorem C, `research/calc/fold_proof.py` (session 3): Krawczyk + signs + positive Jacobi field + embeddedness; rests on Arb and on the implementation (checked 4 ways, see notes/09) |
| 35 | $\Sigma_0$ has no common MOTS anywhere for $d > d_c$ | C | shooting in a window + literature agreement |
| 36 | For $N = \psi^{-2}$ the new common horizon is marginally trapped ($\theta_- < 0$); for $N = 1$ it is not | C | sign of $p$ on $S_c$ (leading order) |
| 37 | Under predictability, the BL holes' black-hole regions are joined by $t_*$ | G | needs claim 3 for marginally trapped surfaces (closure of $\mathcal B$) and predictability of the development |
| 38 | Computer-assisted proof of (F) | N (computer-assisted) | achieved in session 3, see claim 34. History: a session-2 "verified" run was **invalid** (python-flint series cap of 10), and the session-2 diagnosis of a structural axis barrier was **wrong** (the cause was low-order steps at large step sizes) |
| 39 | $d_c = 0.766197428240 \pm 1.1\times10^{-13}\,M$ (1.532394856480 ± 2e-13 bare masses) | N (computer-assisted) | Theorem C; contains Hui & Lin's values |

## Event horizons of EMR mergers (notes/11, session 4)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 40 | Theorem D: in Hintz's Λ > 0 EMR spacetimes, with $\mathcal B$ = events that cannot signal the region beyond the cosmological horizon, for small ε: both holes' interiors lie in $\mathcal B$ with closed trapped surfaces; the horizons are localised to within any margins $c$, $\eta$; the holes lie in different components of $\mathcal B\cap\Sigma_t$ before the plunge and in the same one after it; $\mathcal B$ is path-connected | N | proof in notes/11 §5. Inputs: Hintz Thm 1.1, Def. 4.1, eqs. (1.2), (2.19), (6.2) (**read**); uniformity on the tube and over $X$ is our reading; elementary causal theory. Potentially new: Hintz lists the event horizon and the connectedness of $\mathcal B$ as open |
| 41 | Spheres $\{t_* = c, r = r_0\}$ inside a Kerr–de Sitter black hole are trapped for every spin and every Λ ≥ 0 ($\partial_rE > 0$ term by term) | R | exact identity (`tests/test_eh_merger.py`); probably known |
| 42 | Kerr–de Sitter (Carter form) solves Ric = Λg | R | exact at rational parameters (slow test); SdS symbolic (fast); Λ = 0 reduces to Kerr |
| 43 | Adapted slicing (Lemma 4): a time function equal to Fermi time near $C$ and to $t_*$ away from it exists for small ε, using a logarithmic cut-off | N (small) | elementary proof in notes/11 §4 |
| 44 | *(superseded by 45, 46 in session 5)* Theorem D₀ (Λ = 0, Hintz Rem. 6.2): merger statements 1–4 of Theorem D with $\mathscr I^+$ | G | conditional on (Ex) existence via KS23, (K-esc) far points of the initial slice escape, (K-trap) points with $r \le r_H - c_K$ do not. KS23 §3.8 proves (K-esc), (K-trap) in their own coordinates (**read**); the anchoring of KS's $r$ to Kerr's $r$ on the initial layer is *not checked* |
| 45 | Theorem D for Λ = 0: Hintz's gluing up to $t_* = T_{\rm late}$ followed by Hintz 2026's stability development. With $\mathcal B$ = events that cannot signal arbitrarily large distances, conclusions 1–5 hold as for Λ > 0 | N | notes/11 §5. Inputs: Hintz III (as in 40) and Hintz 2026 (arXiv:2606.28253) Thms 1.1, 13.1, consequence (3) (**read**); smallness of its decay constants for small data is our reading. Removes the KS anchoring problem of claim 44 |
| 46 | Theorem E (Λ = 0, textbook $\mathcal B_{\mathscr I} = M \setminus J^-(\mathscr I^+)$): both holes in $\mathcal B_{\mathscr I}$ with trapped surfaces; different components before the plunge, the same after: an event-horizon merger exactly as defined in notes/07 | N | notes/11 §6, given Theorem D's inputs and (N), which is Lemma 6 (claim 49) |
| 47 | In Kerr, $\tilde u = t - r_* + a\arctan(r/a)$ is a time function outside the horizon; outgoing principal null geodesics from far points of a slice avoid the future of the near region | R / N (small) | exact identity (`tests/test_eh_merger.py`) + elementary argument (notes/11, Lemma 5) |
| 48 | $\partial_t + \frac12\partial_r$ is uniformly timelike ($g \le -\frac14$) in Kerr for $r \ge 10M$, $|a| \le M$ | R | exact polynomial check (test) |
| 49 | Lemma 6: (N) holds. Points far out near the late initial slice reach $\mathscr I^+$, because in the 2026 paper's gauge the outgoing null component $g(\partial_0, \partial_0)$ is $O(r^{-1-\delta})$ | N (reading) | Kerr part checked exactly ($-a^2\sin^2\theta/4r^2$). The boost (Lemma 3.19) and perturbation (Def. 5.3: $h_{00} \in \pi^C h$) parts are **read**. The support of the Nash–Moser correction in $t_* \ge 1$ (Def. 12.1, Thm 12.3) was confirmed in session 6. Uniformity up to the corner of $\mathscr I^+$ is our reading; an independent expert check is advised |
| 50 | Corollary D1: on the slow scale the event-horizon merger time converges to the time $C$ crosses $r_H$ ($\epsilon \to 0$) | N (small) | notes/11, from Theorem D(4) |
| 62 | Lemma 7: (N) holds at all late times. Every point of the late development with $r \ge R_N$ reaches $\mathscr I^+$ | N (reading) | escape curve of Lemma 6 started at any late time. For $t_* \ge 1$: $h$ from Hintz 2026 (13.2), with $h_{00} = O(\rho_I^{1+2\gamma^\Upsilon})$ (proof of its Lemma 11.5), $\chi_K = 0$ near $\mathscr I^+$ (Def. 3.17), boost term (3.44). For $t_* < 1$: Lemma 6, which uses the support of the Nash–Moser correction (**read**). Frame reading confirmed against (3.29)–(3.36) and Def. 5.2 in the session-9 review, which also corrected the case split. Outside expert check still advised |
| 63 | Corollary E1: for Λ = 0, $\mathcal B_{\mathscr I} = \mathcal B$. So localisation and connectedness (Theorem D, 3 and 5) hold for the textbook black-hole region | N | Lemma 1(iv) + Lemma 7 + Theorem D (notes/11 §6); mass-mismatch identity checked in `tests/test_eh_merger.py` |

## Gap 3: no common apparent horizon at t = 0 (sessions 6–7; notes/12)

| # | Claim | Cat. | Basis |
| --- | --- | --- | --- |
| 51 | Lemma G1 (BL, d = 0.78): spheres of radius $s \le 0.155$ about a puncture have $H < 0$, so closed minimal surfaces enclosing that puncture avoid them | N (computer-assisted) | `research/calc/bl_global.py`, interval boxes including $s = 0$; maximum principle |
| 52 | Lemma G2 (BL, d = 0.78): spheres $\lvert x\rvert = R \ge 1.05$ have $H > 0$, so every closed minimal surface lies in $B(0, 1.05)$ | N (computer-assisted) | interval boxes for $R \le 2$, analytic bound beyond |
| 53 | A common minimal surface, if one exists, has a G-invariant outermost representative: an embedded sphere whose profile runs from $(0, z_0)$ to $z = 0$, meeting the plane once at a right angle, with $z_0 \in [0.545, 1.05]$ | N (small) | outermost MOTS (Andersson–Metzger) is unique, hence symmetric; stable outermost components are spheres; G1, G2 |
| 54 | Shooting landscape at d = 0.78: the curves dive into a throat for $z_0 < z_{\rm top} = 0.5896816$, the top of the individual horizon. Just above it they hug the individual horizon, pass a neck of radius $\approx 4(z_0 - z_{\rm top})$ at $z \approx 0.2126$ (its bottom) and flare out of $B(0, 1.05)$. For $z_0 \in [0.605, 0.87]$ they reach the plane with $F \ge 0.017 > 0$; beyond that they leave the ball | C | floating point (`bl_global.landscape`, rtol $10^{-13}$) |
| 55 | Obstruction: necks of every size near $z_{\rm top}$. In $\sigma$-time ($ds = \rho\,d\sigma$) the axis is a line of hyperbolic fixed points; a proof needs an isolating-block passage argument there. The stability Jacobi field changes sign only after the neck, and the Penrose bound $A \le 16\pi$ is not violated before it | G → resolved | analysis in idea_log; desingularised system checked in `tests/test_bl_global.py`. *Resolved in session 7* by claims 57–58 (thin necks are unstable), not by isolating blocks |
| 56 | Lohner's method in Arb (`research/calc/lohner.py`): parallelepiped enclosures with QR frames, dual-number Jacobians, Lagrange remainders over a Picard box | R (method) | tests against rotation, a hyperbolic saddle and a Riccati equation (exact solutions); standard method (Lohner 1987), our implementation |
| 57 | Thin-neck lemma (BL, d = 0.78): for a waist of radius $\lambda \le 0.00488$ at height $z_w \in [0.1797, 0.2422]$, the $\lambda$-variation Jacobi field changes sign within $\lvert\sigma\rvert \le 1.5$ on both sides, so the surface is unstable | N (computer-assisted) | `research/calc/gap3_neck.py`: 640 pieces, margin $\ge 0.527$; catenoid limit and Jacobi equation by SymPy; fails as it must at $\sigma = 1.0$; a second integrator (`gap3_neck_lohner.py`) agrees at test points; stability of the outermost surface (claim 4) |
| 58 | Approach lemma: a curve with $u < 0$, $\rho \le \rho_e$ either closes on the axis or reaches a waist with $\rho_w \le \lvert v_e\rvert + K\rho_e^2/2$, $\lvert z_w - z_e\rvert \le (2\rho_e/e + K\rho_e^2/2)/\sqrt{1 - K\rho_e - K^2\rho_e^2/4}$ | N (small) | elementary ODE estimate (notes/12 §3); inequalities in `tests/test_gap3.py`; $K = 16.53$ by interval arithmetic |
| 59 | **Theorem F**: equal-mass time-symmetric BL data at $d = 0.78$ ($M = 1$) contain no closed embedded minimal surface enclosing both punctures; the apparent horizon at $t = 0$ is the two individual horizons | N (computer-assisted) | claims 51–53, 57, 58 and the validated sweep: 2,452 $z_0$-intervals covering $[0.543, 1.051]$, each ending in E1/E2/E3/N (`research/calc/gap3_sweep.py`, log in `research/calc/data/`); 2.8 h on 4 cores. Potentially new: Dain–Weinstein–Yamada 2011 prove non-existence only for sufficiently large separation |
| 60 | Falsification of the sweep: at $d = 0.76 < d_c$ the $2^{-14}$ intervals containing the two real common minimal surfaces are not excluded, while an interval between them is | R (test) | `tests/test_gap3.py` |
| 61 | Floating point: the thin-neck sign pattern persists for necks up to radius 0.04 (with $\sigma = \pm2$), so medium necks are unstable too | C | point evaluations of the certified integrator; region certificate not achieved (parameter wrapping) |

## Gap 3 over an interval; the head-on example (sessions 9–10; notes/12 §7, 13)

| # | Claim | Status | Evidence |
| --- | --- | --- | --- |
| 64 | Small-neck lemma over $d \in [0.765625, 0.78125]$: $\lambda \le 3\cdot2^{-8}$, $\zeta_w \in [-120, -68]\cdot2^{-9}$, both directions | N (computer-assisted) | `gap3_neck.certify_region_rel`, nine runs, 2,664 pieces, margin $\ge 0.0029$; record `research/calc/data/gap3_neck_interval.txt` (run once, replay command there) |
| 65 | Local fold lemma: $F > 0$ on $z_0 \in [0.62891, 0.63696]$ for every $d \in (d_c, 0.767578]$ | N (computer-assisted) | Theorem C + `gap3_fold.py`: part 1 $F_{zz} \ge 125.1$, $F_d \ge 0.591$ (33 pieces); part 2 $F_d \ge 0.177$ (198 pieces); records in `research/calc/data/gap3_fold_*.txt` (run once); thin-run cross-check against Theorem C (`tests/test_gap3_fold.py`, slow) |
| 66 | **Theorem F′**: for every $d \in (d_c, 0.767578125]$ no closed embedded minimal surface encloses both punctures | N (computer-assisted) | claims 64, 65, interval barriers (`tests/test_bl_global.py`), band certificate and sweep (`gap3_interval.py`, verify of the stored log); see notes/12 §7 for the replay status |
| 67 | **Theorem H**: near-critical BL head-on development: two individual horizons and no common one at $t = 0$; a common MOTS from $t_* = c(d - d_c) + O((d - d_c)^2)$, $c \in [2.3768, 2.3936]$ | N (computer-assisted + analytic), $\delta_0$ non-explicit | claims 33, 38, 66; notes/13 |

## Toward explicit binaries (session 10; notes/14, research/route_map.md)

| # | Claim | Status | Evidence |
| --- | --- | --- | --- |
| 68 | Proposition P: the BL development is Petrov type I on an open dense subset of the initial slice for every $d > 0$ (type D on the axis) | N (small), analytically proved | closed-form eigenvalue differences on $z = 0$; exact rational check at $d = 0.78$; `research/calc/petrov_bl.py`, `tests/test_petrov_bl.py` |
| 69 | P ⇒ no algebraically special ansatz (single KS, KSK, RT, Kundt, twisting II/D) reproduces the BL development on an open set meeting the type-I region | N (small) | P + Goldberg–Sachs |
| 70 | Axisymmetric locally AF radiative vacuum with a second Killing field is boost-rotation symmetric; its sources need nodal singularities or negative masses | E | Pravda–Pravdová, gr-qc/0003067 (section read); Bičák gr-qc/0004016 §11 (read); Bičák–Schmidt 1984 (*not read*) |
| 71 | Proposition S: the geodesic-slicing development of BL data is a convergent series in $t$ near each point, with explicit recursive coefficients; it solves the full vacuum equations | R (standard theorem applied) | Cauchy–Kovalevskaya + constraint propagation; code validated against exact Schwarzschild free fall and the Hamiltonian constraint (`tests/test_bl_series.py`) |
| 72 | Near the critical surface ($d = d_c$, $(\rho, z) = (0.25, 0)$) the series' convergence radius is $T \approx 2.7$ | C (numerical) | ratio estimates through $t^{14}$, `research/calc/data/bl_series_dc.json` |
| 73 | For conformally flat time-symmetric vacuum data, $B_{ij} = -t\,C_{ij}(\gamma_0) + O(t^3) = O(t^3)$; for BL the $t^3$ coefficient is nonzero | R (identity) + C (numerical) | `bl_series.magnetic_weyl`; Schwarzschild control $B \equiv 0$ |
| 74 | Route J: for $d - d_c$ small, the apparent-horizon merger region near $S_c$ lies inside the convergence domain of S | N (small), qualitative | analytic dependence on $d$ + claims 67, 71; explicit bounds missing |
| 75 | An explicit (category A) dynamical binary metric | **not found** | routes A–J screened, `research/route_map.md` |
