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
| 49 | Lemma 6: (N) holds. Points far out near the late initial slice reach $\mathscr I^+$, because in the 2026 paper's gauge the outgoing null component $g(\partial_0, \partial_0)$ is $O(r^{-1-\delta})$ | N (reading) | Kerr part checked exactly ($-a^2\sin^2\theta/4r^2$). The boost (Lemma 3.19) and perturbation (Def. 5.3: $h_{00} \in \pi^C h$) parts, and the support properties in the proof of Thm 13.1, are **read** but rest on our reading; independent expert check advised |
