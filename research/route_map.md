# Route map: explicit dynamical binary-black-hole spacetimes

*Created session 10 (2026-10-31). Target: an explicit dynamical binary-black-hole metric of 4D,
asymptotically flat (AF) vacuum GR, $R_{\mu\nu} = 0$, $\Lambda = 0$, $G = c = 1$. First restriction:
nonspinning, axisymmetric, head-on. Result categories (A–E) and levels (1–5) as defined in
[`notes/14-explicit-binaries.md`](../notes/14-explicit-binaries.md) §1.*

**Status legend.**

- **rejected**: a precise obstruction applies, either from the literature or proved here.
- **tested**: a decisive calculation was run.
- **promising**: an executable plan with an identified next calculation.
- **untested**: only screened.

No route has produced a category-A metric. **No new exact solution was found.** *Session 12:* route E now gives an
**evaluable** category-B metric with certified errors (`research/calc/bl_metric.py`, notes/16). Three further
closed-form classes are excluded for the BL head-on spacetime: algebraically special (B), conformally flat
geodesic slices, and silent flow (notes/16 §5).

## Summary

| Route | Target category | Status | Decisive result so far |
| --- | --- | --- | --- |
| A. Axisymmetric reduction with radiation | A | screened, not integrable | 1 Killing field: 2+1 gravity + a wave map, no known integrable structure (§A) |
| B. Kerr–Schild-type decompositions | A | **rejected** (single KS, any algebraically special form) | BL head-on spacetime is Petrov type I on an open dense set (proved, exact) |
| C. Interaction corrections to two-hole seeds | A/D | untested beyond literature | superposed KS data are not vacuum (literature); first correction is an unrestricted PDE |
| D. Matched asymptotic zones | D | untested | finite-order matching is not exact; the transition-region residual is $O(\epsilon^n)$, not zero |
| E. Time series from exact data | **B** (local) | **tested, works** | BL development in geodesic gauge = convergent series (Cauchy–Kovalevskaya), coefficients computed to $t^{14}$, validated; certified convergence time $T_{\min} = 5.556\times10^{-4}$ near $S_c$ (sessions 11–12); **evaluator** with certified enclosures of $g$, $g^{-1}$, $\partial g$, $\partial^2g$ (session 12) |
| F. Characteristic / tetrad formulations | B/C | screened | exact null data exist (Wessels 1998), but they are data, not a spacetime |
| G. Solution-generating methods (Ernst, inverse scattering, Bäcklund) | A | **rejected** for regular AF radiating binaries (literature) | need 2 commuting Killing fields; AF radiative ⇒ boost-rotation symmetric ⇒ struts or negative mass |
| H. Resummation of expansions | B | untested | no convergence theorem for PN/PM; formal |
| I. Numerical/symbolic discovery | E | **tested** (session 12) | exact universal $t^4$ coefficient; conformally flat slices and silent flow excluded (certified); Padé exploratory (notes/16 §5) |
| J. Near-critical BL as the explicit example | B + C | **done at an explicit separation** (session 16) | At $d = 0.7662$: no common MOTS at $t = 0$ (F′); none in the shell $\Omega'$ for $t \le 5\times10^{-6}$ (Theorem N); one at $t = 6.5\times10^{-6}$ (Theorem M), and at $t = 1.2\times10^{-5}$ for every $d \in [d_c^+, 0.7662]$ (Theorem M′). All inside the certified series metric. Barriers and foliations replace the spectral gap (notes/20) |

## A · Symmetry reduction with retained radiation

- **Question.** Does the axisymmetric, twist-free vacuum system (the head-on case) admit a closed
  functional ansatz with genuine radiation?
- **Representation.** Write the metric as
  $g = e^{2\Phi}\,d\varphi^2 + e^{-2\Phi} h_{ab}\,dx^a dx^b$, with $\partial_\varphi$
  hypersurface orthogonal and $h$ a Lorentzian 3-metric on the orbit space. The vacuum equations
  become $\Box_h \Phi = 0$ and $R_{ab}(h) = 2\,\partial_a\Phi\,\partial_b\Phi$ (Geroch
  reduction). Axis regularity requires $e^{2\Phi} \sim \rho^2$ with no conical defect. 2+1 gravity has no
  local degrees of freedom, so the single radiative degree of freedom is $\Phi$ (the "+"
  polarization).
- **Known.** Integrable structure (Ernst, Belinski–Zakharov) appears only with a *second*
  commuting Killing field. Route G shows that this forces boost-rotation symmetry.
- **Smallest decisive calculation.** A symmetry (Lie point) analysis of the reduced system, to find
  any further reduction compatible with two separated sources and a regular axis. Not run. The
  plausible outcome is "only the known reductions", which would be a restricted no-go, not a metric.
- **Success / falsification.** Success means a closed family with two sources and a regular axis.
  Falsified if every symmetry reduction adds a Killing field (then route G applies).
- **Cost.** Moderate (Lie symmetry software; the PDE system is small).

## B · Kerr–Schild and related decompositions

- **Question.** Can $g = \eta + H\,k\otimes k$, or another algebraically special form, describe a
  head-on binary?
- **Known (literature).**
  - Vacuum KS metrics have a geodesic, shear-free $k$ and are algebraically special
    (Goldberg–Sachs). The vacuum KS solutions are classified (Kerr–Schild 1965; Debney–Kerr–Schild
    1969; Stephani et al. 2003, ch. 32). *Not re-checked in the book; consistent statements in the
    sources listed in `research/literature.md`.*
  - Superposed KS data are not vacuum (e.g. arXiv:2405.03289).
- **Decisive calculation (done, `research/calc/petrov_bl.py`).** The vacuum development of
  time-symmetric BL data is of Petrov type I on an open dense subset of the initial slice, for
  every separation $d > 0$, hence on a spacetime neighbourhood of it.
  - On the plane $z = 0$, the eigenvalue differences of $E_{ij} = R_{ij}(\gamma)$ are given in
    closed form: $\lambda_\varphi - \lambda_z = 3d^2/(2r^4(2r+1)) \ne 0$.
  - At $d = 0.78$ the discriminant is an exactly positive rational.
- **Consequence.** No algebraically special ansatz reproduces this spacetime on any open set
  meeting the type-I region: single-KS, Kerr–Schild–Kundt, Robinson–Trautman, Kundt, twisting
  type II/D.
- **Scope.** This is about the BL development. It does not rule out other head-on binaries in
  general, though generic configurations are expected to be type I.
- **Status.** Rejected for single-vector KS. **Untested:** two-vector KS
  $\eta + H_1 k_1k_1 + H_2 k_2k_2$, which need not be algebraically special. The decisive test is
  whether the $O(m_1m_2)$ cross term of $\mathrm{Ric}$ lies in the image of the linearized operator
  restricted to $\{q_1k_1k_1 + q_2k_2k_2\}$: 10 equations for 2 functions, expected inconsistent.

## C · Nonlinear interaction corrections to two-hole seeds

- **Question.** Does a finite family of correction functions close under the full equations?
- **Known.** Superposed KS / BL data are used as approximate data in numerical relativity, and the
  constraint corrections are generic elliptic solutions (Matzner–Huq–Shoemaker; Moreno–Núñez–Sarbach,
  gr-qc/0205060).
- **Assessment.** Without a structural reason, the first correction is an unrestricted solution of
  the linearized equations with a quadratic source. That is the original problem in perturbative
  form, which is route D or H, not an explicit metric. **Untested**, deprioritised.

## D · Matched inner / near / radiation zones

- **Assessment.** Matched asymptotics produce controlled approximations (category D) at finite
  order. The residual of an interpolated metric in the transition zone is $O(\epsilon^n)$, not zero.
  Hintz's gluing (note 11) is the rigorous version for the extreme mass ratio, and it is an
  existence theorem (C), not a formula. **Untested** here.

## E · Time series from exact initial data (executed)

- **Question.** Is the vacuum development of BL data representable by an explicit convergent
  series, and what does that representation cover?
- **Representation.** Geodesic slicing ($N = 1$, $\beta = 0$), comoving coordinates,
  $\gamma = \sum_k g_k t^{2k}$, $K = t\sum_k \kappa_k t^{2k}$, with an explicit recursion
  (`research/calc/bl_series.py`): $g_1 = -\mathrm{Ric}(\gamma_0)$, and so on.
- **Convergence (theorem, standard).** In second-order form the system is of Kovalevskaya type.
  The data are analytic away from the punctures, so the series converges for $|t| < T(x)$, and the
  constraints propagate. This makes it **category B** (exact, convergent) at **level 2** (a local
  dynamical spacetime). *Session 11:* $T \ge T_{\min} = 2.988\times10^{-4}$ on the tube of radius
  $0.02$ about $S_c$, certified (notes/15, Proposition T). Elsewhere $T(x)$ is not computed.
- **Validation.**
  - Single centre: matches the exact Schwarzschild free fall from rest, $(r/R_0)^2$, through
    $t^{10}$ to $10^{-12}$.
  - Hamiltonian constraint satisfied order by order to rounding.
  - $g_1$ agrees with the independent Cartesian formula.
- **What it is not.** It is not a closed form, and it does not cover a substantial pre-merger
  interval. Geodesic slicing reaches the crunch at finite $t$: for Schwarzschild,
  $T = \pi\sqrt{R_0^3/8M}$, which is $\pi M$ at the throat.
  - *Corrected in session 13 (note 17 §7):* at that time the comoving observer reaches the $r = 0$
    curvature singularity, with no caustic before it. The *event* is physical.
  - What is gauge-dependent is that this slicing reaches it, and when.
  - For BL the meaning of the numerical $T \approx 2.1$–$2.8$ is undetermined: it could be a
    curvature singularity, a caustic, or a complex singularity.
- **Next.** *Done near $S_c$ (session 11):* an abstract Cauchy–Kovalevskaya bound with constants
  certified in Arb. A sharper scale would be needed to reach the numerical $T \approx 2.7$.

## F · Characteristic, double-null, tetrad formulations

- **Known.** Exact analytic *characteristic data* for axisymmetric non-rotating vacuum spacetimes
  exist (Wessels 1998, CQG 15, 2509, gr-qc/9811074). They are data on a null hypersurface, not a
  spacetime. An exact far-zone expansion alone is not an exact two-black-hole spacetime.
- **Possible use.** Characteristic data posed on the outgoing null cone of the near-critical BL
  region would let the radiation content be specified, an advantage over BL, whose radiation content
  is fixed implicitly. **Untested.**

## G · Solution-generating methods (literature obstruction)

- **Hypotheses of the methods.** Ernst / Belinski–Zakharov inverse scattering / Bäcklund
  transformations need two commuting Killing fields (stationary axisymmetric, or two spacelike).
- **Obstruction.**
  - An axisymmetric, locally AF (a piece of $\mathscr I^+$), radiative vacuum or electrovacuum
    spacetime with a second Killing field forming a 2D Lie algebra with the axial one is
    **boost-rotation symmetric**. This is stated in Pravda–Pravdová, *Boost-rotation symmetric
    spacetimes – review*, Czech. J. Phys. 50 (2000) 333, gr-qc/0003067, relevant section read,
    citing Bičák–Schmidt, J. Math. Phys. 25 (1984) 600 (vacuum) and Bičák–Pravdová 1998
    (electrovacuum). *Bičák–Schmidt 1984 itself not read.*
  - In these spacetimes the objects are uniformly accelerated. "The acceleration can be caused by
    nodal singularities or by mutual gravitational interaction. However, in the case with no nodal
    singularities, negative masses occur" (same review). The same is stated in Bičák,
    gr-qc/0004016 §11, read.
  - Stationary or static two-black-hole equilibria are excluded (Neugebauer–Hennig 2009,
    Hennig–Neugebauer 2011; Bach–Weyl and Israel–Khan need a strut).
- **Consequence.** Any exact binary from these methods that is AF and radiative has struts or
  strings, or negative masses, or uniformly accelerated (not mutually infalling) components. A
  regular head-on collision of two positive masses has a single (axial) Killing field, where these
  methods do not apply. **Rejected** for the stated target. Not a no-go for other methods.

## H · Resummation

PN/PM expansions of the two-body metric are asymptotic, with no convergence theorem in the
strong-field region. An effective-one-body metric is an effective Hamiltonian's metric, not the
physical spacetime. A Padé or Borel prescription is not a proof. **Untested;** it needs a concrete
summability hypothesis first.

## I · Numerical/symbolic discovery

Route E's coefficients serve as discovery data. Ratio estimates $|g_k/g_{k+1}|$ estimate $T^2$
(numerical evidence only; for Schwarzschild they converge slowly to $T^2$ from above). No pattern
promoted to a claim.

## J · Near-critical BL: the explicit example (partly done, session 11)

- **Combination.**
  - Theorem H (note 13; Proposition B + Theorem F′): an apparent-horizon merger at
    $t_*(d) = c(d - d_c) + O((d - d_c)^2)$, $c \in [2.3768, 2.3936]$.
  - Route E: near the critical surface $S_c$ the spacetime *is* the convergent series for
    $|t| < T$.
- **Done (notes/15).**
  - *Proposition T* (computer-assisted). For every $d \in [d_c - 10^{-13}, 0.767578125]$, the
    development equals its convergent series on $U \times (-T_{\min}, T_{\min})$. Here $U$ is the
    closed tube of radius $0.02$ about $S_c$ and $T_{\min} = 2.988\times10^{-4}$. The proof is the
    abstract CK theorem with constants re-derived from Baouendi–Goulaouic; $C$ and $K_0$ are
    certified on 2,815 complex cells. The note-14 argument "$T$ is analytic in $d$, hence bounded
    below" was not a proof and is replaced by this.
  - *Corollary J*. For $d \in (d_c, d_c + \delta_1)$ the local common-MOTS formation lies inside
    that region: category B for the metric, C for the merger. $\delta_1$ is not explicit;
    $\approx 1.25\times10^{-4}$ by the linear law (heuristic).
- **Missing.**
  - An explicit $\delta$ (blocking estimates, notes/15 §7.3):
    - a certified spectral gap $\lambda_2(L_{S_c})$ over all Fourier modes;
    - certified second derivatives of the true shooting function on $\Sigma_t$, from a validated
      MOTS solver for $(\gamma(t), K(t))$ fed by the Proposition T enclosures.
  - Coverage of the whole Theorem F′ interval. It needs $T \gtrsim 0.0033$, eleven times
    $T_{\min}$.
  - Control of common MOTSs outside $\mathcal U$ for $0 < t < t_*$. Without it, the local branch is
    not shown to be the first common apparent horizon.
- **Next decisive calculation.** The $m \ge 0$ spectral gap of $L_{S_c}$ (validated ODE
  eigenvalue enclosures on the profile), then the validated shooting solver on $\Sigma_t$.

## E → A · Using the evaluable metric to search for a closed form (session 12)

The B-class representation is now evaluable, and its coefficients are computed exactly (rational arithmetic
at rational-distance points) or in balls. Four experiments were run (notes/16 §5):

| Experiment | Candidate class | Result |
| --- | --- | --- |
| A | the $t^4$ coefficient is a universal covariant expression in $\mathrm{Ric}$, $\Delta\mathrm{Ric}$ | holds (proved): $\tfrac56\mathrm{Ric}\cdot\mathrm{Ric} - \tfrac16\lvert\mathrm{Ric}\rvert^2\gamma_0 - \tfrac1{12}\Delta\mathrm{Ric}$; structural, not a closed form |
| C1 | geodesic slices stay conformally flat ($-dt^2 + \Phi^4\cdot$ flat) | **excluded** at order $t^2$ (Cotton tensor $\ne 0$, certified) |
| C2 | silent flow ($B_{ij} = 0$; ODEs along worldlines; contains every local-free-fall ansatz) | **excluded** at order $t^3$ (certified) |
| D | Padé in $t^2$ | exploratory; see notes/16 §5 |

**Next for A.** Untested classes, each with a bounded test:

- conformally flat or silent *maximal* slicings (a different foliation);
- two-centre Kerr–Schild with two null vectors (route B);
- a variable change built from the Weyl invariants (Experiment B).

**Next for B.** A validated energy-estimate continuation step: the abstract-CK route is capped
near $T \approx 10^{-3}$ (notes/16 §2).

## Session 13 · Correction, extension, geometric variables

- **E (the representation).**
  - The certified time of sessions 11–12 used a misread abstract-CK constant. It is corrected, and
    extended with the proved Theorem CK′.
  - Now $T = 8.683\times10^{-4}$ on the whole tube and $1.545\times10^{-3}$ on the waist region (notes/17 §3).
  - The CK route's apparent "ceiling" of session 12 was conditional, not a theorem.
- **Toward an energy method** (the next route for B): Einstein–Bianchi in Gauss gauge, which is
  symmetric hyperbolic in the present chart, with an a posteriori $H^s$ estimate. Inspected and not
  executed; the requirements are in notes/17 §3.
- **A (closed form).** Two more bounded classes are excluded for BL:
  - pointwise algebraic dependence on the tidal tensor (at $t^6$);
  - closed isotropic tidal laws $E' = f(E)$ (at $t^4$).

  The first orders are algebraic and universal for conformally flat time-symmetric data (notes/17
  §§5–6). Still untested: two-vector Kerr–Schild; maximal-slicing analogues of C1, C2; a recurrence
  in the order (needs more exact orders).

## Session 14 · Bottlenecks, local/nonlocal split, foliations

- **Map.** notes/18 §3 gives six bottlenecks with evidence, and §4 audits the scope of every
  earlier exclusion.
- **E → A (H1).** The geodesic development is exactly silent through $t^4$. The first nonlocal
  source $(\mathrm{curl}\,H)_3$ has a covariant law. Near the holes the local (silent) collapse
  dominates; elsewhere the nonlocal part does (exploratory).
  - *Status:* tested. The silent part is kept as a local model; a global silent factorisation is
    rejected (exploratory).
- **Gauge (H2).** Conformally flat leaves are excluded for every bounded maximal-type lapse at three
  separations (certified, Proposition L). Pointwise, a conformally flat lapse exists, because there
  are only two conditions per point.
  - *Status:* tested, rejected for that family.
  - *Open:* arbitrary lapses, an overdetermined 2×1 fifth-order system.
- **Not executed, ranked** (notes/18 §5):
  - H8, the $m_+m_-$ interaction sector;
  - H6, the Geroch-reduced 2+1 test;
  - H4, two-vector Kerr–Schild;
  - H7, Misner data;
  - H5, the energy step.
- **Next decisive experiment.** The singular set $t = \tau(x)$ of the series, from the coefficients
  to $t^{24}$ along a ray. A real, smoothly varying $\tau$ favours a resummation in $s/\tau(x)^2$;
  a complex one needs nonlocal continuation.

## Session 15 · Spatial coupling in Einstein–Bianchi form

- **Formulation (route E, toward energy methods).**
  - The Gauss-gauge E–H system is derived and verified exactly against the series.
  - Only the frame form is symmetric hyperbolic.
  - The $L^2$ curvature energy identity is proved.
  - *Open:* commuted energies, the domain of dependence of the tube, restart errors.
- **E → A / B (coefficient hierarchy).**
  - The $t^6$ law is a theorem for all conformally flat vacuum data.
  - The curl-degree grading, with its proved lag lemma, organises all orders into linear transport
    problems along the silent flow (main route M1).
  - The first experiment did not falsify M1.
  - *Next decisive experiment:* the $j = 1$ sector as an ODE in $t$ up to the silent collapse.

## Session 16 · Route J completed at an explicit separation (notes/20)

- **Method.** Existence of a MOTS by Andersson–Metzger barriers; exclusion in a shell by a foliation
  with $\theta_+ > 0$ (convex ridges allowed). Both are built from level sets of the first-order
  expansion, on which the true expansion is $h + R$, with $R$ a certified mean-value remainder.
  The quantitative Lyapunov–Schmidt programme (spectral gap, $C^{2,\alpha}$ bounds, MOTS solver on
  $\Sigma_t$) is no longer needed.
- **Result.** At $d = 0.7662$:
  - no common MOTS at $t = 0$ (F′);
  - one at $t = 6.5\times10^{-6}$ (Theorem M), $1.06$ times the linear-law onset.

  - none in the shell $\Omega'$ for $t \le 5\times10^{-6}$ (Theorem N).

  The local birth time is bracketed in $[0.81, 1.06]\,t_*^{\rm lin}$.
- **Next decisive calculation.** Theorem F′ for $0 < t \le t_0$: a CK certificate on
  $B(0, 1.05)$ minus the throat balls, needed only for $t \lesssim 10^{-5}$, plus the F′ sweep for the
  perturbed (level-$h$, $\lvert h\rvert \le \sup\lvert R\rvert$) shooting curves. That makes the birth time of
  the first common apparent horizon of the whole slice explicit.
- **Then.** A $d$-interval version (Lohner brackets in $(z_0, d)$), and the time range of the CK
  certificate (larger $d - d_c$).
- **Unchanged.**
  - Category A: no candidate.
  - Event horizons: need the global development.
