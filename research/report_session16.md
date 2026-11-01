# Report: an explicit apparent-horizon merger, session 16

*2026-11-01. Branch `claude/hopeful-euler-w9ghbg`, based on `main` at `90f023b` (PRs #1–#6).
Mathematics: [`notes/20-explicit-merger.md`](../notes/20-explicit-merger.md).*

**The session closes the gap "no certified merger at any explicit $d$".** For equal-mass
Brill–Lindquist data at $d = 0.7662$:

- the time-symmetric slice has no common apparent horizon (Theorem F′, inherited);
- the geodesic-slicing slice at $t = 6.5\times10^{-6}M$ (and at $1.2\times10^{-5}M$) contains one
  (**Theorem M**);
- the whole event takes place where the spacetime metric is the certified convergent series of route E.

The complementary exclusion (no MOTS in the shell before $t_0 = 5\times10^{-6}$, Theorem N) is **not
completed**. The foliation certificate works on test pieces except for the ridge condition over
$z_0$-pieces near the fold, and the fix is identified (notes/20 §4). So the birth time has an
explicit upper bound, $1.06$ times the linear law, but no certified lower bound.

No closed-form (category-A) metric was found, and nothing is claimed about event horizons.

## Obstacle investigated

Theorem H (notes/13) proved the apparent-horizon merger only for $d \in (d_c, d_c + \delta_0)$ with
$\delta_0$ not computed. The route planned in notes/15 §7.3 made $\delta_0$ explicit through a
quantitative Lyapunov–Schmidt reduction. It needed:

- a certified spectral gap of the Jacobi operator of $S_c$ over all Fourier modes;
- $C^{2,\alpha}$ bounds on the nonlinear expansion map;
- a validated MOTS solver on the evolved slices.

None of these existed.

## The idea

**Existence and exclusion of MOTSs need only $C^0$ information**:

- *existence*: the true $\theta_+$ on two surfaces (the Andersson–Metzger barrier theorem);
- *exclusion*: the true $\theta_+$ on a foliation (the comparison principle).

Neither needs the spectrum.

The difficulty is that $\theta_+ \approx 10^{-6}$ arises from cancellation between terms of size $1$.
It is removed by choosing the surfaces as **level sets of the first-order expansion**,
$\theta_+^{(1)} = H[\gamma_0] + t\,\mathrm{tr}_S\kappa_0 = h$. They are profiles of an ODE, so validated
shooting produces them exactly. The true expansion on them is then $h + R$, where $R$ is a
mean-value remainder with crude enclosures, bounded with the CK tails of notes/16.

For the foliation, leaves whose shooting curve meets the plane at a non-right angle are allowed. With
$F > 0$ the equatorial ridge is convex, and a smooth MOTS cannot touch it. So no closing condition is
needed.

## What was obtained

| Item | Status | Evidence |
| --- | --- | --- |
| Lemma 1: $\theta_+$ of a surface of revolution in any axisymmetric metric | derived | SymPy level-set check, non-conformally-flat metric with $\gamma_{xz} \ne 0$, $K \ne 0$ |
| Lemma 2: $\theta_+ = h + R$ on level surfaces, $R$ by the mean-value form | proved | conventions checked: $\Theta(p_0) - h \approx 10^{-17}$ on a floating-point level surface; true $\theta_+ - h \approx 10^{-11}$ from the certified evaluator |
| Lemma 3: foliation exclusion with convex ridges | proved (elementary) | notes/20 §2.3 |
| barrier theorem (AM 2009, Thms 3.1 and 4.1) | literature, **read in the PDF** | hypotheses recorded in notes/20 §2.3 and `research/literature.md` |
| **Theorem M** at $t_1 = 1.2\times10^{-5}$ and $t_1' = 6.5\times10^{-6}$ | computer-assisted | `bl_barrier.json`, `bl_barrier_near.json` |
| Theorem N, $t \in [0, 5\times10^{-6}]$ | **not completed**: Lemma 3's hypotheses certified on test pieces, except $F > 0$ over $z_0$-pieces near the fold | `bl_foliate.run`; diagnosis in notes/20 §4 |
| near $S_c$, small $t$: outer trapped but not trapped | numerical observation | $p = \mathrm{Ric}(\nu,\nu) \in [-0.58, 0.11]$ on $S_c$ |

## Tools

- `research/calc/bl_barrier.py`:
  - the level-$h$ shooting (fold_proof's integrator, copied with the source as a parameter;
    bit-identical regression);
  - the general expansion and the remainder;
  - validated secant refinement;
  - barrier certification;
  - polar nesting and the domain check.
- `research/calc/bl_foliate.py`:
  - gap3_fold's Lohner machinery with the level-$h$ source, $t$ as a set parameter, and
    **sub-step enclosures** (centre polynomial plus the Picard flow-Jacobian applied to the
    spread);
  - the foliation certificate and the extreme-leaf comparison (not completed: notes/20 §4).
- `tests/test_bl_barrier.py`: 8 fast tests and 1 slow one. The fast ones cover the formulas, the
  conventions, the Lohner sub-steps on an exact flow, and the stored records. The slow one checks
  that the integrator copy is bit-identical.

## Lessons

- **Moore brackets wrap by $\sim 3\times10^7$ along these curves.** $\delta = 10^{-8}$ blows up in the
  last segment. $\delta = 10^{-10}$ survives, but reaches the plane with $\rho$-widths of $6\times10^{-3}$:
  fine for crude remainder bounds, useless for comparing surfaces $3\times10^{-4}$ apart. Lohner sets
  over the same bracket have widths $\sim 10^{-10}$.
- **The floating-point roots are good only to $3\times10^{-10}$**, even at `rtol = 1e-13`. Brackets
  of $10^{-10}$ need validated secant steps. $F$ is linear across them to $10^{-17}$.
- **Midpoint–radius squaring.** On a Lohner step hull near the axis, $u/\rho \in [0.45, 2.2]$ and arb
  computes its square as $[-1.31, 4.83]$, so $\gamma_{TT} \ni 0$. Use $\cos\alpha = \sqrt{1 - \sin^2\alpha}$
  there, and short steps near the axis.
- **Interleaved logs.** `tail -f` on two logs printed results under the wrong run once, and two
  secant seeds were swapped. The secant step corrected it, but read logs per file.

## Proved versus computed versus open

| Status | Items |
| --- | --- |
| proved | Lemmas 1–3; the topological step (the MOTS separates the shell) |
| computer-assisted (new) | Theorem M ($t_1$, $t_1'$) |
| inherited | Theorems C, F′; Proposition T, CK′ certificates; AM barrier theorem |
| open | Theorem N (local lower bound); the first common apparent horizon of the whole slice; event horizons; a $d$-interval version; a closed form |

## Next single decisive task

**First, Theorem N** by the convexity argument of notes/20 §4: second $z_0$-jets on Lohner pieces plus
point-$z_0$ runs, about 15 runs. **Then Theorem F′ for $0 < t \le t_0$.** This turns the local bracket into a global statement: the common
apparent horizon of the whole slice is born in $(t_0, t_1']$. It needs two things:

1. a CK certificate on $B(0, 1.05)$ minus the throat balls, needed only for $t \lesssim 10^{-5}$;
2. the F′ sweep rerun for the level-$h$ curves with $\lvert h\rvert \le \sup\lvert R\rvert$ and $t$ as a parameter.
   Its exclusion events are open conditions, and the outermost common MOTS is unique, hence
   symmetric. So the reduction of notes/12 §1 carries over unchanged.

## Reproduce

```bash
python -m pytest tests/test_bl_barrier.py
```

```bash
python -m research.calc.bl_barrier
```
