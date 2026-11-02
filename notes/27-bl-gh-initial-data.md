# 27 — Exact Brill–Lindquist initial data in harmonic gauge

**2026-11-02. Scope:** Exact initial 50-field GH construction off the punctures. Not a finite-time 3+1 Einstein validation, coordinate-overlap certificate, or merger solution.

## Exact vacuum constraints and harmonic initialization

Let
\[
p=1+\frac1{4r_+}+\frac1{4r_-},\qquad
r_\pm=\sqrt{x^2+y^2+(z\mp d/2)^2},
\]
and $\gamma_{ij}=p^4\delta_{ij}$, $K_{ij}=0$ at $t=0$. Off punctures, $\Delta_{\delta}p=0$. The vacuum Hamiltonian constraint is
\[
{}^{(3)}R+K^2-K_{ij}K^{ij}=-8p^{-5}\Delta p=0;
\]
the momentum constraints vanish identically.

Choose the initial ADM gauge $N=p^{-2}$, $\beta^i=0$, $\partial_tN=\partial_t\beta^i=0$ and a prescribed GH source $H_\mu=0$. Then $\partial_tg_{\mu\nu}=0$ and
\[
g_{\mu\nu}=\mathrm{diag}(-p^{-4},p^4,p^4,p^4),\quad
\Pi_{\mu\nu}=0,\quad \Phi_{i\mu\nu}=\partial_i g_{\mu\nu}.
\]
In particular,
\[
\Phi_{i00}=4p^{-5}\partial_i p,\quad
\Phi_{ijk}=4p^3(\partial_i p)\delta_{jk},\quad
\Phi_{i0j}=0.
\]
Thus all $10+10+30=50$ first-order GH fields are explicitly defined.

Because $\sqrt{-g}=Np^6=p^4$, the exact densitized inverse is
\[
\boxed{\sqrt{-g}\,g^{\mu\nu}
=\mathrm{diag}(-p^8,1,1,1).}
\]
Using the zero initial time derivatives and spatial derivatives of the constant diagonal entries gives
$\partial_\mu(\sqrt{-g}g^{\mu\nu})|_{t=0}=0$ identically. Therefore $\Gamma^\nu=0$, so $\mathcal C_\mu=H_\mu+\Gamma_\mu=0$ **on the initial slice**. Reduction constraints $\partial_i g-\Phi_i=0$ and $\partial_i\Phi_j-\partial_j\Phi_i=0$ hold exactly by definition and commutation of spatial derivatives. Symbolic code tests these continuum identities for the two Newtonian puncture potentials.

## Critical limitations

The initial normal derivative of the GH gauge constraint has not been separately computed from the **full nonlinear** first-order GH evolution; the standard constraint-propagation theorem must still be applied with exact hypotheses. The punctures are excluded and the coordinate data are not global periodic $H^3$ fields. The existing CK-series BL development uses geodesic lapse $1$, whereas this GH construction chooses $p^{-2}$: **a validated positive-time coordinate transition and common-chart overlap have not been constructed**. No finite time slab, moving boundary, Schwarzschild nonzero-error metric, or full merger is certified.

## Reproduce

    python -m research.calc.bl_gh_initial_data
    python -m pytest -q tests/test_bl_gh_initial_data.py

See Friedrich–Rendall (2000), https://arxiv.org/abs/gr-qc/0002074; Lindblom et al. (2006), https://arxiv.org/abs/gr-qc/0512093; and existing repo notes 17, 19 and 22.

**Novelty:** new-to-repository explicit initial-data handoff, not a new mathematical-relativity theorem.
