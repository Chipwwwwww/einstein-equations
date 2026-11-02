# 29 — Variable-background H³ inverse-metric functional calculus and failure tests

**2026-11-02. Outcome: a verified, standard Banach-algebra lemma for a spatially nonconstant reference matrix, with exact-rational certificate and negative controls. No new result relative to established Sobolev functional calculus.** This is deliberately separated from the nonflat Kasner full 50-field finite-slab certificate of note28. Neither result certifies the Brill–Lindquist merger domain.

## Theorem (variable-background inverse / Lorentzian tube)

Let \(M=(\mathbb R/2\pi\mathbb Z)^3\) and equip \(H^3(M;\mathrm{Mat}_{4\times4})\) with the normalized Fourier-Sobolev Frobenius norm. Let \(\bar g(x)\) be a smooth Lorentzian 4×4 metric field with positive spatial 3×3 block and with certified numerical bounds

\[
\|\bar g^{-1}\|_{H^3}\le B,\qquad
\sup_x\|\bar g^{-1}(x)\|_{\mathrm{op}}\le J,\qquad
\lambda_{\min}(\bar g_{ij}(x))\ge \sigma>0.
\]

Let \(R>0\) satisfy

\[
q=48^2 B R<1,\qquad 6JR<1,\qquad 6R<\sigma.
\]

Then **every symmetric 4×4** \(g\in H^3\) satisfying \(\|g-\bar g\|_{H^3}\le R\) is Lorentzian with positive spatial metric and

\[
\boxed{\|g^{-1}\|_{H^3}\le \frac{B}{1-q}},\qquad
\boxed{\|g^{-1}-h^{-1}\|_{H^3}\le48^2\Bigl(\frac{B}{1-q}\Bigr)^2\|g-h\|_{H^3}}
\]

for every \(g,h\) in the same ball. All constants are computable directly from a **spatially varying reference**; no assumption that \(\bar g\) is close to Minkowski in \(H^3\) is made beyond the separate smallness of \(R\).

**Proof.** The matrix-field Fourier Sobolev norm has the Banach-algebra bound \(\|AB\|_3\le48\|A\|_3\|B\|_3\) by the same convolution argument as note25, using Frobenius submultiplicativity. Write \(\delta=g-\bar g\) and \(E=\bar g^{-1}\delta\). Then \(\|E\|_3\le48BR\). For all \(n\ge1\), \(\|E^n\bar g^{-1}\|_3\le 48^n\|E\|_3^nB\le q^nB\). Consequently

\[
(\bar g+\delta)^{-1}=\sum_{n=0}^{\infty}(-E)^n\bar g^{-1}
\]

converges absolutely in \(H^3\) with norm at most \(B/(1-q)\). The inverse identity \(g^{-1}-h^{-1}=g^{-1}(h-g)h^{-1}\), multiplied twice in \(H^3\), gives the displayed Lipschitz bound. By the embedding \(\|\delta\|_{L^\infty,op}\le6\|\delta\|_3\le6R\), the pointwise matrices \(\bar g+\theta\delta=\bar g(I+\theta\bar g^{-1}\delta)\) are invertible for all \(\theta\in[0,1]\); since they are symmetric their Lorentzian inertia cannot change. Weyl's inequality gives \(\lambda_{\min}(g_{ij})\ge\sigma-6R>0\). QED.

**Explicit nonconstant example.** Take \(\bar g=\mathrm{diag}(-f,f,1,1)\), \(f=1+(1/1000)\cos x_1\). The exact Fourier norm is \(\|\cos x_1\|_3=2\), so \(\|f-1\|_3=1/500\). The scalar inverse Neumann estimate gives \(\|f^{-1}-1\|_3\le (1/500)/(1-48/500)\), so \(\|\bar g^{-1}\|_{H^3}\le4+2[(1/500)/(1-48/500)]<5\). Pointwise \(J=2\), \(\sigma=999/1000\) are certified. With \(R=10^{-6}\), \(q=36/3125\), the full-ball inverse bound is \(15625/3089\), with Lipschitz bound \(562500000000/9541921\). Both are independent of grid sampling. Test programs deliberately reject: a Neumann nonconvergence, an uncertified Lorentzian homotopy, a spatial-positivity failure, and invalid signed radii.

**Limitations.** This is an established functional-calculus technique with very pessimistic factors 2304; it does not by itself bound the complete GH algebraic source on an **arbitrary** reference metric (one still needs the lapse branch, coefficient Jacobians and domain constants), and does not turn a punctured black-hole domain into a torus. On a curved physical spatial domain, substitute **proved** Sobolev algebra and embedding constants (not the periodic 48 and 6 without justification).

Reproduce:

    python -m research.calc.gh_reference_neumann
    python -m pytest -q tests/test_gh_reference_neumann.py

This lemma isolates a potential future improvement: sharper reference-metric-adapted matrix products and localized Sobolev constants can reduce the current overconservative \(48^2\) growth; no unproved sharp factor is used in the certificate.