# What the earlier XCTS 0.00747 diagnostic measures

**DOCUMENTED:** At pinned commit `ec2a35874d4a5105521a563498fdc8f33d784a82`,
`ObserveNorms` implements `L2Norm` as the **unweighted nodal RMS**:

```math
\operatorname{L2Norm}(H)=\sqrt{N^{-1}\sum_{p=1}^{N} H_p^2}.
```

There are no quadrature, Jacobian or proper-volume weights in this reduction.
The reported `Volume` is a separate column, not its denominator. DG face nodes
are counted as stored, and the distribution of nodes changes the statistic.
This is not a continuum geometric $L^2$ norm. Even `L2IntegralNorm` is a distinct
operation; its availability does not change the executed `L2Norm` choice.

**DOCUMENTED:** The XCTS computed field called `HamiltonianConstraint` is

```math
H_s=4\widetilde\Delta\psi-\tfrac12\left[
\psi\widetilde R+\psi^5(K^2-K_{ij}K^{ij}-16\pi\rho)\right].
```

The conformal scalar-curvature identity gives the **continuum normalization**

```math
H_s=-\tfrac12\psi^5 C_{\rm ADM},\qquad
C_{\rm ADM}=R+K^2-K_{ij}K^{ij}-16\pi\rho.
```

Thus the saved field is conformally weighted, with the factor of one half
explicitly noted in the source as matching SpEC. It is not the raw unweighted
ADM scalar, a relative error, or a residual divided by a local curvature scale.
Here the matter density is zero. Dimensional quantities use the calculation's
input mass normalization; measured total horizon mass differs from its target.

**NUMERICAL:** For the old converged observation (AMR label 1, $N=20304$):

| Calculation from saved data | Value |
| --- | ---: |
| Reported `L2Norm(HamiltonianConstraint)` | 0.007465810823616376 |
| Independently recomputed RMS of the saved field | 0.007465810823616378 |
| Maximum absolute saved field | 0.10525351344973821 |
| RMS after continuum-normalization rescaling $-2H_s/\psi^5$ | 0.008022960868554517 |
| RMS of $R+K^2-K_{ij}K^{ij}$ using saved SpatialRicci/metric/K | 0.00802192191930164 |
| Maximum pointwise difference between the latter two constructions | 0.0011294279877320512 |

The last two constructions do **not** agree to roundoff on this coarse output.
The continuum identity alone does not force separately discretized derived
quantities to obey it. The source of this numerical disagreement was not
resolved here; it must not be silently removed or described as a rigorous
constraint-error bound. Neither the 0.00747 statistic nor the Newton algebraic
residual certifies continuum initial constraints.

Recompute with `python check_hamiltonian.py` in the retained workspace; the
script reads the preserved old XCTS HDF5 files. Exact results are in
`evidence/hamiltonian-definition.json`.

Primary source: [norm definition and reduction](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/ParallelAlgorithms/Events/ObserveNorms.hpp#L65),
[local sum of squares](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/ParallelAlgorithms/Events/ObserveNorms.cpp#L69),
[XCTS constraint field](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/PointwiseFunctions/Xcts/SpacetimeQuantities.cpp#L362).
