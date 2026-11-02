# SpECTRE output pilot: XCTS smoke completed; GH evolution blocked

**2026-11-02. EXECUTED / NUMERICAL. No physical evolution interval, common-horizon
formation, GH evolution residual, continuum enclosure or exact vacuum certificate.**

This is a new binary initial-data calculation, not SXS:BBH:1132. The independent
branch starts at PR #12 head `f235eff91a862afe74d369b219f7fd151771d9c4`.
Inherited mathematical certificates were neither changed nor independently reviewed.

## Resource decision and first blocker

The host exposes five logical CPUs but cgroup `cpu.max=400000 100000`, hence a
four-core quota; `memory.max=17179869184` (16 GiB), no swap, and approximately
30 GiB initially free disk. No scheduler allocation or hard session walltime was
exposed; shell CPU time was unlimited. The smoke container was capped at four
cores, 10 GiB RAM and 256 PIDs; the solve had a 300-second timeout. Sampled memory
use during the solve was about 1.44 GiB, **not a measured peak**.

The published merger experiment [1] cost 58,000--117,000 core-hours on 560 cores.
These are not cost predictions or lower bounds for a short coarse run. They do
not justify launching a production merger on this machine. No evolution with a
demonstrated resource-fit configuration was available.

**Actual blocking dependency:** the official deployment distribution includes
`SolveXcts` and Python support, but no `EvolveGhBinaryBlackHole` or
`EvolveGhSingleBlackHole`. The attempted `bbh start-inspiral` stopped with
`ValueError: Executable not found: EvolveGhBinaryBlackHole` (exit 1). A GH C++
parser test and GH initial-data import therefore did not run. A full source
build was not attempted; this is not a proof that a build or any merger is
impossible with four cores. The early-stop path was a bounded actual XCTS test.

## Exact software and inputs

- Release **2026.10.28**, source commit
  `ec2a35874d4a5105521a563498fdc8f33d784a82`, verified against the release tag and
  the executable/HDF5 build information.
- Official Linux/amd64 OCI manifest:
  `sha256:2906ee3720535453985fb1935ded92407412e42ed4fb2b6423cc429cb7fc4fd1`.
  This pins the packaged binaries and dependency filesystem, not just a mutable tag.
- `evidence/provenance.json` records the actual solver binary hash, resource
  limits and local runtime preparation. `evidence/build-libraries.txt`,
  `python-freeze.txt` and `dpkg-versions.txt` record numerical/build dependencies.
  Among them: Charm++ 8.0.0, Boost 1.90.0, Blaze 3.9, GSL 2.8, C++ HDF5 1.14.6,
  Clang 21.1.8, Python 3.14.4, NumPy 2.5.3 and h5py 3.16.0.
- `InitialData.generated.yaml` is the unchanged pipeline-rendered input;
  `InitialData.yaml` is the **actual parser-validated and executed input**.
  The only solver-input changes were observing every elliptic AMR iteration,
  adding Hamiltonian/momentum fields and replacing ADM-integral observation
  with inexpensive lapse/constraint norms. Physical parameters were unchanged.

The official image is 5,227,878,922 compressed bytes across 29 layers. The local
Docker daemon uses `vfs`. To avoid repeated parent-layer copies, official layers
were streamed, checksum-verified, filtered to a runtime subset, and flattened.
Build-only files and selected optional Python packages were omitted; numerical
sources and binaries were not edited. A first extraction stopped at its 5 GiB
guard, then build objects/static archives were removed and the last layer was
re-extracted and verified. An omitted dynamic-loader path was repaired with the
identical loader already extracted from the official image. Layer verification
and the derived image ID are recorded. Use the full pinned image for portable
reproduction; allow substantially more Docker storage with `vfs`.

## Commands that actually ran

Inside the prepared runtime, with `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`:

```sh
spectre bbh generate-id -q 1 --chi-A 0 0 0 --chi-B 0 0 0 \
  -D 12 -e 0 -L 0 -P 5 --no-control --no-submit \
  -o /work/id-generated --num-procs 4
spectre validate /work/id-generated/InitialData.yaml
# Apply the observation-only changes described above; checked-in actual input follows.
spectre validate /work/id-smoke/InitialData.yaml
cd /work/id-smoke
timeout --signal=TERM 300 SolveXcts --input-file InitialData.yaml +p4
timeout --signal=TERM 120 spectre bbh postprocess-id InitialData.yaml \
  --control False --horizon-l-max 8
python /work/inspect_output.py /work/id-smoke --spectre-maps
# Failed before any evolution:
spectre bbh start-inspiral /work/id-smoke/InitialData.yaml \
  --lev 0 --no-submit -o /work/gh-not-run
```

Parser, generation, solve and postprocessing logs, plus a labeled failed-handoff
excerpt, are in `evidence/`. Trailing whitespace is removed; raw log hashes and
external locations are recorded in provenance.
`SolveXcts --version` was also tried and rejected: this executable reports build
information at startup, while `spectre --version` is the Python CLI option.

For reproduction on a machine with suitable Docker storage:

```sh
bash research/nr_output_pilot/reproduce.sh /absolute/path/to/new-empty-output
```

The wrapper uses the full pinned official image, the checked-in input, the same
thread count and timeouts. The full-image reproduction wrapper was syntax-checked but was not replayed
as a second solve on this vfs machine. Floating-point/file hashes need not be bitwise identical on
another architecture or reduction schedule. No HDF5 outputs are committed.

## Binary and numerical initial-data result

**EXECUTED:** target masses 0.5 each, zero target spins; centers at `(+-6,0,0)`;
conformal masses 0.41; input orbital frequency `0.02179820375911994` and expansion
`-0.00012832805047143658`. The domain has two excision spheres of grid radius
0.7626, surrounding shells of radius 3.2, envelope radius 48, and outer boundary
radius 100000 in input coordinates. Separation 12 was a smaller-than-the-D=16
documentation example chosen **only for output smoke testing**. Its suitability
for a scientifically useful merger has not been demonstrated.

The PN orbital-parameter generator emitted `ODEintWarning: Excess work done`.
No eccentricity reduction, physical mass/spin control loop, resolution study or
junk-radiation measurement was performed. The conformal-background spins being
zero does not make the measured horizons exactly nonspinning.

The solver exited 0 after **82 seconds**, at 05:34:01 UTC. Newton converged in
six iterations to algebraic residual `8.460949365323763e-10`. The final numerical
constraint diagnostics saved under `ConstraintDiagnostics.dat` are:

- min lapse: `0.6052010356096055`;
- solver-reported `L2Norm(HamiltonianConstraint)`: `0.007465810823616376`;
- momentum components: `0.0014476238328703652`, `0.0014789225379495824`,
  `0.001466714760087738`;
- maximum sampled absolute Hamiltonian constraint: `0.10525351344973821`.

The small algebraic residual is **not** a bound on continuum Einstein constraints.
The coarse physical-constraint diagnostics are substantial.

Postprocessing exited 0 and found the two **individual** apparent horizons:

| Horizon | Measured Christodoulou mass | Measured dimensionless spin z |
| --- | ---: | ---: |
| AhA | 0.5306984682689151 | -0.0026973695369417175 |
| AhB | 0.5306976998090253 | -0.0027587647148993185 |

Measured total mass differs from the target normalization; do not identify
input separation 12 with exactly twelve measured total masses. **No common
horizon was sought or found, and no formation time exists for this experiment.**

## Actual HDF5 output and the two observations

Host outputs: `/workspace/scratch/spectre-pilot/id-smoke/` (container `/work/id-smoke/`).
Checksums and a complete dataset inventory are in `evidence/output.json`.

| File | Bytes | SHA256 |
| --- | ---: | --- |
| BbhVolume0.h5 | 20351695 | `337097246b72d30672e2c137ccaeee41b9493f4cc32e1d280db5a73d186ee03d` |
| BbhReductions.h5 | 12077512 | `de2bb4f3944300def7e434879e206783a850e1bda00e79b1f99b8b7f8ed8f97f` |
| Horizons.h5 | 124904 | `ca3248a6d9700d0fd490a0e397983448bf568499b28aacbe61cee08acd31cb90` |

The saved volume observations have values **0 and 1**: these are **elliptic AMR
iteration labels on one initial slice**, not physical times. Observation 0 is
the initial guess; observation 1 contains the converged solve. Both have 84
elements, 20304 grid points, 44 blocks, tensor-product Legendre bases and
Gauss--Lobatto quadrature. Extents are `(6,6,6)`, `(6,6,7)`, `(6,6,10)`.
No positive-time spacetime interval was saved.

Paths below are relative to `/VolumeData.vol/ObservationId.../`; the exact IDs
and every component name are in the JSON inventory. All stated shapes occur at
**both** observations.

| Actual dataset names/families | Shape per component | Type |
| --- | --- | --- |
| `Lapse`, `ConformalFactor`, `HamiltonianConstraint` | `(20304,)` | float64 |
| `Shift_x/y/z`, `ShiftExcess_x/y/z`, `MomentumConstraint_x/y/z` | `(20304,)` | float64 |
| `SpatialMetric_xx/yx/yy/zx/zy/zz` | `(20304,)` | float64 |
| `InverseSpatialMetric_*`, `ExtrinsicCurvature_*`, `SpatialRicci_*` (six each) | `(20304,)` | float64 |
| `SpatialChristoffelSecondKind_*` (18) | `(20304,)` | float64 |
| `InertialCoordinates_x/y/z`, `RadiallyCompressedCoordinates_x/y/z` | `(20304,)` | float64 |
| `total_extents` | `(252,)` | uint64 |
| `bases`, `quadratures` | `(252,)` | int32 |
| `grid_names` | `(2173,)` | int8 (serialized names) |
| `functions_of_time` | `(42742,)` | int8 (serialized object) |
| `connectivity` | `(108000,)` | int32 |
| `BlockId`, `ElementId` | `(12000,)` | uint64 |

At the parent volume group: `domain` has shape `(35475,)`, int8;
`global_functions_of_time` has shape `(42742,)`, int8. Reduction datasets include
`ConstraintDiagnostics.dat` `(2,8)` and `NewtonRaphsonResiduals.dat` `(7,5)`, both
float64. Horizon diagnostics `AhA.dat`, `AhB.dat` are `(1,9)`, float64, accompanied
by surface coefficients and coordinates.

## Maps, boundaries and diagnostics

**EXECUTED / NUMERICAL:** SpECTRE bindings deserialized the saved domain and all
four functions `ShapeA`, `SizeA`, `ShapeB`, `SizeB`. Element IDs, meshes, block
maps and functions reconstructed every saved inertial coordinate at both
observations with maximum sampled difference **0**. Static-map finite-difference
velocity and function first derivatives were also zero.

The observation-specific functions cover only `[0,2.220446049250313e-14]` and
`[1,1.0000000000000222]`. **They do not cover the gap.** The separate global
objects cover `[0,infinity)` in this stationary-shape-map initial-data case and
were successfully evaluated at intermediate arguments. These arguments are
not a simulated physical time interval. Complete histories of dynamically
updated maps, control-system cleaning, restarts and remeshing remain untested.

Element and boundary shapes are recoverable from the deserialized block/element
maps, their faces and the recorded domain input. Initial excision and outer
geometry are therefore available. This experiment did not explicitly enclose
boundary normals/fluxes or excision characteristic speeds, nor test moving
excision/outer boundaries. A found apparent horizon is not an outflow proof.

Numerical checks on all saved nodes found positive lapse, positive-definite
spatial metric and Lorentzian signature for the four-metric reconstructed using
`ShiftExcess`, as the pipeline does for its inertial shift. This checks algebraic
initial-slice geometry only. **GH reduction constraints and an evolution
residual could not be evaluated:** no evolved Pi/Phi or physical time samples
were produced.

## Residual-readiness gate and next experiment

| Requirement | Actual result |
| --- | --- |
| XCTS parser, solve, double output and individual horizons | EXECUTED |
| GH input import and C++ parser | BLOCKED: missing executable |
| SpacetimeMetric, Pi, Phi, GaugeH, SpacetimeDerivGaugeH | Absent from actual output |
| Current supported gauge-output names | DOCUMENTED in pinned executable source [2]; not parser-tested |
| Definite initial spatial domain/coordinates | EXECUTED reconstruction |
| Moving mesh velocity and complete evolving map histories | UNVERIFIED |
| Merger-spanning output | NOT RUN |
| Sampled GH evolution residual | NOT COMPUTABLE FROM THESE OUTPUTS |
| Rigorous continuum residual | NOT ESTABLISHED |

`proposed_gh_observation.yaml` is explicitly an **unparsed proposal**, not an
executed configuration. It requests the 50 GH fields, H and spacetime derivative
of H, lapse, shift and spatial metric every synchronized slab in double precision.
The rendered evolution input must also preserve actual damping prescriptions,
boundary conditions, domain/maps, all transition/restart information and horizon
diagnostics. A separate ringdown event is required after the domain transition.

Current output stores 60 float64 scalar components, about 9.75 MB raw tensor
data per observation at this mesh, plus metadata/connectivity. The full GH
state alone costs `400 * number_of_grid_points * number_of_observations` bytes.
At one million points and 401 observations it costs 160.4 GB before gauge,
coordinates and geometry. Frequent merger-window output requires its own budget.

The smallest next falsification experiment is **not another merger run**:
obtain/build the pinned GH executable on a suitable machine, parser-test the
event, import this initial-data format, and save three closely spaced physical
times on a fixed mesh. Independently reconstruct coordinates, lapse/shift and
Pi/Phi, evaluate the actual damped-GH operator, include face jumps, and compare
time differences and constraints. Reject the approach immediately if global
map histories cannot evaluate every intermediate time, or these reconstructions
do not converge under time/mesh refinement. These tests remain numerical.

For a mathematically defined continuum reconstruction, select explicit spatial
and temporal polynomials and a regularity-preserving treatment of DG faces and
time/domain transitions. Enclose maps/inverses and all GH residual derivatives
required by the chosen theorem, including boundary and initial-constraint
defects, using validated arithmetic over entire cells/slabs. Samples, elliptic
solver tolerance and spectral-tail fits alone cannot supply those bounds.

## Primary sources inspected

1. [Lovelace, Nelli et al., merger demonstration and costs](https://arxiv.org/abs/2410.00265),
   with [actual ancillary inputs](https://arxiv.org/src/2410.00265v2/anc).
2. [Pinned binary GH executable and observable fields](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/Evolution/Executables/GeneralizedHarmonic/EvolveGhBinaryBlackHole.hpp).
3. [Installation and packaged-executable limitations](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/docs/Installation/Installation.md),
   [build requirements](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/docs/Installation/BuildSystem.md),
   [container recipe](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/containers/Dockerfile.buildenv).
4. [XCTS pipeline](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/support/Pipelines/Bbh/InitialData.py),
   [initial template](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/support/Pipelines/Bbh/InitialData.yaml),
   [inspiral template](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/support/Pipelines/Bbh/Inspiral.yaml).
5. [HDF5 volume format and map-history distinction](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/IO/H5/VolumeData.hpp),
   [actual GH operator](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/Evolution/Systems/GeneralizedHarmonic/TimeDerivative.tpp).
