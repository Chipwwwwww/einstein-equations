# Binary GH executable/parser pass; startup blocks three physical observations

**2026-11-02 continuation of PR #13. EXECUTED:** the actual binary-black-hole GH
target was built and its C++ input parser passed. The attempted execution then
segfaulted during global-cache construction, before a successful numerical
initial-data import was demonstrated. **No GH HDF5 files or physical-time
observations were produced.** No continuum residual or exact-vacuum claim.

The research branch started clean at remote/local head
`2f7ef4d56974a64a006fc7e902be5cb5d2c67bfb`. The previous three XCTS output-file
checksums were rechecked and remain identical; see
[`evidence/preserved-xcts-checksums.json`](evidence/preserved-xcts-checksums.json).
The old elliptic observations are still not physical-time observations.

## Gate results

| Requirement | Result | Evidence / limitation |
| --- | --- | --- |
| GH executable obtained | PASS, EXECUTED | Actual `EvolveGhBinaryBlackHole`, built, linked, help launched; evolution startup fails below |
| Actual C++ input parsed | PASS, EXECUTED | `gh-parser.log`; run also prints `Option parsing completed` |
| XCTS initial data imported | UNKNOWN / not reached | Crash in GlobalCache construction; no successful import or missing-field diagnosis |
| Three physical-time GH observations saved | FAIL, EXECUTED directory inspection | No HDF5 files; actual time list empty |
| All 50 GH fields present | UNKNOWN / not observed | Requested in parsed input; no datasets exist |
| Gauge and damping history recoverable | UNKNOWN | Prescriptions retained in parsed input, but no evolving fields/FoT history |
| Maps valid between observations | UNKNOWN / not reached | No new observations; old stationary XCTS check does not answer this |
| Boundary/excision geometry recoverable | UNKNOWN for GH output | Defined in parsed input; no serialized GH domain or boundary data produced |
| Sampled residual computable | FAIL from these outputs | No GH time samples; no sampled residual evaluated |
| Rigorous continuum residual certified | FAIL / not established | Not attempted or implied by this output pilot |

## Pinned executable and resource decision

- SpECTRE **2026.10.28**, source commit
  `ec2a35874d4a5105521a563498fdc8f33d784a82`.
- Original official amd64 image manifest:
  `sha256:2906ee3720535453985fb1935ded92407412e42ed4fb2b6423cc429cb7fc4fd1`.
  Existing image/source were reused, not repeatedly unpacked in full.
- Compiler: **Ubuntu Clang 21.1.8**, C++20, Release `-O3 -march=native
  -mno-avx512f`; exact definitions/includes/options and link command in
  `evidence/build-support-check.log` and `target-dry-run.log`.
- Prior image-pinned dependency/package inventories remain in the parent
  `evidence/` directory. Build tools/headers/static libraries missing from the
  earlier trimmed runtime were selectively restored from digest-verified
  layers; selected extraction total 1,707,759,832 bytes plus small linker support.
  Layer records and restoration scripts are included. The layer JSON summaries
  were reconstructed from their checksum-verifying extraction logs because
  the linker-support restoration reused the original metadata filename.
- Host: **AMD EPYC 9V74**, five visible CPUs, four-core cgroup quota,
  **16 GiB RAM**, initially approximately **13 GiB free disk** this session.
  No external compute allocation/access was configured.
- Before building: **20-minute absolute build deadline**, one compiler process,
  **12 GiB container limit**, four-core quota, **4 GiB minimum free disk**.
  ccache disabled to avoid its documented multi-GB cache. Targeted reuse was
  estimated as five missing libraries plus one large translation unit; fit
  within the 20-minute ceiling was uncertain and tested, not assumed.
- Actual successful build interval including corrected failures:
  **06:33:07--06:45:53 UTC (766 seconds)**. Recorded cgroup memory peak
  **8,686,776,320 bytes (8.09 GiB)**, no OOM events. Disk remained above 8 GiB
  through parser/startup tests. The run had a separate 300-second timeout with
  10-second kill grace, four threads, the same memory cap and disk guard.

Binary: **434,484,680 bytes**, SHA256
`468efeedf0e0c7ec074a938ba1d7ba2e69979c398fb3cd90eda18e7ac9ccb2e6`.
Host copy: `/workspace/scratch/spectre-pilot/gh-executable-gate/bin/EvolveGhBinaryBlackHole`.
Container: `spectre-output-pilot:/sxscollaboration/spectre/build/bin/EvolveGhBinaryBlackHole`.
It is a native build for this host, not a claim of portability to another CPU.

The pinned release has no binary assets; the inspected nearby releases publish
CCE archives, not the BBH GH executable. The existing container offered the
least costly available local path. No other physical system was substituted.

## Build and parser commands actually executed

Within the reused container build tree `/sxscollaboration/spectre/build`:

```sh
export CCACHE_DISABLE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
make -j1 CharmModuleInit/fast AmrEvents/fast GeneralRelativityAnalyticData/fast \
  ApparentHorizonFinderCriteria/fast ControlSystem/fast EvolutionTriggers/fast \
  EvolveGhBinaryBlackHole/fast
```

The first compile rejected the packaged PCH: target CPU `graniterapids` versus
this host's `znver4`. The corrected attempt removed only that generated PCH,
rebuilt `SpectrePch/fast`, then reran the targeted build (see `build_driver.sh`).
It compiled the large BBH object successfully. Link metadata then encountered
Git's dubious-ownership guard. The final **link-only** retry used:

```sh
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory
export GIT_CONFIG_VALUE_0=/sxscollaboration/spectre
make -j1 EvolveGhBinaryBlackHole/fast
EvolveGhBinaryBlackHole --help +p1
sha256sum /sxscollaboration/spectre/build/bin/EvolveGhBinaryBlackHole
```

Each build retry used only the remaining time to the original deadline. Logs
distinguish the resolved PCH and ownership failures from the unresolved startup
failure. The CMake `fast` targets reused existing pinned static archives rather
than triggering unrelated full-library rebuilds. Native PCH, main-module object
and the five missing libraries were rebuilt. **Scientific source was not
edited.** The final source checkout reports only a missing container-recipe file
already absent in the trimmed packaging, not a modified numerical source.

After the executable was obtained:

```sh
spectre bbh start-inspiral /work/id-smoke/InitialData.yaml \
  --refinement-level 0 --polynomial-order 5 --no-submit --num-procs 4 \
  --param FinalTime=0.0002 -o /work/gh-generated
python /work/gh-executable-gate/prepare_input.py \
  /work/gh-generated/Inspiral.yaml /work/gh-three-time/Inspiral.yaml
spectre validate /work/gh-three-time/Inspiral.yaml
```

`--no-submit` generated the input; its log's `Run ...` message is not evidence
of execution. The actual later run is logged separately.

`Inspiral.generated.yaml` is preserved. `Inspiral.executed.yaml` is the input
that **passed the C++ parser and was supplied to the execution attempt**.
`prepare_input.py` records every changed YAML path and old/new value in
`evidence/input-changes.json`. It changed slab size to `1e-4`, initial step to
`1e-5`, output names, output events, termination at `2e-4`, and disabled the
user CCE extraction event. Replacing the slab events removes history cleaning
to retain the short FoT history. Domain, gauge, damping, spatial discretization,
initial data and control-system parameters were unchanged.

The parser accepted double-precision output for `SpacetimeMetric`, `Pi`, `Phi`,
`GaugeH`, `SpacetimeDerivGaugeH`, `Lapse`, `Shift`, `SpatialMetric` and
`GaugeConstraint`. This establishes parsing, **not successful export**.

## First unresolved failing command

In container working directory `/work/gh-three-time`:

```sh
timeout --signal=TERM --kill-after=10 300 \
  EvolveGhBinaryBlackHole --input-file Inspiral.yaml +p4
```

It returned **exit 1 immediately** at **06:49:11 UTC**. SpECTRE caught a
segmentation fault; the relevant stack is:

```text
ah::Storage::LockedPreviousSurface<Frame::Distorted>::LockedPreviousSurface(...&&)
Parallel::GlobalCache_detail::make_mutable_cache_tag_storage<...>
Parallel::GlobalCache<EvolutionMetavars>::GlobalCache(...)
Parallel::Main<EvolutionMetavars>::Main(...)
```

The complete backtrace is in `evidence/run.log`. Source location to investigate:
[Storage.cpp, move constructor](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/ParallelAlgorithms/ApparentHorizonFinder/Storage.cpp#L150).
The cause is **UNVERIFIED**. The backtrace does not establish an upstream bug,
an incompatible reused artifact, or a missing initial-data field. In particular,
no successful import, excision-speed test, boundary test or evolution step was
demonstrated. Dependent gates stopped; no substitute executable or physics
workaround was attempted.

## Outputs and smallest resolving test

Actual run directory: `/workspace/scratch/spectre-pilot/gh-three-time/`.
It contains input, recorded changes and `run.log`, **zero HDF5 files**.
The planned times `[0, 0.0001, 0.0002]` never became observations. Consequently
dataset names/shapes/types, evolving mesh/domain compatibility, maps between
times, reduction constraints and sampled GH residuals have **no new results**.
The explicit empty inventory is `evidence/gh-output-inventory.json`.
File/binary/log hashes are in `SHA256SUMS` and `evidence/result.json`.

The smallest resolving experiment is to reproduce/localize the same global-cache
startup failure with a fully consistent target build (and useful debug symbols),
before any numerical evolution. A clean build is needed to distinguish source
behavior from the targeted reuse of packaged archives. Once startup succeeds,
reuse this parser-validated three-time input; do not extend it to a merger run.

The separate XCTS normalization investigation is
[`HAMILTONIAN.md`](HAMILTONIAN.md); it does not replace the failed three-time gate.

## Validation and primary sources

Executed: actual target build/link/help; C++ parser; bounded failed startup;
HDF5 absence inspection; numerical recomputation of the old XCTS diagnostic;
preservation checksums. Scripts were syntax-checked, the preparation script's
output was checked against the parsed input, and artifact checksums checked.
No old mathematical certificate or historical repository test suite was rerun.

Pinned primary source: [installation/build requirements](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/docs/Installation/Installation.md),
[build options](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/docs/Installation/BuildSystem.md),
[BBH target/dependencies](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/src/Evolution/Executables/GeneralizedHarmonic/CMakeLists.txt),
[pipeline](https://github.com/sxs-collaboration/spectre/blob/ec2a35874d4a5105521a563498fdc8f33d784a82/support/Pipelines/Bbh/Inspiral.py).
