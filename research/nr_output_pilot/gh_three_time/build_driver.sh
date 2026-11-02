#!/bin/bash
set -euo pipefail
cd /sxscollaboration/spectre/build
export CCACHE_DISABLE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
rm -f CMakeFiles/SpectrePch.dir/cmake_pch.hxx.pch
make -j1 SpectrePch/fast
# CMake fast targets reuse existing, hash-pinned libraries and avoid unrelated rebuilds.
make -j1 CharmModuleInit/fast AmrEvents/fast GeneralRelativityAnalyticData/fast ApparentHorizonFinderCriteria/fast ControlSystem/fast EvolutionTriggers/fast EvolveGhBinaryBlackHole/fast
