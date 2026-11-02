#!/usr/bin/env bash
# Requires sufficient image storage; vfs layer copying can exceed small disks.
# Reproduces XCTS output only. This script does not launch GH evolution.
set -euo pipefail
pilot_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
output_dir=${1:?Usage: reproduce.sh NEW_OUTPUT_DIRECTORY}
mkdir -p -- "$output_dir"
output_dir=$(cd -- "$output_dir" && pwd)
if [[ -n $(ls -A -- "$output_dir") ]]; then
  echo 'Refusing to overwrite a nonempty output directory.' >&2
  exit 1
fi
image='sxscollaboration/spectre@sha256:2906ee3720535453985fb1935ded92407412e42ed4fb2b6423cc429cb7fc4fd1'
docker run --rm --platform linux/amd64 --cpus=4 --memory=10g \
  --mount "type=bind,src=$pilot_dir,dst=/pilot,readonly" \
  --mount "type=bind,src=$output_dir,dst=/work" \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -w /work \
  --entrypoint /bin/bash "$image" -c '
    set -euo pipefail
    cp /pilot/InitialData.yaml .
    spectre validate InitialData.yaml > parser.log 2>&1
    timeout --signal=TERM 300 SolveXcts --input-file InitialData.yaml +p4 > solve.log 2>&1
    timeout --signal=TERM 120 spectre bbh postprocess-id InitialData.yaml \
      --control False --horizon-l-max 8 > horizons.log 2>&1
    python /pilot/inspect_output.py /work --spectre-maps > evidence.json
  '
