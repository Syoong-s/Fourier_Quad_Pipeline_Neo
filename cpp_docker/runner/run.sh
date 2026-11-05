#!/usr/bin/env bash
set -euo pipefail
FQ_RUNNER_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$FQ_RUNNER_DIR/common.sh"
if [[ "${1:-}" == --check ]]; then
    fq_load "${2:?Usage: bash run.sh --check /path/cpppipeline.env}"
    fq_check
    exit
fi
fq_load "${1:?Usage: bash run.sh /path/cpppipeline.env}"
# Slurm supplies task count and placement from the user's job allocation.
exec srun "${FQ_SRUN_ARGS[@]}" --export=ALL --mpi=pmi2 --kill-on-bad-exit=1 \
    bash "$FQ_RUNNER_DIR/runtime.sh" "$FQ_ENV_FILE" run
