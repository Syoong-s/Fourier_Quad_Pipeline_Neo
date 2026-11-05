#!/usr/bin/env bash
# Container entry. run.sh invokes this once per MPI task.
set -euo pipefail
FQ_RUNNER_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$FQ_RUNNER_DIR/common.sh"
fq_load "$1"
action="$2"
# Keep the container environment explicit, including after site module loads.
while IFS= read -r name; do
    case "$name" in
        APPTAINERENV_*|SINGULARITYENV_*|APPTAINER_BIND|APPTAINER_BINDPATH|SINGULARITY_BIND|SINGULARITY_BINDPATH)
            unset "$name" ;;
    esac
done < <(compgen -e)
unset LD_PRELOAD
fq_export() { export "APPTAINERENV_$1=$2" "SINGULARITYENV_$1=$2"; }
if [[ "$action" == run ]]; then
    # These are per-task MPI connection values, not a Slurm allocation check.
    while IFS= read -r name; do
        case "$name" in SLURM_*|PMI_*|PMI2_*) fq_export "$name" "${!name}" ;; esac
    done < <(compgen -e)
    fq_export OMPI_MCA_ess pmi
fi
for assignment in "${FQ_CONTAINER_ENV[@]}"; do
    fq_export "${assignment%%=*}" "${assignment#*=}"
done
fq_export HOME /tmp
args=(exec --cleanenv --no-eval --no-home --no-mount hostfs,bind-paths,cwd
      --pwd /workspace/src_pipe)
source_mode=ro
[[ "$action" != build ]] || source_mode=rw
args+=(--bind "$FQ_SOURCE_HOST:/workspace/src_pipe:$source_mode"
       --bind "$FQ_CONFIG_HOST:/fq-config/pipeline.ini:ro")
for ((i=0; i<${#FQ_BIND_SOURCES[@]}; ++i)); do
    args+=(--bind "${FQ_BIND_SOURCES[i]}:${FQ_BIND_TARGETS[i]}:${FQ_BIND_MODES[i]}")
done
if [[ "$action" == run ]]; then
    exec "$FQ_RUNTIME" "${args[@]}" "$FQ_SIF" \
        /workspace/src_pipe/Fourier_Quad_Pipe --config /fq-config/pipeline.ini
elif [[ "$action" == build ]]; then
    exec "$FQ_RUNTIME" "${args[@]}" "$FQ_SIF" \
        bash -c 'set -euo pipefail; make clean; exec make -j"$1" "STACK_PREFIX=${CPP_STACK_PREFIX}"' \
        fq-build "$FQ_BUILD_JOBS"
else
    printf 'Unknown runner action: %s\n' "$action" >&2
    exit 2
fi
