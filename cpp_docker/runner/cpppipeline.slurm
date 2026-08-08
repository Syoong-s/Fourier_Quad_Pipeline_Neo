#!/usr/bin/env bash

#SBATCH --job-name=cpppipeline
#SBATCH --partition=cpu
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=cpppipeline-%j.out
#SBATCH --error=cpppipeline-%j.err

set -euo pipefail

RUNNER_DIR="${CPP_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${RUNNER_DIR}/cpppipeline.env}"
RUNNER="${RUNNER_DIR}/run-apptainer.sh"

# ==========================================
# Function: Stop a pipeline job with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Load configured host modules
# Method: Validate the Bash module array before loading the requested toolchain
# ==========================================
load_host_modules() {
    local module_declaration

    module_declaration="$(declare -p HPC_MODULES 2>/dev/null)" ||
        die "HPC_MODULES must be a Bash array in ${ENV_FILE}"
    [[ "${module_declaration}" == "declare -a "* ]] ||
        die "HPC_MODULES must be a Bash indexed array in ${ENV_FILE}"

    if ((${#HPC_MODULES[@]} == 0)); then
        return
    fi
    command -v module >/dev/null 2>&1 ||
        die "module is unavailable but HPC_MODULES is not empty"
    module load "${HPC_MODULES[@]}"
}

# ==========================================
# Function: Load paths, modules, and the Slurm PMI2 launcher contract
# Method: Source trusted configuration and validate scalar and array settings
# ==========================================
load_job_environment() {
    local srun_declaration

    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    [[ -x "${RUNNER}" ]] || die "Apptainer runner is not executable: ${RUNNER}"
    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    if [[ -n "${SITE_ENV_SCRIPT:-}" ]]; then
        [[ -r "${SITE_ENV_SCRIPT}" ]] ||
            die "SITE_ENV_SCRIPT is not readable: ${SITE_ENV_SCRIPT}"
        # shellcheck disable=SC1090
        source "${SITE_ENV_SCRIPT}"
    fi

    load_host_modules
    srun_declaration="$(declare -p SRUN_ARGS 2>/dev/null)" ||
        die "SRUN_ARGS must be a Bash array in ${ENV_FILE}"
    [[ "${srun_declaration}" == "declare -a "* ]] ||
        die "SRUN_ARGS must be a Bash indexed array in ${ENV_FILE}"
    [[ "${MPI_LAUNCH_MODE:-}" == "srun" ]] ||
        die "MPI_LAUNCH_MODE must be srun"
    [[ "${SLURM_MPI_TYPE:-}" == "pmi2" ]] ||
        die "SLURM_MPI_TYPE must be pmi2"
    [[ "${CPP_SOURCE_CONTAINER:-}" == /* ]] ||
        die "CPP_SOURCE_CONTAINER must be an absolute path"
    [[ "${CPP_EXPO_LIST_CONTAINER:-}" == /* ]] ||
        die "CPP_EXPO_LIST_CONTAINER must be an absolute path"
    [[ "${CPP_EXECUTABLE:-}" == /* ]] ||
        die "CPP_EXECUTABLE must be an absolute path"
}

# ==========================================
# Function: Launch the compiled pipeline across the Slurm allocation
# Method: Let srun create PMI2 state before each rank enters the same SIF
# ==========================================
launch_pipeline() {
    local pipeline_arguments=("$@")
    local container_command

    [[ -n "${SLURM_NTASKS:-}" ]] ||
        die "SLURM_NTASKS is unavailable; submit this file with sbatch"
    command -v srun >/dev/null 2>&1 || die "srun is unavailable"
    srun --mpi=list 2>&1 | grep -qw "${SLURM_MPI_TYPE}" ||
        die "Slurm does not advertise --mpi=${SLURM_MPI_TYPE}"

    if ((${#pipeline_arguments[@]} == 0)); then
        pipeline_arguments=("${CPP_EXPO_LIST_CONTAINER}")
    fi

    container_command=(
        "${RUNNER}"
        --env-file "${ENV_FILE}"
        --pwd "${PROCESS_DATA_CONTAINER}"
        --
        "${CPP_EXECUTABLE}"
        "${pipeline_arguments[@]}"
    )

    srun "${SRUN_ARGS[@]}" \
        --mpi="${SLURM_MPI_TYPE}" \
        --ntasks="${SLURM_NTASKS}" \
        "${container_command[@]}"
}

# ==========================================
# Function: Validate the immutable runtime and launch a compiled cpp_Standard
# Method: Check the SIF once, verify the executable, then start PMI2 ranks
# ==========================================
main() {
    load_job_environment
    printf 'Slurm job %s: %s nodes, %s ranks\n' \
        "${SLURM_JOB_ID:-unknown}" \
        "${SLURM_JOB_NUM_NODES:-unknown}" \
        "${SLURM_NTASKS:-unknown}"

    CPP_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" --check
    CPP_RUNNER_ENV_FILE="${ENV_FILE}" \
        "${RUNNER}" --pwd "${CPP_SOURCE_CONTAINER}" -- \
        test -x "${CPP_EXECUTABLE}"
    launch_pipeline "$@"
}

main "$@"
