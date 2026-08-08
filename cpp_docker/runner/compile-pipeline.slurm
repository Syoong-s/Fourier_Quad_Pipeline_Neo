#!/usr/bin/env bash

#SBATCH --job-name=cpp-compile
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=cpp-compile-%j.out
#SBATCH --error=cpp-compile-%j.err

set -euo pipefail

RUNNER_DIR="${CPP_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${RUNNER_DIR}/cpppipeline.env}"
RUNNER="${RUNNER_DIR}/run-apptainer.sh"

# ==========================================
# Function: Stop compilation with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Load configured host modules
# Method: Validate the Bash module array before loading runtime prerequisites
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
# Function: Load the trusted runner configuration
# Method: Source optional site setup, load runtime modules, and validate build inputs
# ==========================================
load_job_environment() {
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

    [[ "${CPP_SOURCE_CONTAINER:-}" == /* ]] ||
        die "CPP_SOURCE_CONTAINER must be an absolute path"
    [[ "${CPP_EXECUTABLE:-}" == /* ]] ||
        die "CPP_EXECUTABLE must be an absolute path"
    [[ "${CPP_BUILD_JOBS:-}" =~ ^[1-9][0-9]*$ ]] ||
        die "CPP_BUILD_JOBS must be a positive integer"
    [[ "${CPP_MAKE_CLEAN:-}" =~ ^[01]$ ]] ||
        die "CPP_MAKE_CLEAN must be 0 or 1"
    [[ "${SLURM_CPUS_PER_TASK:-1}" =~ ^[1-9][0-9]*$ ]] ||
        die "SLURM_CPUS_PER_TASK must be a positive integer"
    ((CPP_BUILD_JOBS <= SLURM_CPUS_PER_TASK)) ||
        die "CPP_BUILD_JOBS exceeds SLURM_CPUS_PER_TASK"
}

# ==========================================
# Function: Compile the mounted cpp_Standard source once
# Method: Optionally clean and invoke its portable Makefile inside the SIF
# ==========================================
compile_pipeline() {
    if [[ "${CPP_MAKE_CLEAN}" == "1" ]]; then
        CPP_RUNNER_ENV_FILE="${ENV_FILE}" \
            "${RUNNER}" --pwd "${CPP_SOURCE_CONTAINER}" -- \
            make -C "${CPP_SOURCE_CONTAINER}" clean
    fi

    CPP_RUNNER_ENV_FILE="${ENV_FILE}" \
        "${RUNNER}" --pwd "${CPP_SOURCE_CONTAINER}" -- \
        make -C "${CPP_SOURCE_CONTAINER}" -j"${CPP_BUILD_JOBS}"
}

# ==========================================
# Function: Validate the SIF and compile cpp_Standard on one compute node
# Method: Check the immutable stack, build once, and require the final executable
# ==========================================
main() {
    [[ -n "${SLURM_JOB_ID:-}" ]] || die "submit this file with sbatch"
    load_job_environment
    CPP_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" --check
    compile_pipeline
    CPP_RUNNER_ENV_FILE="${ENV_FILE}" \
        "${RUNNER}" --pwd "${CPP_SOURCE_CONTAINER}" -- \
        test -x "${CPP_EXECUTABLE}"
    printf 'Pipeline compilation passed: %s\n' "${CPP_EXECUTABLE}"
}

main "$@"
