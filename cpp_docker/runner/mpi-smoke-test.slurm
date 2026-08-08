#!/usr/bin/env bash

#SBATCH --job-name=cpp-mpi-smoke
#SBATCH --partition=cpu
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00
#SBATCH --output=cpp-mpi-smoke-%j.out
#SBATCH --error=cpp-mpi-smoke-%j.err

set -euo pipefail

RUNNER_DIR="${CPP_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${RUNNER_DIR}/cpppipeline.env}"
RUNNER="${RUNNER_DIR}/run-apptainer.sh"
TEST_SOURCE_DIR="${RUNNER_DIR}/tests"
SMOKE_DIR=

# ==========================================
# Function: Stop a smoke test with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Load the host modules requested by the HPC configuration
# Method: Validate the Bash module array and load it before MPI or Apptainer use
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
    module load "${HPC_MODULES[@]}" ||
        die "failed to load one or more modules from HPC_MODULES"
}

# ==========================================
# Function: Remove only the unique shared smoke directory created by this job
# Method: Validate the mktemp result and recursively remove that exact path
# ==========================================
cleanup_smoke_directory() {
    if [[ -n "${SMOKE_DIR}" && -d "${SMOKE_DIR}" ]]; then
        rm -rf -- "${SMOKE_DIR}"
    fi
}

# ==========================================
# Function: Load paths, optional site setup, modules, and srun arguments
# Method: Source trusted configuration and validate every array-valued setting
# ==========================================
load_job_environment() {
    local srun_declaration

    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    [[ -x "${RUNNER}" ]] || die "Apptainer runner is not executable: ${RUNNER}"
    [[ -r "${TEST_SOURCE_DIR}/mpi_identity.cpp" ]] ||
        die "MPI identity source is missing"
    [[ -r "${TEST_SOURCE_DIR}/science_stack_smoke.cpp" ]] ||
        die "science-stack source is missing"

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
}

# ==========================================
# Function: Compile MPI and scientific-stack tests into shared temporary storage
# Method: Bind read-only test sources and one unique writable job directory
# ==========================================
compile_smoke_programs() {
    [[ -n "${HPC_SHARED_SCRATCH_HOST:-}" ]] ||
        die "HPC_SHARED_SCRATCH_HOST is empty"
    [[ -d "${HPC_SHARED_SCRATCH_HOST}" && -w "${HPC_SHARED_SCRATCH_HOST}" ]] ||
        die "shared smoke parent is not writable: ${HPC_SHARED_SCRATCH_HOST}"

    SMOKE_DIR="$(mktemp -d \
        "${HPC_SHARED_SCRATCH_HOST%/}/.cpppipeline-smoke.${SLURM_JOB_ID}.XXXXXX")"
    trap cleanup_smoke_directory EXIT

    CPP_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" \
        --bind "${TEST_SOURCE_DIR}:/tests:ro" \
        --bind "${SMOKE_DIR}:/smoke:rw" \
        --pwd /smoke -- \
        mpicxx -O2 -std=c++17 /tests/mpi_identity.cpp \
        -o /smoke/mpi_identity

    CPP_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" \
        --bind "${TEST_SOURCE_DIR}:/tests:ro" \
        --bind "${SMOKE_DIR}:/smoke:rw" \
        --pwd /smoke -- \
        mpicxx -O2 -std=c++17 /tests/science_stack_smoke.cpp \
        -o /smoke/science_stack_smoke \
        -lcfitsio -lfftw3 -llapack -lblas -lm
}

# ==========================================
# Function: Launch one shared smoke executable across the Slurm allocation
# Method: Let srun establish PMI2 before every rank enters the same SIF
# ==========================================
launch_smoke_program() {
    local executable_name="$1"
    local container_command=(
        "${RUNNER}"
        --env-file "${ENV_FILE}"
        --bind "${SMOKE_DIR}:/smoke:rw"
        --pwd /smoke
        --
        "/smoke/${executable_name}"
    )

    command -v srun >/dev/null 2>&1 || die "srun is unavailable"
    srun "${SRUN_ARGS[@]}" \
        --mpi="${SLURM_MPI_TYPE}" \
        --ntasks="${SLURM_NTASKS}" \
        "${container_command[@]}"
}

# ==========================================
# Function: Build and run multi-node C++ MPI validation programs
# Method: Check the SIF, compile once, report placement, and verify scientific libraries
# ==========================================
main() {
    [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_NTASKS:-}" ]] ||
        die "submit this file with sbatch"
    load_job_environment
    command -v srun >/dev/null 2>&1 || die "srun is unavailable"
    srun --mpi=list 2>&1 | grep -qw "${SLURM_MPI_TYPE}" ||
        die "Slurm does not advertise --mpi=${SLURM_MPI_TYPE}"
    CPP_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" --check
    compile_smoke_programs
    launch_smoke_program mpi_identity
    launch_smoke_program science_stack_smoke
    printf 'Multi-node C++ MPI smoke test completed successfully.\n'
}

main "$@"
