#!/usr/bin/env bash

#SBATCH --job-name=cpp-sif-build
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=cpp-sif-build-%j.out
#SBATCH --error=cpp-sif-build-%j.err

set -euo pipefail

RUNNER_DIR="${CPP_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${RUNNER_DIR}/cpppipeline.env}"
TEMPORARY_SIF=
TEMPORARY_CHECKSUM=

# ==========================================
# Function: Stop SIF construction with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Remove only unfinished artifacts created by this job
# Method: Validate nonempty exact paths and unlink files left before atomic rename
# ==========================================
cleanup_temporary_artifacts() {
    if [[ -n "${TEMPORARY_SIF}" && -f "${TEMPORARY_SIF}" ]]; then
        rm -f -- "${TEMPORARY_SIF}"
    fi
    if [[ -n "${TEMPORARY_CHECKSUM}" && -f "${TEMPORARY_CHECKSUM}" ]]; then
        rm -f -- "${TEMPORARY_CHECKSUM}"
    fi
}

# ==========================================
# Function: Build a SIF from one reviewed Docker archive
# Method: Isolate caches, build a temporary file, checksum it, and rename atomically
# ==========================================
main() {
    local output_directory
    local output_name

    [[ -n "${SLURM_JOB_ID:-}" ]] || die "submit this file with sbatch"
    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    [[ -n "${CPP_DOCKER_ARCHIVE:-}" && "${CPP_DOCKER_ARCHIVE}" == /* ]] ||
        die "CPP_DOCKER_ARCHIVE must be an absolute path"
    [[ -r "${CPP_DOCKER_ARCHIVE}" ]] ||
        die "Docker archive is not readable: ${CPP_DOCKER_ARCHIVE}"
    [[ -n "${CPP_SIF:-}" && "${CPP_SIF}" == /* ]] ||
        die "CPP_SIF must be an absolute path"
    [[ ! -e "${CPP_SIF}" ]] ||
        die "refusing to overwrite existing SIF: ${CPP_SIF}"
    [[ ! -e "${CPP_SIF}.sha256" ]] ||
        die "refusing to overwrite existing checksum: ${CPP_SIF}.sha256"
    [[ -d "${APPTAINER_CACHE_DIR:-}" && -w "${APPTAINER_CACHE_DIR}" ]] ||
        die "APPTAINER_CACHE_DIR must be an existing writable directory"
    [[ -d "${APPTAINER_TMP_DIR:-}" && -w "${APPTAINER_TMP_DIR}" ]] ||
        die "APPTAINER_TMP_DIR must be an existing writable directory"

    export APPTAINER_CACHEDIR="${APPTAINER_CACHE_DIR}"
    export APPTAINER_TMPDIR="${APPTAINER_TMP_DIR}"

    command -v apptainer >/dev/null 2>&1 || die "apptainer is unavailable"
    output_directory="$(dirname "${CPP_SIF}")"
    output_name="$(basename "${CPP_SIF}")"
    [[ -d "${output_directory}" && -w "${output_directory}" ]] ||
        die "SIF output directory is not writable: ${output_directory}"
    TEMPORARY_SIF="${output_directory}/.${output_name}.partial.${SLURM_JOB_ID}.sif"
    TEMPORARY_CHECKSUM="${output_directory}/.${output_name}.sha256.partial.${SLURM_JOB_ID}"
    [[ ! -e "${TEMPORARY_SIF}" && ! -e "${TEMPORARY_CHECKSUM}" ]] ||
        die "temporary output already exists for job ${SLURM_JOB_ID}"
    trap cleanup_temporary_artifacts EXIT

    umask 002
    apptainer build "${TEMPORARY_SIF}" \
        "docker-archive:${CPP_DOCKER_ARCHIVE}"
    [[ -s "${TEMPORARY_SIF}" ]] ||
        die "SIF build did not create a nonempty image"
    mv -- "${TEMPORARY_SIF}" "${CPP_SIF}"
    TEMPORARY_SIF=
    sha256sum "${CPP_SIF}" >"${TEMPORARY_CHECKSUM}"
    mv -- "${TEMPORARY_CHECKSUM}" "${CPP_SIF}.sha256"
    TEMPORARY_CHECKSUM=
    printf 'Created SIF: %s\n' "${CPP_SIF}"
    cat "${CPP_SIF}.sha256"
}

main "$@"
