#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${SCRIPT_DIR}/cpppipeline.env}"
TEMPORARY_SIF=
TEMPORARY_CHECKSUM=

# ==========================================
# Function: Stop SIF acquisition with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Locate the user-space HPC container runtime
# Method: Honor APPTAINER_BIN, then try Apptainer and Singularity
# ==========================================
select_runtime() {
    if [[ -n "${APPTAINER_BIN:-}" ]]; then
        command -v "${APPTAINER_BIN}" >/dev/null 2>&1 ||
            die "container runtime is unavailable: ${APPTAINER_BIN}"
        printf '%s\n' "${APPTAINER_BIN}"
    elif command -v apptainer >/dev/null 2>&1; then
        printf '%s\n' apptainer
    elif command -v singularity >/dev/null 2>&1; then
        printf '%s\n' singularity
    else
        die "neither apptainer nor singularity is available"
    fi
}

# ==========================================
# Function: Remove only unfinished artifacts created by this process
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
# Function: Pull one immutable SIF without overwriting an existing image
# Method: Pull to a temporary file, checksum it, and rename both outputs atomically
# ==========================================
main() {
    local output_directory
    local output_name
    local runtime_command

    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    [[ -n "${OCI_IMAGE_URI:-}" ]] || die "OCI_IMAGE_URI is empty in ${ENV_FILE}"
    [[ -n "${CPP_SIF:-}" && "${CPP_SIF}" == /* ]] ||
        die "CPP_SIF must be an absolute path"
    [[ ! -e "${CPP_SIF}" ]] ||
        die "refusing to overwrite existing SIF: ${CPP_SIF}"
    [[ ! -e "${CPP_SIF}.sha256" ]] ||
        die "refusing to overwrite existing checksum: ${CPP_SIF}.sha256"

    output_directory="$(dirname "${CPP_SIF}")"
    [[ -d "${output_directory}" && -w "${output_directory}" ]] ||
        die "SIF output directory is not writable: ${output_directory}"

    runtime_command="$(select_runtime)"
    output_name="$(basename "${CPP_SIF}")"
    TEMPORARY_SIF="${output_directory}/.${output_name}.partial.$$.sif"
    TEMPORARY_CHECKSUM="${output_directory}/.${output_name}.sha256.partial.$$"
    [[ ! -e "${TEMPORARY_SIF}" && ! -e "${TEMPORARY_CHECKSUM}" ]] ||
        die "temporary SIF output already exists"
    trap cleanup_temporary_artifacts EXIT

    umask 002
    "${runtime_command}" pull "${TEMPORARY_SIF}" "docker://${OCI_IMAGE_URI}"
    [[ -s "${TEMPORARY_SIF}" ]] ||
        die "runtime returned without creating a nonempty SIF: ${CPP_SIF}"
    mv -- "${TEMPORARY_SIF}" "${CPP_SIF}"
    TEMPORARY_SIF=
    sha256sum "${CPP_SIF}" >"${TEMPORARY_CHECKSUM}"
    mv -- "${TEMPORARY_CHECKSUM}" "${CPP_SIF}.sha256"
    TEMPORARY_CHECKSUM=
    printf 'Created SIF: %s\n' "${CPP_SIF}"
    cat "${CPP_SIF}.sha256"
}

main "$@"
