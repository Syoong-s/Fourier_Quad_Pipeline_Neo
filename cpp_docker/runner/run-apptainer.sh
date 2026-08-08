#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${CPP_RUNNER_ENV_FILE:-${SCRIPT_DIR}/cpppipeline.env}"
CONTAINER_PWD=
CHECK_ONLY=0
EXTRA_BINDS=()
RUN_COMMAND=()
RUNTIME_ARGUMENTS=()

# ==========================================
# Function: Display the Apptainer runner command-line contract
# Method: Print supported configuration, bind, working-directory, and check options
# ==========================================
usage() {
    printf '%s\n' \
        "Usage: $0 [--env-file FILE] [--pwd CONTAINER_DIR]" \
        "          [--bind SRC:DEST[:ro|rw]] [--check] [--] [COMMAND [ARG...]]"
}

# ==========================================
# Function: Stop runner execution with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Require one nonempty shell variable
# Method: Resolve the variable indirectly and reject missing configuration
# ==========================================
require_variable() {
    local variable_name="$1"
    [[ -n "${!variable_name:-}" ]] ||
        die "required variable ${variable_name} is empty in ${ENV_FILE}"
}

# ==========================================
# Function: Validate one host bind source directory
# Method: Require an absolute path to an existing directory visible on this node
# ==========================================
require_host_directory() {
    local variable_name="$1"
    local directory_path="${!variable_name}"

    [[ "${directory_path}" == /* ]] ||
        die "${variable_name} must be an absolute host path"
    [[ -d "${directory_path}" ]] ||
        die "${variable_name} does not exist on $(hostname): ${directory_path}"
}

# ==========================================
# Function: Validate one container destination path
# Method: Require an absolute path and reject bind-delimiter characters
# ==========================================
require_container_path() {
    local variable_name="$1"
    local container_path="${!variable_name}"

    [[ "${container_path}" == /* ]] ||
        die "${variable_name} must be an absolute container path"
    [[ "${container_path}" != *:* && "${container_path}" != *,* ]] ||
        die "${variable_name} contains an unsupported bind delimiter"
}

# ==========================================
# Function: Require one configuration value to be a Bash indexed array
# Method: Inspect the sourced declaration and reject scalars or associative arrays
# ==========================================
require_indexed_array() {
    local variable_name="$1"
    local variable_declaration

    variable_declaration="$(declare -p "${variable_name}" 2>/dev/null)" ||
        die "${variable_name} must be a Bash indexed array in ${ENV_FILE}"
    [[ "${variable_declaration}" == "declare -a "* ]] ||
        die "${variable_name} must be a Bash indexed array in ${ENV_FILE}"
}

# ==========================================
# Function: Locate the user-space HPC container runtime
# Method: Honor APPTAINER_BIN, then try Apptainer and Singularity in order
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
# Function: Parse runner options without consuming the application command
# Method: Shift known flags and retain all trailing arguments for container exec
# ==========================================
parse_arguments() {
    while (($#)); do
        case "$1" in
            --env-file)
                (($# >= 2)) || die "--env-file requires a path"
                ENV_FILE="$2"
                shift 2
                ;;
            --pwd)
                (($# >= 2)) || die "--pwd requires a container path"
                CONTAINER_PWD="$2"
                shift 2
                ;;
            --bind)
                (($# >= 2)) || die "--bind requires SRC:DEST[:ro|rw]"
                EXTRA_BINDS+=("$2")
                shift 2
                ;;
            --check)
                CHECK_ONLY=1
                shift
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                die "unknown option: $1"
                ;;
            *)
                break
                ;;
        esac
    done

    RUN_COMMAND=("$@")
}

# ==========================================
# Function: Load and validate all shared-filesystem and container paths
# Method: Source the trusted env file and check every required runtime input
# ==========================================
load_configuration() {
    local variable_name
    local required_variables=(
        CPP_SIF
        CPP_SOURCE_HOST
        CPP_SOURCE_CONTAINER
        ASTROMETRY_CAT_HOST
        ASTROMETRY_CAT_CONTAINER
        SOURCE_CAT_HOST
        SOURCE_CAT_CONTAINER
        FLAT_PATH_HOST
        FLAT_PATH_CONTAINER
        PROCESS_DATA_HOST
        PROCESS_DATA_CONTAINER
        CPP_IMAGE_ID_EXPECTED
        CPP_GXX_VERSION_EXPECTED
        CPP_OPENMPI_VERSION_EXPECTED
        MPI_LAUNCH_MODE
        SLURM_MPI_TYPE
        HPC_SCRUB_MPI_ENV
    )

    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    for variable_name in "${required_variables[@]}"; do
        require_variable "${variable_name}"
    done

    [[ "${CPP_SIF}" == /* ]] || die "CPP_SIF must be an absolute path"
    [[ -r "${CPP_SIF}" ]] || die "SIF image is not readable: ${CPP_SIF}"
    [[ -n "${CONTAINER_PWD}" ]] || CONTAINER_PWD="${CPP_SOURCE_CONTAINER}"
    [[ "${MPI_LAUNCH_MODE}" == "srun" ]] ||
        die "MPI_LAUNCH_MODE must be srun"
    [[ "${SLURM_MPI_TYPE}" == "pmi2" ]] ||
        die "SLURM_MPI_TYPE must be pmi2"
    [[ "${HPC_SCRUB_MPI_ENV}" == "1" ]] ||
        die "HPC_SCRUB_MPI_ENV must be 1"
    require_indexed_array HPC_EXTRA_BINDS
    require_indexed_array HPC_PASSTHROUGH_ENV
    require_indexed_array HPC_CONTAINER_ENV

    for variable_name in \
        CPP_SOURCE_HOST \
        ASTROMETRY_CAT_HOST \
        SOURCE_CAT_HOST \
        FLAT_PATH_HOST \
        PROCESS_DATA_HOST; do
        require_host_directory "${variable_name}"
    done

    for variable_name in \
        CPP_SOURCE_CONTAINER \
        ASTROMETRY_CAT_CONTAINER \
        SOURCE_CAT_CONTAINER \
        FLAT_PATH_CONTAINER \
        PROCESS_DATA_CONTAINER; do
        require_container_path "${variable_name}"
    done

    # Optional pipeline output mounts (empty = use default under PROCESS_DATA)
    for variable_name in \
        SCIENCE_ROOT_HOST SCIENCE_ROOT_CONTAINER \
        DQ_ROOT_HOST DQ_ROOT_CONTAINER \
        EXTCAT_INPUT_HOST EXTCAT_INPUT_CONTAINER \
        REARR_OUTPUT_HOST REARR_OUTPUT_CONTAINER \
        EXPOLIST_DIR_HOST EXPOLIST_DIR_CONTAINER \
        FD_OUTPUT_HOST FD_OUTPUT_CONTAINER; do
        if [[ -n "${!variable_name:-}" ]]; then
            case "${variable_name}" in
                *_HOST)      require_host_directory "${variable_name}" ;;
                *_CONTAINER) require_container_path "${variable_name}" ;;
            esac
        fi
    done
}

# ==========================================
# Function: Verify the optional immutable SIF checksum
# Method: Validate the configured SHA256 syntax and compare it on the host once
# ==========================================
verify_sif_checksum() {
    local actual_checksum

    if [[ -z "${CPP_SIF_SHA256_EXPECTED:-}" ]]; then
        return
    fi
    [[ "${CPP_SIF_SHA256_EXPECTED}" =~ ^[[:xdigit:]]{64}$ ]] ||
        die "CPP_SIF_SHA256_EXPECTED must contain exactly 64 hexadecimal digits"
    command -v sha256sum >/dev/null 2>&1 ||
        die "sha256sum is required for SIF verification"
    actual_checksum="$(sha256sum "${CPP_SIF}" | awk '{print $1}')"
    [[ "${actual_checksum,,}" == "${CPP_SIF_SHA256_EXPECTED,,}" ]] ||
        die "SIF checksum mismatch: ${CPP_SIF}"
}

# ==========================================
# Function: Build the clean container environment for Slurm PMI2 execution
# Method: Forward scheduler/PMI state and only explicitly allowlisted extra values
# ==========================================
append_runtime_environment() {
    local environment_assignment
    local environment_name

    while IFS= read -r environment_name; do
        case "${environment_name}" in
            SLURM_*|PMI_*|PMI2_*)
                RUNTIME_ARGUMENTS+=(
                    --env "${environment_name}=${!environment_name}"
                )
                ;;
        esac
    done < <(compgen -e | LC_ALL=C sort)

    for environment_name in "${HPC_PASSTHROUGH_ENV[@]}"; do
        [[ "${environment_name}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] ||
            die "invalid name in HPC_PASSTHROUGH_ENV: ${environment_name}"
        if [[ -v "${environment_name}" ]]; then
            RUNTIME_ARGUMENTS+=(
                --env "${environment_name}=${!environment_name}"
            )
        fi
    done

    RUNTIME_ARGUMENTS+=(--env "OMPI_MCA_ess=pmi")
    for environment_assignment in "${HPC_CONTAINER_ENV[@]}"; do
        [[ "${environment_assignment}" == *=* ]] ||
            die "HPC_CONTAINER_ENV entries must use NAME=value"
        environment_name="${environment_assignment%%=*}"
        [[ "${environment_name}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] ||
            die "invalid name in HPC_CONTAINER_ENV: ${environment_name}"
        RUNTIME_ARGUMENTS+=(--env "${environment_assignment}")
    done
}

# ==========================================
# Function: Execute a command in the immutable SIF with Compose-equivalent binds
# Method: Build clean PMI2-aware arguments and replace the wrapper with the runtime
# ==========================================
run_container() {
    local runtime_command="$1"
    local bind_specification
    RUNTIME_ARGUMENTS=(
        exec
        --cleanenv
        --no-home
        --pwd "${CONTAINER_PWD}"
        --bind "${CPP_SOURCE_HOST}:${CPP_SOURCE_CONTAINER}:rw"
        --bind "${ASTROMETRY_CAT_HOST}:${ASTROMETRY_CAT_CONTAINER}:ro"
        --bind "${SOURCE_CAT_HOST}:${SOURCE_CAT_CONTAINER}:ro"
        --bind "${FLAT_PATH_HOST}:${FLAT_PATH_CONTAINER}:ro"
        --bind "${PROCESS_DATA_HOST}:${PROCESS_DATA_CONTAINER}:rw"
        --env "CPP_SOURCE_CONTAINER=${CPP_SOURCE_CONTAINER}"
        --env "SCIENCE_ROOT_CONTAINER=${SCIENCE_ROOT_CONTAINER:-}"
        --env "DQ_ROOT_CONTAINER=${DQ_ROOT_CONTAINER:-}"
        --env "PROCESS_DATA_CONTAINER=${PROCESS_DATA_CONTAINER}"
        --env "EXTCAT_INPUT_CONTAINER=${EXTCAT_INPUT_CONTAINER:-}"
        --env "REARR_OUTPUT_CONTAINER=${REARR_OUTPUT_CONTAINER:-}"
        --env "EXPOLIST_DIR_CONTAINER=${EXPOLIST_DIR_CONTAINER:-}"
        --env "FD_OUTPUT_CONTAINER=${FD_OUTPUT_CONTAINER:-}"
    )

    # Optional pipeline mounts (bound only when HOST is set)
    if [[ -n "${SCIENCE_ROOT_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${SCIENCE_ROOT_HOST}:${SCIENCE_ROOT_CONTAINER}:ro")
    fi
    if [[ -n "${DQ_ROOT_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${DQ_ROOT_HOST}:${DQ_ROOT_CONTAINER}:ro")
    fi
    if [[ -n "${EXTCAT_INPUT_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${EXTCAT_INPUT_HOST}:${EXTCAT_INPUT_CONTAINER}:ro")
    fi
    if [[ -n "${REARR_OUTPUT_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${REARR_OUTPUT_HOST}:${REARR_OUTPUT_CONTAINER}:rw")
    fi
    if [[ -n "${EXPOLIST_DIR_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${EXPOLIST_DIR_HOST}:${EXPOLIST_DIR_CONTAINER}:rw")
    fi
    if [[ -n "${FD_OUTPUT_HOST:-}" ]]; then
        RUNTIME_ARGUMENTS+=(--bind "${FD_OUTPUT_HOST}:${FD_OUTPUT_CONTAINER}:rw")
    fi

    for bind_specification in "${EXTRA_BINDS[@]}"; do
        RUNTIME_ARGUMENTS+=(--bind "${bind_specification}")
    done
    for bind_specification in "${HPC_EXTRA_BINDS[@]}"; do
        RUNTIME_ARGUMENTS+=(--bind "${bind_specification}")
    done
    append_runtime_environment

    if ((CHECK_ONLY)); then
        RUN_COMMAND=(
            bash -c
            'set -euo pipefail
             test "${CPP_IMAGE_ID}" = "${CPP_IMAGE_ID_EXPECTED}"
             test "$(g++ -dumpfullversion)" = "${CPP_GXX_VERSION_EXPECTED}"
             mpicxx --showme:version |
                 grep -F "Open MPI ${CPP_OPENMPI_VERSION_EXPECTED}"
             ompi_info --all | grep -F "MCA ess: pmi"
             test -r "${PMI2_PREFIX}/include/slurm/pmi2.h"
             test -r "${PMI2_PREFIX}/lib/libpmi2.so"
             pkg-config --exact-version=4.6.4 cfitsio
             pkg-config --exact-version=3.3.11 fftw3
             test -r "${CPP_STACK_PREFIX}/include/eigen3/Eigen/Dense"
             test -d "${CPP_SOURCE_CONTAINER}"
             [[ -z "${SCIENCE_ROOT_CONTAINER:-}" || -d "${SCIENCE_ROOT_CONTAINER}" ]]
             [[ -z "${DQ_ROOT_CONTAINER:-}" || -d "${DQ_ROOT_CONTAINER}" ]]
             test -d "${PROCESS_DATA_CONTAINER}"
             [[ -z "${EXTCAT_INPUT_CONTAINER:-}" || -d "${EXTCAT_INPUT_CONTAINER}" ]]
             [[ -z "${REARR_OUTPUT_CONTAINER:-}" || -d "${REARR_OUTPUT_CONTAINER}" ]]
             [[ -z "${EXPOLIST_DIR_CONTAINER:-}" || -d "${EXPOLIST_DIR_CONTAINER}" ]]
             [[ -z "${FD_OUTPUT_CONTAINER:-}" || -d "${FD_OUTPUT_CONTAINER}" ]]
             printf "Apptainer image and bind checks passed on %s.\n" "$(hostname)"'
        )
    elif ((${#RUN_COMMAND[@]} == 0)); then
        RUN_COMMAND=(bash)
    fi

    exec "${runtime_command}" "${RUNTIME_ARGUMENTS[@]}" \
        --env "CPP_IMAGE_ID_EXPECTED=${CPP_IMAGE_ID_EXPECTED}" \
        --env "CPP_GXX_VERSION_EXPECTED=${CPP_GXX_VERSION_EXPECTED}" \
        --env "CPP_OPENMPI_VERSION_EXPECTED=${CPP_OPENMPI_VERSION_EXPECTED}" \
        "${CPP_SIF}" "${RUN_COMMAND[@]}"
}

# ==========================================
# Function: Coordinate argument parsing, configuration checks, and container exec
# Method: Resolve all inputs before replacing the wrapper process
# ==========================================
main() {
    local runtime_command

    parse_arguments "$@"
    load_configuration
    require_container_path CONTAINER_PWD
    if ((CHECK_ONLY)); then
        verify_sif_checksum
    fi
    runtime_command="$(select_runtime)"
    run_container "${runtime_command}"
}

main "$@"
