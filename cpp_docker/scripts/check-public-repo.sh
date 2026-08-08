#!/usr/bin/env bash

set -euo pipefail

REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_MAKEFILE="${REPOSITORY_DIR}/../cpp_Standard/Makefile"

# ==========================================
# Function: Reject archives, oversized payloads, and private path markers
# Method: Inspect the public project without following external directories
# ==========================================
check_public_payload() {
    local private_home_marker="/home/""alatrion"
    local private_cluster_marker="/lustre/home/""acct-"

    if find "${REPOSITORY_DIR}" -type f \
        \( -name '*.tar' -o -name '*.tar.gz' -o -name '*.tar.bz2' \
           -o -name '*.tar.xz' -o -name '*.sif' \) \
        -print -quit | grep -q .; then
        printf 'Archive or SIF found in public project.\n' >&2
        return 1
    fi

    if find "${REPOSITORY_DIR}" -type f -size +95M -print -quit | grep -q .; then
        printf 'A public project file exceeds 95 MiB.\n' >&2
        return 1
    fi

    if grep -R -n -F \
        -e "${private_home_marker}" \
        -e "${private_cluster_marker}" \
        --exclude=.env \
        --exclude=cpppipeline.env \
        --exclude=site-env.sh \
        --exclude-dir=.git "${REPOSITORY_DIR}" "${PIPELINE_MAKEFILE}"; then
        printf 'Private path marker found in public assets.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Validate the portable Makefile and locked single-image inputs
# Method: Require overridable stack paths, one runtime, PMI2, and a current lock
# ==========================================
check_build_contract() {
    [[ -s "${REPOSITORY_DIR}/pixi.lock" ]] ||
        { printf 'pixi.lock is missing or empty.\n' >&2; return 1; }
    grep -F 'STACK_PREFIX ?=' "${PIPELINE_MAKEFILE}" >/dev/null
    grep -F 'FROM ${RUNTIME_IMAGE} AS runtime' \
        "${REPOSITORY_DIR}/Dockerfile" >/dev/null
    grep -F 'cpppipeline.gxx.version="12.3.0"' \
        "${REPOSITORY_DIR}/Dockerfile" >/dev/null
    grep -F 'cpppipeline.openmpi.version="4.1.8"' \
        "${REPOSITORY_DIR}/Dockerfile" >/dev/null
    grep -F 'cpppipeline.pmi2.client="slurm-25.11.2"' \
        "${REPOSITORY_DIR}/Dockerfile" >/dev/null
    grep -F 'cpppipeline.slurm.direct-launch="true"' \
        "${REPOSITORY_DIR}/Dockerfile" >/dev/null
    grep -F '  default:' "${REPOSITORY_DIR}/pixi.lock" >/dev/null
    if grep -F \
        -e 'name: openmpi' \
        -e 'name: libpmix' \
        "${REPOSITORY_DIR}/pixi.lock" >/dev/null; then
        printf 'Pixi lock unexpectedly contains MPI or PMIx.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Validate the example env as the runner's single configuration source
# Method: Source it and assert version, path, PMI2, and array contracts
# ==========================================
check_environment_contract() {
    bash -eu -o pipefail -c '
        source "$1"
        [[ "$(declare -p HPC_MODULES)" == "declare -a "* ]]
        [[ "$(declare -p HPC_EXTRA_BINDS)" == "declare -a "* ]]
        [[ "$(declare -p HPC_PASSTHROUGH_ENV)" == "declare -a "* ]]
        [[ "$(declare -p HPC_CONTAINER_ENV)" == "declare -a "* ]]
        [[ "$(declare -p SRUN_ARGS)" == "declare -a "* ]]
        [[ "${CPP_SOURCE_CONTAINER}" == "/workspace/src_pipe" ]]
        [[ "${CPP_EXPO_LIST_CONTAINER}" == "${PROCESS_DATA_CONTAINER%/}/expo_list.list" ]]
        [[ "${CPP_EXECUTABLE}" == "${CPP_SOURCE_CONTAINER%/}/Fourier_Quad_Pipe" ]]
        [[ "${CPP_BUILD_JOBS}" =~ ^[1-9][0-9]*$ ]]
        [[ "${CPP_MAKE_CLEAN}" =~ ^[01]$ ]]
        [[ "${CPP_IMAGE_ID_EXPECTED}" == "gxx12.3-openmpi4.1.8-pmi2" ]]
        [[ "${CPP_GXX_VERSION_EXPECTED}" == "12.3.0" ]]
        [[ "${CPP_OPENMPI_VERSION_EXPECTED}" == "4.1.8" ]]
        [[ "${MPI_LAUNCH_MODE}" == "srun" ]]
        [[ "${SLURM_MPI_TYPE}" == "pmi2" ]]
        [[ "${HPC_SCRUB_MPI_ENV}" == "1" ]]
    ' _ "${REPOSITORY_DIR}/runner/cpppipeline.env.example"
}

# ==========================================
# Function: Validate Docker Compose archive bind variables
# Method: Require fixed Science/DQ host inputs, container destinations, and
#         read-only Compose mounts in the public local configuration template.
# ==========================================
check_compose_environment_contract() {
    bash -eu -o pipefail -c '
        source "$1"
    ' _ "${REPOSITORY_DIR}/.env.example"
}

# ==========================================
# Function: Validate shell syntax, LF endings, and executable metadata
# Method: Check every operational runner and script deterministically
# ==========================================
check_shell_assets() {
    local executable_file
    local executable_files=(
        "${REPOSITORY_DIR}/scripts/check-public-repo.sh"
        "${REPOSITORY_DIR}/scripts/verify-image.sh"
        "${REPOSITORY_DIR}/runner/build-sif.slurm"
        "${REPOSITORY_DIR}/runner/compile-pipeline.slurm"
        "${REPOSITORY_DIR}/runner/cpppipeline.slurm"
        "${REPOSITORY_DIR}/runner/inspect-cluster-mpi.sh"
        "${REPOSITORY_DIR}/runner/mpi-smoke-test.slurm"
        "${REPOSITORY_DIR}/runner/pull-sif.sh"
        "${REPOSITORY_DIR}/runner/run-apptainer.sh"
    )

    if grep -R -I -l $'\r$' --exclude-dir=.git "${REPOSITORY_DIR}" | grep -q .; then
        printf 'CRLF line ending found in public project.\n' >&2
        return 1
    fi

    for executable_file in "${executable_files[@]}"; do
        bash -n "${executable_file}"
        [[ -x "${executable_file}" ]] ||
            { printf 'Runner is not executable: %s\n' "${executable_file}" >&2; return 1; }
    done

    if find "${REPOSITORY_DIR}/runner" -maxdepth 1 -type f \
        -name '*pilogin*' -print -quit | grep -q .; then
        printf 'Site-specific runner file found.\n' >&2
        return 1
    fi
    if grep -F -e 'MPI_LAUNCHER' -e 'mpirun' \
        "${REPOSITORY_DIR}/runner/cpppipeline.slurm" \
        "${REPOSITORY_DIR}/runner/mpi-smoke-test.slurm" \
        "${REPOSITORY_DIR}/runner/run-apptainer.sh" >/dev/null; then
        printf 'Host MPI launcher reference found in an execution path.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Run all deterministic publication checks
# Method: Execute payload, build, env, and shell checks in order
# ==========================================
main() {
    check_public_payload
    check_build_contract
    check_environment_contract
    check_compose_environment_contract
    check_shell_assets
    printf 'Public repository checks passed.\n'
}

main "$@"
