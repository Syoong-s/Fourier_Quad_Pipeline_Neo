#!/usr/bin/env bash

set -euo pipefail

# ==========================================
# Function: Print one command if available without changing cluster state
# Method: Resolve the executable and execute its version or configuration query
# ==========================================
run_if_available() {
    local command_name="$1"
    shift

    if command -v "${command_name}" >/dev/null 2>&1; then
        printf '\n[%s]\n' "${command_name}"
        "${command_name}" "$@" || true
    fi
}

# ==========================================
# Function: Audit compiler, MPI, Slurm launch plugins, and Apptainer
# Method: Run read-only version and configuration commands for compatibility review
# ==========================================
main() {
    run_if_available gcc --version
    run_if_available g++ --version
    run_if_available mpicxx --showme:version
    run_if_available mpicxx --showme:command
    run_if_available mpirun --version
    run_if_available ompi_info --config
    run_if_available srun --mpi=list
    run_if_available apptainer --version
    run_if_available singularity --version
}

main "$@"
