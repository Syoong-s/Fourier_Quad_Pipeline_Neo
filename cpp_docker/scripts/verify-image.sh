#!/usr/bin/env bash

set -euo pipefail

IMAGE_NAME="${1:-cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2}"
REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_IMAGE_ID=gxx12.3-openmpi4.1.8-pmi2
EXPECTED_GXX=12.3.0
EXPECTED_OPENMPI=4.1.8

# ==========================================
# Function: Stop image verification with one concise diagnostic
# Method: Print an error to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Verify image metadata, versions, libraries, and mount points
# Method: Run exact assertions in a disposable container
# ==========================================
verify_runtime() {
    docker run --rm "${IMAGE_NAME}" bash -lc "
        test \"\${CPP_IMAGE_ID}\" = \"${EXPECTED_IMAGE_ID}\"
        test \"\$(g++ -dumpfullversion)\" = \"${EXPECTED_GXX}\"
        mpicxx --showme:version | grep -F \"Open MPI ${EXPECTED_OPENMPI}\"
        mpirun --version | grep -F \"${EXPECTED_OPENMPI}\"
        ompi_info --all | grep -F 'MCA ess: pmi'
        ompi_info --all | grep -F 'MCA schizo: slurm'
        pkg-config --exact-version=4.6.4 cfitsio
        pkg-config --exact-version=3.3.11 fftw3
        pkg-config --exact-version=3.3.11 fftw3f
        test -r \"\${PMI2_PREFIX}/include/slurm/pmi2.h\"
        test -r \"\${PMI2_PREFIX}/lib/libpmi2.so\"
        test -r \"\${CPP_STACK_PREFIX}/include/eigen3/Eigen/Dense\"
        test -e \"\${CPP_STACK_PREFIX}/lib/liblapack.so\"
        test -e \"\${CPP_STACK_PREFIX}/lib/libblas.so\"
        test -d /workspace/src_pipe
        test -d /data/archive/science
        test -d /data/archive/dqmask
        test -d /data/catalogs/AstroDir
        test -d /data/catalogs/ExtSrcDir
        test -d /data/calib/FlatDir
        test -d /data/DataProcess
    "
}

# ==========================================
# Function: Compile and execute the bundled MPI and scientific-stack tests
# Method: Mount sources read-only and build only in a disposable container directory
# ==========================================
verify_compiler_linkage() {
    docker run --rm \
        --volume "${REPOSITORY_DIR}/runner/tests:/tests:ro" \
        "${IMAGE_NAME}" bash -lc '
        mkdir -p /tmp/cpppipeline-smoke
        mpicxx -O2 -std=c++17 /tests/mpi_identity.cpp \
            -o /tmp/cpppipeline-smoke/mpi_identity
        mpicxx -O2 -std=c++17 /tests/science_stack_smoke.cpp \
            -o /tmp/cpppipeline-smoke/science_stack_smoke \
            -lcfitsio -lfftw3 -llapack -lblas -lm
        mpirun --mca plm isolated --oversubscribe \
            -n 2 /tmp/cpppipeline-smoke/mpi_identity
        mpirun --mca plm isolated --oversubscribe \
            -n 2 /tmp/cpppipeline-smoke/science_stack_smoke
    '
}

# ==========================================
# Function: Verify that the image is a source-external development runtime
# Method: Reject pipeline source files and compiled pipeline artifacts in the image
# ==========================================
verify_source_absence() {
    docker run --rm "${IMAGE_NAME}" bash -lc '
        test ! -e /workspace/src_pipe/main.cpp
        test ! -e /workspace/src_pipe/Fourier_Quad_Pipe
        test -z "$(find /workspace/src_pipe -xdev -type f -print -quit)"
    '
}

# ==========================================
# Function: Run the complete image verification suite
# Method: Check Docker access, versions, linkage, MPI execution, and source exclusion
# ==========================================
main() {
    command -v docker >/dev/null 2>&1 || die "docker is unavailable"
    docker info >/dev/null
    verify_runtime
    verify_compiler_linkage
    verify_source_absence
    printf 'Portable image verification passed: %s\n' "${IMAGE_NAME}"
}

main "$@"
