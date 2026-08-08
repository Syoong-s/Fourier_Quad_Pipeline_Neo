# cpp_Standard portable container environment

This directory builds one x86_64 Linux runtime for `cpp_Standard`. The same
OCI image is converted once to an Apptainer SIF and used on every Slurm site
that provides Apptainer or Singularity and the Slurm PMI2 launch plugin.
Pipeline source, compressed Science/DQ archives, catalogues, calibration data,
processing data, and outputs remain on bind-mounted host storage.

## Runtime contract

| Component | Version or interface |
| --- | --- |
| G++ | 12.3.0 |
| OpenMPI | 4.1.8 |
| Slurm client interface | PMI2 from Slurm 25.11.2 |
| OpenMPI Slurm integration | direct-launch environment components enabled |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / BLAS | 3.11.0 |
| OpenBLAS | 0.3.33 |

There is one Docker target (`runtime`), one Pixi environment (`default`), one
image ID (`gxx12.3-openmpi4.1.8-pmi2`), and one SIF. The cluster never loads a
host compiler or host OpenMPI for the application. Slurm starts each rank with
`srun --mpi=pmi2`; the application links only the MPI and scientific libraries
inside the SIF. OpenMPI is configured with both PMI and Slurm direct-launch
support; this does not link the application to a host Slurm library.

This contract is portable across x86_64 Slurm clusters with PMI2. It is not a
claim of compatibility with ARM systems, non-Slurm schedulers, or Slurm sites
that expose only PMIx. The baseline image uses shared memory within a node and
routable TCP across nodes. Vendor-fabric acceleration requires separate
site-qualified testing.

## Build and verify

From this directory:

```text
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh \
  cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
bash scripts/check-public-repo.sh
```

The verification compiles and runs two MPI ranks inside Docker, exercises the
scientific stack, checks exact component versions and PMI support, and confirms
that pipeline source is absent from the image.

## Local compilation

Copy `.env.example` to `.env`, fill every host path, and run:

```text
docker compose run --rm FourierQuad-CPP
make -C /workspace/src_pipe clean
make -C /workspace/src_pipe -j4
```

The Makefile uses `mpicxx`, C++17, and the image search paths. Its optional
`STACK_PREFIX` remains available outside the container. Catalogue and
flat-field destinations must match the compile-time strings in
`cpp_Standard/include/process_main/LensingConfig.hpp`.

The fixed local bind interface is:

- source: `CPP_SOURCE_HOST` → `/workspace/src_pipe` (read/write);
- Science archive: `SCIENCE_ROOT_HOST` → `SCIENCE_ROOT_CONTAINER` (read-only);
- DQMask archive: `DQ_ROOT_HOST` → `DQ_ROOT_CONTAINER` (read-only);
- astrometry/source catalogues and flat calibration (read-only);
- processing data and generated initializer output (read/write).

Pass `SCIENCE_ROOT_CONTAINER` and `DQ_ROOT_CONTAINER` to `--science-root` and
`--dq-root`. Use `PROCESS_DATA_CONTAINER` for `--output-root`. Pipeline arguments
need no container-interface change: repeat `--dataset TARGET:PREFIX` for a batch
and repeat `--contains TOKEN` for OR matching.

## HPC deployment

The generic runner supports two acquisition paths:

- save the reviewed image with `docker save`, then submit
  `runner/build-sif.slurm`;
- set a digest-pinned `OCI_IMAGE_URI`, then run `runner/pull-sif.sh`.

Both paths create a SHA256 sidecar and refuse to overwrite an existing SIF.
Copy `runner/cpppipeline.env.example` to `runner/cpppipeline.env`, edit the
shared paths, then validate in this order:

1. `run-apptainer.sh --check`;
2. `compile-pipeline.slurm`;
3. `mpi-smoke-test.slurm` on one node and then multiple nodes;
4. `cpppipeline.slurm` only after science paths and the exposure list are
   reviewed.

See [the English runner guide](runner/README.md) and
[the Chinese runner guide](runner/README-CN.md). Dependency provenance and
license boundaries are recorded in [SOURCES.md](SOURCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
