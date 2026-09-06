# Fourier_Quad C++ container

This directory builds an x86_64 Linux toolchain image for `cpp_Standard` or
`cpp_Lite`. Source, configuration, catalogs, observation data, and outputs
stay outside the image on bind-mounted storage.

> 中文版：[README-CN.md](README-CN.md)

## Runtime

The Rocky Linux 8.10 image contains G++ 12.3.0, OpenMPI 4.1.8 with PMI2,
CFITSIO 4.6.4, FFTW 3.3.11, Eigen 3.4.0, LAPACK 3.11.0, and OpenBLAS 0.3.33.

Its portable HPC baseline is x86_64, Slurm `pmi2`,
Apptainer/Singularity, a shared filesystem, and routable TCP. Other
architectures, PMIx-only sites, schedulers, or vendor fabrics require separate
validation.

## Build and verify

### Pull the GHCR image

```bash
docker pull ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest
```

### Download the source and build

Download `cpp-docker-<tag>.zip` from
[GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases),
extract it, and run:

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

The verifier checks component versions, the scientific stack, two MPI ranks,
PMI support, and the absence of pipeline source from the image.

## Local use

```bash
cp .env.example .env
# Set CPP_SOURCE_HOST and the host paths required by the selected phases.
docker compose run --rm FourierQuad-CPP
```

### Common `.env` parameters

Copy `.env.example`, then update these values for the host layout. A
`*_HOST` value is a host path; its paired `*_CONTAINER` value is the
absolute path used by `pipeline.ini` or CLI inside the container.

| Parameter | Typical setting | Constraint |
|---|---|---|
| `IMAGE_NAME` | Image tag to run or build locally. | Must match the pulled image or `docker build -t` value. |
| `BUILD_JOBS` | Image-build parallelism. | Size for available CPU and memory. |
| `HOST_UID`, `HOST_GID` | Current host-user UID/GID. | Set when outputs must be directly writable by the host user. |
| `CPP_SOURCE_HOST` | `cpp_Standard` or `cpp_Lite` source directory. | Mounted read/write at `/workspace/src_pipe` so configuration and build products persist. |
| `SCIENCE_ROOT_*`, `DQ_ROOT_*` | Science/DQ archive paths. | Enable through `compose.optional.yaml` when required; Lite always needs DQ for main processing. |
| `ASTROMETRY_CAT_*`, `SOURCE_CAT_*`, `FLAT_PATH_*` | Gaia, external-source, and flat directories. | Container paths must match effective INI/CLI values; flat is needed only by its Standard branch. |
| `PROCESS_DATA_*` | Writable processing directory. | Holds exposure lists, intermediates, and results; default container path is `/data/DataProcess`. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Optional independent phase paths. | Add `compose.optional.yaml` when used; otherwise prefer locations below processing data. |

Inside the container:

```bash
make -C /workspace/src_pipe -j4
cp /workspace/src_pipe/pipeline.example.ini \
   /workspace/src_pipe/pipeline.ini
# Edit pipeline.ini with container paths.
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --config /workspace/src_pipe/pipeline.ini
```

Core binds are source, astrometry/source catalogs, flat calibration, and
writable processing data. Science/DQ archives and independent phase paths are
optional and should be mounted only when the selected phases need them.
Run-time paths belong in `pipeline.ini` or supported CLI options; fixed
`config/*.hpp` values require rebuilding.

For `process_astrocat`, bind the raw Gaia directory read-only and write tiles
to a writable path, normally below the processing-data bind. Its producer
output is independent of `[lensing].astrometry_cat`; configure that consumer
path separately. Standard must also set `[lensing].astrometry_cat_type = 2`;
Lite already uses the one-degree layout.

## Slurm

Convert the same image to one SIF and follow the
[runner guide](runner/README.md). The supported launch boundary is:

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

Dependency sources and licenses are recorded in [SOURCES.md](SOURCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
