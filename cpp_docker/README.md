# Fourier_Quad C++ container

This directory builds one x86_64 Linux toolchain image for `cpp_Standard` or
`cpp_Lite`. Pipeline source, configuration, catalogs, observation data, and
outputs are not copied into the image; they remain on bind-mounted storage.

> 中文版：[README-CN.md](README-CN.md)

## Runtime

| Component | Version |
|---|---|
| Rocky Linux | 8.10 |
| G++ | 12.3.0 |
| OpenMPI | 4.1.8 with PMI2 support |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / OpenBLAS | 3.11.0 / 0.3.33 |

The portable HPC baseline is x86_64, Slurm `pmi2`, Apptainer/Singularity, a
shared filesystem, and routable TCP. PMIx-only sites, ARM, other schedulers,
and vendor-fabric acceleration need separate qualification.

## Build and verify

### Pull the GHCR image

```bash
docker pull ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest
```

### Download the source and build

Download `cpp-docker-<tag>.zip` for the selected tag from
[GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases), extract it,
then run:

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

The verification checks component versions, the scientific stack, two MPI
ranks, PMI support, and the absence of pipeline source from the image.

## Local use

```bash
cp .env.example .env
# Set CPP_SOURCE_HOST and every host path used by the selected phases.
docker compose run --rm FourierQuad-CPP
```

### Common `.env` parameters

Copy `.env.example`, then adapt the table to the host directories and selected phases.
`*_HOST` values are host paths; `*_CONTAINER` values are the absolute paths seen by the
pipeline and `pipeline.ini` inside the container.

| Parameter | Typical change | Constraint |
|---|---|---|
| `IMAGE_NAME` | Set the image tag to run or build locally. | It must match `docker build -t` or the pulled image. |
| `BUILD_JOBS` | Set image-build parallelism. | Size it for available host CPU and memory. |
| `HOST_UID`, `HOST_GID` | Set the current host user's UID/GID. | Change when host users must directly own and edit outputs. |
| `CPP_SOURCE_HOST` | Point to a `cpp_Standard` or `cpp_Lite` source directory. | It is mounted read/write at `/workspace/src_pipe` for `pipeline.ini` and build products. |
| `SCIENCE_ROOT_HOST/CONTAINER` | Set the Science-image archive and its container path. | Used by `process_init`; match `[init].science_root` or CLI and enable it through `compose.optional.yaml`. |
| `DQ_ROOT_HOST/CONTAINER` | Set the DQ-mask archive and its container path. | Required whenever DQ is read and always for Lite; match `[init].dq_root` or CLI. |
| `ASTROMETRY_CAT_HOST/CONTAINER` | Point to the Gaia catalog directory. | The container path must match `[lensing].astrometry_cat`; Type 1/2 must also match `[lensing].astrometry_cat_type`. |
| `SOURCE_CAT_HOST/CONTAINER` | Point to normalized external-source tiles. | Match the effective `[extcat].output_directory`, `[lensing].source_cat`, or `--extcat-output` path. |
| `FLAT_PATH_HOST/CONTAINER` | Point to flat calibration data. | Needed only for Standard with `[lensing].include_flat=1`; match `[lensing].flat_path`. |
| `PROCESS_DATA_HOST/CONTAINER` | Point to writable processing storage. | Holds exposure lists, intermediates, and results; container default is `/data/DataProcess`. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Set only for phases that need independent mounts. | Add `compose.optional.yaml` and use the matching container paths in INI/CLI; otherwise prefer locations below processing data. |

Inside the container:

```bash
make -C /workspace/src_pipe -j4
cp /workspace/src_pipe/pipeline.example.ini \
   /workspace/src_pipe/pipeline.ini
# Edit pipeline.ini with container paths, not host paths.
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --config /workspace/src_pipe/pipeline.ini
```

Core binds are source (`/workspace/src_pipe`), astrometry/source catalogs,
flat calibration, and writable processing data. Science/DQ archives and the
extcat/rearr/exposure-list/FD mounts are optional; enable the latter group with
`compose.optional.yaml` only when needed.

The current program can set astrometry/source paths and other run-selectable
values in `pipeline.ini`. Fixed numerical constants still come from
`config/*.hpp` and require rebuilding.

For `process_astrocat`, expose the raw Gaia directory through a suitable
read-only bind and configure its container path as `[astrocat].input_directory`
or `--astrocat-input`. Catalog binds are read-only, so configure
`[astrocat].output_directory` or `--astrocat-output` as a writable location,
normally below the processing-data bind. This producer output is independent
of `[lensing].astrometry_cat`; a later consumer run must separately point that
setting at the generated directory and set
`[lensing].astrometry_cat_type = 2`.

## Slurm

The same image can be converted to one SIF. Copy
`runner/cpppipeline.env.example`, then follow the
[runner guide](runner/README.md). The supported launch boundary is:

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

Dependency sources and licenses are recorded in [SOURCES.md](SOURCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
