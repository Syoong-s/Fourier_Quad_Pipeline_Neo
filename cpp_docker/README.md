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
