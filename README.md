# Fourier_Quad C++ Pipeline

Modern C++17 and MPI implementation of the Fourier_Quad weak-lensing shear
measurement pipeline. It can prepare exposure lists, process CCD images through
the nine-stage shear pipeline, rearrange output catalogs, and perform the
field-distortion shear test from one executable.

> 中文使用说明见 [CPP_GUIDE_CN.md](CPP_GUIDE_CN.md). The complete English
> reference is [CPP_GUIDE.md](CPP_GUIDE.md), and every configuration option is
> listed in [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md).

## Choose a version

| Version | Recommended for |
|---|---|
| `cpp_Standard/` | New users and the complete feature set, including optional PCA PSF reconstruction |
| `cpp_Lite/` | A smaller frozen production variant with unused feature branches removed |

Both versions build `Fourier_Quad_Pipe` and use the same command-line options.
Start with `cpp_Standard` unless you specifically need the Lite configuration.

Clone the repository with the URL shown by GitHub's **Code** button, or download
the selected source and Docker archives from the repository's **Releases**
page. The Docker quick start requires both the source tree and `cpp_docker/`.

## Quick start with Docker

The published image provides a reproducible G++/OpenMPI scientific stack. It
does not contain the pipeline source or observation data; both are mounted from
the host when the container starts.

Requirements: Docker with Compose support and an x86-64 Linux host.

### 1. Prepare an exposure list

For a first run, use an isolated writable processing tree rather than the raw
archive. A top-level exposure list contains one per-exposure list and its chip
count on each line:

```text
"/data/DataProcess/g2019/stamps/exposure_001.list" 5
```

The referenced per-exposure list contains one Science FITS path per line. All
paths must be valid inside the container. The integrated initializer can create
this tree and both list types automatically; see
[Starting from compressed archives](#starting-from-compressed-archives).

### 2. Configure catalog paths and stages

In `cpp_Standard/config/` or `cpp_Lite/config/`:

- set the astrometry catalog, source catalog, and flat-field paths in
  `LensingConfig.hpp` to their container paths;
- review the nine-stage `PROCESS_stage` switch and science parameters;
- review the default process switches in `ProcessConfig.hpp`.

Command-line options override the process switches for each run. Use
`./Fourier_Quad_Pipe --help` after compiling to see the effective interface.

### 3. Configure the container mounts

From the repository root:

```bash
cd cpp_docker
cp .env.example .env
```

Edit `.env` and set:

- `IMAGE_NAME=ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest`;
- `CPP_SOURCE_HOST` to the absolute host path of `cpp_Standard/` or
  `cpp_Lite/`;
- the catalog, flat-field, and writable processing-data host paths;
- `HOST_UID` and `HOST_GID` to your numeric user and group IDs when necessary.

The catalog `*_CONTAINER` destinations must match the paths compiled into
`LensingConfig.hpp`.

### 4. Compile and run the numerical pipeline

Pull the image and enter the container:

```bash
docker compose pull
docker compose run --rm FourierQuad-CPP
```

Inside the container:

```bash
make clean
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-extcat false \
  --run-init false \
  --run-main true \
  --run-rearr false \
  --run-fd false \
  --expo-list /data/DataProcess/expo_g2019.list
```

This explicit command runs only the numerical nine-stage pipeline, regardless
of the defaults in `ProcessConfig.hpp`. Increase or reduce the MPI rank count to
match the number of exposures and available resources.

### 5. Find the results

The primary products are written below the dataset root derived from the FITS
paths:

- `result/<EXPOSURE>_all.cat`: combined per-exposure shear catalog;
- `stamps/`: type-specific intermediate products;
- `expo_info.dat`: exposure-level diagnostics.

If catalog rearrangement is enabled, spatially partitioned catalogs and
`catalog_summary.txt` are also written below the configured rearrangement
output directory.

## Starting from compressed archives

The initializer reads Science and DQ `.fits.fz` archives, creates an isolated
processing tree, extracts chip images, and generates `expo_<target>.list`.
With the Science and DQ archives mounted, a chained run is:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-extcat false \
  --run-init true \
  --run-main true \
  --run-rearr true \
  --run-fd false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dqmask \
  --output-root /data/DataProcess \
  --dataset g2019:c4d_19 \
  --contains v1 \
  --existing fail
```

Repeat `--dataset TARGET:PREFIX` to process multiple datasets sequentially and
repeat `--contains TOKEN` to accept any matching archive basename. Use
`--existing resume` only after reviewing an interrupted output tree.

For Docker, define the optional host paths in `cpp_docker/.env` and include
`compose.optional.yaml` when starting the container. The optional Compose file
also exposes mounts for external-catalog input, rearranged catalogs, exposure
lists, and field-distortion output:

```bash
docker compose -f compose.yaml -f compose.optional.yaml \
  run --rm FourierQuad-CPP
```

Every host variable referenced by `compose.optional.yaml` must name an existing
directory. Point unused writable outputs to dedicated empty directories rather
than to scientific input data.

See [the initializer output contract](CPP_GUIDE.md#initializer-output-contract)
before processing a new archive layout.

## Build from source

Required software:

- a C++17 compiler and MPI with `mpicxx`;
- CFITSIO, FFTW3/FFTW3F, Eigen3, LAPACK, and BLAS;
- GNU Make.

When the headers and libraries share one prefix:

```bash
export SCIENCE_PREFIX=/path/to/scientific-stack
make -C cpp_Standard -j4 STACK_PREFIX="$SCIENCE_PREFIX"
mpirun -np 4 ./cpp_Standard/Fourier_Quad_Pipe \
  --run-extcat false --run-init false --run-main true \
  --run-rearr false --run-fd false \
  --expo-list /path/to/expo_list.list
```

Replace `cpp_Standard` with `cpp_Lite` for the Lite variant. If Eigen is in a
different prefix, also pass `EIGEN_INCLUDE=/path/to/eigen3`.

## Run on a Slurm cluster

The included Apptainer/Singularity runner uses the same image and executable:

```bash
cd cpp_docker/runner
cp cpppipeline.env.example cpppipeline.env
```

Edit `cpppipeline.env` with the GHCR image, SIF destination, source/data mounts,
site modules, and Slurm settings. Then validate in order:

*Note:* GHCR image is available as `ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest`.

```bash
bash pull-sif.sh
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch mpi-smoke-test.slurm
sbatch cpppipeline.slurm \
  --run-extcat false --run-init false --run-main true \
  --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_g2019.list
```

Pin `OCI_IMAGE_URI` by digest for production. The supplied runtime expects a
Slurm site with PMI2 support; verify the cluster interface before a multi-node
science run. See [the runner guide](cpp_docker/runner/README.md) for resource
templates and site-specific configuration.

## Manual for AI

Here also provides a manual wraped as a codex/claude code plugin, which gives agent abilities to help you build environment, switch parameters and run pipeline. For more details, please refer to the [Codex/Claude Code Plugin](https://github.com/Syoong-s/FQNeoAIManual).


## Documentation

- [C++ source, run modes, Docker, and runner guide](CPP_GUIDE.md)
- [中文完整指南](CPP_GUIDE_CN.md)
- [Complete parameter reference](CPP_PIPELINE_PARAMETERS.md)
- [Container environment](cpp_docker/README.md)
- [HPC runner](cpp_docker/runner/README.md)

## License

Distributed under the [MIT License](LICENSE).
