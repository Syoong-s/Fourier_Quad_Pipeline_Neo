# Fourier_Quad C++ Pipeline

> 中文版：[README_CN.md](README_CN.md)

## Overview

MPI-parallel C++17 implementation of the Fourier_Quad weak-lensing shear
measurement pipeline for DECam data. The unified executable
`Fourier_Quad_Pipe` provides six top-level phases, from catalog preparation
and archive initialization through shear catalogs and the field-distortion
(FD) test.

## Pipeline variants

| Pipeline | Source directory | Purpose |
|---|---|---|
| C++ Lite | [`cpp_Lite`](cpp_Lite/) | Fixed production path with unused alternate branches removed. |
| C++ Standard | [`cpp_Standard`](cpp_Standard/) | Full branch set, including optional scientific paths. |

Both variants use the same executable name, INI format, and command-line
interface. Lite physically omits Standard-only branches; they cannot be enabled
with configuration values.

## Quick start

### 1. Download a Release source package

Open [GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)
and download the source package for the variant you plan to run:

| Pipeline | Release source package |
|---|---|
| C++ Lite | `cpp-lite-<tag>.zip` |
| C++ Standard | `cpp-standard-<tag>.zip` |

Use a fixed Release so the software version used for an analysis can be
recorded and reproduced.

### 2. Prepare the environment

Install the dependencies in [Environment prerequisites](#environment-prerequisites).

### 3. Prepare input data

Prepare the four input classes in
[Input data requirements](#input-data-requirements). Science images, a Gaia
catalog, and an external source catalog are required. DQ masks depend on the
selected variant and configuration.

### 4. Configure the Pipeline

From the extracted source directory, copy the complete run-time template and
edit the paths, datasets, phases, and run-time science choices:

```bash
cp pipeline.example.ini pipeline.ini
```

Configuration precedence is:

```text
compiled defaults < --config INI < command-line options
```

Fixed settings that are not represented in the INI remain in `config/*.hpp`
and require rebuilding after a change. See
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) for the complete
run-time and compile-time reference.

*For `cpp_Lite`, first set ordinary compiled defaults in `Initialize.hpp`;
use the INI or CLI for per-run overrides.*

### 5. Build and run

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Docker and HPC documentation entrances:

| Environment | Guide |
|---|---|
| Local Docker | [C++ Docker guide](cpp_docker/README.md) |
| HPC / Slurm + Apptainer | [C++ runner guide](cpp_docker/runner/README.md) |

See the [C++ guide](CPP_GUIDE.md) for complete build variables, phase
selection, inputs, outputs, run modes, and failure rules.

## Environment prerequisites

| Environment | Prerequisite | Reference stack |
|---|---|---|
| Regular Linux | 64-bit Linux; MPI C++ compiler with C++17 support; OpenMPI or compatible MPI; CFITSIO; FFTW3 (double and single precision); Eigen3; BLAS / LAPACK | GCC 12.3.0; OpenMPI 4.1.8; CFITSIO 4.6.4; FFTW 3.3.11; Eigen 3.4.0 |
| HPC / Slurm | Regular Linux prerequisites; shared filesystem; Slurm; PMI2-compatible launch; Apptainer / Singularity | Same versions as the published x86_64 container image |

## Input data requirements

| Input | Purpose | Minimum requirement |
|---|---|---|
| Science images | Exposures used for source detection, shape measurement, and weak-lensing processing. | Supported Science FITS/FZ data matching the configured archive, exposure, and CCD naming rules. |
| Gaia catalog | Sky-coordinate reference for source matching and astrometric calibration. | Covers the Science footprint and follows the configured Type 1 or Type 2 tile layout. |
| External source catalog | Supplies sky positions, `zp`, and photometry for matching and downstream processing. | Canonical one-degree tiles containing `ra`, `dec`, `zp`, and at least one observed-band magnitude. |
| DQ masks (*configuration-dependent*) | Marks bad, saturated, defective, or otherwise invalid pixels. | Required by Lite and by every Standard configuration that reads DQ data. |

The exact catalog filenames, schemas, producer/consumer path relationships,
and archive conventions are documented in the
[C++ guide](CPP_GUIDE.md). Configure run-time values in `pipeline.ini`; use
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) to identify settings
that remain compile-time only.

## Detailed guides

- [C++ Pipeline guide](CPP_GUIDE.md) /
  [中文指南](CPP_GUIDE_CN.md)
- [C++ parameter reference](CPP_PIPELINE_PARAMETERS.md)
- [C++ Docker guide](cpp_docker/README.md) /
  [中文指南](cpp_docker/README-CN.md)
- [C++ Slurm/Apptainer runner](cpp_docker/runner/README.md) /
  [中文指南](cpp_docker/runner/README-CN.md)

## License

Repository-authored code is distributed under the [MIT License](LICENSE).
Container dependencies retain their upstream licenses; see the
[third-party notices](cpp_docker/THIRD_PARTY_NOTICES.md).
