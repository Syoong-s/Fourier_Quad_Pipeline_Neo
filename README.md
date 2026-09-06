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
edit the paths, datasets, phases, and run-time science choices。

Available configuration methods and their priorities are as follows:

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

### Science images

Science images are the primary scientific exposures, not calibration catalogs or
an output directory. They must use a FITS/FZ format supported by the selected
variant, be readable through its FITS logic, and follow the exposure/CCD naming
and list conventions described in the detailed guide.

### Gaia catalog

The Gaia catalog supplies accurate RA/Dec reference positions for object
matching and astrometric calibration. It must cover the actual Science-image
footprint and be stored directly under the configured `ASTROMETRY_CAT`
directory. Every consumed tile starts with one header line; subsequent rows use
the first two numeric fields as RA and Dec. Rows may be comma- or
whitespace-separated and additional fields are ignored.

`LensingConfig::AstroCatType` selects one of two filename layouts. Type `1` is
the legacy large-tile layout below. Type `2` uses the one-degree layout and may
be generated from raw flat-directory files by `process_astrocat`.

**Filename convention:**

- `|Dec| < 80°`: `gaia_<p|m><D>_<RR>.cat`, where
  `D = floor(|Dec| / 10) + 1` (1-8) and `RR = floor(RA / 10)` (00-35,
  zero-padded).
- `|Dec| >= 80°`: `gaia_<p|m>9.cat`, without an RA suffix.
- `p` denotes nonnegative Dec; `m` denotes negative Dec.
- Each file must contain one header line.

> Examples:
> 1. gaia_p1_00.cat covers `0° <= RA < 10°` and `0° <= Dec < 10°`
> 2. gaia_m3_12.cat covers `120° <= RA < 130°` and `-30° < Dec <= -20°`
> 3. gaia_p9.cat covers `0° <= RA < 360°` and `80° <= Dec <= 90°`

*For Type 1, to ensure that stars can still be selected for exposures located at the edges of the 10°×10° grid, it is recommended to expand the upper and lower Dec coverage limits of a single star catalog by 2° based on the aforementioned limits, and expand the upper and lower RA limits by 2°, 4°, and 6° within the 0°, 30°, and 60° ranges, respectively.*

Type 2 files use
`<ASTROMETRY_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`. The default
prefix in `config/pathconfig.hpp` is `astra_`; it does not include `RA_`.
RA bounds use three digits and Dec bounds use signed `p`/`m` two-digit values.
For example, `astra_RA_123_124_Dec_m05_m04.dat` covers
`123° <= RA < 124°`, `-5° <= Dec < -4°`.

`process_astrocat` reads only direct regular files from `--astrocat-input`,
produces these Type-2 tiles, and removes exact/one-ULP coordinate duplicates.
Its `--astrocat-output` option controls only publication. To consume the result,
set `ASTROMETRY_CAT` to that directory and set `AstroCatType=2` separately.

### External source catalog

The minimum schema is intentionally survey- and band-independent:

| Field | Meaning |
|---|---|
| `ra` | Right Ascension. |
| `dec` | Declination. |
| `zp` | The catalog `zp` quantity consumed by the selected Pipeline configuration. |
| One observed-band magnitude | A magnitude in any one band used by the selected analysis. |

Additional colors, redshifts, object classes, shapes, and flags may be retained,
but they are not part of this minimum input contract. Configure actual column
positions, names, delimiter, header handling, and projection in
`config/ExtCatConfig.hpp` or with its documented CLI overrides.

**Filename convention:**

- 1° × 1° tile:
  `<SOURCE_CAT_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`.
- The default `SOURCE_CAT_TILE_PREFIX` in `config/pathconfig.hpp` is `extern_`;
  the prefix does not include `RA_`.
- RA boundaries use three digits. Dec boundaries use `p` or `m` plus a two-digit
  absolute value. Each upper boundary is one degree above its lower boundary.
- Each file must contain one header line.

> Example:
> `extern_RA_123_124_Dec_m05_m04.dat` covers `123° <= RA < 124°` and
> `-5° <= Dec < -4°`.

Use the C++ `process_extcat` tool to convert raw downloads into these canonical
pipeline tiles. Raw downloader filenames are source artifacts and are not a
second C++ naming convention.

### DQ masks (optional)

DQ masks identify pixels that must not participate in scientific measurements,
including bad pixels, saturation, detector defects, and other invalid regions.
They are optional only when the chosen configuration disables DQ access. Standard
users can select the relevant mask mode with `include_Mask`; current C++ Lite is
fixed to per-chip DQ masks, so Lite runs must provide them. If DQ masks are
omitted, verify that no configured path or active branch still reads them.

The exact catalog filenames, schemas, producer/consumer path relationships,
and archive conventions are documented in the
[C++ guide](CPP_GUIDE.md). Configure run-time values in `pipeline.ini`; use
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) to identify settings
that remain compile-time only.

## Codex/Claude Code Plugin

Here also provides a manual wraped as a codex/claude code plugin, which gives agent abilities to help you build environment, switch parameters and run pipeline. For more details, please refer to the [Codex/Claude Code Plugin](https://github.com/Syoong-s/FQNeoAIManual).

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
