# Fourier_Quad C++ Pipeline

> 中文版：[README_CN.md](README_CN.md)

MPI-parallel C++17 implementation of the Fourier_Quad weak-lensing pipeline for
DECam data. One executable, `Fourier_Quad_Pipe`, can repartition an external
catalog, initialize compressed Science/DQ archives, run the nine numerical
stages, rearrange the resulting catalogs, and perform the field-distortion (FD)
shear test.

Chinese users can start with [CPP_GUIDE_CN.md](CPP_GUIDE_CN.md).

## Choose a variant

| Variant | Use it when |
|---|---|
| [`cpp_Standard`](cpp_Standard/) | You need optional flat, mask, astrometry, external-PSF, hybrid-PSF, or PCA PSF branches. |
| [`cpp_Lite`](cpp_Lite/) | You use the fixed production path: Gaia astrometry, DQ masks, external source catalogs, local-polynomial PSF, and no PCA. |

Both variants share the executable name and command-line interface. Lite omits
the unused branches from its source; it is not Standard with different defaults.

## Quick start

### 1. Download a Release source package

Open [GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)
and download only the source package for the variant you plan to run:

| Variant | Release source package |
|---|---|
| C++ Lite | `cpp-lite-<tag>.zip` |
| C++ Standard | `cpp-standard-<tag>.zip` |

Use a fixed Release so the software version used for an analysis can be
recorded and reproduced.

### 2. Prepare the environment

Install the dependencies in [Environment prerequisites](#environment-prerequisites).

### 3. Prepare input data

Prepare the four input classes described in
[Input data requirements](#input-data-requirements). Science images, the Gaia
catalog, and the external source catalog are required. DQ masks depend on the
selected variant and configuration.

### 4. Configure the pipeline

From the extracted `cpp_Lite` or `cpp_Standard` directory, copy the complete
run-time template and edit paths, datasets, phases, and run-time science choices:

```bash
cp pipeline.example.ini pipeline.ini
```

The configuration precedence is:

```text
compiled defaults < --config INI < command-line options
```

Fixed numerical settings that are not represented in the INI remain in
`config/*.hpp` and require rebuilding after changes. The complete run-time and
compile-time reference is [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md).

### 5. Build and run

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Docker and HPC execution guides:

| Implementation | Local Docker | HPC / Slurm + Apptainer |
|---|---|---|
| C++ | [C++ Docker guide](cpp_docker/README.md) | [C++ runner guide](cpp_docker/runner/README.md) |

See [CPP_GUIDE.md](CPP_GUIDE.md) for complete build variables, phase selection,
run modes, inputs, outputs, and failure rules.

## Environment prerequisites

| Category | Prerequisite | Regular Linux | HPC / Slurm | Notes |
|---|---|---:|---:|---|
| Operating system | 64-bit Linux | Required | Required | A modern Linux distribution is recommended. |
| C++ compiler | MPI C++ wrapper with C++17 support | Required | Required | The pinned container uses GCC 12.3.0. |
| MPI | OpenMPI or a compatible MPI implementation | Required | Required | The pinned container uses OpenMPI 4.1.8. |
| CFITSIO | CFITSIO development library | Required | Required | FITS/FZ I/O; container version 4.6.4. |
| FFTW3 | Double- and single-precision FFTW3 libraries | Required | Required | Fourier transforms; container version 3.3.11. |
| Eigen3 | Eigen3 headers | Required | Required | C++ numerical operations; container version 3.4.0. |
| BLAS / LAPACK | BLAS and LAPACK libraries | Required | Required | The container uses LAPACK 3.11.0 with OpenBLAS 0.3.33. |
| Build tools | `make`, shell, standard GNU tools | Required | Required | Builds and helper scripts. |
| Shared filesystem | All ranks/nodes see the same inputs and outputs | No | Required | Required for multi-node jobs. |
| Slurm | Site scheduler and `srun` | No | Required for Slurm jobs | Batch allocation and launch. |
| PMI2-compatible launch | Slurm `pmi2` support and compatible MPI | No | Required for the repository HPC container path | Used for direct Slurm MPI launch. |
| Apptainer / Singularity | Rootless container runtime | No | Required when using the repository HPC container path | Must be available on compute nodes. |
| Writable processing directory | Shared output/work directory writable by all ranks | Required | Required | Source data and outputs should not be mixed. |

## Input data requirements

| Input | Required | Purpose | Minimum requirement |
|---|---:|---|---|
| Science images | Yes | Exposures on which the pipeline performs source detection, shape measurement, and weak-lensing processing. | Supported Science FITS/FZ data whose file organization matches the configured exposure and CCD recognition rules. |
| Gaia catalog | Yes | High-precision sky-coordinate reference for source matching and astrometric calibration. | Covers the Science-image footprint; each file has one header line followed by rows whose first two fields are numeric `ra` and `dec`. |
| External source catalog | Yes | Supplies sky positions, `zp`, and photometry for external-source matching and downstream processing. | Contains `ra`, `dec`, `zp`, and a magnitude in at least one observed band. |
| DQ masks | Configuration-dependent (optional input class) | Marks bad, saturated, defective, or otherwise invalid pixels. | May be omitted only when the selected variant and configuration do not read DQ masks. |

### Science images

Science images are the primary scientific exposures, not calibration catalogs or
an output directory. They must use a FITS/FZ format supported by the selected
variant, be readable through its FITS logic, and follow the exposure/CCD naming
and list conventions described in the detailed guide.

### Gaia catalog

The Gaia catalog supplies accurate RA/Dec reference positions for object
matching and astrometric calibration. It must cover the actual Science-image
footprint and be stored directly under `[lensing].astrometry_cat` (whose
compiled default is `ASTROMETRY_CAT`). Every tile consumed by the pipeline
must have one header line, followed by rows whose first two numeric fields are
RA and Dec. Rows may be comma- or whitespace-separated; additional fields are
ignored. `[lensing].astrometry_cat_type` selects one of two layouts:

The RA and Dec used in the following lookup formula are the WCS reference
coordinates read from the current Science CCD header.

**Type 1 (legacy large tiles, the default):**

- `|Dec| < 80°`: `gaia_<p|m><D>_<RR>.cat`, where
  `D = floor(|Dec| / 10) + 1` (1-8) and `RR = floor(RA / 10)` (00-35,
  zero-padded).
- `|Dec| >= 80°`: `gaia_<p|m>9.cat`, without an RA suffix.
- `p` denotes nonnegative Dec; `m` denotes negative Dec.
- Each file must contain one header line.

> Examples:
>
> 1. `gaia_p1_00.cat` covers `0° <= RA < 10°` and `0° <= Dec < 10°`.
> 2. `gaia_m3_12.cat` covers `120° <= RA < 130°` and `-30° < Dec <= -20°`.
> 3. `gaia_p9.cat` covers `0° <= RA < 360°` and `80° <= Dec <= 90°`.

*To ensure that stars can still be selected for exposures located at the edges
of the 10° × 10° grid, it is recommended to expand a single catalog's upper and
lower Dec limits by 2°. Around absolute-Dec ranges beginning at 0°, 30°, and
60°, respectively, expand its upper and lower RA limits by 2°, 4°, and 6°.*

**Type 2 (one-degree tiles):**

- `<ASTROMETRY_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`. The default
  prefix in `config/pathconfig.hpp` is `astra_` and does not include `RA_`.
  RA boundaries use three digits and Dec boundaries use signed two-digit values.
- Example: `astra_RA_123_124_Dec_m05_m04.dat` covers
  `123° <= RA < 124°` and `-5° <= Dec < -4°`.

The optional `process_astrocat` phase converts direct regular files from a raw
Gaia directory to Type 2 tiles, removes coordinate pairs that are exact or
within one ULP in both RA and Dec (including tile boundaries), and writes the
required `RA    DEC` header. Its `[astrocat].output_directory` or
`--astrocat-output` setting controls only the producer output directory; it is
not checked against and does not update `[lensing].astrometry_cat`. To consume
the result later, set `[lensing].astrometry_cat` to that directory and set
`[lensing].astrometry_cat_type = 2` separately.

### External source catalog

The minimum schema is survey- and band-independent:

| Field | Meaning |
|---|---|
| `ra` | Right Ascension. |
| `dec` | Declination. |
| `zp` | The photometric-redshift quantity consumed by the selected pipeline configuration. |
| One observed-band magnitude | A magnitude in any one band used by the selected analysis. |

Additional colors, redshifts, object classes, shapes, and flags may be retained,
but they are not part of this minimum input contract. Configure the actual
column positions, delimiter, header handling, and projection in the `[extcat]`
INI section, through the corresponding command-line options, or with their
compiled defaults in `config/ExtCatConfig.hpp`.

**Filename convention:**

- 1° × 1° tile:
  `<SOURCE_CAT_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`.
- The default `SOURCE_CAT_TILE_PREFIX` in `config/pathconfig.hpp` is `extern_`;
  the prefix does not include `RA_`.
- RA boundaries use three digits. Dec boundaries use `p` or `m` plus a two-digit
  absolute value. Each upper boundary is one degree above its lower boundary.
- Each file must contain one header line.

> Example: `extern_RA_123_124_Dec_m05_m04.dat` covers
> `123° <= RA < 124°` and `-5° <= Dec < -4°`.

### DQ masks (optional)

DQ masks identify pixels that must not participate in scientific measurements,
including bad pixels, saturation, detector defects, and other invalid regions.
They are optional only when the chosen variant and configuration disable DQ
access. Standard users can select the relevant mask mode with
`[lensing].include_mask`;
Lite is fixed to per-chip DQ masks, so Lite runs must provide them. If DQ masks
are omitted, verify that no active branch or configured path still reads them.

## Pipeline

Top-level phases always run in this order:

```text
process_astrocat -> process_extcat -> process_init -> process_main -> process_rearr -> process_fd
```

`process_astrocat` publishes deduplicated one-degree Gaia tiles. Its output
directory is independent of the Gaia directory later consumed through
`[lensing].astrometry_cat`.

`process_main` contains nine prime-gated stages from preprocessing through
catalog combination. The full stage and data contracts are in
[CPP_GUIDE.md](CPP_GUIDE.md); run-time and compile-time settings are separated
in [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md).

Principal outputs are generated exposure lists, per-exposure `*_all.cat`
catalogs, rearranged `subcat_*.cat` catalogs, and `fdout/FD_test_comb.dat`.
Input archives and catalogs are read in place.

## Containers and HPC

[`cpp_docker`](cpp_docker/) provides an x86_64 Docker toolchain image. The
image contains compilers and libraries, not pipeline source or observation
data; both remain bind-mounted. For Slurm, use the production
[Apptainer runner](cpp_docker/runner/README.md) after confirming that the site
provides the PMI2 launch plugin.

## Documentation

| Document | Purpose |
|---|---|
| [C++ guide](CPP_GUIDE.md) / [中文指南](CPP_GUIDE_CN.md) | Build, configure, run, inputs, outputs, and failure rules |
| [Parameter reference](CPP_PIPELINE_PARAMETERS.md) | One complete Standard/Lite table per configuration header, including INI/CLI overrides |
| [Container guide](cpp_docker/README.md) / [中文](cpp_docker/README-CN.md) | Docker image and local container workflow |
| [Slurm runner](cpp_docker/runner/README.md) / [中文](cpp_docker/runner/README-CN.md) | Apptainer/Singularity deployment |

## License

Repository-authored code is distributed under the [MIT License](LICENSE).
Container dependencies retain their upstream licenses; see
[third-party notices](cpp_docker/THIRD_PARTY_NOTICES.md).
