# Fourier_Quad C++ Pipeline

MPI-parallel C++17 implementation of the Fourier_Quad weak-lensing shear
measurement pipeline for DECam data. The `Fourier_Quad_Pipe` executable prepares
catalogs, initializes Science/DQ archives, measures shear, rearranges catalogs,
and performs field-distortion analysis.

Build and run the pipeline directly on Linux or an HPC system. Docker and
Apptainer are optional ways to provide its software dependencies.

## Choose a variant

| Variant | Source | Use case |
| --- | --- | --- |
| **Standard** | [`cpp_Standard/`](cpp_Standard/) | Full set of scientific branches, including optional flat correction, external PSF, hybrid PSF and PCA paths. |
| **Lite** | [`cpp_Lite/`](cpp_Lite/) | Fixed production branches with unused alternatives removed. |

Both build `Fourier_Quad_Pipe`. Their configuration files share the same INI
structure, but Lite rejects runtime keys for branches it does not implement.
See the [variant comparison](CPP_GUIDE.md#choose-a-variant).

## Build and run without containers

Use a C++17 MPI compiler and install CFITSIO, FFTW3 (double and single precision),
Eigen3, BLAS and LAPACK. Download or clone the source, then:

```bash
cd cpp_Standard                 # Or cpp_Lite
cp pipeline.example.ini pipeline.ini
# Edit pipeline.ini: input paths, datasets, enabled phases and scientific settings.
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Before running, follow [data preparation](DATA_PREPARATION.md). The example INI
contains commented entries, not a ready-to-run survey configuration. Uncomment
and fill the values you need to override. Otherwise the program retains the
selected variant's compiled defaults, including its enabled phases and original
site-specific paths.

```text
Prepare data → Choose Standard or Lite → Configure → Build → Run → Inspect outputs
```

Only `--config PATH` (or `--config=PATH`) is accepted by the executable. With **no
arguments**, it uses the compiled `config/*.hpp` values and does not look for an
INI automatically. Editing a header requires rebuilding; editing an INI does not.

## Documentation

| Guide | What it covers |
| --- | --- |
| [C++ user guide](CPP_GUIDE.md) | Native build, variant selection, six phases, INI configuration, run recipes and outputs. |
| [Data preparation](DATA_PREPARATION.md) | Science images, Gaia tile layouts, external source catalogs, DQ masks and exposure lists. |
| [Complete parameter reference](CPP_PIPELINE_PARAMETERS.md) | Every configuration-header parameter, Standard/Lite defaults, INI mappings, constraints and rebuild requirements. |
| [Docker](DOCKER.md) | Optional workstation deployment, JSON configuration and one-command launch. |
| [HPC containers](HPC.md) | Optional Apptainer deployment called from your own job script. |
| [Container paths](PATHS.md) | Every path interface, shared mounts, related paths, additional directories and files. |
| [Upgrade from `9925f4f`](UPGRADE.md) | File replacement, removals and migration from the older command-line options. |

## Processing overview

| Phase | Input | Main result |
| --- | --- | --- |
| `process_astrocat` | Raw two-column Gaia catalogs | One-degree Gaia tiles |
| `process_extcat` | Raw external source catalogs | Canonical source tiles |
| `process_init` | Science/DQ FITS/FZ archives | Extracted CCD images and exposure lists |
| `process_main` | Images, catalogs and selected calibrations | Per-exposure shear catalogs |
| `process_rearr` | Per-exposure catalogs | Spatially partitioned subcatalogs |
| `process_fd` | Rearranged catalogs and metadata | Field-distortion measurements |

Enable phases in the INI `[process]` section. Standard defaults to initialization,
main processing, rearrangement and FD; Lite defaults to initialization and main
processing. Catalog preparation is disabled by default in both.

## Optional container deployment

Configure and compile once using the relevant guide. Subsequent launches are:

```bash
# Docker, from the repository root
python3 cpp_docker/fqctl.py run
```

```bash
# Apptainer, called once inside your own job script
bash /absolute/repo/cpp_docker/runner/run.sh /absolute/run/cpppipeline.env
```

Containers supply the toolchain, not observation data or scientific choices.
Both launchers pass an explicit `--config` pointing to a generated INI outside
the source directory.

## Versioning and license

Use a fixed [release](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)
or commit for an analysis. Record the source version, local header changes,
INI and software environment with your results. Existing releases may contain
separate Standard/Lite archives; preserve the directory structure of the package
you download.

An optional [Codex/Claude Code manual](https://github.com/Syoong-s/FQNeoAIManual)
is maintained separately.

Repository-authored code uses the [MIT License](LICENSE). Container dependencies
retain their [upstream licenses](THIRD_PARTY_NOTICES.md).
