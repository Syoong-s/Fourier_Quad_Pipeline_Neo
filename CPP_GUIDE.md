# C++ user guide

[Home](README.md) · [Data preparation](DATA_PREPARATION.md) · [Complete parameter reference](CPP_PIPELINE_PARAMETERS.md)

This guide covers **native Linux and HPC execution** of `Fourier_Quad_Pipe`.
Docker is not required. The optional [Docker](DOCKER.md) and
[Apptainer](HPC.md) guides cover container-specific setup only.

## Choose a variant

| Feature | Standard | Lite |
| --- | --- | --- |
| Source directory | `cpp_Standard/` | `cpp_Lite/` |
| Executable | `Fourier_Quad_Pipe` | `Fourier_Quad_Pipe` |
| Astrometry | Gaia or identity mapping | Gaia only |
| Super-flat correction | Optional | Disabled |
| DQ handling | Configurable mask mode | Per-chip DQ masks |
| Source selection | External-catalog or no-catalog branch | External catalog |
| PSF input | Exposure stars or external PSF | Exposure stars |
| Deblending | Configurable | Enabled |
| PSF model | Local polynomial or hybrid | Local polynomial |
| PCA/multi-scale reconstruction | Optional | Removed |
| Default enabled phases | Init, main, rearrangement, FD | Init, main |

Lite's removed branches cannot be restored by adding INI keys. Ordinary Lite
camera/workflow defaults are centralized in `Initialize.hpp`; advanced settings
remain in `config/*.hpp`. Standard uses its `config/*.hpp` headers directly.

Lite always consumes one-degree Gaia tiles. Standard defaults to the legacy
large-tile layout and can select one-degree tiles with
`[lensing].astrometry_cat_type = 2`. See [data preparation](DATA_PREPARATION.md).

## Build

### Dependencies

| Requirement | Purpose |
| --- | --- |
| C++17 MPI compiler wrapper, such as `mpicxx` | Build and launch the parallel program |
| CFITSIO | Read and write FITS data |
| FFTW3, double and single precision | Fourier transforms |
| Eigen3 | Matrix calculations |
| BLAS and LAPACK | Numerical linear algebra |
| GNU Make | Build the selected variant |

The repository's optional reference toolchain uses GCC 12.3.0, OpenMPI 4.1.8,
CFITSIO 4.6.4, FFTW 3.3.11 and Eigen 3.4.0. Native builds can use a compatible
site-provided stack; compile and run against compatible MPI installations.

```bash
cd cpp_Standard                 # Or cpp_Lite
make -j4
```

For dependencies under a common prefix:

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack -j4
```

For a separate Eigen installation:

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack \
     EIGEN_INCLUDE=/opt/eigen/include/eigen3 -j4
```

| Make variable | Use |
| --- | --- |
| `CXX` | MPI C++ compiler wrapper; defaults to `mpicxx` |
| `STACK_PREFIX` | Adds the prefix's include, Eigen and library directories |
| `EIGEN_INCLUDE` | Directory containing `Eigen/Dense` |

The build targets are `all` and `clean`. After changing configuration headers or
the toolchain, run `make clean` and rebuild. Do not compile or clean the same
source tree concurrently with another build or a running job.

## Configure

```text
Compiled defaults in config/*.hpp (and Lite Initialize.hpp)
                         ↓ overridden only by explicit INI entries
                      --config PATH
```

From the selected source directory:

```bash
cp pipeline.example.ini pipeline.ini
```

The template is an annotated catalog of settings. Its commented blank entries
retain the compiled defaults; **it does not disable phases or replace the
original site's data paths automatically**. Review paths, datasets, detector
geometry, source-catalog schema and enabled phases before running.

| INI section | What to configure |
| --- | --- |
| `[process]` | Top-level phase switches, exposure lists and downstream outputs |
| `[astrocat]` | Raw Gaia input, tile output and header/replacement policy |
| `[extcat]` | Raw source catalogs, column schema, projection and tile output |
| `[init]` | Science/DQ roots, dataset names, archive filters and rerun policy |
| `[lensing]` | Main-stage selection, surviving scientific branches, catalogs and geometry |

### Omitted is not the same as empty

```ini
[lensing]
# pixel_size =                  # Keep the compiled value.
pixel_size = 0.2628              # Explicitly override it.
```

An omitted or commented key preserves its compiled value. An active empty string
clears a string setting; an empty numeric or boolean value is an error. Empty
lists are accepted only where that setting allows them. Use quotes around paths
containing `#` or `;`; those characters otherwise start an inline comment.
Booleans accept `true/false`, `1/0`, `yes/no` and `on/off`.

For native runs, put paths directly in your INI and keep it wherever convenient.
The native program does not require JSON, env files or the container launcher.

### Settings that require rebuilding

The parameter reference marks all settings with an INI override. Change those
in the INI for the next run. For all other settings, edit the relevant header and
rebuild. Change source values rather than derived sizes, array lengths or catalog
indices. In Lite, edit `Initialize.hpp` when the corresponding header value is an
alias of `Initialize::*`.

## Run

```bash
mpirun -np 4 ./Fourier_Quad_Pipe --config /absolute/path/pipeline.ini
```

Use at least two MPI processes for the normal parallel workflow. On an HPC
system, request resources in your own job script and use the site's supported
MPI launcher. Native execution does not require Apptainer.

| Invocation | Result |
| --- | --- |
| `Fourier_Quad_Pipe` | Uses compiled defaults only; no automatic INI lookup |
| `Fourier_Quad_Pipe --config file.ini` | Applies this INI over compiled defaults |
| `Fourier_Quad_Pipe --config=file.ini` | Equivalent explicit INI selection |

All other application options, including the former `--help` and positional
exposure-list argument, have been removed. MPI launch options such as
`mpirun -np` belong to the MPI launcher and are unchanged. Supply a single
`--config`; repeated or incomplete selections are errors.

## Phase order

`process_astrocat` and `process_extcat` run once. Each enabled dataset phase then
runs for each dataset in sequence; a collective failure stops subsequent work.

| Order | `[process]` key | Work |
| --- | --- | --- |
| 1 | `run_process_astrocat` | Publish one-degree Gaia tiles |
| 2 | `run_process_extcat` | Publish canonical external source tiles |
| 3 | `run_process_init` | Extract archives and publish exposure lists |
| 4 | `run_process_main` | Run the numerical shear pipeline |
| 5 | `run_process_rearr` | Spatially partition exposure catalogs |
| 6 | `run_process_fd` | Measure shear in field-distortion bins |

At least one phase must be enabled. An initializer run hands its generated
absolute `expo_<target>.list` to the later dataset phases automatically.

### Main numerical stages

Set `[lensing].process_stage` to the product of the primes for the stages to run.
The full product, `223092870`, enables all nine stages. A partial run requires
valid products from its prerequisite stages; **Stage 9 requires Stage 8**.

| Stage | Prime | Work |
| --- | --- | --- |
| 1 | 2 | Background/noise preprocessing and Gaia matching |
| 2 | 3 | Astrometric solution |
| 3 | 5 | Source detection, deblending and star candidates |
| 4 | 7 | Star-candidate power spectra |
| 5 | 11 | PSF selection and modeling |
| 6 | 13 | Galaxy power spectra |
| 7 | 17 | Fourier_Quad estimators and morphology |
| 8 | 19 | Exposure diagnostics |
| 9 | 23 | Catalog combination and calibration |

## Run recipes

Each block below is an INI fragment to merge into your selected variant's
configuration, not an additional command-line option list. Keep the required
paths and scientific settings for your own data.

### Initialize archives only

```ini
[process]
run_process_astrocat = false
run_process_extcat = false
run_process_init = true
run_process_main = false
run_process_rearr = false
run_process_fd = false

[init]
science_root = /data/archive/science
dq_root = /data/archive/dq
output_root = /data/work
datasets = g2019:c4d_19
contains = v1
existing = fail
```

Choose `existing = resume` to validate and reuse compatible outputs, or
`overwrite` to replace outputs intentionally. Neither policy makes an old output
schema compatible with a different program version.

### Main processing only

```ini
[process]
run_process_astrocat = false
run_process_extcat = false
run_process_init = false
run_process_main = true
run_process_rearr = false
run_process_fd = false
expo_list = /data/work/expo_g2019.list

[init]
datasets = g2019:c4d_19
```

Configure Gaia/source catalogs and detector geometry as described in
[data preparation](DATA_PREPARATION.md). Existing list contents must point to
readable image paths with the established dataset layout.

### Prepare Gaia tiles only

```ini
[process]
run_process_astrocat = true
run_process_extcat = false
run_process_init = false
run_process_main = false
run_process_rearr = false
run_process_fd = false

[astrocat]
input_directory = /data/gaia/raw
output_directory = /data/gaia/tiles
add_header = true
existing_policy = fail
```

Publication does not change the main consumer's Gaia directory or filename
layout. Match both when subsequently running the main phase.

### Prepare external source tiles only

```ini
[process]
run_process_astrocat = false
run_process_extcat = true
run_process_init = false
run_process_main = false
run_process_rearr = false
run_process_fd = false

[extcat]
input_directory = /data/sources/raw
output_directory = /data/sources/tiles
```

Set the raw schema, header/delimiter handling, optional projection and magnitude
columns in `[extcat]`. Lite's later main run also needs
`[lensing].source_cat` to point to the generated tiles.

Launch any of these configurations with the same command:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe --config /absolute/path/pipeline.ini
```

## Outputs

Initialization reads compressed archives in place and publishes:

```text
output_root/
├── expo_<target>.list
├── fits_<target>.list
└── <target>/
    ├── science/
    ├── dqmask/
    ├── expolists/
    ├── stamps/
    ├── astrometry/
    └── result/
```

An initializer manifest accompanies the published products. Per-exposure chip
lists are under `<target>/expolists/`; numerical products are under `stamps/`.

| Product | Purpose |
| --- | --- |
| `<dataset>/result/<exposure>_all.cat` | Per-exposure shear catalog |
| `<dataset>/baked/subcat_*.cat` | Default rearranged catalogs, spatially grouped |
| `<dataset>/fdout/FD_test_comb.dat` | Default field-distortion measurement table |

The default full catalog has **49 columns**: 18 external source fields,
`EXPO_NUM`, `ccD_NUM`, 28 Stage-7 fields and exposure `Chi2`. Its width changes
with the configured external prefix. Standard's no-external-catalog branch has
31 columns. Regenerate downstream products when their schema is incompatible;
do not mix earlier 48/30-column products with the current layout.

Absolute rearrangement/FD output directories override their base. Relative
output names are resolved under the configured base or the dataset root. The
parameter reference documents the exposure-list and output fallback rules.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| MPI compiler or library not found | `CXX`, `STACK_PREFIX`, `EIGEN_INCLUDE` and the loaded site environment |
| Unknown option | Use only `--config`; move scientific values into the INI |
| Unknown INI key | Key spelling and whether the selected variant supports it |
| A commented setting appears to have no effect | Comments preserve compiled defaults; uncomment and supply a valid value |
| Gaia/source tiles cannot be found | Directory, filename layout, prefix and sky coverage |
| Wrong catalog fields | One-based raw columns, projection flags, retained magnitude bands and complete row width |
| Rerun refuses existing products | Select `fail`, `resume` or `overwrite` intentionally; keep schemas compatible |
| Stage 9 rejected | Include Stage 8, with factor 19, in `process_stage` |
| Outputs appear in an unexpected dataset | Inspect paths inside the exposure/chip lists, not just `output_root` |

Keep the exact source version, header changes, INI, software versions and job
resource settings with each analysis.
