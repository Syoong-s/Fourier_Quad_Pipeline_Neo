# 🛠️ C++ Pipeline User Guide

[🏠 Home](README.md) · [🗂️ Data Preparation](DATA_PREPARATION.md) · [📋 Parameter Reference](CPP_PIPELINE_PARAMETERS.md) · [🐳 Docker](DOCKER.md) · [🏔️ HPC](HPC.md)

This guide covers **native Linux and HPC deployment** of `Fourier_Quad_Pipe`, including environment prerequisites, variant selection, INI configuration, compilation, MPI execution, and output inspection.

---

## Choose a Variant

Both source trees compile the same executable (`Fourier_Quad_Pipe`) and share the same INI configuration structure, but target different operational needs:

| Feature / Branch | `cpp_Standard/` (Research & Full) | `cpp_Lite/` (Production) |
| :--- | :--- | :--- |
| **Executable** | `Fourier_Quad_Pipe` | `Fourier_Quad_Pipe` |
| **Primary Header Entry** | `config/*.hpp` directly | `Initialize.hpp` + `config/*.hpp` |
| **Astrometry Reference** | Gaia (Type 1 or Type 2) or identity mapping | Gaia (fixed Type 2 one-degree tiles) |
| **Super-Flat Correction** | Optional (`[lensing].include_flat`) | Disabled (physically removed) |
| **DQ Mask Handling** | Configurable (`[lensing].include_mask = 0..3`) | Fixed per-chip DQ masks |
| **Source Selection** | External catalog or standalone detection | External catalog required |
| **PSF Input** | Exposure stars or external PSF images | Exposure stars only |
| **Source Deblending** | Configurable (`[lensing].deblending`) | Always enabled |
| **PSF Spatial Model** | Local polynomial or hybrid model | Local polynomial only |
| **PCA / Multi-Scale PSF** | Optional (`[lensing].psf_ms = 1`) | Removed |
| **Default Enabled Phases** | `init`, `main`, `rearr`, `fd` | `init`, `main` |

> [!IMPORTANT]
> Branches physically removed from `cpp_Lite/` cannot be restored at runtime; supplying Standard-only INI keys to the Lite executable results in a configuration validation error.

---

## Environment and Dependencies

### Required Toolchain & Libraries

| Dependency | Purpose | Reference Stack Version |
| :--- | :--- | :--- |
| **C++17 MPI Wrapper** (`mpicxx`) | Parallel C++17 compilation and MPI dispatch | GCC `12.3.0` + OpenMPI `4.1.8` |
| **CFITSIO** | Reading and writing FITS / `.fits.fz` images | `4.6.4` |
| **FFTW3** (`fftw3` & `fftw3f`) | Double- and single-precision Fast Fourier Transforms | `3.3.11` |
| **Eigen3** | Header-only linear algebra and matrix operations | `3.4.0` |
| **BLAS & LAPACK** | Dense numerical linear algebra and eigensolvers | OpenBLAS / Reference LAPACK |
| **GNU Make** | Build automation (`all` and `clean` targets) | `4.x+` |

> [!TIP]
> If you prefer not to build the C++ scientific stack on your host machine, use the prebuilt toolchain container via [Docker](DOCKER.md) or [HPC Apptainer](HPC.md).

---

## Configuration

The pipeline uses a **two-layer configuration architecture** that separates compile-time constants from runtime overrides:

```mermaid
flowchart LR
    H["📄 Compiled Headers<br/>config/*.hpp & Initialize.hpp"] -->|"Baseline defaults"| R["⚙️ Runtime Config Parser"]
    I["📝 Explicit INI File<br/>--config pipeline.ini"] -->|"Uncommented keys override"| R
    R -->|"MPI Broadcast (Rank 0 → All)"| E["🚀 Fourier_Quad_Pipe Execution"]
```

### Setting Up `pipeline.ini`

Copy the annotated template inside your chosen variant directory:

```bash
cd cpp_Standard                 # Or cpp_Lite
cp pipeline.example.ini pipeline.ini
```

> [!WARNING]
> `pipeline.example.ini` contains **commented blank entries** (`# key =`), not a ready-to-run survey setup. Commented keys retain the compiled header defaults (including original site-specific paths and default enabled phases). Always uncomment and populate your paths, datasets, and active phases before launching.

| INI Section | Scope | Key Responsibilities |
| :--- | :--- | :--- |
| `[process]` | Workflow Control | Top-level phase switches (`run_process_*`), explicit `expo_list`, and `rearr`/`fd` output paths. |
| `[astrocat]` | Gaia Tiling | Raw Gaia directory, output tile directory, header handling, and overwrite policy. |
| `[extcat]` | Source Tiling | Raw external catalog directory, 1-based column indices, magnitude mapping, and chunk size. |
| `[init]` | Archive Extraction | Science/DQ archive roots, `output_root`, `datasets` (`TARGET:PREFIX`), `contains` filters, and rerun policy. |
| `[lensing]` | Shear Pipeline | Prime stage selector (`process_stage`), camera geometry (`pixel_size`, `chipnx`, `chipny`), and catalog paths. |

### INI Syntax Rules: Omitted vs. Empty

```ini
[lensing]
# pixel_size =                  # Commented: retains compiled header value (0.2628).
pixel_size = 0.2628             # Active: explicitly overrides runtime value.
```

| INI Entry State | Runtime Behavior |
| :--- | :--- |
| **Omitted or commented (`# key =`)** | Retains the compiled default from `config/*.hpp` (or `Initialize.hpp`). |
| **Active empty string (`key =`)** | Clears a string setting (enabling downstream path derivation where supported). |
| **Active empty number or boolean** | Rejected with a configuration error (does not fall back to defaults). |
| **Active empty list** | Clears token/dataset filters; column projection lists must remain non-empty. |

- **Quoting & Comments:** Wrap paths containing `#`, `;`, or commas in quotes (`"..."`); unquoted `#` or `;` starts an inline comment.
- **Booleans:** Accept `true`/`false`, `1`/`0`, `yes`/`no`, and `on`/`off` (case-insensitive).
- **Compile-Time Only Parameters:** Settings marked `INI override: No` in the [Parameter Reference](CPP_PIPELINE_PARAMETERS.md) must be edited in `config/*.hpp` (or Lite's `Initialize.hpp`) followed by `make clean && make -j4`.

---

## Build and Run

### 1. Compile the Executable

Build directly with `make` when libraries are on standard system paths, or pass path overrides:

```bash
cd cpp_Standard                 # Or cpp_Lite

# Standard system paths
make -j4

# Custom scientific stack prefix
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack -j4

# Separate Eigen3 include path
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack \
     EIGEN_INCLUDE=/opt/eigen/include/eigen3 -j4
```

| Make Variable | Default | Description |
| :--- | :--- | :--- |
| `CXX` | `mpicxx` | MPI C++17 compiler wrapper. |
| `STACK_PREFIX` | *(empty)* | Adds `<prefix>/include`, `<prefix>/include/eigen3`, and `<prefix>/lib`. |
| `EIGEN_INCLUDE` | *(empty)* | Explicit directory containing `Eigen/Dense`. |

> [!CAUTION]
> Always run `make clean` before rebuilding after modifying any `.hpp` header or switching compilers. Never run concurrent builds or clean operations in the same source directory while a job is executing.

### 2. Launch with MPI

Use at least **2 MPI ranks** for normal distributed processing:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe --config /absolute/path/to/pipeline.ini
```

| Command-Line Invocation | Execution Behavior |
| :--- | :--- |
| `./Fourier_Quad_Pipe` | Uses compiled header defaults only (no automatic `pipeline.ini` lookup). |
| `./Fourier_Quad_Pipe --config file.ini` | Applies `file.ini` overrides on top of compiled defaults. |
| `./Fourier_Quad_Pipe --config=file.ini` | Equivalent single-argument syntax. |

*(All other application CLI flags are rejected; MPI launcher flags such as `mpirun -np` belong to the MPI runtime.)*

---

## Workflow Phases & Main Stages

### Top-Level Phase Execution Order

Catalog preparation phases (`astrocat`, `extcat`) execute once per run. Enabled dataset phases (`init`, `main`, `rearr`, `fd`) then execute sequentially for each configured `TARGET:PREFIX` dataset:

| Order | Phase Name | `[process]` Switch | Input | Primary Output |
| :---: | :--- | :--- | :--- | :--- |
| **1** | `process_astrocat` | `run_process_astrocat` | Raw 2-column Gaia files | Canonical $1^\circ \times 1^\circ$ Gaia tiles (`astra_*.dat`) |
| **2** | `process_extcat` | `run_process_extcat` | Raw external source catalogs | Canonical $1^\circ \times 1^\circ$ source tiles (`extern_*.dat`) |
| **3** | `process_init` | `run_process_init` | Science & DQ `.fits.fz` archives | Extracted CCD FITS tree + `expo_<target>.list` |
| **4** | `process_main` | `run_process_main` | Extracted CCDs + reference tiles | Per-exposure shear catalogs (`<exposure>_all.cat`) |
| **5** | `process_rearr` | `run_process_rearr` | Per-exposure shear catalogs | Spatially partitioned subcatalogs (`baked/subcat_*.cat`) |
| **6** | `process_fd` | `run_process_fd` | Spatial subcatalogs + metadata | Field-Distortion shear recovery (`fdout/FD_test_comb.dat`) |

### 9-Stage Scientific Pipeline (`process_main`)

Inside `process_main`, exposures are dynamically load-balanced across MPI ranks while chips within each exposure are processed sequentially. Active stages are controlled via **prime factor decomposition** in `[lensing].process_stage`:

```mermaid
flowchart TD
    S1["Stage 1 (Prime 2)<br/>Background & Noise Normalization + Gaia Match"] --> S2["Stage 2 (Prime 3)<br/>PU Polynomial Astrometric Solution"]
    S2 --> S3["Stage 3 (Prime 5)<br/>SNR Source Detection, Deblending & Star Candidates"]
    S3 --> S4["Stage 4 (Prime 7)<br/>FFT Power Spectra of Star Candidates"]
    S4 --> S5["Stage 5 (Prime 11)<br/>PSF Star Selection (χ²) & Spatial Polynomial / PCA Model"]
    S5 --> S6["Stage 6 (Prime 13)<br/>FFT Power Spectra & Noise Subtraction of Galaxies"]
    S6 --> S7["Stage 7 (Prime 17)<br/>Fourier_Quad Deconvolution & 5 Shear Estimators (g1, g2, de, h1, h2)"]
    S7 --> S8["Stage 8 (Prime 19)<br/>Exposure-Level Diagnostics Aggregation"]
    S8 --> S9["Stage 9 (Prime 23)<br/>Catalog Assembly, Quality Cuts & Distortion Calibration"]
```

- **Full 9-Stage Product:** `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 = 223092870` (default).
- **Partial Runs:** Multiply only the primes of the stages you wish to execute (prerequisite stage products must already exist on disk; **Stage 9 (`23`) strictly requires Stage 8 (`19`)**).

---

## Common Run Recipes

Merge the relevant snippet below into your `pipeline.ini` and launch with `mpirun -np 4 ./Fourier_Quad_Pipe --config /path/to/pipeline.ini`.

<details>
<summary><b>Recipe 1: Prepare Gaia & External Source Tiles Only</b></summary>

```ini
[process]
run_process_astrocat = true
run_process_extcat   = true
run_process_init     = false
run_process_main     = false
run_process_rearr    = false
run_process_fd       = false

[astrocat]
input_directory  = /data/gaia/raw
output_directory = /data/gaia/tiles
add_header       = true
existing_policy  = fail

[extcat]
input_directory  = /data/sources/raw
output_directory = /data/sources/tiles
```

> Match `[lensing].astrometry_cat` (and `[lensing].astrometry_cat_type = 2` in Standard) and `[lensing].source_cat` to these output directories in subsequent measurement runs.
</details>

<details>
<summary><b>Recipe 2: Initialize Science/DQ Archives Only</b></summary>

```ini
[process]
run_process_astrocat = false
run_process_extcat   = false
run_process_init     = true
run_process_main     = false
run_process_rearr    = false
run_process_fd       = false

[init]
science_root = /data/archive/science
dq_root      = /data/archive/dq
output_root  = /data/work
datasets     = g2019:c4d_19
contains     = v1
existing     = fail             # fail | resume | overwrite
```
</details>

<details>
<summary><b>Recipe 3: Run Main Shear Measurement on an Existing Exposure List</b></summary>

```ini
[process]
run_process_astrocat = false
run_process_extcat   = false
run_process_init     = false
run_process_main     = true
run_process_rearr    = false
run_process_fd       = false
expo_list            = /data/work/expo_g2019.list

[init]
datasets = g2019:c4d_19

[lensing]
astrometry_cat = /data/gaia/tiles
source_cat     = /data/sources/tiles
process_stage  = 223092870
```
</details>

---

## Output Directory Structure & Products

```text
<output_root>/
├── expo_<target>.list                  # Top-level exposure list (per-exposure chip lists)
├── fits_<target>.list                  # Extracted Science FITS manifest
└── <target>/                           # Dataset root (e.g., g2019/)
    ├── science/                        # Extracted per-CCD Science FITS
    ├── dqmask/                         # Extracted per-CCD DQ mask FITS
    ├── expolists/                      # Per-exposure chip lists
    ├── astrometry/                     # WCS headers and astrometric residuals
    ├── stamps/                         # Intermediate stamps, power spectra, PSF models
    ├── result/
    │   └── <exposure>_all.cat          # Per-exposure combined shear catalog (49 cols default)
    ├── baked/
    │   └── subcat_*.cat                # Spatially partitioned subcatalogs (process_rearr)
    └── fdout/
        └── FD_test_comb.dat            # Field-distortion calibration table (process_fd)
```

- **Catalog Width:** By default (`ext_cat = 1` with 18 external columns), `<exposure>_all.cat` has **49 columns**: `18` external source fields + `EXPO_NUM` + `ccD_NUM` + `28` Stage-7 shear/morphology fields + exposure `Chi2`. Standard's `ext_cat = 0` branch outputs **31 columns**.

---

## Troubleshooting

| Symptom / Error | Cause & Resolution |
| :--- | :--- |
| **MPI compiler or library not found** | Check `CXX`, `STACK_PREFIX`, `EIGEN_INCLUDE`, and loaded `module` environment. |
| **Unknown CLI option** | Only `--config PATH` is accepted; move all phase/path options into `pipeline.ini`. |
| **Unknown INI key** | Verify key spelling and confirm the key is not a Standard-only branch used with `cpp_Lite`. |
| **Commented INI entry has no effect** | `# key =` retains compiled defaults; remove `#` and provide an explicit value. |
| **Missing Gaia or source tiles** | Check directory path, Type-1 vs. Type-2 naming (`astrometry_cat_type`), prefix, and sky footprint. |
| **Stage 9 (`23`) rejected** | Include Stage 8 (`19`) in `[lensing].process_stage`. |
| **Rerun refuses existing outputs** | Set `[init].existing` (or catalog `existing_policy`) to `resume` or `overwrite`. |
| **Outputs written to unexpected path** | Downstream phases derive dataset roots from the FITS paths inside `expo_<target>.list`, not solely from `[init].output_root`. |
