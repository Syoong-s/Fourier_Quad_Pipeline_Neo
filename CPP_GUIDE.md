# Fourier_Quad C++ Pipeline Guide

Comprehensive guide for the C++17 (`cpp_Standard` / `cpp_Lite`) pipeline: source
structure, pipeline stages, configuration, external catalog, building, run modes,
initializer output layout, Docker environment, and HPC runner. For the project
overview and the Fortran pipeline see [`README.md`](README.md) and
[`F77_GUIDE.md`](F77_GUIDE.md). The full parameter reference is
[`CPP_PIPELINE_PARAMETERS.md`](CPP_PIPELINE_PARAMETERS.md).

`cpp_Standard` is the full build (includes PCA `PSFRecons`); `cpp_Lite` is the
frozen-branch simplified build with `PSFRecons` removed. See
`cpp_Lite/REFACTOR_NOTES.md` for the Lite change log.

> **中文文档**：请参阅 [CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)

---

## Source Structure

### Source directories

#### `cpp_Standard/` — Full C++17 pipeline

| File | Description |
|---|---|
| `main.cpp` | MPI entry point, workflow option parsing, and four-phase ordering. |
| `include/ProcessConfig.hpp` | Shared workflow defaults and phase switches. |
| `src/process_init/`, `include/process_init/` | Archive initializer implementation and headers. |
| `src/process_main/process_main.cpp`, `include/process_main/process_main.hpp` | Exposure-list loading and Stage 1–9 orchestration. |
| `src/process_rearr/`, `include/process_rearr/` | Self-contained `_all.cat` sky partitioning, MPI redistribution, sorted subcatalogs, and summary output. |
| `include/process_rearr/ProcessRearrConfig.hpp` | Rearrangement-only parameters and derived `external columns + 1 + ichi2` row width. |
| `include/process_main/LensingConfig.hpp` | Configuration constants (equivalent to `para.inc` + `cust_para.inc` + `sig_para.inc`). |
| `src/process_main/PreProcess.cpp`, `include/process_main/PreProcess.hpp` | **Stage 1**: pre-processing. |
| `src/process_main/Astrometry.cpp`, `include/process_main/Astrometry.hpp` | **Stage 2**: astrometric calibration. |
| `src/process_main/SourceExtractor.cpp`, `include/process_main/SourceExtractor.hpp` | **Stage 3**: source detection and extraction. |
| `src/process_main/FourierTransformSt1.cpp`, `include/process_main/FourierTransformSt1.hpp` | **Stage 4**: first-stage Fourier transform. |
| `src/process_main/PSFModel.cpp`, `include/process_main/PSFModel.hpp` | **Stage 5**: PSF modeling. |
| `src/process_main/PSFRecons.cpp`, `include/process_main/PSFRecons.hpp` | PSF PCA reconstruction (`PSF_Ms=1` only). |
| `src/process_main/FourierTransformSt2.cpp`, `include/process_main/FourierTransformSt2.hpp` | **Stage 6**: second-stage Fourier transform. |
| `src/process_main/ShearMeasurement.cpp`, `include/process_main/ShearMeasurement.hpp` | **Stage 7**: Fourier\_Quad shear estimation. |
| `src/process_main/ExposureInfo.cpp`, `include/process_main/ExposureInfo.hpp` | **Stage 8**: per-exposure statistics. |
| `src/process_main/CatalogCombiner.cpp`, `include/process_main/CatalogCombiner.hpp` | **Stage 9**: catalog combination and calibration. |
| `src/process_main/` and `include/process_main/` support modules | FITS I/O, linear algebra, image processing, MPI scheduling, and shared numerical utilities. |
| `src/process_main/NumericalRecipes.cpp`, `include/process_main/NumericalRecipes.hpp` | Numerical Recipes port (RNG, sorting, interpolation). |
| `src/process_main/MPIScheduler.cpp`, `include/process_main/MPIScheduler.hpp` | MPI initialization and task distribution. |
| `src/process_main/ExStar.cpp`, `include/process_main/ExStar.hpp` | Star extraction and classification. |
| `Makefile` | Build file. Uses `mpicxx`, C++17, links against CFITSIO, FFTW, LAPACK. |

#### `cpp_Lite/` — Simplified C++17 pipeline

Uses the same integrated `process_extcat` / `process_init` / `process_main` /
`process_rearr` directory layout and runtime option contract as
`cpp_Standard/`, but its scientific modules retain
the frozen Lite branches and `PSFRecons.cpp/.hpp` is absent. See
`cpp_Lite/REFACTOR_NOTES.md` for the detailed change log.



## Source layout

- `include/ProcessConfig.hpp`: workflow defaults and default phase switches.
- `include/process_extcat/`, `src/process_extcat/`: external-catalog schema,
  parsing, MPI byte-range partitioning, and deterministic tile publication.
- `include/process_init/`, `src/process_init/`: initializer wrapper plus the
  preserved `Initializer` and `FitsExtractor` modules.
- `include/process_main/`, `src/process_main/`: `LensingConfig`, all numerical
  modules, exposure-list loading, and the complete Stage 1–9 orchestration.
- `include/process_rearr/`, `src/process_rearr/`: self-contained `_all.cat`
  schema validation, full-sky weighted k-d partitioning, MPI redistribution,
  sorted subcatalog publication, and summary output.
- Each `cpp_Standard` / `cpp_Lite` root contains only the executable entry point, build file,
  documentation, and phase implementation trees.


## Pipeline Stages

The pipeline consists of 9 stages. Stage execution is controlled by the
`PROCESS_stage` parameter (defined in `para.inc` / `LensingConfig.hpp`), which is
the product of prime factors. A stage runs when `PROCESS_stage` is divisible by
its prime. The default value `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` enables
all stages.

| Stage | Prime | Function | Description |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | Read FITS images, apply flat-field and mask corrections, estimate background noise (F6 mode-bar estimator). |
| 2 | 3 | `proc_astrometry` / `Astrometry` | Astrometric calibration using Gaia reference catalog; WCS fitting. |
| 3 | 5 | `proc_source` / `SourceExtractor` | Source detection, deblending, and stamp extraction. |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | First-stage Fourier transform of galaxy stamps. |
| 5 | 11 | `proc_PSF` / `PSFModel` | PSF modeling from stellar stamps via local polynomial fitting. Optional PCA reconstruction (`PSF_Ms=1`). |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | Second-stage Fourier transform. |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | Fourier\_Quad shear estimation from fourth-order Fourier moments. |
| 8 | 19 | `proc_info` / `ExposureInfo` | Collect per-exposure statistics (PSF FWHM, star count, etc.). |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | Combine shear catalogs across exposures and apply calibration corrections. |

To disable a stage, divide `PROCESS_stage` by its prime factor. For example,
setting `PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19` (omit 23) skips
catalog combination.

---



## Configuration

### C++ (`LensingConfig.hpp`)

`cpp_Standard` consolidates all parameters from the three Fortran include files
into `include/process_main/LensingConfig.hpp`. `cpp_Lite` uses the same relative
path but retains only its frozen parameter subset; see
`cpp_Lite/REFACTOR_NOTES.md`. Catalog column indices are shifted to 0-based for
C++. See the complete
[`C++ Pipeline External Inputs and Parameter Reference`](CPP_PIPELINE_PARAMETERS.md)
for every runtime option, external input, `ProcessConfig.hpp` default,
`LensingConfig.hpp` parameter, valid value, and Standard/Lite difference.

---



## External Source Catalog

When `ext_cat = 1`, the C++ and Fortran pipelines read one-degree DES Y6 GOLD
tiles from the directory configured as `SOURCE_CAT`. The
[`gen_src_cat`](gen_src_cat/README.md) utilities download or repartition those
tiles with the filename convention expected by the catalog readers. The Python
downloader writes the DES Y6 GOLD 18-column schema. The C++ MPI repartitioner
preserves every raw column by default or selects any ordered list when explicit
projection is enabled; it is also integrated into `cpp_Standard` and `cpp_Lite`
as `process_extcat`, the optional first phase before `process_init` and
`process_main`.

The generated exposure `_all.cat` files can then be repartitioned by celestial
region with the self-contained fourth phase `process_rearr`. Its pass-through
width uses `EXTCAT_TOTAL_COLUMNS` (default 18); its dedicated config derives the
complete default width as `18 + 1 + ichi2(25) = 44`. RA/Dec remain configured
raw one-based positions when explicit projection is disabled and are converted
automatically to projection positions when it is enabled. Outputs default to
each dataset's `rearranged_catalog/` directory.

For C++ external catalogs, set the one-based raw positions
`EXTCAT_RA_COLUMN_ONE_BASED`, `EXTCAT_DEC_COLUMN_ONE_BASED`, and
`EXTCAT_ZP_COLUMN_ONE_BASED` in the selected `ProcessConfig.hpp`. Their defaults
are `5`, `6`, and `17`, matching DES Y6 GOLD `ra`, `dec`, and `dnf_z`.
`process_main` converts only these three fields; other external catalog columns may be
arbitrary strings and no fixed 18-column layout is required. When explicit
projection is enabled, all three raw fields must be selected and their output
positions are derived automatically from the projection order. Runtime jobs may
override the positions with `--extcat-ra-column`, `--extcat-dec-column`, and
`--extcat-zp-column`. `process_rearr` uses the same RA/Dec mapping rule but
requires every complete `_all.cat` field to be finite numeric data.

```bash
cd gen_src_cat
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pyvo
python query_y6gold_sync_mp_v2.py
```

Review the sky bounds, row limit, output directory, and query concurrency in the
script before running it. For C++, set the primary `SOURCE_CAT` path in
`LensingConfig.hpp`, or pass `--extcat-output`; `EXTCAT_OUTPUT_DIRECTORY`
follows that value so generation and numerical processing use the same path.
Raw local catalogs can instead be tiled inside the pipeline with
`--run-extcat true --extcat-input PATH`. Fortran still uses `SOURCE_CAT` in
`para.inc`. See
[`gen_src_cat/README.md`](gen_src_cat/README.md) for the Python schema, C++
projection modes, and catalog-generation behavior.

---



## Compiler and libraries

Local verification used GCC/G++ 15.2.0, Open MPI 5.0.10, CFITSIO 4.6.3, FFTW
3.3.10, Eigen 3.4.0, and OpenBLAS/LAPACK 0.3.33 on Linux in WSL2.

The existing portable HPC target is GCC/G++ 12.3.0, Open MPI 4.1.8, CFITSIO
4.6.4, FFTW 3.3.11, Eigen 3.4.0, and LAPACK 3.11.0. A cluster may use equivalent
site modules as long as one matching MPI C++ wrapper compiles and launches the
program.


## Build

When all headers and libraries are visible through the compiler's normal search
paths:

```bash
make -j4
```

For one consolidated scientific-stack prefix containing headers, Eigen under
`include/eigen3`, and libraries under `lib`:

```bash
make STACK_PREFIX="${STACK_PREFIX}" -j4
```

When Eigen and the linked libraries use different prefixes, pass portable local
overrides explicitly:

```bash
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" -j4
```

No Windows-native compiler or wrapper is required. Cluster builds use the same
Makefile after loading the site's compiler, MPI, and scientific-library modules.

Run the focused external-catalog column reader test with:

```bash
make test-extcat-reader
```

Run the self-contained rearrangement unit and MPI integration tests with:

```bash
make test-rearr
```


## Defaults and option syntax

Edit `include/ProcessConfig.hpp` to set the normal dataset and default execution
mode. `RUN_PROCESS_EXTCAT`, `RUN_PROCESS_INIT`, and `RUN_PROCESS_REARR` default
to `false`; `RUN_PROCESS_MAIN` defaults to `true`. Every command-line option is
optional and overrides its configured default. Both `--name value` and
`--name=value` are accepted in any order.

The `EXTCAT_*` values configure raw-catalog discovery, the output directory,
parsing policies, MPI task size, optional ordered column selection, and raw
RA/Dec/ZP columns. With explicit selection disabled, every raw input field is
preserved in place. `process_main` reads only these three configured fields;
unselected catalog fields need not be numeric. `EXTCAT_TOTAL_COLUMNS` gives
`process_rearr` the pass-through external width; explicit projection instead
uses the projection-list length because that is the emitted external schema.
Set the primary catalog path with `SOURCE_CAT` in
`include/process_main/LensingConfig.hpp`. `EXTCAT_OUTPUT_DIRECTORY` is a
read-only reference to that value, so `process_extcat` writes where
`process_main` reads. `--extcat-output` overrides both for one invocation.

`DATASETS` stores paired target/prefix values, and `CONTAINS` stores the archive
basename tokens accepted with OR semantics. For example:

```cpp
inline const std::vector<DatasetSpec> DATASETS = {
    {"g2013", "c4d_13"},
    {"g2014", "c4d_14"},
    {"g2019", "c4d_19"},
};
inline const std::vector<std::string> CONTAINS = {"v1", "v2"};
```

| Option | Purpose |
|:---|:---|
| `--run-extcat BOOL` | Enable or disable external-catalog repartitioning. |
| `--run-init BOOL` | Enable or disable archive initialization at runtime. |
| `--run-main BOOL` | Enable or disable the numerical pipeline at runtime. |
| `--run-rearr BOOL` | Enable or disable `_all.cat` spatial rearrangement; it follows `process_main` when both run. |
| `--extcat-input PATH` | Directory containing any number of raw text catalogs. |
| `--extcat-output PATH` | Destination tile directory and effective C++ `SOURCE_CAT`. |
| `--extcat-contains TEXT` | Repeatable case-sensitive basename token; repeats use OR. |
| `--extcat-recursive BOOL` | Enable or disable recursive input discovery. |
| `--extcat-delimiter MODE` | `auto`, `whitespace`, `comma`, or `tab`. |
| `--extcat-header MODE` | `auto`, `present`, or `absent`. |
| `--extcat-columns LIST` | One or more one-based raw indices; output fields follow this exact order and width. |
| `--extcat-ra-column N` | Raw one-based RA index; overrides `ra` header discovery. |
| `--extcat-dec-column N` | Raw one-based Dec index; overrides `dec` header discovery. |
| `--extcat-zp-column N` | Raw one-based photometric-redshift (`dnf_z`) index used by `process_main`. |
| `--extcat-chunk-mib N` | Positive MPI byte-range task size in MiB. |
| `--extcat-malformed POLICY` | `fail` or `skip` malformed rows. |
| `--extcat-existing POLICY` | `fail` or transactionally `overwrite` generated tiles. |
| `--science-root PATH` | Original read-only Science `.fits.fz` repository. |
| `--dq-root PATH` | Original read-only DQ `.fits.fz` repository. |
| `--output-root PATH` | Parent of the target directory and generated lists. |
| `--dataset TARGET:PREFIX` | One paired dataset; repeat the option for a batch. |
| `--target NAME` | Legacy single-dataset target; cannot be mixed with `--dataset`. |
| `--prefix TEXT` | Legacy single-dataset prefix; cannot be mixed with `--dataset`. |
| `--contains TEXT` | Accepted basename token; repeat for OR matching. |
| `--existing MODE` | `fail`, `resume`, or `overwrite`; default is `fail`. |
| `--f77-max-path N` | Initializer-only maximum generated path length; `0` disables the check. |
| `--expo-list PATH` | Exposure list used in main/rearr-only mode. |
| `--help` | Print the effective command contract. |

Boolean values accept `true`/`false`, `1`/`0`, and `on`/`off`. One legacy
positional exposure-list path is retained as a compatibility alias, but new jobs
should use `--expo-list`. The first explicit `--dataset` replaces configured
`DATASETS`, and subsequent occurrences append. `--contains` follows the same
replacement/append rule for `CONTAINS`. Other duplicate scalar options use the
last value. The first explicit `--extcat-contains` similarly replaces
`EXTCAT_FILENAME_TOKENS`; later occurrences append. Dataset target names must be
unique within one invocation.

`--f77-max-path` belongs to `process_init`: its default value of 150 protects
paths intended to remain compatible with the legacy Fortran workflow. It does
not truncate or reject paths in `process_main`. Main-process paths are
`std::string` values and are instead subject to the selected filesystem and I/O
library limits (including CFITSIO's filename capacity for FITS products).


## Run modes

External-catalog-only execution accepts any number of matching input files and
does not require a configured dataset:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --extcat-input /data/raw_catalogs \
  --extcat-output /data/catalogs/des_y6_chunks \
  --extcat-contains .csv --extcat-contains y6_gold
```

`process_extcat` runs collectively once before the dataset loop. It writes the
same one-degree filenames as `gen_src_cat/query_y6gold_sync_mp_v2.py`. Its
default output preserves the complete raw schema; `--extcat-columns 5,3,4,1`,
for example, writes raw columns 5, 3, 4, and 1 as output columns 1–4. If it
fails, no later phase starts. Output passed onward to `process_main` must
include the raw columns configured by `EXTCAT_RA_COLUMN_ONE_BASED`,
`EXTCAT_DEC_COLUMN_ONE_BASED`, and `EXTCAT_ZP_COLUMN_ONE_BASED`; a
rearrangement-only run requires RA and Dec but does not require ZP. In
pass-through mode the required fields are read at their configured positions.
With explicit projection, the reader automatically maps each raw index to its
position in the ordered projection list.

Main-only local execution:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true \
  --expo-list /data/work/expo_g2019.list
```

Rearrangement-only local execution consumes existing per-exposure `_all.cat`
files referenced by the same exposure list:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main false --run-rearr true \
  --expo-list /data/work/expo_g2019.list
```

All rearrangement-specific parameters are in
`include/process_rearr/ProcessRearrConfig.hpp`. With the default 18-field
external catalog, `ichi2=25` and the complete row width is calculated there as
`18 + 1 + 25 = 44`. The default 0.1-degree grid targets about 500,000 rows per
weighted k-d partition. Outputs are written below each dataset root in
`rearranged_catalog/` as `subcat_NNNNNN.cat` plus
`catalog_summary.txt`. Every data row must have the exact numeric width and
finite values; configured missing catalogs and malformed rows are skipped and
reported.

Initializer-only local execution:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work --dataset g2019:c4d_19
```

Batch initialization pairs every target with its own prefix. Multiple contains
tokens select an archive when any token occurs in its basename:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work \
  --dataset g2013:c4d_13 --dataset g2014:c4d_14 \
  --dataset g2019:c4d_19 \
  --contains v1 --contains v2
```

Chained local execution uses the same initializer options with downstream phase
switches enabled. After successful initialization, `process_main` and
`process_rearr` receive
the normalized absolute `output_root/expo_<target>.list` path returned by
`process_init`. That generated path overrides `--expo-list`, the legacy positional
argument, and every configured exposure-list default.

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main true --run-rearr true \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work --dataset g2019:c4d_19 \
  --existing resume
```

Enable all four switches to build the external catalog first, then initialize,
process, and rearrange every configured dataset. The effective
`--extcat-output` path is also used by the numerical source extractor.

Datasets execute sequentially on the same MPI communicator and stop at the first
failure. In main/rearr-only batch mode, omit `--expo-list`; the driver derives one
`output_root/expo_<target>.list` path per dataset. A single external exposure list
is accepted only for a single downstream-only dataset. In chained batch mode, every
initializer-generated absolute list overrides external exposure-list input for
its corresponding dataset.

On a Slurm cluster, use the same executable arguments with the site launcher,
for example `srun -n 40 ./Fourier_Quad_Pipe ...`. The initializer failure status
is collective; the numerical phase is never entered after initialization fails.


## Initializer output contract

For each `--output-root OUTPUT --dataset TARGET:PREFIX`, initialization builds
the pipeline directory tree below `OUTPUT/TARGET`. The exact order is: create
the fixed type directories idempotently; extract Science/DQ chip images; write
each successful Science exposure's `stamps/<EXPOSURE>.list`; have rank zero
publish the two top-level lists; create chip-product exposure subdirectories
from the published expo list; and finally publish the manifest. Source
`.fits.fz` archives are read in place and are never copied or removed. The
complete fixed directory contract and its chip-product subset are centralized
in `include/OutputLayout.hpp` for both variants.

### Output layout under `OUTPUT/TARGET`

- `science/<EXPOSURE>/<EXPOSURE>_<N>.fits` - Science chip images, sharded one
  subdirectory per exposure; `<N>` is the sequential two-dimensional HDU
  occurrence index (1, 2, 3, ...).
- `dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits` - DQ mask chip images, sharded
  one subdirectory per exposure; `<CCDNUM>` is the FITS `CCDNUM` header value.
- `stamps/` - per-exposure chip lists (`<EXPOSURE>.list`) written during
  extraction, plus all `process_main` intermediate products in type-specific
  subdirectories (below).
- `astrometry/dat_Astro/`, `astrometry/Head/`, `astrometry/dat_Chk/` -
  astrometry solutions (`<P>_astro.dat`), WCS `.head` files, and check data.
- `result/` - final per-exposure products, including `<EXPOSURE>_all.cat`
  consumed by `process_rearr`.

### process_main intermediate products (under `stamps/`)

Type-specific subdirectories replace the former flat `stamps/`, `rescale/`,
`starxy/`, `fits_psfresi/`, `dat_pcs/`, and `dat_starcomp/` directories:

`Norm/`, `cat_Orig/`, `dat_StarInfo/`, `dat_StarCanInfo/`, `dat_SrcInfo/`,
`dat_PsfFit/`, `dat_Shear/`, `dat_ExpoInfo/`, `dat_StarComp/`, `dat_StarCompV2/`,
`dat_Rescale/`, `dat_StarXY/`, `dat_Pcs/`, `fits_StarCan/`, `fits_StarCanN/`,
`fits_StarCanP/`, `fits_StarP/`, `fits_Src/`, `fits_Noise/`, `fits_SrcP/`,
`fits_PsfLocal/`, `fits_PsfSrc/`, `fits_PsfResi/`.

Every chip-scoped product is sharded one level further by exposure:
`<TYPE>/<EXPOSURE>/<CHIP><SUFFIX>`. For example, a normalized chip is written
as `stamps/Norm/<EXPOSURE>/<CHIP>_norm.fits`, and its astrometry solution is
`astrometry/dat_Astro/<EXPOSURE>/<CHIP>_astro.dat`. Exposure-scoped products
such as `.head`, `_star_info_expo.dat`, `_star_power_expo.fits`,
`_PSF_source.fits`, `_expo_info.dat`, and `_all.cat` remain directly in their
existing type directories. Rank zero creates these chip-product
`<EXPOSURE>/` directories idempotently only after publishing and re-reading
`expo_TARGET.list`.

### process_main path and output-failure contract

The Stage 1--9 producer/consumer chain has been audited against this layout.
Chip products are constructed on both write and read paths with
`OutputLayout::chipPath`; exposure products remain directly in their declared
type directories. The checked chain is: astrometry/normalization, WCS/check,
source and star-candidate extraction, star power, PSF products, source power,
shear, exposure information, and final catalog combination. DQ input is read
from the same `dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits` contract published
by initialization. No path-layer or suffix mismatch was found.

All `process_main` text outputs use the checked `MainIO::OutputFile` stream.
FITS outputs use the checked `FitsIO` creation/write/close paths. A create,
write, flush, or close failure emits one `Output creation failed` diagnostic
containing MPI rank, operation, output path, and the OS/CFITSIO reason, then
terminates the MPI job with `MPI_Abort`. This avoids leaving the master or
another worker blocked in the dynamic scheduler.

### Exposure-list generation

During extraction each rank writes `stamps/<EXPOSURE>.list` (the chip image
paths it produced) directly from the extraction result - no post-hoc disk
re-scan is performed. After all ranks finish, rank zero scans `stamps/`, sorts
the per-exposure lists, and atomically publishes the first two files below:

- `OUTPUT/expo_TARGET.list` - top-level list; each line is
  `"<OUTPUT/TARGET/stamps/<EXPOSURE>.list>"  <chip count>`.
- `OUTPUT/fits_TARGET.list` - flat list of every Science chip image path.

Rank zero then reads the published expo list and creates each chip-product
exposure subdirectory. Only after that step succeeds does it atomically publish
`OUTPUT/init_TARGET_manifest.json`. Manifest schema version 2 records the
`exposure_directories_created` completion flag as well as all active basename
filters in the `filename_tokens` array.

Science chips are numbered by two-dimensional HDU occurrence; DQ chips are
numbered by `CCDNUM`. Downstream stages derive the dataset root from a Science
chip path via `getDir(image, 3)` (three levels up:
`science/<EXPOSURE>/<file>` -> `science` -> `OUTPUT/TARGET`), and per-chip DQ
masks are read from `dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits`.


## Docker Environment


The Docker environments provide a reproducible build toolchain without requiring
manual installation of compilers and scientific libraries.


### Docker directories

#### `cpp_docker/`

Builds a portable HPC runtime image:

| Component | Version |
|---|---|
| Base image | Rocky Linux 8.10 |
| G++ | 12.3.0 (conda-forge) |
| OpenMPI | 4.1.8 (Slurm PMI2 direct-launch) |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / OpenBLAS | 3.11.0 / 0.3.33 |

Key files: `Dockerfile`, `compose.yaml`, `pixi.toml`, `.env.example`,
`scripts/verify-image.sh`, `scripts/check-public-repo.sh`, `runner/` (HPC
deployment scripts), `SOURCES.md`, `THIRD_PARTY_NOTICES.md`.

---


### Quick start (C++)

```bash
cd cpp_docker
cp .env.example .env
# Edit .env: set CPP_SOURCE_HOST, catalog paths, and PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-CPP
```

Inside the container:

```bash
cd /workspace/src_pipe
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```


### Verifying an image

Each Docker directory includes a verification script:

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
bash cpp_docker/scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

For detailed Docker environment documentation, see:
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)
- [`cpp_docker/README.md`](cpp_docker/README.md) / [`cpp_docker/README-CN.md`](cpp_docker/README-CN.md)

---



## HPC Deployment

Both Docker environments include `runner/` directories with Slurm/Apptainer
deployment scripts. The typical workflow is:

1. Pull the GHCR image and convert to a SIF:
   ```bash
   bash f77_docker/runner/pull-sif.sh    # or cpp_docker/runner/pull-sif.sh
   ```

2. Configure the environment:
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # Edit paths for your cluster
   ```

3. Audit MPI compatibility (read-only):
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. Submit the pipeline:
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

For detailed HPC runner documentation, see:
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)
- [`cpp_docker/runner/README.md`](cpp_docker/runner/README.md) / [`cpp_docker/runner/README-CN.md`](cpp_docker/runner/README-CN.md)

---
