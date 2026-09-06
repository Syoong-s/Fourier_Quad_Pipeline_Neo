# Fourier_Quad C++ Pipeline Guide

This guide explains how to build, configure, and run `Fourier_Quad_Cpp`.

> 中文版：[CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)

## Program layout

Both [`cpp_Standard`](cpp_Standard/) and [`cpp_Lite`](cpp_Lite/) build
`Fourier_Quad_Pipe`. Each variant contains:

- `main.cpp`: MPI lifecycle and six-phase dispatch;
- `pipeline.example.ini`: complete run-time configuration template;
- `config/`: compiled defaults and fixed numerical settings;
- `include/` and `src/`: phase modules, shared catalog layout, and run-time
  configuration;
- `include/general/` and `src/general/`: exposure-list, path, MPI, scheduler,
  output-layout, and numerical utilities;
- `Makefile`: portable C++17/MPI build.

Standard retains optional scientific branches. Lite fixes and physically
removes eight branches:

| Setting | Lite behavior |
|---|---|
| astrometry | Gaia only (`ASTROMETRY_trivial=0`) |
| flat | disabled (`include_FLAT=0`) |
| mask | per-chip DQ mask (`include_Mask=2`) |
| sources | external catalog (`ext_cat=1`) |
| PSF input | stars in the exposure (`ext_PSF=0`) |
| deblending | enabled |
| PSF model | local polynomial (`PSF_type=1`) |
| PCA/multi-scale | disabled (`PSF_Ms=0`) |

## Processing model

The driver invokes six phases in a fixed order. `process_astrocat` and
`process_extcat` run once, in that order; the other enabled phases run once per
dataset. Datasets are processed sequentially and the first collective failure
stops the run.

| Phase | CLI switch | Purpose |
|---|---|---|
| `process_astrocat` | `--run-astrocat` | Repartition raw two-column Gaia catalogs into deduplicated one-degree tiles. |
| `process_extcat` | `--run-extcat` | Repartition raw External source catalog files into sky tiles. |
| `process_init` | `--run-init` | Discover `.fits.fz` archives, extract Science/DQ chips, and publish exposure lists. |
| `process_main` | `--run-main` | Run the nine-stage numerical shear pipeline. |
| `process_rearr` | `--run-rearr` | Partition `*_all.cat` rows into spatial subcatalogs. |
| `process_fd` | `--run-fd` | Recover mean shear in field-distortion bins. |

`process_main` uses a prime-product stage selector:

| Stage | Prime | Work |
|---:|---:|---|
| 1 | 2 | background/noise preprocessing and Gaia matching |
| 2 | 3 | astrometric solution |
| 3 | 5 | source detection, deblending, and star candidates |
| 4 | 7 | star-candidate power spectra |
| 5 | 11 | PSF selection and modeling |
| 6 | 13 | galaxy power spectra |
| 7 | 17 | Fourier_Quad estimators and morphology measurements |
| 8 | 19 | exposure diagnostics |
| 9 | 23 | catalog combination and calibration |

A stage runs when `process_stage` is divisible by its prime. Stage 9 requires
Stage 8. The full product `223092870` enables all stages.

## Build

Required libraries are CFITSIO, FFTW3 (double and float), Eigen3, LAPACK, and
BLAS. Use an MPI C++ wrapper with C++17 support.

```bash
cd cpp_Standard                 # or cpp_Lite
make -j4
./Fourier_Quad_Pipe --help
```

For libraries installed below one prefix:

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack -j4
```

For a separate Eigen installation:

```bash
make STACK_PREFIX=/opt/science-stack \
     EIGEN_INCLUDE=/opt/eigen/include/eigen3 -j4
```

The live Standard and Lite Makefiles expose only `all` and `clean`. `make` is
equivalent to `make all`; use `make clean` before rebuilding after a compiled
configuration or toolchain change.


## Configure a run

Supports three configuration methods:
1. Configuration file pipeline.ini;
Copy the template belonging to the selected variant:

```bash
cp pipeline.example.ini pipeline.ini
```
2. Command-line arguments;
3. Setting compiler constants in the header files under the config directory.

Configuration is resolved as:

```text
compiled defaults < INI file < CLI options
```

The INI has five sections:

- `[process]`: phase switches and downstream output/list paths;
- `[astrocat]`: raw two-column Gaia input, independent tile output, header mode, and existing-output policy;
- `[extcat]`: raw catalog parsing, schema, and tile publication;
- `[init]`: archive roots, dataset pairs, and existing-output policy;
- `[lensing]`: run-selectable stage, branch, catalog, and camera settings.

Lite rejects Standard-only lensing keys. Unknown sections, keys, malformed
values, and inconsistent stage/schema choices are errors. The INI is applied
transactionally on every rank before any phase runs.

CLI options accept `--name value` and `--name=value`. Booleans accept
`true/false`, `1/0`, `yes/no`, and `on/off`. The first explicit `--dataset`,
`--contains`, or `--extcat-contains` replaces its configured list; later
occurrences append values. A single bare argument remains a legacy alias for
`--expo-list`.

Use `./Fourier_Quad_Pipe --help` for the exhaustive current option list and
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) for INI keys.

For a run-time setting, prefer the selected variant's INI template and use CLI
only for per-run overrides. For a header-only setting, edit the matching file
below `cpp_Standard/config/` or `cpp_Lite/config/`, change source parameters
rather than derived constants, then run `make clean` and rebuild. The parameter
reference has one complete Standard/Lite comparison table for each configuration
domain and marks every value that can avoid a rebuild.

### Centralized path configuration

The selected variant's `config/pathconfig.hpp` is the sole physical definition
point for compiled input/output paths, workflow output/list names,
rearrangement filenames, and initializer/product relative directories. Existing
namespaces and RuntimeConfig keys remain unchanged.

- Keep `AstroCatConfig::ASTROCAT_OUTPUT_DIRECTORY` initialized from
  `LensingConfig::ASTROMETRY_CAT` as two coupled symbols. After RuntimeConfig is
  built, `[astrocat].output_directory` and `[lensing].astrometry_cat` are
  independent producer/consumer fields and must both be set when they should
  name the same newly published tile directory.
- Keep `ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY` as a `const std::string&` to
  `LensingConfig::SOURCE_CAT`. `[lensing].source_cat` and
  `[extcat].output_directory` update the same effective RuntimeConfig directory;
  `--extcat-output` has final precedence.
- `FLAT_PATH` and `PSF_PATH` exist only in Standard because Lite physically
  removed those optional branches.
- `ProcessRearrConfig`'s fixed filenames and the two `OutputLayout` directory
  arrays have no RuntimeConfig/INI/CLI override. Editing them—or any compiled
  default directly in `pathconfig.hpp`—requires `make clean && make`.

INI and CLI values override RuntimeConfig copies only; they do not modify the
header or its compile-time relationships.

### Common and data-source-dependent parameters

Review the following entry points before each run. Put persistent runtime settings in the
INI and use the listed CLI options only for one-off overrides. Only fixed settings marked
as build time require editing the selected variant's header and running
`make clean && make`. Do not edit derived sizes or runtime-resolved catalog offsets alone.

| Category | Parameters (current defaults) | How to change | When to change / constraints |
|---|---|---|---|
| Top-level phases | `[process].run_process_astrocat/run_process_extcat/run_process_init/run_process_main/run_process_rearr/run_process_fd` | INI; runtime `--run-astrocat`, `--run-extcat`, `--run-init`, `--run-main`, `--run-rearr`, `--run-fd` | Select phases for this invocation. Standard defaults to `false/false/true/true/true/true`; Lite defaults to `false/false/true/true/false/false`. |
| Science/DQ archives and datasets | `[init].science_root`, `dq_root`, `output_root`, `datasets`, `contains` | INI; runtime `--science-root`, `--dq-root`, `--output-root`, `--dataset`, `--contains` | Change for another observing archive, basename prefix, filter token, or output root. Lite requires per-CCD DQ masks. |
| Exposure lists and phase outputs | `[process].expo_list`, `rearr_output_directory`, `rearr_output_base_directory`, `rearranged_expo_list_filename`, `rearranged_expo_list_directory`, `fd_expo_list`, `fd_output_directory`, `fd_output_base_directory` | INI; corresponding CLI: `--expo-list`, `--rearr-output-dir`, `--rearr-output-base`, `--rearr-list-name`, `--rearr-list-dir`, `--fd-expo-list`, `--fd-output-dir`, `--fd-output-base` | Change for downstream-only execution or a different rearrangement/FD output or exposure-list location. |
| Fixed generated layout | `SKIP_DIRECTORY_NAME`, `SUBCAT_PREFIX`, `SUBCAT_EXTENSION`, `SUMMARY_FILENAME`, `NON_CHIP_BASE_DIRECTORIES`, `CHIP_PRODUCT_DIRECTORIES` | `config/pathconfig.hpp`, build time | Change only when the published catalog naming or relative output-directory contract changes; rebuild and regenerate affected products. |
| Gaia catalog tiling | `[astrocat].input_directory`, `output_directory`, `add_header=true`, `existing_policy=fail` | INI; runtime `--astrocat-input`, `--astrocat-output`, `--astrocat-add-header`, `--astrocat-existing` | Change for another raw Gaia catalog or rerun policy. Output belongs only to `process_astrocat` and is not propagated to `[lensing].astrometry_cat`. |
| Gaia catalog layout | `[lensing].astrometry_cat_type=1`, `astrometry_cat`; `ASTROMETRY_TILE_PREFIX="astra_"` | layout/path by INI at runtime; prefix in `config/pathconfig.hpp` at build time | `1` reads legacy large `gaia_*.cat` tiles; `2` accumulates one-degree `<prefix>RA_*.dat` tiles from `process_astrocat`. Change the consumer directory with the layout; rebuild after changing the prefix. |
| External-catalog discovery and publication | `[extcat].input_directory`, `output_directory`; `SOURCE_CAT_TILE_PREFIX="extern_"` | paths by INI/CLI at runtime; prefix in `config/pathconfig.hpp` at build time | Change the raw external-catalog directory or normalized tile directory. Output must not equal or sit below input; it is also the effective `SOURCE_CAT`. The prefix is shared by producer and consumer. |
| External-catalog schema | `[extcat].total_columns`, `use_explicit_columns`, `input_columns`, `use_explicit_coordinate_columns`, `ra_column`, `dec_column`, `zp_column` | INI; projection and RA/Dec/ZP can be overridden with `--extcat-columns`, `--extcat-ra-column`, `--extcat-dec-column`, `--extcat-zp-column` | Change for another survey or column order. Explicit projection must retain RA, Dec, photo-z, and fields consumed by enabled phases; `CatalogLayout` resolves complete widths and downstream offsets. |
| FD magnitude mapping | `[extcat].mag_g_column`, `mag_r_column`, `mag_i_column`, `mag_z_column`, `mag_y_column` | INI, runtime; `0` marks an absent band | Change for another survey/band schema. Explicit projection should retain the selected magnitude; `process_fd` requires at least one band and chooses the first available in `i -> z -> r -> g -> y` order, without hand-edited FD offsets. |
| Calibration paths and Standard branches | `[lensing].flat_path`, `psf_path`, `astrometry_trivial=0`, `include_flat=0`, `include_mask=2`, `ext_cat=1`, `ext_psf=0`, `psf_type=1`, `psf_ms=0` | Standard INI runtime settings; Lite rejects keys for deleted branches | Change for another flat/external-PSF source or an alternate science branch. Lite is fixed to Gaia, no flat, per-CCD DQ, external sources, frame PSF, local polynomial, and no PCA. |
| Image and detector geometry | `[lensing].ccd_split=2`, `pixel_size=0.2628`, `nmax_chip=62`, `chipnx=2046`, `chipny=4094`; build-time `ns=64`, `chip_margin=8` | First five are runtime INI fields; edit the last two in `config/LensingConfig.hpp` and rebuild | Review together for another camera, amplifier layout, pixel scale, PSF-map geometry, or stamp size. Stage 1 reads Science FITS axes dynamically; successful Standard Hybrid PSF maps and FD bounds use effective runtime `chipnx/chipny`. This repository has no fixed `npx/npy` settings. |
| Numerical stages | `[lensing].process_stage=223092870` | INI, runtime | Select the nine main stages by prime factors; Stage 9 (23) requires Stage 8 (19). |
| Source detection and pixel threshold | `saturation_thresh=25000` | `config/LensingConfig.hpp`, build time | Recalibrate on representative data and rebuild when the image source, gain, or saturation definition changes. |
| FD detector rules | `bad_ccds={}`, `chip_mask_edge=50` | `config/FDConfig.hpp`, build time | Inclusive bounds derive from effective runtime geometry as `[edge, chipnx-edge]` and `[edge, chipny-edge]`. Change the edge or bad-CCD list for another camera; `n_bad_ccds` is derived. |

See [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) for every individual
Standard/Lite default, legal value, INI/CLI override, and rebuild rule. Preserve a baseline
configuration and validate coupled changes on the smallest representative dataset first.

## Common run modes

Initializer and numerical stages from one INI:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini
```

Main-only, using an existing exposure list:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

Archive initialization with CLI overrides:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini \
  --science-root /data/archive/science \
  --dq-root /data/archive/dqmask \
  --output-root /data/work \
  --dataset g2019:c4d_19 --existing resume
```

Gaia-catalog-only, with CLI overrides:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-astrocat true --run-extcat false --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --astrocat-input /data/raw_gaia --astrocat-output /data/gaia/tiles \
  --astrocat-add-header true --astrocat-existing fail
```

`--astrocat-output` controls only `process_astrocat`. It is not compared with,
propagated to, or required to match `[lensing].astrometry_cat`.

External-catalog-only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs \
  --extcat-output /data/catalogs/tiles
```

At least one phase must be enabled. If initialization runs successfully, its
generated absolute `expo_<target>.list` is used by later phases. In
downstream-only mode, omit `--expo-list` only when `output_root` and the dataset
name can derive it unambiguously.

On Slurm, replace `mpirun` with the site-supported launcher. The supplied
container runner requires `srun --mpi=pmi2`; see
[cpp_docker/runner/README.md](cpp_docker/runner/README.md).

## Inputs and outputs

The canonical definitions and minimum requirements for Science images, the
Gaia catalog, the External source catalog, and DQ masks are in the root
[Input data requirements](README.md#input-data-requirements). DQ masks are a
configuration-dependent input class: Lite always reads per-chip DQ masks,
whereas Standard may omit them only when the selected
`[lensing].include_mask` mode and active code path do not access DQ data.

An exposure list contains one chip-list path per nonblank record. A trailing
legacy chip count is accepted. Quoted paths are supported by the current C++
reader.

Initialization reads compressed archives in place and creates a dataset tree
below `output_root/<target>/`, including `science/`, `dqmask/`, `stamps/`, and
`result/`. It also publishes `expo_<target>.list`, `fits_<target>.list`, and an
initializer manifest.

The most important downstream products are:

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

`_all.cat` is a per-exposure shear catalog. By default, each row contains the
external source catalog fields, the original 1-based `EXPO_NUM`, one CCD number,
and 25 pipeline fields (45 columns total). Regenerate Stage 9, rearrangement,
and FD products after this schema change; legacy 44-column rows are incompatible.
`subcat_*.cat` spatially repartitions those rows by Dec and RA; repeated
measurements at identical coordinates are adjacent, which simplifies
deduplication. `FD_test_comb.dat` is the field-distortion test table written by
`process_fd` and is used to calibrate shear measurements.

## Frequent errors

- Do not enable Stage 9 without Stage 8.
- The external-catalog output must not equal or be nested below its input.
- With explicit column projection, include raw RA, Dec, and photo-z columns.
- Do not use Standard-only `[lensing]` keys with Lite.
- Do not run concurrent builds that clean the same source tree.
- Container paths, not host paths, belong in INI/CLI arguments executed inside
  a container.

## Container deployment

For local Docker use, see [cpp_docker/README.md](cpp_docker/README.md). For
Slurm/Apptainer, see [cpp_docker/runner/README.md](cpp_docker/runner/README.md).
The image is a toolchain/runtime image: it never supplies the pipeline source,
configuration, catalogs, or observation data.
