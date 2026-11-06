# Quick Start

[🏠 Home](../README.md) · [🗂️ Data Preparation](DATA_PREPARATION.md) · [📋 C++ User Guide](CPP_GUIDE.md) · [📋 Full Parameter Reference](CPP_PIPELINE_PARAMETERS.md)

This guide covers **site-, camera-, and dataset-specific settings** that a new user should review before running pipeline. It does not teach numerical tuning, alternative PSF algorithms, estimator choices, or internal science branches.

## 0. Install files
Click [here](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases/download/v1.5.9/cpp_Standard-v1.5.9.zip) to download the latest main program.
For users who utilize the docker image, click [here](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases/download/v1.5.9/cpp_docker-v1.5.9.zip) to download the latest launcher.

## 1. Create a run configuration

The pipeline combines compiled defaults in `cpp_Standard/config/*.hpp` with optional INI overrides.

```bash
cd cpp_Standard
cp pipeline.example.ini pipeline.ini
```
> The values in the examples are illustrative, **not a validated configuration for a new instrument**.

| Entry | Meaning |
| --- | --- |
| Omitted key or `# key =` | Keep the compiled default from `config/*.hpp`. |
| `key = value` | Override the compiled default. (empty value will clear the default !)|

## 2. Identify the input data and output locations

Prepare the source files following the [Data Preparation Guide](DATA_PREPARATION.md).

| INI setting | What to provide |
| --- | --- |
| `[init].science_root` | Directory containing raw Science FITS/FZ archives. |
| `[init].dq_root` | (Optional) Directory containing corresponding DQ FITS/FZ archives, when the selected mask mode needs them. |
| `[init].output_root` | Writable parent directory for extracted datasets and all generated science products. |
| `[init].datasets` | Comma-separated `TARGET:PREFIX` pairs, e.g. `gband:c4d_`. The target names the dataset output folder; the prefix matches archive basenames. |
| `[init].contains` | (Optional) Case-sensitive archive basename filter, OR-matched. **Review for every dataset**. |
| `[lensing].astrometry_cat` | Directory of prepared Gaia reference tiles consumed by the main pipeline. |
| `[lensing].source_cat` | Directory of prepared canonical external source-catalog tiles. |

For example:

```ini
[init]
science_root = /archive/science
dq_root = /archive/dqmask
output_root = /data/processed
datasets = gband:c4d_
contains = v1

[lensing]
astrometry_cat = /data/catalogs/gaia_tiles
source_cat = /data/catalogs/source_tiles
```

## 3. Review camera-specific geometry

The compiled defaults assume DECam. Verify the following values against the physical detector and extracted FITS images:

| INI setting | Compiled default | When to change |
| --- | --- | --- |
| `[lensing].pixel_size` | `0.2628` arcsec/pixel | Different angular pixel scale. |
| `[lensing].nmax_chip` | `62` | Different physical CCD count. |
| `[lensing].chipnx` | `2046` pixels | Different CCD width. |
| `[lensing].chipny` | `4094` pixels | Different CCD height. |
| `[lensing].ccd_split` | `2` | Whole CCD (`1`) rather than two amplifier halves (`2`). |

For illustration **only**, a verified 16-CCD camera might need different values:

```ini
[lensing]
pixel_size = 0.30
nmax_chip = 16
chipnx = 4096
chipny = 4096
ccd_split = 1
```

**Do not use these illustrative values without checking your instrument.** A different detector may require code-level validation beyond changes to these five settings.

Input-format conventions that are **compile-time only** also need review when changing surveys:

| Header | Setting | Current convention |
| --- | --- | --- |
| `config/InitConfig.hpp` | `ARCHIVE_SUFFIX` | `.fits.fz` mark archives for processing |
| `config/InitConfig.hpp` | `CCDNUM_KEYWORD` | `CCDNUM` identifies physical CCDs |
| `config/InitConfig.hpp` | `DQ_STEM_REPLACE_FROM` / `DQ_STEM_REPLACE_TO` | `ood` → `ooi` basename tokens differ between DQ and Science. |

## 4. Match your reference catalogs

> The default expects `astra_` / `extern_` for external catalog tile files, if not, must change the **compile-time only** variables:
> `ASTROMETRY_TILE_PREFIX` / `SOURCE_CAT_TILE_PREFIX`


**Gaia tiles:** For `[lensing].astrometry_cat_type = 1` users, make sure each tile contains a margin of half-width in each direction of the camera. In practice, for a camera coverage of $1.0\times1.5$ arcsec, a margin must be no less than $0.5\times0.75$ arcsec on both sides.

**External source catalog:** canonical tiles must include RA, Dec, photometric redshift (`zp`), and at least one observed-band magnitude. If the source survey has a different column schema, check the following one-based positions:

| INI key in `[extcat]` | What to verify |
| --- | --- |
| `total_columns` | Width of each pass-through catalog row. |
| `ra_column`, `dec_column`, `zp_column` | Raw one-based field indices. |
| `mag_g_column`, `mag_r_column`, `mag_i_column`, `mag_z_column`, `mag_y_column` | Raw magnitude indices; use `0` for a missing band. |
| `use_explicit_columns`, `input_columns` | Set together only for an ordered raw-column projection. |
| `use_explicit_coordinate_columns` | Enable if explicit RA/Dec indices should replace header discovery. |

> Do not copy column positions from another source catalog!

**If original catalogs are unavailable in format**, enable `[process].run_process_astrocat` and/or `run_process_extcat` and configure the corresponding `[astrocat].input_directory` / `output_directory` and `[extcat].input_directory` / `output_directory`. Raw input and tile output must be separate directories.

## 5. Select only the workflow phases you need

For a first pilot that initializes exposures and runs the main measurement using already prepared catalogs, set:

```ini
[process]
run_process_astrocat = false
run_process_extcat = false
run_process_init = true
run_process_main = true
run_process_rearr = false
run_process_fd = false
```

## 6. Build and validate a small run

```bash
cd cpp_Standard
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

If a compile-time header was changed:

```bash
make clean
make -j4
```

## 7. Check List
Before running the pipeline, check:

- [ ] Science/DQ archive names and physical `CCDNUM` values agree.
- [ ] `datasets` and `contains` select the intended exposures.
- [ ] The output directory is writable.
- [ ] Gaia layout, filenames, and sky coverage agree.
- [ ] External source columns and available magnitude bands match the catalog.
- [ ] Pixel scale, CCD geometry, and amplifier partition were verified.
- [ ] The pilot enables only the intended workflow phases.

For every advanced parameter, see the [Full Parameter Reference](CPP_PIPELINE_PARAMETERS.md).