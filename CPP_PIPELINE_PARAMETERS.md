# C++ Pipeline Parameters Reference

This file is the authoritative reference for **all adjustable parameters** in the
`cpp_Standard` and `cpp_Lite` executables. Parameters are organized into four
tables, one per pipeline function: `process_extcat`, `process_init`,
`process_main`, and `process_rearr`.

**Column meanings**

| Column | Meaning |
|:---|:---|
| Parameter file name | The constant name as it appears in the configuration header (`ProcessConfig.hpp`, `LensingConfig.hpp`, or `ProcessRearrConfig.hpp`). |
| CLI parameter | The command-line option that overrides the compiled default at runtime. `—` marks a compile-time-only parameter (no CLI override; rebuild required). |
| Options | All accepted values; the default is suffixed with `*`. For single-valued constants the sole value is shown with `*`. |
| Function description | What the parameter controls and what each option value does. |

**Configuration layers and precedence**

1. **`ProcessConfig.hpp`** — workflow defaults. A matching CLI option overrides the
   default for one invocation; rebuilding is not required.
2. **`LensingConfig.hpp`** — scientific and numerical constants. Changes require
   rebuilding the selected C++ variant.
3. **`ProcessRearrConfig.hpp`** — parameters used only by `process_rearr`.
   Changes require rebuilding.

CLI options accept both `--name value` and `--name=value`. Boolean values accept
`true`, `false`, `1`, `0`, `on`, and `off`. Duplicate scalar options use the last
value. The first explicit `--dataset`, `--contains`, or `--extcat-contains`
clears its configured list; later occurrences append.

The root driver calls the four functions in order: `process_extcat` (once),
then `process_init` → `process_main` → `process_rearr` (per dataset). At least
one phase must be enabled. Within a dataset, an enabled `process_rearr` always
follows an enabled `process_main`.

---

## Table 1 — process_extcat: External Catalog Repartitioning

Splits raw external catalogs into 0.1-degree sky tiles. All parameters below
are from `ProcessConfig.hpp` unless noted.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_EXTCAT` | `--run-extcat` | `true`, `false*` | Phase switch. `true` runs `process_extcat` as the first pipeline phase (requires non-empty input/output directories); `false` skips it. |
| `EXTCAT_INPUT_DIRECTORY` | `--extcat-input` | Path string (default `""*`) | Root directory containing raw catalog files. Required when `--run-extcat=true`. |
| `EXTCAT_OUTPUT_DIRECTORY` | `--extcat-output` | Path string (default = `SOURCE_CAT`) | Output tile directory. Defaults to `LensingConfig::SOURCE_CAT`; a CLI override also updates the effective `SOURCE_CAT` for `process_main`. Must not equal or be nested below the input directory. |
| `EXTCAT_FILENAME_TOKENS` | `--extcat-contains` | List of strings (default empty*; repeatable) | Case-sensitive substring filters matched against file basenames. A file matches if any token matches; an empty list accepts all regular files. The first CLI use clears the compiled list. |
| `EXTCAT_RECURSIVE` | `--extcat-recursive` | `true*`, `false` | `true` recurses into subdirectories when discovering raw catalogs; `false` scans only the top level. |
| `EXTCAT_DELIMITER` | `--extcat-delimiter` | `auto*`, `whitespace`, `comma`, `tab` | Selects the raw table delimiter. `auto` chooses comma when present, otherwise tab when present, otherwise whitespace. |
| `EXTCAT_HEADER_MODE` | `--extcat-header` | `auto*`, `present`, `absent` | Controls header handling. `auto` recognizes case-insensitive `ra`/`dec` column names and can classify the leading record as a header; `present` requires a header; `absent` treats all records as data. Leading blank and `#` comment lines are always skipped. |
| `EXTCAT_MALFORMED_POLICY` | `--extcat-malformed` | `fail*`, `skip` | `fail` stops collectively on the first malformed data row; `skip` skips malformed rows and reports the count. |
| `EXTCAT_EXISTING_POLICY` | `--extcat-existing` | `fail*`, `overwrite` | `fail` rejects existing generated tiles; `overwrite` transactionally replaces the complete generated tile set. |
| `EXTCAT_CHUNK_MIB` | `--extcat-chunk-mib` | Positive integer (default `64*`) | Approximate newline-aligned MPI byte-range task size in MiB. Controls task granularity, not the final tile size. |
| `EXTCAT_TOTAL_COLUMNS` | — | `18*` | Canonical column count for the DES Y6 GOLD reference schema. Compile-time constant; sets the pass-through output width when explicit projection is disabled. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | `false` preserves all raw fields in place (pass-through); `true` enables explicit column projection. Setting `--extcat-columns` enables this. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated 1-based indices (default `1–18*`) | Ordered output column projection when `use_explicit_columns=true`. Output column 1 uses the first listed raw field, and so on. Repeated indices are allowed; list length sets output width. |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` | `--extcat-ra-column`, `--extcat-dec-column` | `false*`, `true` | `false` uses named `ra`/`dec` header columns or compiled defaults for sky tiling; `true` enables explicit coordinate column indexing. Setting either `--extcat-ra-column` or `--extcat-dec-column` enables this. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | One-based raw RA column index used for sky tiling. Must be distinct from Dec and ZP. Also selects the RA consumed by `process_main`. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | One-based raw Dec column index used for sky tiling. Must be distinct from RA and ZP. Also selects the Dec consumed by `process_main`. |

> **Note:** `EXTCAT_ZP_COLUMN_ONE_BASED` is configured alongside these CLI options
> but consumed only by `process_main`; it appears in Table 3.

---

## Table 2 — process_init: Archive Discovery and Extraction

Discovers Science/DQ FITS archives, extracts per-chip images, and publishes
exposure lists. All parameters below are from `ProcessConfig.hpp`.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_INIT` | `--run-init` | `true`, `false*` | Phase switch. `true` runs archive discovery/extraction and publishes exposure lists; `false` skips it. |
| `SCIENCE_ROOT` | `--science-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/g"*`) | Read-only multi-HDU Science FITS/FZ archive repository. Archives are selected by dataset prefix and filename tokens. |
| `DQ_ROOT` | `--dq-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"*`) | Read-only multi-HDU DQ FITS/FZ archive repository paired with Science archives. |
| `OUTPUT_ROOT` | `--output-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/g_band_v1"*`) | Parent of dataset directories and `expo_<target>.list`. Contains extracted chip images and generated intermediate products. |
| `DATASETS` | `--dataset` | List of `TARGET:PREFIX` pairs (Std default `{"gband","c4d_"}*`; Lite default `{"g2019","c4d_19"}*`) | One or more target/prefix pairs. Repeatable; the first `--dataset` clears the compiled list. Legacy `--target`/`--prefix` set a single pair and cannot mix with `--dataset`. |
| `CONTAINS` | `--contains` | List of strings (default `{"v1"}*`; repeatable) | Case-sensitive basename substring filters for archive file discovery. A file matches if any token matches. The first CLI use clears the compiled list. |
| `EXISTING` | `--existing` | `fail*`, `resume`, `overwrite` | `fail` rejects existing output; `resume` keeps existing files and continues; `overwrite` replaces all existing output. |
| `F77_MAX_PATH` | `--f77-max-path` | Non-negative integer (Std `150*`; Lite `149*`) | Maximum path length for generated exposure and chip-list files. `0` disables the limit. |
| `EXPO_LIST` | `--expo-list` (or positional `LEGACY_EXPO_LIST`) | Path string (default `""*`) | Single exposure list for main/rearr-only mode. When `--run-init=true`, the initializer output takes precedence. Cannot serve multiple datasets in downstream-only mode. |

---

## Table 3 — process_main: Numerical Pipeline Stages

Executes the Stage 1–9 Fourier_Quad shear pipeline. CLI-overridable parameters
are from `ProcessConfig.hpp`; all others are compile-time constants from
`LensingConfig.hpp` (rebuild required to change).

### 3a. Runtime (CLI-overridable) parameters

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_MAIN` | `--run-main` | `true*`, `false` | Phase switch. `true` runs the numerical Stage 1–9 pipeline; `false` skips it. |
| `SOURCE_CAT` / `EXTCAT_OUTPUT_DIRECTORY` | `--extcat-output` | Path string (default `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"*`) | External-catalog tile directory read by `process_main`. `--extcat-output` updates the effective `SOURCE_CAT` before processing. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | Controls external-catalog column resolution. `false` uses pass-through positions; `true` resolves RA/Dec/ZP through the explicit projection. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated 1-based indices (default `1–18*`) | Projection order used to resolve RA/Dec/ZP tile positions when `use_explicit_columns=true`. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | Raw RA field consumed by `process_main` for source-catalog matching. Resolved through projection when enabled. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | Raw Dec field consumed by `process_main`. Resolved through projection when enabled. |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `--extcat-zp-column` | Positive integer (default `17*`) | Raw photometric-redshift (`dnf_z`) field consumed by `process_main`. Not needed for sky tiling. Resolved through projection when enabled. |
| `EXPO_LIST` | `--expo-list` (or positional, or init output) | Path string (default `""*`) | Exposure-list file. Each non-empty line identifies a per-exposure chip-list file. An initializer-generated list takes precedence in chained execution. |

### 3b. Stage control (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `PROCESS_stage` | — | Product of primes (default `2·3·5·7·11·13·17·19·23 = 223092870*`) | Stage control bitmask. Each prime factor enables one stage: `2`→Pre-process, `3`→Astrometry, `5`→Source extraction, `7`→FFT Stage 1, `11`→PSF model, `13`→FFT Stage 2, `17`→Shear measurement, `19`→Exposure info, `23`→Catalog combination. Stage 9 (`23`) requires Stage 8 (`19`); the pipeline rejects `23` without `19`. |
| `ASTROMETRY_trivial` | — | `0*`, `1` (Std only; Lite frozen to `0`) | `0` uses Gaia-based astrometric solution; `1` uses trivial (identity) astrometry. Lite implements only `0`. |
| `include_FLAT` | — | `0*`, `1` (Std only; Lite frozen to `0`) | `0` disables super-flat multiplication; `1` enables it using `FLAT_PATH`. Lite implements only `0`. |
| `include_Mask` | — | `0`, `1`, `2*`, `3` (Std only; Lite frozen to `2`) | `0` no mask; `1` legacy mask branch; `2` per-chip DQ mask from `dirOutput/dqmask/<exposure>/`; `3` combines legacy and DQ mask. Lite implements only `2`. |

### 3c. Image / CCD size (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `npx` | — | `3000*` | Nominal CCD width in pixels. |
| `npy` | — | `5000*` | Nominal CCD height in pixels. |
| `strl` | — | `150*` | Maximum string length for internal path/buffer allocations. |

### 3d. Split and background parameters (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `ext_cat` | — | `0`, `1*` (Std only; Lite frozen to `1`) | `0` disables external source catalog; `1` enables `SOURCE_CAT` usage. Lite implements only `1`. |
| `ext_PSF` | — | `0*`, `1` (Std only; Lite frozen to `0`) | `0` measures PSF from frame stars; `1` uses external PSF image from `PSF_PATH`. Lite implements only `0`. |
| `CCD_split` | — | `1`, `2*` | `1` whole-chip amplifier region for background/noise fits; `2` two-amplifier split. Other values are unsafe. |
| `blocksize` | — | `200*` | Pixel width/height of blocks sampled for background estimation. |
| `nct` | — | `12*` | Number of rectangular-monomial terms in the background surface fit. Must not exceed available stable background samples. |
| `ncx` | — | `3*` | Number of successive x powers in the background basis before y power increments. With `nct=12`, the basis spans `x^0..2` for `y^0..3`. |

### 3e. PSF selection and configuration (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `psf_order` | — | `8*` | Reserved legacy PSF-order constant; currently unused. |
| `npo` | — | `64*` | Reserved legacy PSF selection size; currently used only to derive `nstar_min`. |
| `npox` | — | `8*` | Reserved legacy PSF basis-width constant; currently unused. |
| `nstar_min` | — | `npo·3/2 = 96*` | Base exposure-wide star threshold. PSF star selection rejects an exposure when total candidate count is below `2·nstar_min` (192 by default). |
| `npl` | — | `10*` | Number of ordered 2D polynomial terms fitted per PSF Fourier pixel; `10` includes terms through total degree 3. |
| `nplx` | — | `2*` | Compatibility value copied into PSF routines; currently not consumed downstream. |
| `nstar_min_local` | — | `16*` | Minimum finite stars required for one local chip PSF fit. |
| `step_psf` | — | `100*` (Std only) | Standard very-local PSF map grid spacing in pixels; used only with `PSF_type=2`. Absent in Lite. |
| `n_neighbor` | — | `5*` (Std only) | Standard number of nearest stars used by the very-local PSF branch. Absent in Lite. |
| `deblending` | — | `0`, `1*` (Std only; Lite frozen to `1`) | `0` disables source deblending; `1` always applies deblending. Lite implements only `1`. |
| `PSF_type` | — | `1*`, `2` (Std only; Lite frozen to `1`) | `1` uses local polynomial PSF fit; `2` uses very-local PSF map with nearest-neighbor interpolation. Lite implements only `1`. |
| `PSF_Ms` | — | `0*`, `1` (Std only; Lite frozen to `0`) | `0` disables multi-scale/PCA PSF reconstruction; `1` enables it (Standard only). When `1`, all PCA parameters in §3j become active. |

### 3f. Stamp geometry and detection (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `ns` | — | `64*` | Square source/star stamp width in pixels. FFT routines assume the configured dimensions consistently. |
| `nsns` | — | `ns·ns = 4096*` | Derived stamp pixel count. |
| `chip_margin` | — | `8*` | Extra pixels around the half-stamp extraction region. |
| `ns_2` | — | `ns/2 = 32*` | Derived half-stamp size. |
| `nl_2` | — | `ns_2 + chip_margin = 40*` | Derived half-width of the larger extraction region. |
| `nl` | — | `2·nl_2 = 80*` | Derived larger extraction width. |
| `flag_thresh` | — | `3*` | Maximum reserved/flagged edge distance used when accepting a source stamp. |
| `chip_edge_margin` | — | `chip_margin = 8*` | Derived alias used to reject sources too near chip edges. |
| `dz_thresh` | — | `0.1*` | Maximum redshift difference for keeping neighboring entries during external-catalog deblending. |
| `source_thresh` | — | `2.0*` | Detection/defect connected-pixel threshold in noise-sigma units. |
| `core_thresh` | — | `4.0*` | Detection-core/peak threshold in noise-sigma units. |
| `flat_thresh` | — | `0.01*` (Std only) | Reserved Standard legacy flat threshold; currently unused. Absent in Lite. |
| `area_max` | — | `ns·ns = 4096*` | Derived maximum connected-region workspace/area. |
| `area_thresh` | — | `6*` | Minimum connected-region pixel count. |

### 3g. Smoothing and detection limits (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `gal_smooth` | — | `0*` Std / `2*` Lite | Galaxy/noise FFT power smoothing. `0` none; `1` linear 5×5 hole-aware; `2` log-power 5×5 hole-aware. This intentional default difference changes main-galaxy processing. |
| `star_smooth` | — | `0`, `1`, `2*` | Star FFT power smoothing with the same modes. In star power normalization, `0` uses four neighboring central pixels while values `≥1` use the central pixel. |
| `SNR_PSF` | — | `100.0*` | S/N threshold used to select PSF-star candidates; some preliminary selection uses half this value. |
| `saturation_thresh` | — | `25000.0*` | Raw-pixel saturation cutoff and normalized peak rejection reference. |
| `pixel_size` | — | `0.2628` arcsec* | Detector pixel scale used to convert PSF sizes to angular units. |

### 3h. Catalog and memory dimensions (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `len_g` | — | `40*` | Galaxy stamps per FITS layout row/block used by stamp I/O. |
| `len_s` | — | `15*` | Star stamps per FITS layout row/block. |
| `ngal_max` | — | `4000*` | Maximum galaxies stored per chip. Larger detected catalogs are truncated at this limit. |
| `nstar_max` | — | `2000*` | Maximum stars stored per chip. Memory use includes arrays scaling as `nstar_max²`. |
| `npara` | — | `25*` | Number of per-source/per-star parameter slots in internal tables. Must cover all configured indices. |
| `len_sam` | — | `50*` | Number of exposure-wide selected-star stamps placed in one FITS layout row/block. |
| `npd` | — | `33*` | Number of astrometric PU distortion coefficients per coordinate. Must match astrometry file serialization and fitting code. |
| `NMAX_EXPO` | — | `25000*` | Maximum exposure records allocated for aggregation. |
| `NMAX_CHIP` | — | `62*` | Maximum chips represented in fixed-capacity PSF/exposure arrays. |

### 3i. Mode-bar noise-plane estimator (compile-time)

These parameters act together and should normally remain synchronized with the
validated estimator convention rather than be tuned independently.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `sig_blocksize` | — | `200*` | Pixel side length of one robust seed block. |
| `sig_block_max` | — | `sig_blocksize² = 40000*` | Derived maximum pixels stored per block. |
| `sig_max_blocks` | — | `2048*` | Maximum block seeds retained per amplifier. |
| `sig_min_block_pixels` | — | `1000*` | Minimum valid pixels required for one block seed. |
| `sig_min_block_triples` | — | `1000*` | Minimum valid pixel triples required for a block estimate. |
| `sig_min_blocks` | — | `4*` | Minimum accepted blocks required for the estimator. |
| `sig_hist_nbin` | — | `256*` | Histogram bin count for locating the sky mode. |
| `sig_hist_range` | — | `6.0*` | Histogram half/range scale in robust-width units. |
| `sig_min_mode_count` | — | `500*` | Minimum sample count supporting mode estimation. |
| `sig_min_lower_count` | — | `1000*` | Minimum lower-side sample count for width estimation. |
| `sig_lower_quantile` | — | `0.3173105*` | Lower-side quantile used by the mode-width estimator. |
| `sig_clip_k` | — | `3.0*` | Symmetric clipping threshold in sigma units. |
| `sig_rdil` | — | `2*` | Radius used to dilate the estimator's private brightness mask. |
| `sig_clip_niter` | — | `2*` | Number of clipped plane-fit iterations. |
| `sig_min_fit_triples` | — | `1000*` | Minimum triples entering a noise-plane fit. |
| `sig_min_fit_frac` | — | `0.20*` | Minimum surviving fraction of possible fit triples. |
| `sig_median_ratio` | — | `1.2678405*` | Calibration from robust lower-side width to sigma convention. |
| `sig_plane_min` | — | `1.0e-8*` | Minimum positive fitted noise-plane value. |
| `sig_max_plane_ratio` | — | `4.0*` | Maximum allowed variation ratio across a fitted plane. |
| `sig_pivot_min` | — | `1.0e-8*` | Minimum numerical pivot accepted by the small plane solve. |
| `sig_scale_s1` | — | `0.673475*` | Retained Stage-1 calibration candidate. |
| `sig_scale_s2` | — | `1.027786*` | Stage-2 calibration candidate used by the current pipeline. |
| `sig_scale` | — | `sig_scale_s2 = 1.027786*` | Active derived selector converting the fitted plane to the published `2·sigma²` convention. |

### 3j. Standard-only multi-scale/PCA PSF parameters (compile-time)

These values are compiled only by `cpp_Standard` and are active only when
`PSF_Ms=1`. `cpp_Lite` removes the PCA implementation and all of these
parameters.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `rescale_size` | — | `1.2*` | Target PSF size used to derive the residual rescaling factor. |
| `procs_pn` | — | `40*` | Process-group size passed to PCA reconstruction scheduling. |
| `work_pn` | — | `10*` | Concurrent worker count within the PCA scheduler grouping. |
| `nblocks` | — | `2*` | Number of spatial blocks per CCD axis; default creates a 2×2 block grid. |
| `n_pcs` | — | `100*` | Maximum number of residual principal components stored and fitted. |
| `npp6th` | — | `28*` | Number of ordered 2D sixth-degree polynomial terms used for PCA coefficient surfaces. |
| `pca_negative_eigenvalue_threshold` | — | `-1.0e-5*` | Eigenvalue below which a PCA covariance result is classified as invalid. |
| `nmax_star_pchip` | — | `1000000*` | Reserved legacy per-chip PCA star capacity; currently unused. |

### 3k. File-system paths (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `SOURCE_CAT` | `--extcat-output` | Path string (default `"/lustre/.../des_y6_cat"*`) | Primary external source-catalog tile directory. Overridden at runtime by `--extcat-output`; also used as the extcat output default. |
| `ASTROMETRY_CAT` | — | Path string (default `"/lustre/.../gaia_cat_sorted"*`) | Directory of Gaia reference tiles expected by `generateGaiaFileName`. |
| `FLAT_PATH` | — | Path string (default `"/lustre/.../DES_super_flat/i2014"*`; Std only) | Per-chip flat FITS files used for super-flat multiplication when `include_FLAT=1`. Absent in Lite. |
| `PSF_PATH` | — | Path string (default `"hahahaha"*`; Std only) | Directory containing `PSF.fits` used when `ext_PSF=1`. Absent in Lite. |

### 3l. Internal catalog column indices (compile-time, 0-based)

These are zero-based positions in the C++ per-source result rows. They are not
the 18 raw external-catalog projection indices. Changing them changes the
internal/output layout and requires coordinated reader/writer changes.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `isig` | — | `3*` | Reserved; currently unused. |
| `istar` | — | `4*` | Star indicator / classification field. |
| `ipeak` | — | `4*` | Reserved alias; currently unused. |
| `i_imax` | — | `5*` | Peak x index. |
| `i_jmax` | — | `6*` | Peak y index. |
| `ih_flux` | — | `7*` | Half-light/source flux field. |
| `ih_area` | — | `8*` | Source area field. |
| `iflag` | — | `9*` | Quality flag. |
| `iPSF` | — | `10*` | PSF/SNR-related stored field. |
| `iSNR_F` | — | `11*` | Fourier S/N field. |
| `ira` | — | `12*` | Right ascension. |
| `idec` | — | `13*` | Declination. |
| `igf1` | — | `14*` | Field-distortion component 1. |
| `igf2` | — | `15*` | Field-distortion component 2. |
| `ig1` | — | `16*` | Fourier_Quad shear estimator 1. |
| `ig2` | — | `17*` | Fourier_Quad shear estimator 2. |
| `ide` | — | `18*` | Fourier_Quad normalization estimator. |
| `ih1` | — | `19*` | Higher-order estimator 1. |
| `ih2` | — | `20*` | Higher-order estimator 2. |
| `icos2` | — | `21*` | Spin-2 cosine term. |
| `isin2` | — | `22*` | Spin-2 sine term. |
| `iparity` | — | `23*` | WCS parity. |
| `ichi2` | — | `24*` | Zero-based index of the 25th field (exposure chi2); also consumed by `process_rearr` as `ProcessRearrConfig::ichi2 = 25`. |

### 3m. Calibration and camera geometry (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `g1_c` | — | `-0.001*` | Additive correction applied to field-distortion component 1 during final catalog combination. |
| `g2_c` | — | `-0.0003*` | Additive correction applied to field-distortion component 2 during final catalog combination. |
| `chi2_thresh` | — | `0.01*` | Maximum accepted exposure/PSF chi-square diagnostic during catalog combination. |
| `chipnx` | — | `2046*` | Nominal science CCD width used to normalize PSF coordinates. Must match the instrument geometry. |
| `chipny` | — | `4094*` | Nominal science CCD height used to normalize PSF coordinates. Must match the instrument geometry. |
| `Camera_ccd_num` | — | `62*` | Instrument CCD count used by exposure and PCA allocations. |

---

## Table 4 — process_rearr: Catalog Rearrangement

Redistributes per-exposure `_all.cat` rows into spatially sorted subcatalogs.
Compile-time parameters are from `ProcessRearrConfig.hpp`; CLI-overridable
parameters are from `ProcessConfig.hpp`. All rearrangement parameters are
identical in Standard and Lite.

### 4a. Runtime (CLI-overridable) parameters

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_REARR` | `--run-rearr` | `true`, `false*` | Phase switch. `true` rearranges generated `_all.cat` rows by sky region; runs after `process_main` when both are enabled, or independently on existing results. `false` skips it. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | Controls the effective external-field width via `externalCatalogColumns()`. `false` uses `EXTCAT_TOTAL_COLUMNS`; `true` uses the explicit projection length. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated 1-based indices (default `1–18*`) | Projection list whose length sets the runtime-effective external width when `use_explicit_columns=true`. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | Raw RA field position used by rearr to derive the effective RA column in the `_all.cat` row. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | Raw Dec field position used by rearr to derive the effective Dec column in the `_all.cat` row. |
| `EXPO_LIST` | `--expo-list` (or positional, or init output) | Path string (default `""*`) | Exposure-list file. Each per-exposure chip list is used to derive `result/<PREFIX>_all.cat` paths. |

### 4b. Derived column layout (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `ichi2` | — | `LensingConfig::ichi2 + 1 = 25*` | Derived count of the 24 shear fields plus exposure Chi2 appended after CCD_NUM. It is a count, not the zero-based last index. |
| `CCD_COLUMN_COUNT` | — | `1*` | Derived fixed CCD_NUM field count. |
| `ALL_CAT_TOTAL_COLUMNS` | — | `EXTCAT_TOTAL_COLUMNS + 1 + ichi2 = 44*` | Compile-time pass-through `_all.cat` width. Compile-time validation permits a different positive external width; focused tests verify the shipped 44-column default. |
| `externalCatalogColumns(options)` | — | `18*` (pass-through) or projection length | Runtime-effective external width. Explicit projection output contains exactly its selected fields. |
| `allCatalogColumns(options)` | — | Effective external width `+ 1 + ichi2*` | Runtime-effective exact row width used by the parser and MPI transfers. |

### 4c. Spatial partitioning and output (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `SKY_GRID_DEGREES` | — | `0.1*` | Full-sky bin width in both RA and Dec. |
| `RA_BIN_COUNT` | — | `3600*` | Full-sky RA grid dimension. |
| `DEC_BIN_COUNT` | — | `1800*` | Full-sky Dec grid dimension. |
| `SKY_TILE_COUNT` | — | `RA_BIN_COUNT × DEC_BIN_COUNT = 6480000*` | Derived flattened full-sky grid length. |
| `TARGET_SUBCAT_ROWS` | — | `500000*` | Target rows per weighted k-d partition; actual counts follow indivisible 0.1-degree tiles. |
| `OUTPUT_DIRECTORY` | — | `"rearranged_catalog"*` | Absolute destination when absolute; otherwise a child of the dataset root. Empty selects the dataset root itself. |
| `SUBCAT_PREFIX` | — | `"subcat_"*` | Generated partition filename prefix. |
| `SUBCAT_EXTENSION` | — | `".cat"*` | Generated partition filename extension. |
| `SUBCAT_ID_WIDTH` | — | `6*` | Minimum zero-padded partition-ID width. Larger IDs are not truncated. |
| `SUMMARY_FILENAME` | — | `"catalog_summary.txt"*` | Partition count and raw RA/Dec bound report written beside subcatalogs. |
| `OUTPUT_PRECISION` | — | `10*` | Significant digits for catalog values. |
| `SUMMARY_PRECISION` | — | `4*` | Fixed decimals for summary bounds. |
| `SKIP_MISSING_CATALOGS` | — | `true*`, `false` | `true` counts and skips absent per-exposure `_all.cat` files (legacy behavior); `false` fails on missing catalogs. |
| `SKIP_MALFORMED_ROWS` | — | `true*`, `false` | `true` skips and counts rows with wrong width or any non-finite/non-numeric field; `false` fails collectively on the first malformed row. |

`process_rearr` copies the shared input header, sorts every emitted partition by
Dec then RA, and uses source exposure/row only as deterministic tie breakers.
It writes `subcat_NNNNNN.cat` files and the summary below the effective output
directory. Existing same-name files are truncated like the legacy pipeline.

---

## Standard and Lite summary

The runtime command line is intentionally identical in both variants. Their
current configuration differences are:

- Standard default dataset: `gband:c4d_`; Lite: `g2019:c4d_19`.
- Standard `F77_MAX_PATH=150`; Lite `149`.
- Standard `gal_smooth=0`; Lite `gal_smooth=2`.
- Lite freezes Gaia astrometry, no flat, DQ-mask mode, external source catalogs,
  frame-derived PSF stars, deblending, local-polynomial PSF, and no PCA. The
  unselected Standard branches and their exclusive constants are absent.
- Standard retains the optional external PSF, very-local PSF, and PCA branches.

When moving a parameter file between variants, preserve these intentional
differences instead of copying one `LensingConfig.hpp` over the other.
