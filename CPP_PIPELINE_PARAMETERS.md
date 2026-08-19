# C++ Pipeline Parameters Reference

This file is the authoritative reference for **all adjustable parameters** in the
`cpp_Standard` and `cpp_Lite` executables. Parameters are organized by pipeline
function; both variants resolve one shared runtime `CatalogLayout` for extcat,
main, rearrangement, and FD catalog consumers.

**Column meanings**

| Column | Meaning |
|:---|:---|
| Parameter file name | The constant name as it appears in the configuration header (`ProcessConfig.hpp`, `LensingConfig.hpp`, or `ProcessRearrConfig.hpp`). |
| CLI parameter | The command-line option when one exists. `INI [section] key` marks a file-only runtime override. `—` means compile-time-only. |
| Options | All accepted values; the default is suffixed with `*`. For single-valued constants the sole value is shown with `*`. |
| Function description | What the parameter controls and what each option value does. |

**Configuration layers and precedence**

1. **Compiled headers** — `ProcessConfig.hpp`, `ExtCatConfig.hpp`,
   `InitConfig.hpp`, and the live runtime defaults in `LensingConfig.hpp` seed
   the selected variant.
2. **INI file** — `--config PATH` applies `[process]`, `[extcat]`, `[init]`, and
   `[lensing]` overrides without rebuilding.
3. **CLI** — matching command-line values override both compiled and file
   values. The effective precedence is compiled defaults < INI < CLI.
4. **Fixed headers** — catalog indices, array dimensions, algorithm constants,
   and `ProcessRearrConfig.hpp` remain compile-time-only.
5. **`CatalogLayout`** — derived runtime catalog schema. It is not a
   precedence layer: `main` resolves it once from normalized runtime options
   and existing `LensingConfig` source-field indices.

CLI options accept both `--name value` and `--name=value`. Boolean values accept
`true`, `false`, `1`, `0`, `on`, and `off`. Duplicate scalar options use the last
value. The first explicit `--dataset`, `--contains`, or `--extcat-contains`
clears its configured list; later occurrences append.

Rank 0 reads the INI once and broadcasts its exact text; every rank parses and
validates before initializing the immutable store. The root driver calls five
functions in order: `process_extcat` (once), then
`process_init` → `process_main` → `process_rearr` → `process_fd` (per dataset).
All phases may be disabled for a parse/validation dry run.

---

## Table 1 — process_extcat: External Catalog Repartitioning

Splits raw external catalogs into 0.1-degree sky tiles. Every field below is
available in `[extcat]` using its lower-case short name; listed CLI options have
the final precedence.

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
| `EXTCAT_TOTAL_COLUMNS` | — | `18*` | Pass-through external-catalog width. Each C++ variant uses it only while constructing `CatalogLayout` when explicit projection is disabled and verifies integrated extcat output matches it. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | `false` preserves all raw fields in place (pass-through); `true` enables explicit column projection. Setting `--extcat-columns` enables this. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated 1-based indices (default `1–18*`) | Ordered output column projection when `use_explicit_columns=true`. Both variants reject zero/duplicate entries and require configured RA/Dec/ZP. A configured positive magnitude omitted from the list is absent. The unique list length is the output width. |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` | `--extcat-ra-column`, `--extcat-dec-column` | `false*`, `true` | `false` uses named `ra`/`dec` header columns or compiled defaults for sky tiling; `true` enables explicit coordinate column indexing. Setting either `--extcat-ra-column` or `--extcat-dec-column` enables this. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | One-based raw RA column index used for sky tiling. Must be distinct from Dec and ZP. Also selects the RA consumed by `process_main`. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | One-based raw Dec column index used for sky tiling. Must be distinct from RA and ZP. Also selects the Dec consumed by `process_main`. |
| `EXTCAT_MAG_G_COLUMN_ONE_BASED` | — | `0` or positive integer (default `7*`) | Optional one-based raw g-band magnitude identity; `0` means absent. FD uses it after i/z/r fallback. |
| `EXTCAT_MAG_R_COLUMN_ONE_BASED` | — | `0` or positive integer (default `9*`) | Optional one-based raw r-band magnitude identity; `0` means absent. FD uses it after i/z fallback. |
| `EXTCAT_MAG_I_COLUMN_ONE_BASED` | — | `0` or positive integer (default `11*`) | Optional one-based raw i-band magnitude identity; `0` means absent. FD prefers i when available. |
| `EXTCAT_MAG_Z_COLUMN_ONE_BASED` | — | `0` or positive integer (default `13*`) | Optional one-based raw z-band magnitude identity; `0` means absent. FD uses it when i is absent. |
| `EXTCAT_MAG_Y_COLUMN_ONE_BASED` | — | `0` or positive integer (default `15*`) | Optional one-based raw y-band magnitude identity; `0` means absent. FD uses it as the final fallback. |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `--extcat-zp-column` | Positive integer (default `17*`) | One-based raw photometric-redshift identity consumed by main/FD. |

These three required fields and five optional magnitudes are the only named
external fields in each variant's `CatalogLayout`. Other physical columns contribute
to `external_columns` but have no layout member. FD requires at least one
magnitude and selects i -> z -> r -> g -> y; the selected band drives both the
magnitude-range cut and size-magnitude star-bar calculation.

---

## Table 2 — process_init: Archive Discovery and Extraction

Discovers Science/DQ FITS archives, extracts per-chip images, and publishes
exposure lists. Every initializer field is available in `[init]`; listed CLI
options have the final precedence.

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_INIT` | `--run-init` | `true*`, `false` | Phase switch. `true` runs archive discovery/extraction and publishes exposure lists; `false` skips it. |
| `SCIENCE_ROOT` | `--science-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/g"*`) | Read-only multi-HDU Science FITS/FZ archive repository. Archives are selected by dataset prefix and filename tokens. |
| `DQ_ROOT` | `--dq-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"*`) | Read-only multi-HDU DQ FITS/FZ archive repository paired with Science archives. |
| `OUTPUT_ROOT` | `--output-root` | Path string (default `"/lustre/home/acct-phyzj/share/DES/g_band_v1"*`) | Parent of dataset directories and `expo_<target>.list`. Contains extracted chip images and generated intermediate products. |
| `DATASETS` | `--dataset` | List of `TARGET:PREFIX` pairs (both variants default `{"gband","c4d_"}*`) | One or more target/prefix pairs. Repeatable; the first `--dataset` clears the compiled list. Legacy `--target`/`--prefix` set a single pair and cannot mix with `--dataset`. |
| `CONTAINS` | `--contains` | List of strings (default `{"v1"}*`; repeatable) | Case-sensitive basename substring filters for archive file discovery. A file matches if any token matches. The first CLI use clears the compiled list. |
| `EXISTING` | `--existing` | `fail*`, `resume`, `overwrite` | `fail` rejects existing output; `resume` keeps existing files and continues; `overwrite` replaces all existing output. |
| `F77_MAX_PATH` | `--f77-max-path` | Non-negative integer (both variants `0*`) | Maximum path length for generated exposure and chip-list files. `0` disables the limit. |
| `EXPO_LIST` | `--expo-list` (or positional `LEGACY_EXPO_LIST`) | Path string (default `""*`) | Single exposure list for main/rearr-only mode. When `--run-init=true`, the initializer output takes precedence. Cannot serve multiple datasets in downstream-only mode. |

---

## Table 3 — process_main: Numerical Pipeline Stages

Executes the Stage 1–9 Fourier_Quad shear pipeline. Run-selectable values are
stored in `RuntimeConfig`; fixed dimensions, indices, thresholds, and algorithm
constants remain in `LensingConfig.hpp`.

### 3a. Runtime (CLI-overridable) parameters

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_MAIN` | `--run-main` | `true*`, `false` | Phase switch. `true` runs the numerical Stage 1–9 pipeline; `false` skips it. |
| `SOURCE_CAT` / `EXTCAT_OUTPUT_DIRECTORY` | `--extcat-output` | Path string (default `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"*`) | External-catalog tile directory read by `process_main`. `--extcat-output` updates the effective `SOURCE_CAT` before processing. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | Controls external-catalog column resolution. In both variants, `false` uses pass-through positions and `true` resolves all mandatory external fields through the shared layout. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated unique 1-based indices (default `1–18*`) | Projection order used by each variant's `CatalogLayout`; it must retain configured RA/Dec/ZP. Retained positive magnitude identities become available to FD. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | Raw RA field consumed by `process_main` for source-catalog matching. Resolved through projection when enabled. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | Raw Dec field consumed by `process_main`. Resolved through projection when enabled. |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `--extcat-zp-column` | Positive integer (default `17*`) | Raw photometric-redshift (`dnf_z`) field consumed by `process_main` and required by the shared layout in both variants. Resolved through projection when enabled. |
| `EXPO_LIST` | `--expo-list` (or positional, or init output) | Path string (default `""*`) | Exposure-list file. Each non-empty line identifies a per-exposure chip-list file. An initializer-generated list takes precedence in chained execution. |

### 3b. Stage and Standard branch controls (runtime INI)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `PROCESS_stage` | `INI [lensing] process_stage` | Product of primes (default `2·3·5·7·11·13·17·19·23 = 223092870*`) | Stage control bitmask. Each prime factor enables one stage: `2`→Pre-process, `3`→Astrometry, `5`→Source extraction, `7`→FFT Stage 1, `11`→PSF model, `13`→FFT Stage 2, `17`→Shear measurement, `19`→Exposure info, `23`→Catalog combination. Stage 9 (`23`) requires Stage 8 (`19`). |
| `ASTROMETRY_trivial` | `INI [lensing] astrometry_trivial` | `0*`, `1` (Standard only; Lite fixed to `0`) | `0` uses Gaia-based astrometric solution; `1` uses trivial astrometry. The key is rejected by Lite because the alternative branch is absent. |
| `include_FLAT` | `INI [lensing] include_flat` | `0*`, `1` (Standard only; Lite fixed to `0`) | Enables super-flat multiplication with `flat_path`. The key is absent in Lite. |
| `include_Mask` | `INI [lensing] include_mask` | `0`, `1`, `2*`, `3` (Standard only; Lite fixed to `2`) | Selects no/legacy/DQ/combined masking. Lite implements only DQ masking and rejects the key. |

### 3c. Removed nominal image settings

`npx` and `npy` were removed. Failed Stage-1 reads return without fabricating a
nominal image; physical PSF geometry uses runtime `chipnx`/`chipny`.

### 3d. Split and background parameters (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `ext_cat` | `INI [lensing] ext_cat` | `0`, `1*` (Standard only; Lite fixed to `1`) | Selects internal or external source catalogs. Lite is external-only and rejects this key. Standard FD requires `1`. |
| `ext_PSF` | `INI [lensing] ext_psf` | `0*`, `1` (Standard only; Lite fixed to `0`) | Selects frame-star or external-PSF input. Lite rejects this key. |
| `CCD_split` | `INI [lensing] ccd_split` | `1`, `2*` | `1` uses one whole-chip background/noise region; `2` uses two amplifier regions. |
| `blocksize` | — | `200*` | Target pixel width/height for balanced background blocks; the implementation rounds the amplifier dimensions to a block count and covers the complete region. |
| `nct` | — | `12*` | Number of rectangular-monomial terms in the background surface fit. Must not exceed available stable background samples. |
| `ncx` | — | `3*` | Number of successive x powers in the background basis before y power increments. With `nct=12`, the basis spans `x^0..2` for `y^0..3`. |
| `bg_rough_grid_x` / `bg_rough_grid_y` | — | `32*` / `32*` | Deterministic rough-fit sampling grid. The rough bilinear fit is a residual preconditioner and is not subtracted from the image. |
| `bg_min_block_pixels` | — | `1000*` | Minimum weight-valid, finite pixels required for one background block point. |
| `bg_min_clipped_pixels` | — | `200*` | Minimum pixels surviving asymmetric pixel-level clipping inside one background block before its clipped mean is accepted. |
| `bg_min_valid_frac` | — | `0.25*` | Minimum valid-pixel fraction required for one background block point. |
| `bg_clip_low` / `bg_clip_high` | — | `4.0*` / `2.5*` | Asymmetric lower/upper MAD clipping thresholds for block residuals. |
| `bg_fit_clip_sigma` | — | `3.0*` | MAD clipping threshold applied to residuals of the final 12-term model. |
| `bg_fit_max_iter` | — | `4*` | Maximum iterative final-model MAD clipping passes. |
| `bg_min_fit_factor` | — | `3*` | Minimum final-fit block count factor relative to `nct`; the absolute minimum is 30. |

### 3e. PSF selection and configuration

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `npo` | — | `64*` | Reserved legacy PSF selection size; currently used only to derive `nstar_min`. |
| `nstar_min` | — | `npo·3/2 = 96*` | Base exposure-wide star threshold. PSF star selection rejects an exposure when total candidate count is below `2·nstar_min` (192 by default). |
| `npl` | — | `10*` | Number of ordered 2D polynomial terms fitted per PSF Fourier pixel; `10` includes terms through total degree 3. |
| `nstar_min_local` | — | `16*` | Minimum finite stars required for one local chip PSF fit. |
| `step_psf` | — | `100*` (Std only) | Standard very-local PSF map grid spacing in pixels; used only with `PSF_type=2`. Absent in Lite. |
| `n_neighbor` | — | `5*` (Std only) | Standard number of nearest stars used by the very-local PSF branch. Absent in Lite. |
| `deblending` | `INI [lensing] deblending` | `0`, `1*` (Standard only; Lite fixed to `1`) | Selects source deblending. Lite's always-on implementation has no key. |
| `PSF_type` | `INI [lensing] psf_type` | `1*`, `2` (Standard only; Lite fixed to `1`) | Selects local polynomial or very-local PSF fitting. Lite has only the local implementation. |
| `PSF_Ms` | `INI [lensing] psf_ms` | `0*`, `1` (Standard only; Lite fixed to `0`) | Enables Standard multi-scale/PCA reconstruction. Lite has no PCA implementation or key. |

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
| `area_max` | — | `ns·ns = 4096*` | Derived maximum connected-region workspace/area. |
| `area_thresh` | — | `6*` | Minimum connected-region pixel count. |

### 3g. Smoothing and detection limits (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `gal_smooth` | `INI [lensing] gal_smooth` | Non-negative integer (`0*` both variants) | Galaxy/noise FFT power smoothing. `0` disables smoothing; the existing smoothing implementation interprets positive radii/modes. |
| `star_smooth` | `INI [lensing] star_smooth` | Non-negative integer (`2*`) | Star FFT power smoothing and regularization control. |
| `size_fit_rmax` | — | `4*` | Radius in Fourier pixels for the equal-weight low-frequency curvature fits that produce `gal_size_T` and `psf_size_T`. |
| `point_stat_beta` | — | `0.10*` | Single survey-wide shape parameter for the Stage 7 extended template `P(k) exp[-beta (k/kmax)^2]`. Compile-time; tune offline and rebuild rather than fitting per source. |
| `point_stat_k_frac` | — | `0.90*` | Fraction of the PSF reliable Fourier radius retained as `kmax` for `delta_chi2` and `orth_ext`. |
| `point_stat_eps` | — | `1e-20*` | Numerical floor used only to reject singular point-source statistics and stabilize scalar normalization. It is not a science cut. |
| `point_stat_min_corr` | — | `1e-6*` | Minimum absolute normalized source–PSF projection needed for stable `orth_ext`; numerical protection only. |
| `SNR_PSF` | — | `100.0*` | S/N threshold used to select PSF-star candidates; some preliminary selection uses half this value. |
| `saturation_thresh` | — | `25000.0*` | Raw-pixel saturation cutoff and normalized peak rejection reference. |
| `pixel_size` | `INI [lensing] pixel_size` | Finite positive double (`0.2628` arcsec*) | Detector pixel scale used to convert PSF sizes to angular units. |

### 3h. Catalog and memory dimensions (compile-time)

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `ngal_max` | — | `4000*` | Initial source metadata reservation hint; vectors grow beyond it. |
| `nstar_max` | — | `2000*` | Initial star metadata reservation hint; live pairwise storage uses actual candidate counts. |
| `n_user_max` | — | `200*` | Maximum number of flux-ranked image-side astrometry detections passed to `patternMatching()`. Detection, sorting, and dynamic catalog storage remain uncapped; Standard and Lite use the same compile-time value. |
| `src_npara` | — | `12*` | Shared source/PSF/FFT2 internal row width. It exactly covers fields through `iSNR_F` and PSF parameter slot 11. |
| `shear_cat_ncols` | — | `iorth_ext + 1 = 28*` | Stage 7 shear-catalog width through `gal_size_T`, `psf_size_T`, `delta_chi2`, and `orth_ext`. The exposure `chi2` appended later is not part of this count. |
| `expo_cat_ncols` | — | `shear_cat_ncols + 1 = 29*` | Exposure-catalog width containing the 28 shear fields plus the appended exposure `chi2`. FD consumes this width directly. |
| `npd` | — | `33*` | Number of astrometric PU distortion coefficients per coordinate. Must match astrometry file serialization and fitting code. |
| `nmax_chip` | `INI [lensing] nmax_chip` | Positive integer (`62*`) | Runtime chip capacity used for initializer checks and PSF allocations. Exposure chip counts above it fail explicitly; there is no instrument-specific upper ceiling. |

`NMAX_EXPO` was removed. Exposure lists and FD per-exposure arrays use checked
runtime sizes; large MPI vectors are transferred in `INT_MAX`-bounded chunks.

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

### 3k. File-system paths

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `SOURCE_CAT` fallback / `extcat.output_directory` | `--extcat-output` or `INI [lensing] source_cat` | Path string (default `"/lustre/.../des_y6_cat"*`) | One authoritative external tile directory for generation and reading. `source_cat` is an INI alias. |
| `ASTROMETRY_CAT` | `INI [lensing] astrometry_cat` | Path string (default `"/lustre/.../gaia_cat_sorted"*`) | Gaia reference-tile directory used by both variants. |
| `FLAT_PATH` | `INI [lensing] flat_path` | Path string (Standard only) | Super-flat directory used when `include_flat=1`; absent from Lite config. |
| `PSF_PATH` | `INI [lensing] psf_path` | Path string (Standard only) | External PSF directory used when `ext_psf=1`; absent from Lite config. |

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
| `igalsizeT` | — | `24*` | Galaxy low-frequency Fourier-power curvature size `T` in pixel². |
| `ipsfsizeT` | — | `25*` | PSF low-frequency Fourier-power curvature size `T` in pixel². |
| `idelta_chi2` | — | `26*` | Normalized PSF-minus-extended weighted residual; smaller/negative values are more PSF-like and positive values favor the fixed extended template. |
| `iorth_ext` | — | `27*` | Source projection along the PSF-orthogonal extension direction; larger positive values are more extended. |
| `ichi2` | — | `shear_cat_ncols = 28*` | Zero-based index of the appended exposure chi2, immediately after the 28-field shear catalog. |

### 3m. Calibration and runtime camera geometry

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `g1_c` | — | `-0.001*` | Additive correction applied to field-distortion component 1 during final catalog combination. |
| `g2_c` | — | `-0.0003*` | Additive correction applied to field-distortion component 2 during final catalog combination. |
| `chi2_thresh` | — | `0.01*` | Maximum accepted exposure/PSF chi-square diagnostic during catalog combination. |
| `chipnx` | `INI [lensing] chipnx` | Positive integer (`2046*`) | Runtime science CCD width used for PSF maps, bounds, and coordinate normalization. |
| `chipny` | `INI [lensing] chipny` | Positive integer (`4094*`) | Runtime science CCD height used for PSF maps, bounds, and coordinate normalization. |

---

## Table 4 — process_rearr: Catalog Rearrangement

Redistributes per-exposure `_all.cat` rows into spatially sorted subcatalogs.
Compile-time algorithm parameters are from `ProcessRearrConfig.hpp`; runtime
paths and switches are seeded by `ProcessConfig.hpp` and may be overridden by
INI or CLI. Both variants consume the startup schema and do not derive width or
coordinates locally.

### 4a. Runtime (CLI-overridable) parameters

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `RUN_PROCESS_REARR` | `--run-rearr` | Standard: `true*`, `false`; Lite: `true`, `false*` | Phase switch. `true` rearranges generated `_all.cat` rows by sky region; runs after `process_main` when both are enabled, or independently on existing results. `false` skips it. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `--extcat-columns` | `false*`, `true` | The shared `CatalogLayout` uses `EXTCAT_TOTAL_COLUMNS` for pass-through or the unique projection length for explicit mode. Rearrangement does not derive width itself. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `--extcat-columns` | Comma-separated unique 1-based indices (default `1–18*`) | Projection list resolved once at startup in either variant; rearrangement consumes its effective width and coordinates. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `--extcat-ra-column` | Positive integer (default `5*`) | Raw RA identity mapped once into `layout.external.ra`. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `--extcat-dec-column` | Positive integer (default `6*`) | Raw Dec identity mapped once into `layout.external.dec`. |
| `EXPO_LIST` | `--expo-list` (or positional, or init output) | Path string (default `""*`) | Exposure-list file. Each per-exposure chip list is used to derive `result/<PREFIX>_all.cat` paths. |

### 4b. Shared runtime column layout

| Parameter file name | CLI parameter | Options | Function description |
|:---|:---|:---|:---|
| `layout.external_columns` | — | `18*` (pass-through) or unique projection length | Effective external prefix emitted/read by catalog phases. |
| `layout.external.ra/dec/zp` | — | Required effective zero-based indices | Mandatory external positions used by main/rearr/FD. |
| `layout.external.mag_g/mag_r/mag_i/mag_z/mag_y` | — | Optional effective zero-based index or `none` | Optional magnitude positions. Raw config `0` or projection omission produces `none`; FD selects i -> z -> r -> g -> y. |
| `layout.ccd` | — | `layout.external_columns` | Absolute zero-based CCD_NUM position. |
| `layout.source_base` | — | `layout.ccd + 1` | Absolute base of the process-main suffix. |
| `layout.source_columns` | — | `LensingConfig::expo_cat_ncols = 29*` | Count of 28 shear fields plus exposure Chi2. |
| `layout.all_columns` | — | `layout.source_base + layout.source_columns` | Exact runtime row width used by rearrangement and FD readers. |

There is no fixed total-width constant. The shipped pass-through layout in either variant
layout has `18 + 1 + 29 = 48` columns. A valid explicit projection with `N`
unique external fields has `N + 30` columns; missing mandatory fields fail at
startup. Extcat/main/rearr may run with no magnitude; FD fails collectively
before catalog I/O unless at least one magnitude position is available.

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

The process/extcat/init CLI is intentionally identical in both variants. Their
runtime lensing domains differ:

- Lite freezes Gaia astrometry, no flat, DQ-mask mode, external source catalogs,
  frame-derived PSF stars, deblending, local-polynomial PSF, and no PCA. The
  unselected Standard branches and their runtime keys are absent.
- Standard retains the optional external PSF, very-local PSF, and PCA branches.
- Both variants resolve one pipeline-level runtime `CatalogLayout`; Lite still
  defaults `process_rearr` and `process_fd` off while Standard defaults them on.

When moving a parameter file between variants, preserve these intentional
differences instead of copying one `LensingConfig.hpp` over the other.
