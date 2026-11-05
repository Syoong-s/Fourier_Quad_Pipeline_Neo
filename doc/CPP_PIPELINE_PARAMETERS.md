# 📋 C++ Pipeline Parameter Reference

[🏠 Home](../README.md) · [🛠️ C++ User Guide](CPP_GUIDE.md) · [🗂️ Data Preparation](DATA_PREPARATION.md) · [🔗 Container Paths](PATHS.md)

Complete reference for all compile-time configuration headers (`cpp_Standard/config/*.hpp`) and runtime `pipeline.ini` keys in `Fourier_Quad_Pipe`.

---

## ⚙️ How to Use This Reference

> [!NOTE]
> **Configuration Precedence:** `Compiled Header Defaults` → `Explicitly Passed --config pipeline.ini`.
> For general INI syntax rules, omitted vs. empty behavior, and build commands, see the **[C++ User Guide (CPP_GUIDE.md)](CPP_GUIDE.md#configuration)**.

| Table Value / Notation | Meaning |
| :--- | :--- |
| **`[section].key`** in *INI override* | Runtime-configurable via `--config pipeline.ini` (no rebuild required). |
| **`No`** in *INI override* | Compile-time constant only; edit the header in `cpp_Standard/config/*.hpp` and run `make clean && make -j4`. |
| **Derived constants** | Shown with their formula (e.g., `ns * ns = 4096`); always modify the source parameter rather than derived expressions or column indices. |

---

## 🧭 Quick Jump

| Category | Sections |
| :--- | :--- |
| **Workflow & Data Ingestion** | [Workflow (`ProcessConfig`)](#workflow-processconfig) · [Archive Initialization (`InitConfig`)](#archive-initialization-initconfig) · [Tile Naming](#tile-naming) · [Gaia Preparation (`AstroCatConfig`)](#gaia-preparation-astrocatconfig) · [External Catalog (`ExtCatConfig`)](#external-catalog-preparation-extcatconfig) |
| **Core Lensing Pipeline** | [Catalog & Calibration Paths](#catalog-and-calibration-paths-lensingconfig) · [Branches & Geometry](#branches-and-geometry-lensingconfig) · [PSF Selection](#psf-selection-lensingconfig) · [Stamps & Detection](#stamps-and-detection-lensingconfig) · [Background](#background-estimation-lensingconfig) · [Noise Products](#noise-products-lensingconfig) · [Sigma Estimation](#sigma-estimation-lensingconfig) · [Catalog Indices](#catalog-indices-lensingconfig) · [Calibration & PCA](#calibration-and-pca-lensingconfig) |
| **Post-Processing & Layout** | [Spatial Rearrangement (`ProcessRearrConfig`)](#spatial-rearrangement-processrearrconfig) · [Output Layout (`OutputLayout`)](#output-layout-outputlayout) · [Field Distortion (`FDConfig`)](#field-distortion-fdconfig) · [Coupled Values Summary](#runtime-selection-and-coupled-values) |

---

## Workflow (`ProcessConfig`)

Edit `config/ProcessConfig.hpp` for phase defaults and `config/pathconfig.hpp` for paths. At least one phase must be enabled. Phases execute in the fixed order shown in the C++ guide.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `RUN_PROCESS_ASTROCAT` (`bool`) | `false` | `[process].run_process_astrocat` | Prepare Gaia tiles once, before other phases. |
| `RUN_PROCESS_EXTCAT` (`bool`) | `false` | `[process].run_process_extcat` | Prepare external source tiles once, before dataset phases. |
| `RUN_PROCESS_INIT` (`bool`) | `true` | `[process].run_process_init` | Extract Science/DQ archives and publish dataset lists. |
| `RUN_PROCESS_MAIN` (`bool`) | `true` | `[process].run_process_main` | Run the numerical measurement stages. |
| `RUN_PROCESS_REARR` (`bool`) | `true` | `[process].run_process_rearr` | Spatially partition per-exposure catalogs. |
| `RUN_PROCESS_FD` (`bool`) | `true` | `[process].run_process_fd` | Run scientific field-distortion analysis; this is not a container diagnostic. |
| `EXPO_LIST` (`const char*`) | `""` | `[process].expo_list` | Top-level exposure-list file for downstream-only processing. Empty derives `output_root/expo_<target>.list`. One explicit list cannot serve multiple main/rearrangement datasets. |
| `REARR_OUTPUT_DIRECTORY` (`const char*`) | `"baked"` | `[process].rearr_output_directory` | Rearranged catalog directory: absolute, or relative to the configured base/dataset root. Empty selects the base itself. |
| `REARR_OUTPUT_BASE_DIRECTORY` (`const char*`) | `""` | `[process].rearr_output_base_directory` | Base for relative rearrangement output. Empty derives the dataset root from input image paths. |
| `REARRANGED_EXPO_LIST_FILENAME` (`const char*`) | `"cat_gband_ori.list"` | `[process].rearranged_expo_list_filename` | Published rearranged-list filename. It is a file name, not a directory. |
| `REARRANGED_EXPO_LIST_DIRECTORY` (`const char*`) | `""` | `[process].rearranged_expo_list_directory` | Directory for the published list. Empty uses the input exposure-list parent. |
| `FD_EXPO_LIST` (`const char*`) | `""` | `[process].fd_expo_list` | Explicit FD input-list file. Empty selects the rearranged-list directory and filename. |
| `FD_OUTPUT_DIRECTORY` (`const char*`) | `"fdout"` | `[process].fd_output_directory` | FD output directory: absolute, or relative to the FD base/dataset root. |
| `FD_OUTPUT_BASE_DIRECTORY` (`const char*`) | `""` | `[process].fd_output_base_directory` | Base for relative FD output. Empty uses the dataset root. |

## Archive initialization (`InitConfig`)

Configured in `config/InitConfig.hpp` and `config/pathconfig.hpp`. Compiled absolute paths are site-specific examples, not portable defaults.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `DatasetSpec::target` (`std::string`) | `"gband"` | `[init].datasets` | One nonempty output-directory component, not `.`, `..` or a path; unique within the dataset list. |
| `DatasetSpec::prefix` (`std::string`) | `"c4d_"` | `[init].datasets` | Nonempty archive basename prefix paired with target. |
| `SCIENCE_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/g"` | `[init].science_root` | Readable Science archive root. |
| `DQ_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"` | `[init].dq_root` | Readable DQ archive root. Check both initialization and selected mask consumers before omitting DQ. |
| `OUTPUT_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/g_band_v1"` | `[init].output_root` | Writable parent for dataset trees and published exposure lists. |
| `DATASETS` (`std::vector<DatasetSpec>`) | `{{"gband", "c4d_"}}` | `[init].datasets` | Comma-separated `TARGET:PREFIX` pairs; enabled dataset phases run sequentially for each pair. |
| `CONTAINS` (`std::vector<std::string>`) | `{"v1"}` | `[init].contains` | Case-sensitive archive basename tokens, OR matched. An active empty list removes this filter. |
| `ARCHIVE_SUFFIX` (`const char*`) | `".fits.fz"` | No | Case-sensitive archive suffix, removed when forming exposure stems. |
| `CCDNUM_KEYWORD` (`const char*`) | `"CCDNUM"` | No | FITS keyword holding physical detector identity; do not substitute HDU sequence numbers. |
| `DQ_STEM_REPLACE_FROM` (`const char*`) | `"ood"` | No | Substring replaced throughout a DQ stem to match the Science exposure. Empty disables replacement. |
| `DQ_STEM_REPLACE_TO` (`const char*`) | `"ooi"` | No | Replacement substring; change together with the source token. |
| `EXISTING` (`const char*`) | `"fail"` | `[init].existing` | `fail`, `resume` or `overwrite`. Choose deliberately when outputs already exist. |
| `F77_MAX_PATH` (`int`) | `0` | `[init].f77_max_path` | Nonnegative generated-path length limit; `0` disables this compatibility guard. |

## Tile naming

Both definitions are in `config/pathconfig.hpp`. The full grammar is `<prefix>RA_DDD_DDD_Dec_[pm]DD_[pm]DD.dat`. Rebuild and regenerate or rename tiles after changing a prefix. Old-prefix files are not treated as current generated products.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `AstroCatConfig::ASTROMETRY_TILE_PREFIX` (`std::string_view`) | `"astra_"` | No | Gaia one-degree tile prefix; excludes the fixed `RA_` token. |
| `ExtCatConfig::SOURCE_CAT_TILE_PREFIX` (`std::string_view`) | `"extern_"` | No | External source one-degree tile prefix; excludes `RA_`. |

## Gaia preparation (`AstroCatConfig`)

Edit `config/AstroCatConfig.hpp` and `config/pathconfig.hpp`. See [Data Preparation: Gaia Catalog](DATA_PREPARATION.md#gaia-astrometric-catalog) for layout conventions and deduplication rules.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `ASTROCAT_INPUT_DIRECTORY` (`const char*`) | `""` | `[astrocat].input_directory` | Readable flat directory of raw RA/Dec files; only direct regular children are read. |
| `ASTROCAT_OUTPUT_DIRECTORY` (`std::string`) | `LensingConfig::ASTROMETRY_CAT` | `[astrocat].output_directory` | Output directory must be disjoint from input: neither may contain the other. Separate from the main-phase consumer at runtime. |
| `ASTROCAT_ADD_HEADER` (`bool`) | `true` | `[astrocat].add_header` | `true`: input starts with data; `false`: skip exactly one line per input. Output always has a header. |
| `ASTROCAT_EXISTING_POLICY` (`const char*`) | `"fail"` | `[astrocat].existing_policy` | `fail` or `overwrite`. Overwrite replaces only files matching the current generated-prefix contract. |

## External catalog preparation (`ExtCatConfig`)

Edit `config/ExtCatConfig.hpp` and `config/pathconfig.hpp`. See [Data Preparation: External Source Catalog](DATA_PREPARATION.md#external-source-catalog) for schema and column-mapping rules.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `EXTCAT_INPUT_DIRECTORY` (`const char*`) | `""` | `[extcat].input_directory` | Readable raw catalog root. |
| `EXTCAT_OUTPUT_DIRECTORY` (`const std::string&`) | `LensingConfig::SOURCE_CAT` | `[extcat].output_directory` | Writable published tile directory. Must not equal or be nested below the raw input. |
| `EXTCAT_FILENAME_TOKENS` (`std::vector<std::string>`) | `{}` | `[extcat].filename_tokens` | Case-sensitive basename tokens with OR matching; empty accepts all basenames. |
| `EXTCAT_RECURSIVE` (`bool`) | `true` | `[extcat].recursive` | Search subdirectories of the raw input. |
| `EXTCAT_DELIMITER` (`const char*`) | `"auto"` | `[extcat].delimiter` | `auto`, `whitespace`, `comma` or `tab`. |
| `EXTCAT_HEADER_MODE` (`const char*`) | `"auto"` | `[extcat].header_mode` | `auto`, `present` or `absent`. |
| `EXTCAT_MALFORMED_POLICY` (`const char*`) | `"fail"` | `[extcat].malformed_policy` | `fail` or `skip`. Skip explicitly accepts loss of malformed rows. |
| `EXTCAT_EXISTING_POLICY` (`const char*`) | `"fail"` | `[extcat].existing_policy` | `fail` or `overwrite` existing generated tiles. |
| `EXTCAT_CHUNK_MIB` (`std::uint64_t`) | `64` | `[extcat].chunk_mib` | Positive MPI byte-range task size in MiB. |
| `EXTCAT_TOTAL_COLUMNS` (`std::size_t`) | `18` | `[extcat].total_columns` | Positive canonical catalog width in pass-through mode. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` (`bool`) | `false` | `[extcat].use_explicit_columns` | Enable ordered projection; must be set explicitly in the INI. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` (`std::vector<std::size_t>`) | `{1, 2, ..., 18}` | `[extcat].input_columns` | Nonempty list of positive raw column indices, in output order. Retain every field needed downstream. |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` (`bool`) | `false` | `[extcat].use_explicit_coordinate_columns` | Use explicit raw RA/Dec positions instead of header discovery; set this flag explicitly. |
| `EXTCAT_RA_COLUMN_ONE_BASED` (`std::size_t`) | `5` | `[extcat].ra_column` | Positive one-based raw RA column. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` (`std::size_t`) | `6` | `[extcat].dec_column` | Positive one-based raw Dec column. |
| `EXTCAT_MAG_G_COLUMN_ONE_BASED` (`std::size_t`) | `7` | `[extcat].mag_g_column` | One-based raw g magnitude column; `0` marks this band absent. |
| `EXTCAT_MAG_R_COLUMN_ONE_BASED` (`std::size_t`) | `9` | `[extcat].mag_r_column` | One-based raw r magnitude column; `0` marks this band absent. |
| `EXTCAT_MAG_I_COLUMN_ONE_BASED` (`std::size_t`) | `11` | `[extcat].mag_i_column` | One-based raw i magnitude column; `0` marks this band absent. |
| `EXTCAT_MAG_Z_COLUMN_ONE_BASED` (`std::size_t`) | `13` | `[extcat].mag_z_column` | One-based raw z magnitude column; `0` marks this band absent. |
| `EXTCAT_MAG_Y_COLUMN_ONE_BASED` (`std::size_t`) | `15` | `[extcat].mag_y_column` | One-based raw y magnitude column; `0` marks this band absent. |
| `EXTCAT_ZP_COLUMN_ONE_BASED` (`std::size_t`) | `17` | `[extcat].zp_column` | Positive one-based raw zp/photo-z column; not a magnitude zero point. |

## Catalog and calibration paths (`LensingConfig`)

Compiled definitions are in `config/pathconfig.hpp`. Runtime overrides do not modify header constants or automatically synchronize independent fields.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `ASTROMETRY_CAT` (`std::string`) | `"/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted"` | `[lensing].astrometry_cat` | Main-phase Gaia tile directory (supports Type 1 and Type 2 layouts). |
| `SOURCE_CAT` (`std::string`) | `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"` | `[lensing].source_cat` / `[extcat].output_directory` | Source-tile input and extcat publication default; the INI aliases update one runtime field. |
| `FLAT_PATH` (`std::string`) | `"/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014"` | `[lensing].flat_path` | Readable super-flat FITS directory when `include_flat = 1`. |
| `PSF_PATH` (`std::string`) | `"hahahaha"` | `[lensing].psf_path` | Placeholder: replace with a readable external PSF root before enabling `ext_psf = 1`. |

## Branches and geometry (`LensingConfig`)

Edit `config/LensingConfig.hpp`.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `DEFAULT_CHIP_COUNT` (`int`) | `62` | `[lensing].nmax_chip` | Positive effective CCD count. |
| `pixel_size` (`double`) | `0.2628` | `[lensing].pixel_size` | Positive angular scale in arcsec/pixel. |
| `chipnx` (`int`) | `2046` | `[lensing].chipnx` | Positive physical CCD width used by geometry-dependent consumers. Stage 1 reads actual FITS axes. |
| `chipny` (`int`) | `4094` | `[lensing].chipny` | Positive physical CCD height; keep consistent with camera geometry. |
| `ASTROMETRY_trivial` (`int`) | `0` | `[lensing].astrometry_trivial` | `0`: Gaia astrometry; `1`: deliberately use identity mapping. |
| `AstroCatType` (`int`) | `1` | `[lensing].astrometry_cat_type` | `1`: legacy large Gaia tiles; `2`: one-degree tiles. |
| `PROCESS_stage` (`int`) | `223092870` | `[lensing].process_stage` | Product of enabled primes `2,3,5,7,11,13,17,19,23`. Stage 9 (`23`) requires Stage 8 (`19`). |
| `include_FLAT` (`int`) | `0` | `[lensing].include_flat` | `0` off; `1` applies super-flat correction. |
| `include_Mask` (`int`) | `2` | `[lensing].include_mask` | Mask-input selector `0–3`. Review all active consumers before omitting DQ. |
| `include_BGsub` (`int`) | `1` | No | `0` off; `1` subtracts the fitted Science background. |
| `PreprocsType` (`int`) | `2` | No | `1`: historical F77 background/sigma estimators; `2`: current C++ estimators. |
| `ext_cat` (`int`) | `1` | `[lensing].ext_cat` | `0` disables external source catalog; `1` enables it. FD requires `1`. |
| `ext_PSF` (`int`) | `0` | `[lensing].ext_psf` | `0`: use exposure stars; `1`: read external PSF images. |
| `CCD_split` (`int`) | `2` | `[lensing].ccd_split` | `1`: whole CCD; `2`: amplifier halves. |
| `deblending` (`int`) | `1` | `[lensing].deblending` | `0` off; `1` enables source deblending. |
| `PSF_type` (`int`) | `1` | `[lensing].psf_type` | `1`: local polynomial; `2`: hybrid PSF. |
| `PSF_Ms` (`int`) | `0` | `[lensing].psf_ms` | `0` off; `1` enables PCA/multi-scale PSF reconstruction. |

## PSF selection (`LensingConfig`)

All entries in this section are compiled settings in `config/LensingConfig.hpp`. Rebuild after editing. Tune scientific cuts on representative data, not merely to make a run pass.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `npo` (`int`) | `64` | No | Exposure PSF sample count; positive. |
| `nstar_min` (`int`) | `npo * 3 / 2 = 96` | No | Derived minimum stars for exposure PSF fitting; change `npo` rather than this expression alone. |
| `npl` (`int`) | `10` | No | Local PSF polynomial coefficient count minus one. |
| `nstar_min_local` (`int`) | `16` | No | Minimum retained stars for a local fit. |
| `PsfGroupingType` (`int`) | `4` | No | `1` F77 selection, `2` threshold graph, `3` mutual KNN, `4` adaptive pair/fraction cuts. |
| `psf_exposure_min_candidates` (`int`) | `60` | No | Positive minimum exposure-wide PSF candidates. |
| `psf_count_pilot_clip_sigma` (`double`) | `3.0` | No | Positive pilot median/MAD clipping multiplier for integer star areas. |
| `psf_count_pilot_clip_iterations` (`int`) | `3` | No | Positive maximum pilot clipping passes. |
| `psf_count_zero_mad_quantile` (`double`) | `0.05` | No | Lower quantile for zero-MAD pilot; upper is `1-q`; require `0 <= q < 0.5`. |
| `psf_count_hist_range_sigma` (`double`) | `5.0` | No | Positive pilot-width multiplier for the local star-area search range; histogram bins remain fixed at two integer counts. |
| `psf_count_locus_sigma` (`double`) | `4.0` | No | Positive lower/upper count-MAD refinement multiplier; outer elbows may only widen final strict bounds. |
| `psf_count_locus_min_samples` (`int`) | `30` | No | Positive minimum samples for an integer star-area locus. |
| `psf_minchi_reference_fraction` (`double`) | `1.0 / 3.0` | No | Exposure top-size reference fraction in `(0,1]`. |
| `psf_minchi_reference_max_per_chip` (`int`) | `5` | No | Positive per-chip cap on reference stars. |
| `psf_minchi_sigma_cut` (`double`) | `4.0` | No | Positive minimum-chi rejection sigma multiplier. |
| `psf_knn_k` (`int`) | `20` | No | Positive number of retained neighbors for KNN grouping. |
| `psf_group_merge_ratio` (`double`) | `0.30` | No | Nonnegative secondary-group relative-size threshold. |
| `psf_group_merge_min_gaia` (`int`) | `1` | No | Nonnegative minimum Gaia matches in a merged group. |
| `psf_gaia_match_radius_pix` (`double`) | `2.0` | No | Positive Gaia match radius in pixels. |
| `psf_gaia_locus_min_matches` (`int`) | `5` | No | Positive minimum Gaia matches supporting the locus. |
| `psf_pair_chi_valid_peak_fraction` (`double`) | `0.3678794411714423216` | No | `exp(-1)`; strict relative-height threshold for a valid pair-chi peak; in `(0,1)`. |
| `psf_bad_fraction_valid_peak_fraction` (`double`) | `0.10` | No | Strict relative-height threshold for a valid bad-pair fraction peak; in `(0,1)`. |
| `psf_adaptive_elbow_search_height_fraction` (`double`) | `0.10` | No | Elbow candidates must lie below this fraction of the smoothed main peak; in `(0,1)`. |
| `psf_press_rejection_enabled` (`bool`) | `false` | No | Enable post-fit PRESS cleanup only after validating its effect on the target data. |
| `psf_press_sigma_cut` (`double`) | `4.0` | No | Positive standardized PRESS rejection sigma multiplier. |
| `psf_press_max_removals` (`int`) | `5` | No | Nonnegative maximum PRESS removals per chip. |
| `psf_loo_min_denom` (`double`) | `1.0e-6` | No | Minimum leave-one-out denominator in `(0,1)`; numerical guard. |
| `step_psf` (`int`) | `100` | No | Positive spatial PSF star sampling step in pixels. |

## Stamps and detection (`LensingConfig`)

All entries without an INI mapping are in `config/LensingConfig.hpp` and require a rebuild. Keep camera geometry, stamp dimensions and detector-specific thresholds consistent.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `ns` (`int`) | `64` | No | Positive even stamp/Fourier-grid side length in pixels; coordinate all stamp/FFT changes. |
| `nsns` (`int`) | `ns * ns = 4096` | No | Derived pixel count in a Science stamp. |
| `chip_margin` (`int`) | `8` | No | Nonnegative extra chip-edge extraction margin in pixels. |
| `ns_2` (`int`) | `ns / 2 = 32` | No | Derived half stamp side. |
| `nl_2` (`int`) | `ns_2 + chip_margin = 40` | No | Derived half expanded extraction side. |
| `nl` (`int`) | `nl_2 * 2 = 80` | No | Derived expanded extraction side. |
| `flag_thresh` (`int`) | `3` | No | Nonnegative maximum accepted source extraction flag. |
| `chip_edge_margin` (`int`) | `chip_margin = 8` | No | Derived alias for chip-edge checks. |
| `dz_thresh` (`double`) | `0.1` | No | Nonnegative redshift tolerance for external-catalog matching. |
| `n_user_max` (`int`) | `500` | No | Maximum flux-ranked detections passed to astrometric pattern matching; a scientific selection limit, not a storage-capacity limit. |
| `ngal_max` (`int`) | `2000` | No | Initial galaxy-vector reservation hint; not a hard source limit. |
| `nstar_max` (`int`) | `1000` | No | Initial star-vector reservation hint; not a hard star limit. |
| `src_npara` (`int`) | `12` | No | Stage-6 source-marker metadata field count; keep consumers consistent. |
| `npd` (`int`) | `33` | No | PU astrometric distortion coefficient count; change only with a coordinated model change. |
| `source_thresh` (`double`) | `2.0` | No | Positive source-detection S/N threshold. |
| `core_thresh` (`double`) | `4.0` | No | Positive source-core threshold. |
| `n_neighbor` (`int`) | `5` | No | Positive neighbor count for deblending. |
| `area_max` (`int`) | `ns * ns = 4096` | No | Derived maximum connected source area. |
| `area_thresh` (`int`) | `6` | No | Positive minimum connected source area in pixels. |
| `size_fit_rmax` (`int`) | `4` | No | Maximum radius for curvature-size fitting; validate numerical/scientific changes. |
| `point_stat_beta` (`double`) | `0.10` | No | Survey-wide fixed point-source Fourier-power template shape; not fitted per source. |
| `point_stat_k_frac` (`double`) | `0.90` | No | Fourier-radius fraction used in point-source statistics. |
| `point_stat_eps` (`double`) | `1.0e-20` | No | Positive denominator floor for point statistics. |
| `point_stat_min_corr` (`double`) | `1.0e-6` | No | Minimum accepted point-template correlation. |
| `SNR_PSF` (`double`) | `100.0` | No | Positive minimum PSF-star signal-to-noise ratio. |
| `saturation_thresh` (`double`) | `25000.0` | No | Saturated-pixel threshold in detector counts; review gain and data calibration. |
| `gal_smooth` (`int`) | `0` | `[lensing].gal_smooth` | Nonnegative galaxy-stamp smoothing selector. |
| `star_smooth` (`int`) | `2` | `[lensing].star_smooth` | Nonnegative star-stamp smoothing selector. |

## Background estimation (`LensingConfig`)

Compiled settings in `config/LensingConfig.hpp`; rebuild after editing.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `nct` (`int`) | `12` | No | Number of background polynomial rectangle coefficients; positive. |
| `ncx` (`int`) | `3` | No | Background rectangle x count; keep consistent with `nct`. |
| `blocksize` (`int`) | `200` | No | Target side length of background blocks in pixels. |
| `bg_rough_grid_x` (`int`) | `32` | No | Positive rough-background grid column count. |
| `bg_rough_grid_y` (`int`) | `32` | No | Positive rough-background grid row count. |
| `bg_min_block_pixels` (`int`) | `1000` | No | Minimum valid pixels per background block. |
| `bg_min_clipped_pixels` (`int`) | `200` | No | Minimum surviving pixels after block clipping. |
| `bg_min_valid_frac` (`double`) | `0.25` | No | Minimum valid fraction per block, in `(0,1]`. |
| `bg_clip_low` (`double`) | `4.0` | No | Positive lower-side background clipping sigma. |
| `bg_clip_high` (`double`) | `2.5` | No | Positive upper-side background clipping sigma. |
| `bg_fit_clip_sigma` (`double`) | `3.0` | No | Positive plane-fit clipping sigma. |
| `bg_fit_max_iter` (`int`) | `4` | No | Maximum clipping iterations; nonnegative. |
| `bg_min_fit_factor` (`int`) | `3` | No | Positive minimum sample-to-coefficient factor for a fit. |

## Noise products (`LensingConfig`)

Compiled settings in `config/LensingConfig.hpp`. Types 1/2 are real-space stamps transformed downstream; Type 3 is already power.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `NstampType` (`int`) | `2` | No | `1` deterministic F77 real-space blank stamp, `2` QC/random physical blank stamp, `3` local covariance Fourier power. |
| `noise_sigma_ratio_min` (`double`) | `0.80` | No | Positive minimum blank-to-source sigma ratio. |
| `noise_sigma_ratio_max` (`double`) | `1.25` | No | Maximum blank-to-source sigma ratio; greater than minimum. |
| `noise_mad_ratio_min` (`double`) | `0.70` | No | Positive minimum blank-to-source MAD ratio. |
| `noise_mad_ratio_max` (`double`) | `1.30` | No | Maximum blank-to-source MAD ratio; greater than minimum. |
| `noise_tail_sigma` (`double`) | `2.5` | No | Positive tail-count sigma threshold. |
| `noise_max_tail_fraction` (`double`) | `0.05` | No | Maximum blank-stamp tail fraction, in `[0,1]`. |
| `noise_max_mask_fraction` (`double`) | `0.02` | No | Maximum masked fraction of a blank stamp, in `[0,1]`. |
| `noise_region_size` (`int`) | `192` | No | Positive even outer local-noise square side; must exceed inner size. |
| `noise_inner_size` (`int`) | `96` | No | Even central exclusion side; at least `nl`; outer-minus-inner must be even. |
| `noise_plane_min_valid_fraction` (`double`) | `0.30` | No | Minimum valid shell fraction for plane fitting, in `(0,1]`. |
| `noise_cov_padding_factor` (`double`) | `2.0` | No | Positive covariance FFT padding multiplier. |
| `noise_cov_fft_size` (`int`) | `ceil(noise_region_size * noise_cov_padding_factor) = 384` | No | Derived FFT side; must be at least `2*noise_region_size-1` for linear autocorrelation. |
| `noise_cov_max_lag` (`int`) | `63` | No | Maximum retained signed covariance lag; `0 <= lag < noise_region_size`. |
| `noise_cov_min_valid_pixels` (`int`) | `4096` | No | Positive minimum covariance-mask pixel count. |
| `noise_cov_min_pair_fraction` (`double`) | `0.50` | No | Minimum valid pair-count fraction for a lag, in `(0,1]`. |
| `noise_cov_sigma_ratio_min` (`double`) | `0.80` | No | Positive lower accepted covariance sigma ratio. |
| `noise_cov_sigma_ratio_max` (`double`) | `1.25` | No | Upper accepted covariance sigma ratio; above the lower bound. |
| `noise_cov_max_negative_fraction` (`double`) | `0.25` | No | Maximum negative-power fraction, in `[0,1]`. |
| `noise_cov_imag_tolerance` (`double`) | `1.0e-10` | No | Nonnegative imaginary FFT residual tolerance. |

## Sigma estimation (`LensingConfig`)

Compiled settings in `config/LensingConfig.hpp`; rebuild after edits. Treat estimator/calibration changes as scientific changes and validate them on representative data.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `sig_blocksize` (`int`) | `200` | No | Positive estimator block side length in pixels. |
| `sig_block_max` (`int`) | `sig_blocksize * sig_blocksize = 40000` | No | Derived maximum pixels per block. |
| `sig_max_blocks` (`int`) | `2048` | No | Positive maximum sampled blocks. |
| `sig_min_block_pixels` (`int`) | `1000` | No | Minimum valid pixels in one block. |
| `sig_min_block_triples` (`int`) | `1000` | No | Minimum valid pixel triples per block. |
| `sig_min_blocks` (`int`) | `4` | No | Minimum blocks for a plane fit. |
| `sig_hist_nbin` (`int`) | `256` | No | Positive mode-finding histogram bin count. |
| `sig_hist_range` (`double`) | `6.0` | No | Positive histogram range in sigma units. |
| `sig_min_mode_count` (`int`) | `500` | No | Minimum samples defining the mode. |
| `sig_min_lower_count` (`int`) | `1000` | No | Minimum lower-side samples. |
| `sig_lower_quantile` (`double`) | `0.3173105` | No | Lower-side width quantile in `(0,1)`; calibration constant. |
| `sig_clip_k` (`double`) | `3.0` | No | Positive symmetric clipping sigma multiplier. |
| `sig_rdil` (`int`) | `2` | No | Positive pixel-sampling stride. |
| `sig_clip_niter` (`int`) | `2` | No | Nonnegative clipping iteration count. |
| `sig_min_fit_triples` (`int`) | `1000` | No | Minimum triples in the final fit. |
| `sig_min_fit_frac` (`double`) | `0.20` | No | Minimum retained fit fraction, in `(0,1]`. |
| `sig_median_ratio` (`double`) | `1.2678405` | No | Positive median-to-sigma calibration factor. |
| `sig_plane_min` (`double`) | `1.0e-8` | No | Positive lower floor for the sigma plane. |
| `sig_max_plane_ratio` (`double`) | `4.0` | No | Maximum plane variation ratio; at least `1`. |
| `sig_pivot_min` (`double`) | `1.0e-8` | No | Positive numerical pivot floor. |
| `sig_scale_s1` (`double`) | `0.673475` | No | Stage-1 noise-calibration candidate. |
| `sig_scale_s2` (`double`) | `1.027786` | No | Stage-2 noise-calibration value. |
| `sig_scale` (`double`) | `sig_scale_s2 = 1.027786` | No | Active calibration selector; produces the published `2*sigma^2` convention. |

## Catalog indices (`LensingConfig`)

Compiled in `config/LensingConfig.hpp`. These derived schema positions must remain consistent with writers and readers. Do not edit them to change raw catalog column mapping; use `[extcat]` instead.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `iid` (`int`) | `1 - 1 = 0` | No | PSF polynomial chi-square; zero-based index, not a selection threshold. |
| `ipixx` (`int`) | `2 - 1 = 1` | No | Source center x; zero-based index, not a selection threshold. |
| `ipixy` (`int`) | `3 - 1 = 2` | No | Source center y; zero-based index, not a selection threshold. |
| `isig` (`int`) | `4 - 1 = 3` | No | Local noise sigma; zero-based index, not a selection threshold. |
| `istar` (`int`) | `5 - 1 = 4` | No | Available PSF-star count; zero-based index, not a selection threshold. |
| `ipeak` (`int`) | `5 - 1 = 4` | No | Historical peak alias; zero-based index, not a selection threshold. |
| `i_imax` (`int`) | `6 - 1 = 5` | No | Peak x; zero-based index, not a selection threshold. |
| `i_jmax` (`int`) | `7 - 1 = 6` | No | Peak y; zero-based index, not a selection threshold. |
| `ih_flux` (`int`) | `8 - 1 = 7` | No | Half-light flux; zero-based index, not a selection threshold. |
| `ih_area` (`int`) | `9 - 1 = 8` | No | Source area; zero-based index, not a selection threshold. |
| `iflag` (`int`) | `10 - 1 = 9` | No | Quality flag; zero-based index, not a selection threshold. |
| `iPSF` (`int`) | `11 - 1 = 10` | No | Local PSF size; zero-based index, not a selection threshold. |
| `iSNR_F` (`int`) | `12 - 1 = 11` | No | Fourier S/N; zero-based index, not a selection threshold. |
| `ira` (`int`) | `13 - 1 = 12` | No | Right Ascension; zero-based index, not a selection threshold. |
| `idec` (`int`) | `14 - 1 = 13` | No | Declination; zero-based index, not a selection threshold. |
| `igf1` (`int`) | `15 - 1 = 14` | No | Field-distortion g1; zero-based index, not a selection threshold. |
| `igf2` (`int`) | `16 - 1 = 15` | No | Field-distortion g2; zero-based index, not a selection threshold. |
| `ig1` (`int`) | `17 - 1 = 16` | No | Fourier_Quad g1 estimator; zero-based index, not a selection threshold. |
| `ig2` (`int`) | `18 - 1 = 17` | No | Fourier_Quad g2 estimator; zero-based index, not a selection threshold. |
| `ide` (`int`) | `19 - 1 = 18` | No | Shear response estimator; zero-based index, not a selection threshold. |
| `ih1` (`int`) | `20 - 1 = 19` | No | Higher-order h1 estimator; zero-based index, not a selection threshold. |
| `ih2` (`int`) | `21 - 1 = 20` | No | Higher-order h2 estimator; zero-based index, not a selection threshold. |
| `icos2` (`int`) | `22 - 1 = 21` | No | Spin-2 cosine; zero-based index, not a selection threshold. |
| `isin2` (`int`) | `23 - 1 = 22` | No | Spin-2 sine; zero-based index, not a selection threshold. |
| `iparity` (`int`) | `24 - 1 = 23` | No | WCS parity; zero-based index, not a selection threshold. |
| `igalsizeT` (`int`) | `25 - 1 = 24` | No | Curvature galaxy size; zero-based index, not a selection threshold. |
| `ipsfsizeT` (`int`) | `26 - 1 = 25` | No | Curvature PSF size; zero-based index, not a selection threshold. |
| `idelta_chi2` (`int`) | `27 - 1 = 26` | No | Point-template delta chi-square; zero-based index, not a selection threshold. |
| `iorth_ext` (`int`) | `28 - 1 = 27` | No | Orthogonal point-template extension; zero-based index, not a selection threshold. |
| `shear_cat_ncols` (`int`) | `iorth_ext + 1 = 28` | No | Derived Stage-7 row width. |
| `expo_cat_ncols` (`int`) | `shear_cat_ncols + 1 = 29` | No | Derived per-exposure suffix width including `Chi2`. |
| `ichi2` (`int`) | `shear_cat_ncols = 28` | No | Zero-based exposure chi-square index within the pipeline suffix. |

## Calibration and PCA (`LensingConfig`)

Compiled in `config/LensingConfig.hpp`; rebuild after any deliberate change.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `g1_c` (`double`) | `0.0` | No | Additive field-distortion g1 correction; recalibrate for another dataset/band. |
| `g2_c` (`double`) | `0.0` | No | Additive field-distortion g2 correction; recalibrate for another dataset/band. |
| `chi2_thresh` (`double`) | `0.1` | No | Nonnegative maximum exposure PSF chi-square; distinct from `FDConfig::chi2_thresh`. |
| `rescale_size` (`double`) | `1.2` | No | Positive target PSF residual rescaling size for `PSF_Ms = 1`. |
| `procs_pn` (`int`) | `40` | No | Positive MPI ranks per PCA scheduling group; coordinate with resource allocation. |
| `work_pn` (`int`) | `10` | No | Positive concurrent PCA workers per group. |
| `nblocks` (`int`) | `2` | No | Positive PCA spatial block count per CCD axis. |
| `n_pcs` (`int`) | `100` | No | Positive maximum PCA principal components. |
| `npp6th` (`int`) | `28` | No | Sixth-degree two-dimensional polynomial term count; derived from the basis. |
| `pca_negative_eigenvalue_threshold` (`double`) | `-1.0e-5` | No | Nonpositive cutoff for invalid PCA eigenvalues. |
| `pi` (`double`) | `3.14159265358979323846` | No | Mathematical pi; do not tune. |
| `arc_convert` (`double`) | `pi / 180.0` | No | Derived degrees-to-radians conversion; do not tune. |

## Spatial rearrangement (`ProcessRearrConfig`)

Compiled in `config/ProcessRearrConfig.hpp`; fixed names are physically defined in `config/pathconfig.hpp`. Every entry requires rebuilding after an edit.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `SKY_GRID_DEGREES` (`double`) | `0.1` | No | Positive sky tile width in degrees; coordinate grid count changes. |
| `RA_BIN_COUNT` (`int`) | `3600` | No | Full-sky RA bin count consistent with `SKY_GRID_DEGREES`. |
| `DEC_BIN_COUNT` (`int`) | `1800` | No | Full-sky Dec bin count consistent with `SKY_GRID_DEGREES`. |
| `SKY_TILE_COUNT` (`std::size_t`) | `RA_BIN_COUNT * DEC_BIN_COUNT = 6480000` | No | Derived total full-sky tile count. |
| `TARGET_SUBCAT_ROWS` (`std::uint64_t`) | `500000` | No | Positive target partition row count; controls file size and memory use. |
| `SKIP_DIRECTORY_NAME` (`std::string_view`) | `"Large_Field"` | No | Directory name excluded from scans. |
| `SUBCAT_PREFIX` (`std::string_view`) | `"subcat_"` | No | Published partition filename prefix. |
| `SUBCAT_EXTENSION` (`std::string_view`) | `".cat"` | No | Published partition filename extension. |
| `SUBCAT_ID_WIDTH` (`int`) | `6` | No | Positive minimum zero-padded partition-ID width. |
| `SUMMARY_FILENAME` (`std::string_view`) | `"catalog_summary.txt"` | No | Published summary filename. |
| `OUTPUT_PRECISION` (`int`) | `10` | No | Positive significant digits in catalog output. |
| `SUMMARY_PRECISION` (`int`) | `4` | No | Nonnegative decimal places in summary bounds. |
| `SKIP_MISSING_CATALOGS` (`bool`) | `true` | No | Continue past missing catalogs; `false` enforces strict completeness. |
| `SKIP_MALFORMED_ROWS` (`bool`) | `true` | No | Continue past malformed rows; `false` enforces strict row validity. |

Catalog widths and offsets are resolved from the effective runtime schema, not configured in this header. Complete rows contain the external prefix, `EXPO_NUM`, `ccD_NUM`, 28 Stage-7 fields and exposure `Chi2`. The default is **18 + 2 + 29 = 49 columns** (`EXPO_NUM=18`, `ccD_NUM=19`, source base 20, `Chi2=48`, all zero-based). With `ext_cat = 0` the catalog has 31 columns. Regenerate older incompatible rows before rearrangement or FD processing.

## Output layout (`OutputLayout`)

Defined in `config/pathconfig.hpp`. They contain relative names below the dataset root, not independent input/output roots. One dataset-root mount normally covers them all. Changing names requires a rebuild and coordinated product regeneration.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `NON_CHIP_BASE_DIRECTORIES` (`std::array<const char*, 16>`) | `16 paths below` | No | Complete fixed dataset directories without a chip suffix. |
| `CHIP_PRODUCT_DIRECTORIES` (`std::array<const char*, 16>`) | `16 paths below` | No | Complete fixed per-chip product-directory prefixes. |

| `NON_CHIP_BASE_DIRECTORIES` | `CHIP_PRODUCT_DIRECTORIES` |
| --- | --- |
| `science` | `stamps/Norm` |
| `dqmask` | `stamps/cat_Orig` |
| `stamps` | `stamps/dat_StarCanInfo` |
| `expolists` | `stamps/fits_StarCan` |
| `result` | `stamps/fits_StarCanN` |
| `stamps/dat_StarInfo` | `stamps/fits_StarCanP` |
| `stamps/svg_StarLocus` | `stamps/dat_SrcInfo` |
| `stamps/fits_StarP` | `stamps/fits_Src` |
| `stamps/fits_PsfSrc` | `stamps/fits_Noise` |
| `stamps/dat_ExpoInfo` | `stamps/fits_SrcP` |
| `stamps/dat_StarComp` | `stamps/dat_PsfFit` |
| `stamps/dat_Rescale` | `stamps/fits_PsfLocal` |
| `stamps/dat_Pcs` | `stamps/dat_Shear` |
| `stamps/dat_StarCompV2` | `stamps/dat_StarXY` |
| `astrometry/Head` | `stamps/fits_PsfResi` |
| `astrometry/dat_Chk` | `astrometry/dat_Astro` |

## Field distortion (`FDConfig`)

All entries are compiled in `config/FDConfig.hpp` and require rebuilding. `StaticMode` is the enum with the three enumerators described in the first row. There is no `[fd]` INI section. Runtime geometry and external-catalog magnitude mapping come from `[lensing]` and `[extcat]`.

| Parameter (type) | Compiled default | INI override | Meaning and constraints |
| --- | --- | --- | --- |
| `FD_STATIC_MODE` (`StaticMode`) | `StaticMode::PDF_SIGMA` | No | `PDF_SIGMA`: PDF chi-square estimate with analytic sigma; `PDF_JACK`: PDF estimate with jackknife sigma; `SWSE_JACK`: ratio estimate with jackknife sigma. |
| `FD_USE_PDF_STATIS` (`bool`) | `mode is PDF_SIGMA or PDF_JACK` | No | Derived PDF sign-test switch; change `FD_STATIC_MODE`, not this expression. |
| `FD_USE_JACKKNIFE` (`bool`) | `mode is PDF_JACK or SWSE_JACK` | No | Derived jackknife-uncertainty switch. |
| `FD_USE_SWSE_DATA` (`bool`) | `mode is SWSE_JACK` | No | Derived SWSE-data-model switch. |
| `FD_PER_EXPOSURE_STAR_BAR` (`bool`) | `false` | No | `false`: single global N-to-1 star bar; `true`: per-exposure N-to-N star bars. |
| `nmax_per_core` (`int`) | `20000000` | No | Maximum sources reserved per MPI rank; review per-rank memory use. |
| `fd_num` (`int`) | `21` | No | Positive field-distortion spatial bin count. |
| `PDF_BINS` (`int`) | `4` | No | Positive equal-probability inner-bin count. |
| `gf_lim` (`float`) | `0.0015` | No | Positive spatial field-distortion half-range, from `-gf_lim` to `+gf_lim`. |
| `NMAX` (`int`) | `200` | No | Positive fine-grid sampling count. |
| `MAX_DUP` (`int`) | `5` | No | Positive maximum duplicate measurements. |
| `N_jack` (`int`) | `50` | No | Positive jackknife region count. |
| `nmax_total` (`int`) | `1000000` | No | Positive maximum sources used in k-means. |
| `Km_iter` (`int`) | `100` | No | Positive k-means iteration count. |
| `snrfcut` (`float`) | `4.0` | No | Minimum Fourier signal-to-noise ratio; nonnegative. |
| `snrlow` (`float`) | `0.0` | No | Optional lower source-S/N bound; zero disables it. |
| `snrhigh` (`float`) | `0.0` | No | Optional upper source-S/N bound; zero disables it; keep above the lower bound when both active. |
| `starcut` (`float`) | `20.0` | No | Point-source size threshold; recalibrate with the data. |
| `chi2_thresh` (`float`) | `0.01` | No | Maximum exposure chi-square; separate from the `LensingConfig` threshold. |
| `flagcut` (`float`) | `0.0` | No | Maximum accepted source quality flag. |
| `imaxcut` (`float`) | `64.0` | No | Upper source-peak x coordinate cut; coordinate with stamp geometry. |
| `jmaxcut` (`float`) | `64.0` | No | Upper source-peak y coordinate cut; coordinate with stamp geometry. |
| `zplow` (`float`) | `0.0` | No | Minimum photometric redshift. |
| `zphigh` (`float`) | `3.0` | No | Maximum photometric redshift; greater than `zplow`. |
| `r_half_thresh` (`float`) | `0.0` | No | Optional half-light-radius threshold; zero disables it. |
| `star_bar_mltp` (`float`) | `3.0` | No | Positive stellar-locus sigma multiplier. |
| `psf_chi2_mltp` (`float`) | `3.0` | No | Positive PSF chi-square sigma multiplier. |
| `n_size_bins` (`int`) | `100` | No | Positive stellar-size histogram bin count. |
| `n_mag_bins` (`int`) | `20` | No | Positive magnitude histogram bin count. |
| `size_min` (`float`) | `-2.0` | No | Lower size-histogram bound. |
| `size_max` (`float`) | `2.0` | No | Upper size-histogram bound, greater than `size_min`. |
| `mag_min_val` (`float`) | `10.0` | No | Lower magnitude-histogram bound. |
| `mag_max_val` (`float`) | `30.0` | No | Upper magnitude-histogram bound, greater than `mag_min_val`. |
| `min_bin_count` (`int`) | `100` | No | Positive minimum sample count for a usable bin. |
| `peak_match_tol` (`float`) | `0.05` | No | Positive stellar-peak matching tolerance. |
| `min_concentration` (`float`) | `0.6` | No | Minimum stellar-locus concentration fraction in `[0,1]`. |
| `star_phy_min` (`float`) | `-0.5` | No | Lower physical stellar-size bound. |
| `star_phy_max` (`float`) | `0.2` | No | Upper physical stellar-size bound, above `star_phy_min`. |
| `stage1_snr` (`float`) | `40.0` | No | Nonnegative S/N threshold for histogram accumulation. |
| `stage2_snr` (`float`) | `0.0` | No | Optional second-pass S/N threshold; zero disables it. |
| `init_win_active` (`float`) | `0.1` | No | Positive initial size window for an active exposure. |
| `init_win_fallback` (`float`) | `0.15` | No | Positive fallback initial size window. |
| `default_s_init` (`float`) | `0.5` | No | Default initial stellar-size center. |
| `clip_nsigma` (`float`) | `3.0` | No | Positive iterative clipping sigma multiplier. |
| `min_clip_limit` (`float`) | `0.015` | No | Positive minimum clipping half-width. |
| `default_s_std` (`float`) | `0.05` | No | Positive default stellar-size scatter. |
| `fallback_scut_default` (`float`) | `0.6` | No | Fallback stellar-size cut. |
| `bad_ccds` (`int[]`) | `{-1}` | No | Backing CCD-ID array. With `n_bad_ccds=0` no CCD is excluded; the `-1` entry is inactive storage, not an active detector ID. |
| `n_bad_ccds` (`int`) | `0` | No | Active count of excluded CCD IDs. Keep this explicit count within the actual `bad_ccds` array length; it is not automatically inferred. |
| `chip_mask_edge` (`int`) | `50` | No | Nonnegative symmetric edge mask in pixels; inclusive bounds derive as `[edge, chipnx-edge]` and `[edge, chipny-edge]`. |

Per-exposure statistics use serialized `EXPO_NUM` identities, not file order or terminal `Chi2`. The runtime allocation follows the MPI-global maximum valid exposure ID; there is no separate configured exposure cap.

## Runtime selection and coupled values

The five INI sections are `[process]`, `[astrocat]`, `[extcat]`, `[init]` and `[lensing]`. All supported keys appear in the tables above. A sample with all keys is provided in `cpp_Standard/pipeline.example.ini`; its entries are commented until you opt in.

| Relationship | User action |
| --- | --- |
| **Main stages** | Select prime factors in `process_stage`; include `19` whenever including `23`. |
| **Gaia producer and consumer** | Set `[astrocat].output_directory` and `[lensing].astrometry_cat` together when reusing newly produced tiles, and set `[lensing].astrometry_cat_type = 2`. |
| **External source tiles** | `[extcat].output_directory` and `[lensing].source_cat` address the same runtime field. Set only one, or make both identical. |
| **Projection** | Set `use_explicit_columns = true` together with `input_columns`; set `use_explicit_coordinate_columns = true` for explicit RA/Dec positions. |
| **Downstream-only lists** | An explicit `expo_list` is an existing file. Without it, the program derives one list per dataset from `output_root`. |
| **Main output** | Dataset roots come from the listed image layout; changing only `output_root` does not relocate an already prepared dataset. |
| **Rebuild** | Any header edit changes compiled settings only after a rebuild. An explicit INI still overrides mapped values. |
