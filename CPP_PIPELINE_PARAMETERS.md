# C++ Pipeline Parameter Reference

[Home](README.md) · [C++ user guide](CPP_GUIDE.md) · [Data preparation](DATA_PREPARATION.md)

This is the complete reference for the user configuration headers of Standard
and Lite, including Lite's `Initialize.hpp`, runtime INI fields, fixed scientific
settings, derived constants, catalog indices and output-layout names. Defaults
below describe the source at `9925f4f`; local header edits take precedence as the
compiled baseline. Function-local implementation constants are not user settings.

## How to use this reference

**Configuration order: compiled header defaults → explicitly selected INI.**
The executable accepts only `--config PATH` or `--config=PATH`. No argument means
compiled configuration only: it does not search for `pipeline.ini` and does not
read a configuration-path environment variable.

| Entry | Effect |
| --- | --- |
| Parameter has an INI mapping | Change that INI key without rebuilding. |
| INI override is **No** | Edit the indicated header and rebuild. |
| INI entry is omitted or commented (`# key =`) | Retain the compiled value. |
| Active string entry is empty (`key =`) | Clear the value; a downstream fallback may exist for that particular field. |
| Active number or Boolean is empty | Configuration error; it does not inherit the compiled value. |
| Active list is empty | Clear token/dataset lists; column projections must remain nonempty. |
| Default is **same** | Same as the Standard default on that row. |
| Default is **Not present** | The symbol/branch does not exist in that variant; do not add its INI key. |

All edits to compiled headers require rebuilding, including edits to defaults
that also have an INI mapping. Derived values are documented for completeness:
change their source setting, not an index or size independently. Boolean values
accept `true/false`, `1/0`, `yes/no` and `on/off`. INI section and key names are
case-insensitive; path values and filename filters are not. Lists are comma
separated; quote strings containing `#`, `;` or delimiters.

## Contents

- [Workflow](#workflow-processconfig), [archive initialization](#archive-initialization-initconfig), [catalog preparation](#gaia-preparation-astrocatconfig), [external catalog schema](#external-catalog-preparation-extcatconfig).
- [Catalog paths](#catalog-and-calibration-paths-lensingconfig), [branches and geometry](#branches-and-geometry-lensingconfig), [PSF selection](#psf-selection-lensingconfig), [stamps and detection](#stamps-and-detection-lensingconfig), [background](#background-estimation-lensingconfig), [noise](#noise-products-lensingconfig), [sigma estimation](#sigma-estimation-lensingconfig), [catalog indices](#catalog-indices-lensingconfig), [calibration and PCA](#calibration-and-pca-lensingconfig).
- [Rearrangement](#spatial-rearrangement-processrearrconfig), [output directories](#output-layout-outputlayout), [FD settings](#field-distortion-fdconfig), [Lite initialization aliases](#lite-initializehpp).

## Workflow (`ProcessConfig`)
Edit `config/ProcessConfig.hpp` for phase defaults and `config/pathconfig.hpp` for paths. Lite phase defaults and `EXPO_LIST` originate in `Initialize.hpp`. At least one phase must be enabled. Phases execute in the fixed order shown in the C++ guide.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `RUN_PROCESS_ASTROCAT` (`bool`) | `false` | same | `[process].run_process_astrocat` | Prepare Gaia tiles once, before other phases. |
| `RUN_PROCESS_EXTCAT` (`bool`) | `false` | same | `[process].run_process_extcat` | Prepare external source tiles once, before dataset phases. |
| `RUN_PROCESS_INIT` (`bool`) | `true` | same | `[process].run_process_init` | Extract Science/DQ archives and publish dataset lists. |
| `RUN_PROCESS_MAIN` (`bool`) | `true` | same | `[process].run_process_main` | Run the numerical measurement stages. |
| `RUN_PROCESS_REARR` (`bool`) | `true` | `false` | `[process].run_process_rearr` | Spatially partition per-exposure catalogs. |
| `RUN_PROCESS_FD` (`bool`) | `true` | `false` | `[process].run_process_fd` | Run scientific field-distortion analysis; this is not a container diagnostic. |
| `EXPO_LIST` (`const char*`) | `""` | same | `[process].expo_list` | Top-level exposure-list file for downstream-only processing. Empty derives output_root/expo_<target>.list. One explicit list cannot serve multiple main/rearrangement datasets. |
| `REARR_OUTPUT_DIRECTORY` (`const char*`) | `"baked"` | same | `[process].rearr_output_directory` | Rearranged catalog directory: absolute, or relative to the configured base/dataset root. Empty selects the base itself. |
| `REARR_OUTPUT_BASE_DIRECTORY` (`const char*`) | `""` | same | `[process].rearr_output_base_directory` | Base for relative rearrangement output. Empty derives the dataset root from input image paths. |
| `REARRANGED_EXPO_LIST_FILENAME` (`const char*`) | `"cat_gband_ori.list"` | same | `[process].rearranged_expo_list_filename` | Published rearranged-list filename. It is a file name, not a directory. |
| `REARRANGED_EXPO_LIST_DIRECTORY` (`const char*`) | `""` | same | `[process].rearranged_expo_list_directory` | Directory for the published list. Empty uses the input exposure-list parent. |
| `FD_EXPO_LIST` (`const char*`) | `""` | same | `[process].fd_expo_list` | Explicit FD input-list file. Empty selects the rearranged-list directory and filename. |
| `FD_OUTPUT_DIRECTORY` (`const char*`) | `"fdout"` | same | `[process].fd_output_directory` | FD output directory: absolute, or relative to the FD base/dataset root. |
| `FD_OUTPUT_BASE_DIRECTORY` (`const char*`) | `""` | same | `[process].fd_output_base_directory` | Base for relative FD output. Empty uses the dataset root. |

## Archive initialization (`InitConfig`)
Standard: `config/InitConfig.hpp` and `config/pathconfig.hpp`. Lite: ordinary dataset, archive and path defaults originate in `Initialize.hpp`; advanced policy settings remain in `config/InitConfig.hpp`. Compiled absolute paths are site-specific examples, not portable defaults.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `DatasetSpec::target` (`std::string`) | `"gband"` | same | `[init].datasets` | One nonempty output-directory component, not ., .. or a path; unique within the dataset list. |
| `DatasetSpec::prefix` (`std::string`) | `"c4d_"` | same | `[init].datasets` | Nonempty archive basename prefix paired with target. |
| `SCIENCE_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/g"` | same | `[init].science_root` | Readable Science archive root. |
| `DQ_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"` | same | `[init].dq_root` | Readable DQ archive root. Check both initialization and selected mask consumers before omitting DQ. |
| `OUTPUT_ROOT` (`const char*`) | `"/lustre/home/acct-phyzj/share/DES/g_band_v1"` | same | `[init].output_root` | Writable parent for dataset trees and published exposure lists. |
| `DATASETS` (`std::vector<DatasetSpec>`) | `{{"gband", "c4d_"}}` | same | `[init].datasets` | Comma-separated TARGET:PREFIX pairs; enabled dataset phases run sequentially for each pair. |
| `CONTAINS` (`std::vector<std::string>`) | `{"v1"}` | same | `[init].contains` | Case-sensitive archive basename tokens, OR matched. An active empty list removes this filter. |
| `ARCHIVE_SUFFIX` (`const char*`) | `".fits.fz"` | same | No | Case-sensitive archive suffix, removed when forming exposure stems. |
| `CCDNUM_KEYWORD` (`const char*`) | `"CCDNUM"` | same | No | FITS keyword holding physical detector identity; do not substitute HDU sequence numbers. |
| `DQ_STEM_REPLACE_FROM` (`const char*`) | `"ood"` | same | No | Substring replaced throughout a DQ stem to match the Science exposure. Empty disables replacement. |
| `DQ_STEM_REPLACE_TO` (`const char*`) | `"ooi"` | same | No | Replacement substring; change together with the source token. |
| `EXISTING` (`const char*`) | `"fail"` | same | `[init].existing` | fail, resume or overwrite. Choose deliberately when outputs already exist. |
| `F77_MAX_PATH` (`int`) | `0` | same | `[init].f77_max_path` | Nonnegative generated-path length limit; 0 disables this compatibility guard. |

## Tile naming
Both definitions are in `config/pathconfig.hpp`; Lite derives their values from `Initialize.hpp`. The full grammar is `<prefix>RA_DDD_DDD_Dec_[pm]DD_[pm]DD.dat`. Rebuild and regenerate or rename tiles after changing a prefix. Old-prefix files are not treated as current generated products.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `AstroCatConfig::ASTROMETRY_TILE_PREFIX` (`std::string_view`) | `"astra_"` | same | No | Gaia one-degree tile prefix; excludes the fixed RA_ token. |
| `ExtCatConfig::SOURCE_CAT_TILE_PREFIX` (`std::string_view`) | `"extern_"` | same | No | External source one-degree tile prefix; excludes RA_. |

## Gaia preparation (`AstroCatConfig`)
Edit `config/AstroCatConfig.hpp` and `config/pathconfig.hpp`; Lite input and header defaults originate in `Initialize.hpp`. Produces Type-2 one-degree tiles. Standard must select `astrometry_cat_type = 2` to consume them; Lite already uses that fixed layout. Set the main-phase Gaia path independently, or use the container generator's linked-path behavior.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `ASTROCAT_INPUT_DIRECTORY` (`const char*`) | `""` | same | `[astrocat].input_directory` | Readable flat directory of raw RA/Dec files; only direct regular children are read. |
| `ASTROCAT_OUTPUT_DIRECTORY` (`std::string (Standard); const char* (Lite)`) | `LensingConfig::ASTROMETRY_CAT` | `""` | `[astrocat].output_directory` | Output directory must be disjoint from input: neither may contain the other. Separate from the main-phase consumer at runtime. |
| `ASTROCAT_ADD_HEADER` (`bool`) | `true` | same | `[astrocat].add_header` | true: input starts with data; false: skip exactly one line per input. Output always has a header. |
| `ASTROCAT_EXISTING_POLICY` (`const char*`) | `"fail"` | same | `[astrocat].existing_policy` | fail or overwrite. Overwrite replaces only files matching the current generated-prefix contract. |

The first two numeric columns are RA and Dec; extra columns are ignored. Rows
without two parseable numbers are skipped. Input coordinates must be finite,
with `0 <= RA <= 360` and `-90 <= Dec <= 90`; RA 360 is stored as zero and Dec 90
belongs to the last tile. Exact/one-ULP coordinate duplicates are removed across
input files and tile boundaries. Output begins with `RA    DEC` and retains
round-trip double precision. See [data preparation](DATA_PREPARATION.md#gaia-catalog)
for both filename conventions and the legacy coverage-overlap recommendation.
## External catalog preparation (`ExtCatConfig`)
Edit `config/ExtCatConfig.hpp` and `config/pathconfig.hpp`. Lite inherits its raw input, base width and RA/Dec/ZP defaults from `Initialize.hpp`. In Standard, `[lensing].source_cat` aliases `[extcat].output_directory`; in Lite they are independent producer/consumer fields.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `EXTCAT_INPUT_DIRECTORY` (`const char*`) | `""` | same | `[extcat].input_directory` | Readable raw catalog root. |
| `EXTCAT_OUTPUT_DIRECTORY` (`const std::string& (Standard); const char* (Lite)`) | `LensingConfig::SOURCE_CAT` | `""` | `[extcat].output_directory` | Writable published tile directory. Must not equal or be nested below the raw input. |
| `EXTCAT_FILENAME_TOKENS` (`std::vector<std::string>`) | `{}` | same | `[extcat].filename_tokens` | Case-sensitive basename tokens with OR matching; empty accepts all basenames. |
| `EXTCAT_RECURSIVE` (`bool`) | `true` | same | `[extcat].recursive` | Search subdirectories of the raw input. |
| `EXTCAT_DELIMITER` (`const char*`) | `"auto"` | same | `[extcat].delimiter` | auto, whitespace, comma or tab. |
| `EXTCAT_HEADER_MODE` (`const char*`) | `"auto"` | same | `[extcat].header_mode` | auto, present or absent. |
| `EXTCAT_MALFORMED_POLICY` (`const char*`) | `"fail"` | same | `[extcat].malformed_policy` | fail or skip. Skip explicitly accepts loss of malformed rows. |
| `EXTCAT_EXISTING_POLICY` (`const char*`) | `"fail"` | same | `[extcat].existing_policy` | fail or overwrite existing generated tiles. |
| `EXTCAT_CHUNK_MIB` (`std::uint64_t`) | `64` | same | `[extcat].chunk_mib` | Positive MPI byte-range task size in MiB. |
| `EXTCAT_TOTAL_COLUMNS` (`std::size_t`) | `18` | same | `[extcat].total_columns` | Positive canonical catalog width in pass-through mode. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` (`bool`) | `false` | same | `[extcat].use_explicit_columns` | Enable ordered projection; must be set explicitly in the INI. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` (`std::vector<std::size_t>`) | `{1, 2, ..., 18}` | same | `[extcat].input_columns` | Nonempty list of positive raw column indices, in output order. Retain every field needed downstream. |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` (`bool`) | `false` | same | `[extcat].use_explicit_coordinate_columns` | Use explicit raw RA/Dec positions instead of header discovery; set this flag explicitly. |
| `EXTCAT_RA_COLUMN_ONE_BASED` (`std::size_t`) | `5` | same | `[extcat].ra_column` | Positive one-based raw RA column. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` (`std::size_t`) | `6` | same | `[extcat].dec_column` | Positive one-based raw Dec column. |
| `EXTCAT_MAG_G_COLUMN_ONE_BASED` (`std::size_t`) | `7` | same | `[extcat].mag_g_column` | One-based raw g magnitude column; 0 marks this band absent. |
| `EXTCAT_MAG_R_COLUMN_ONE_BASED` (`std::size_t`) | `9` | same | `[extcat].mag_r_column` | One-based raw r magnitude column; 0 marks this band absent. |
| `EXTCAT_MAG_I_COLUMN_ONE_BASED` (`std::size_t`) | `11` | same | `[extcat].mag_i_column` | One-based raw i magnitude column; 0 marks this band absent. |
| `EXTCAT_MAG_Z_COLUMN_ONE_BASED` (`std::size_t`) | `13` | same | `[extcat].mag_z_column` | One-based raw z magnitude column; 0 marks this band absent. |
| `EXTCAT_MAG_Y_COLUMN_ONE_BASED` (`std::size_t`) | `15` | same | `[extcat].mag_y_column` | One-based raw y magnitude column; 0 marks this band absent. |
| `EXTCAT_ZP_COLUMN_ONE_BASED` (`std::size_t`) | `17` | same | `[extcat].zp_column` | Positive one-based raw zp/photo-z column; not a magnitude zero point. |

Explicit projection must retain RA, Dec, photo-z and every configured nonzero
magnitude column. FD requires at least one available magnitude and chooses the
first available in **i → z → r → g → y** order. Raw indices are one-based;
downstream column positions are resolved automatically. Do not edit output
indices by hand to compensate for a new input schema.
## Catalog and calibration paths (`LensingConfig`)
Compiled definitions are in `config/pathconfig.hpp`; Lite paths derive from `Initialize.hpp`. Runtime overrides do not modify header constants or automatically synchronize independent fields.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `ASTROMETRY_CAT` (`std::string`) | `"/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted"` | same | `[lensing].astrometry_cat` | Main-phase Gaia tile directory. Standard supports both layouts; Lite always consumes one-degree tiles. |
| `SOURCE_CAT` (`std::string`) | `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"` | Not present | `[lensing].source_cat` / `[extcat].output_directory` | Standard source-tile input and extcat publication default; the INI aliases update one runtime field. |
| `SOURCE_CAT_DEFAULT` (`const char*`) | Not present | `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"` | `[lensing].source_cat` | Lite main-phase source-tile input; independent of extcat publication output. |
| `FLAT_PATH` (`std::string`) | `"/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014"` | Not present | `[lensing].flat_path` | Readable super-flat FITS directory when include_flat=1. |
| `PSF_PATH` (`std::string`) | `"hahahaha"` | Not present | `[lensing].psf_path` | Placeholder: replace with a readable external PSF root before enabling ext_psf=1. |

## Branches and geometry (`LensingConfig`)
Edit `config/LensingConfig.hpp`. Lite camera/stage defaults are aliases of `Initialize.hpp`. Standard-only runtime keys are rejected by Lite, rather than silently ignored.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `DEFAULT_CHIP_COUNT` (`int`) | `62` | same | `[lensing].nmax_chip` | Positive effective CCD count; Lite derives this from Initialize::N_CCD. |
| `pixel_size` (`double`) | `0.2628` | same | `[lensing].pixel_size` | Positive angular scale in arcsec/pixel; Lite derives this from Initialize. |
| `chipnx` (`int`) | `2046` | same | `[lensing].chipnx` | Positive physical CCD width used by geometry-dependent consumers. Stage 1 reads actual FITS axes. |
| `chipny` (`int`) | `4094` | same | `[lensing].chipny` | Positive physical CCD height; keep consistent with camera geometry. |
| `ASTROMETRY_trivial` (`int`) | `0` | Not present | `[lensing].astrometry_trivial` | 0: Gaia astrometry; 1: deliberately use identity mapping. |
| `AstroCatType` (`int`) | `1` | Not present | `[lensing].astrometry_cat_type` | 1: legacy large Gaia tiles; 2: one-degree tiles. Lite has no selector and always uses one-degree tiles. |
| `PROCESS_stage` (`int`) | `223092870` | same | `[lensing].process_stage` | Product of enabled primes 2,3,5,7,11,13,17,19,23. Stage 9 (23) requires Stage 8 (19). |
| `include_FLAT` (`int`) | `0` | Not present | `[lensing].include_flat` | 0 off; 1 applies super-flat correction. |
| `include_Mask` (`int`) | `2` | Not present | `[lensing].include_mask` | Mask-input selector 0–3. Review all active consumers before omitting DQ; Lite always uses per-chip DQ. |
| `include_BGsub` (`int`) | `1` | Not present | No | 0 off; 1 subtracts the fitted Science background. Lite has no switch and follows its fixed preprocessing path. |
| `PreprocsType` (`int`) | `2` | Not present | No | 1: historical F77 background/sigma estimators; 2: current C++ estimators. Scientific comparison setting. |
| `ext_cat` (`int`) | `1` | Not present | `[lensing].ext_cat` | 0 disables external source catalog; 1 enables it. FD currently requires 1. |
| `ext_PSF` (`int`) | `0` | Not present | `[lensing].ext_psf` | 0: use exposure stars; 1: read external PSF images. |
| `CCD_split` (`int`) | `2` | same | `[lensing].ccd_split` | 1: whole CCD; 2: amplifier halves. Lite default comes from Initialize. |
| `deblending` (`int`) | `1` | Not present | `[lensing].deblending` | 0 off; 1 enables source deblending. Lite always deblends. |
| `PSF_type` (`int`) | `1` | Not present | `[lensing].psf_type` | 1: local polynomial; 2: hybrid PSF. Lite keeps local polynomial only. |
| `PSF_Ms` (`int`) | `0` | Not present | `[lensing].psf_ms` | 0 off; 1 enables PCA/multi-scale PSF reconstruction. |

## PSF selection (`LensingConfig`)
All entries in this section are compiled settings in `config/LensingConfig.hpp`. Rebuild after editing. Tune scientific cuts on representative data, not merely to make a run pass.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `npo` (`int`) | `64` | same | No | Exposure PSF sample count; positive. |
| `nstar_min` (`int`) | `npo * 3 / 2 = 96` | same | No | Derived minimum stars for exposure PSF fitting; change npo rather than this expression alone. |
| `npl` (`int`) | `10` | same | No | Local PSF polynomial coefficient count minus one. |
| `nstar_min_local` (`int`) | `16` | same | No | Minimum retained stars for a local fit. |
| `PsfGroupingType` (`int`) | `4` | Not present | No | 1 F77 selection, 2 threshold graph, 3 mutual KNN, 4 adaptive pair/fraction cuts. Lite keeps its fixed adaptive path. |
| `psf_exposure_min_candidates` (`int`) | `60` | same | No | Positive minimum exposure-wide PSF candidates. |
| `psf_count_pilot_clip_sigma` (`double`) | `3.0` | same | No | Positive pilot median/MAD clipping multiplier for integer star areas. |
| `psf_count_pilot_clip_iterations` (`int`) | `3` | same | No | Positive maximum pilot clipping passes. |
| `psf_count_zero_mad_quantile` (`double`) | `0.05` | same | No | Lower quantile for zero-MAD pilot; upper is 1-q; require 0 <= q < 0.5. |
| `psf_count_hist_range_sigma` (`double`) | `5.0` | same | No | Positive pilot-width multiplier for the local star-area search range; histogram bins remain fixed at two integer counts. |
| `psf_count_locus_sigma` (`double`) | `4.0` | same | No | Positive lower/upper count-MAD refinement multiplier; outer elbows may only widen final strict bounds. |
| `psf_count_locus_min_samples` (`int`) | `30` | same | No | Positive minimum samples for an integer star-area locus. |
| `psf_minchi_reference_fraction` (`double`) | `1.0 / 3.0` | same | No | Exposure top-size reference fraction in (0,1]. |
| `psf_minchi_reference_max_per_chip` (`int`) | `5` | same | No | Positive per-chip cap on reference stars. |
| `psf_minchi_sigma_cut` (`double`) | `4.0` | same | No | Positive minimum-chi rejection sigma multiplier. |
| `psf_knn_k` (`int`) | `20` | Not present | No | Positive number of retained neighbors for Standard KNN grouping. |
| `psf_group_merge_ratio` (`double`) | `0.30` | Not present | No | Nonnegative secondary-group relative-size threshold. |
| `psf_group_merge_min_gaia` (`int`) | `1` | Not present | No | Nonnegative minimum Gaia matches in a merged group. |
| `psf_gaia_match_radius_pix` (`double`) | `2.0` | same | No | Positive Gaia match radius in pixels. |
| `psf_gaia_locus_min_matches` (`int`) | `5` | same | No | Positive minimum Gaia matches supporting the locus. |
| `psf_pair_chi_valid_peak_fraction` (`double`) | `0.3678794411714423216` | same | No | exp(-1); strict relative-height threshold for a valid pair-chi peak; in (0,1). |
| `psf_bad_fraction_valid_peak_fraction` (`double`) | `0.10` | same | No | Strict relative-height threshold for a valid bad-pair fraction peak; in (0,1). |
| `psf_adaptive_elbow_search_height_fraction` (`double`) | `0.10` | Not present | No | Standard elbow candidates must lie below this fraction of the smoothed main peak; in (0,1). |
| `psf_type3_elbow_search_height_fraction` (`double`) | Not present | `0.10` | No | Lite elbow search-height fraction; in (0,1). This is the Lite symbol, not a Standard runtime key. |
| `psf_press_rejection_enabled` (`bool`) | `false` | `true` | No | Enable post-fit PRESS cleanup only after validating its effect on the target data. |
| `psf_press_sigma_cut` (`double`) | `4.0` | same | No | Positive standardized PRESS rejection sigma multiplier. |
| `psf_press_max_removals` (`int`) | `5` | same | No | Nonnegative maximum PRESS removals per chip. |
| `psf_loo_min_denom` (`double`) | `1.0e-6` | same | No | Minimum leave-one-out denominator in (0,1); numerical guard. |
| `step_psf` (`int`) | `100` | Not present | No | Positive spatial PSF star sampling step in pixels for the applicable Standard branch. |

## Stamps and detection (`LensingConfig`)
All entries without an INI mapping are in `config/LensingConfig.hpp` and require a rebuild. Keep camera geometry, stamp dimensions and detector-specific thresholds consistent.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `ns` (`int`) | `64` | same | No | Positive even stamp/Fourier-grid side length in pixels; coordinate all stamp/FFT changes. |
| `nsns` (`int`) | `ns * ns = 4096` | same | No | Derived pixel count in a Science stamp. |
| `chip_margin` (`int`) | `8` | same | No | Nonnegative extra chip-edge extraction margin in pixels. |
| `ns_2` (`int`) | `ns / 2 = 32` | same | No | Derived half stamp side. |
| `nl_2` (`int`) | `ns_2 + chip_margin = 40` | same | No | Derived half expanded extraction side. |
| `nl` (`int`) | `nl_2 * 2 = 80` | same | No | Derived expanded extraction side. |
| `flag_thresh` (`int`) | `3` | same | No | Nonnegative maximum accepted source extraction flag. |
| `chip_edge_margin` (`int`) | `chip_margin = 8` | same | No | Derived alias for chip-edge checks. |
| `dz_thresh` (`double`) | `0.1` | same | No | Nonnegative redshift tolerance for external-catalog matching. |
| `n_user_max` (`int`) | `500` | `200` | No | Maximum flux-ranked detections passed to astrometric pattern matching; a scientific selection limit, not a storage-capacity limit. |
| `ngal_max` (`int`) | `2000` | same | No | Initial galaxy-vector reservation hint; not a hard source limit. |
| `nstar_max` (`int`) | `1000` | same | No | Initial star-vector reservation hint; not a hard star limit. |
| `src_npara` (`int`) | `12` | same | No | Stage-6 source-marker metadata field count; keep consumers consistent. |
| `npd` (`int`) | `33` | same | No | PU astrometric distortion coefficient count; change only with a coordinated model change. |
| `source_thresh` (`double`) | `2.0` | same | No | Positive source-detection S/N threshold. |
| `core_thresh` (`double`) | `4.0` | same | No | Positive source-core threshold. |
| `n_neighbor` (`int`) | `5` | Not present | No | Positive neighbor count for Standard deblending. |
| `area_max` (`int`) | `ns * ns = 4096` | same | No | Derived maximum connected source area. |
| `area_thresh` (`int`) | `6` | same | No | Positive minimum connected source area in pixels. |
| `size_fit_rmax` (`int`) | `4` | same | No | Maximum radius for curvature-size fitting; validate numerical/scientific changes. |
| `point_stat_beta` (`double`) | `0.10` | same | No | Survey-wide fixed point-source Fourier-power template shape; not fitted per source. |
| `point_stat_k_frac` (`double`) | `0.90` | same | No | Fourier-radius fraction used in point-source statistics. |
| `point_stat_eps` (`double`) | `1.0e-20` | same | No | Positive denominator floor for point statistics. |
| `point_stat_min_corr` (`double`) | `1.0e-6` | same | No | Minimum accepted point-template correlation. |
| `SNR_PSF` (`double`) | `100.0` | same | No | Positive minimum PSF-star signal-to-noise ratio. |
| `saturation_thresh` (`double`) | `25000.0` | same | No | Saturated-pixel threshold in detector counts; review gain and data calibration. Lite inherits Initialize. |
| `gal_smooth` (`int`) | `0` | same | `[lensing].gal_smooth` | Nonnegative galaxy-stamp smoothing selector. |
| `star_smooth` (`int`) | `2` | same | `[lensing].star_smooth` | Nonnegative star-stamp smoothing selector. |

## Background estimation (`LensingConfig`)
Compiled settings in `config/LensingConfig.hpp`; rebuild after editing.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `nct` (`int`) | `12` | same | No | Number of background polynomial rectangle coefficients; positive. |
| `ncx` (`int`) | `3` | same | No | Background rectangle x count; keep consistent with nct. |
| `blocksize` (`int`) | `200` | same | No | Target side length of background blocks in pixels. |
| `bg_rough_grid_x` (`int`) | `32` | same | No | Positive rough-background grid column count. |
| `bg_rough_grid_y` (`int`) | `32` | same | No | Positive rough-background grid row count. |
| `bg_min_block_pixels` (`int`) | `1000` | same | No | Minimum valid pixels per background block. |
| `bg_min_clipped_pixels` (`int`) | `200` | same | No | Minimum surviving pixels after block clipping. |
| `bg_min_valid_frac` (`double`) | `0.25` | same | No | Minimum valid fraction per block, in (0,1]. |
| `bg_clip_low` (`double`) | `4.0` | same | No | Positive lower-side background clipping sigma. |
| `bg_clip_high` (`double`) | `2.5` | same | No | Positive upper-side background clipping sigma. |
| `bg_fit_clip_sigma` (`double`) | `3.0` | same | No | Positive plane-fit clipping sigma. |
| `bg_fit_max_iter` (`int`) | `4` | same | No | Maximum clipping iterations; nonnegative. |
| `bg_min_fit_factor` (`int`) | `3` | same | No | Positive minimum sample-to-coefficient factor for a fit. |

## Noise products (`LensingConfig`)
Compiled settings in `config/LensingConfig.hpp`. Types 1/2 are real-space stamps transformed downstream; Type 3 is already power. Lite retains the local noise-plane utility constants, but not the covariance-power branch or its `noise_cov_*` parameters.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `NstampType` (`int`) | `2` | Not present | No | Standard: 1 deterministic F77 real-space blank stamp, 2 QC/random physical blank stamp, 3 local covariance Fourier power. Lite keeps physical blank stamps. |
| `noise_sigma_ratio_min` (`double`) | `0.80` | same | No | Positive minimum blank-to-source sigma ratio. |
| `noise_sigma_ratio_max` (`double`) | `1.25` | same | No | Maximum blank-to-source sigma ratio; greater than minimum. |
| `noise_mad_ratio_min` (`double`) | `0.70` | same | No | Positive minimum blank-to-source MAD ratio. |
| `noise_mad_ratio_max` (`double`) | `1.30` | same | No | Maximum blank-to-source MAD ratio; greater than minimum. |
| `noise_tail_sigma` (`double`) | `2.5` | same | No | Positive tail-count sigma threshold. |
| `noise_max_tail_fraction` (`double`) | `0.05` | same | No | Maximum blank-stamp tail fraction, in [0,1]. |
| `noise_max_mask_fraction` (`double`) | `0.02` | same | No | Maximum masked fraction of a blank stamp, in [0,1]. |
| `noise_region_size` (`int`) | `192` | same | No | Positive even outer local-noise square side; must exceed inner size. |
| `noise_inner_size` (`int`) | `96` | same | No | Even central exclusion side; at least nl; outer-minus-inner must be even. |
| `noise_plane_min_valid_fraction` (`double`) | `0.30` | same | No | Minimum valid shell fraction for plane fitting, in (0,1]. |
| `noise_cov_padding_factor` (`double`) | `2.0` | Not present | No | Positive covariance FFT padding multiplier. |
| `noise_cov_fft_size` (`int`) | `ceil(noise_region_size * noise_cov_padding_factor) = 384` | Not present | No | Derived FFT side; must be at least 2*noise_region_size-1 for linear autocorrelation. |
| `noise_cov_max_lag` (`int`) | `63` | Not present | No | Maximum retained signed covariance lag; 0 <= lag < noise_region_size. |
| `noise_cov_min_valid_pixels` (`int`) | `4096` | Not present | No | Positive minimum covariance-mask pixel count. |
| `noise_cov_min_pair_fraction` (`double`) | `0.50` | Not present | No | Minimum valid pair-count fraction for a lag, in (0,1]. |
| `noise_cov_sigma_ratio_min` (`double`) | `0.80` | Not present | No | Positive lower accepted covariance sigma ratio. |
| `noise_cov_sigma_ratio_max` (`double`) | `1.25` | Not present | No | Upper accepted covariance sigma ratio; above the lower bound. |
| `noise_cov_max_negative_fraction` (`double`) | `0.25` | Not present | No | Maximum negative-power fraction, in [0,1]. |
| `noise_cov_imag_tolerance` (`double`) | `1.0e-10` | Not present | No | Nonnegative imaginary FFT residual tolerance. |

## Sigma estimation (`LensingConfig`)
Compiled settings in `config/LensingConfig.hpp`; rebuild after edits. Treat estimator/calibration changes as scientific changes and validate them on representative data.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `sig_blocksize` (`int`) | `200` | same | No | Positive estimator block side length in pixels. |
| `sig_block_max` (`int`) | `sig_blocksize * sig_blocksize = 40000` | same | No | Derived maximum pixels per block. |
| `sig_max_blocks` (`int`) | `2048` | same | No | Positive maximum sampled blocks. |
| `sig_min_block_pixels` (`int`) | `1000` | same | No | Minimum valid pixels in one block. |
| `sig_min_block_triples` (`int`) | `1000` | same | No | Minimum valid pixel triples per block. |
| `sig_min_blocks` (`int`) | `4` | same | No | Minimum blocks for a plane fit. |
| `sig_hist_nbin` (`int`) | `256` | same | No | Positive mode-finding histogram bin count. |
| `sig_hist_range` (`double`) | `6.0` | same | No | Positive histogram range in sigma units. |
| `sig_min_mode_count` (`int`) | `500` | same | No | Minimum samples defining the mode. |
| `sig_min_lower_count` (`int`) | `1000` | same | No | Minimum lower-side samples. |
| `sig_lower_quantile` (`double`) | `0.3173105` | same | No | Lower-side width quantile in (0,1); calibration constant. |
| `sig_clip_k` (`double`) | `3.0` | same | No | Positive symmetric clipping sigma multiplier. |
| `sig_rdil` (`int`) | `2` | same | No | Positive pixel-sampling stride. |
| `sig_clip_niter` (`int`) | `2` | same | No | Nonnegative clipping iteration count. |
| `sig_min_fit_triples` (`int`) | `1000` | same | No | Minimum triples in the final fit. |
| `sig_min_fit_frac` (`double`) | `0.20` | same | No | Minimum retained fit fraction, in (0,1]. |
| `sig_median_ratio` (`double`) | `1.2678405` | same | No | Positive median-to-sigma calibration factor. |
| `sig_plane_min` (`double`) | `1.0e-8` | same | No | Positive lower floor for the sigma plane. |
| `sig_max_plane_ratio` (`double`) | `4.0` | same | No | Maximum plane variation ratio; at least 1. |
| `sig_pivot_min` (`double`) | `1.0e-8` | same | No | Positive numerical pivot floor. |
| `sig_scale_s1` (`double`) | `0.673475` | same | No | Stage-1 noise-calibration candidate. |
| `sig_scale_s2` (`double`) | `1.027786` | same | No | Stage-2 noise-calibration value. |
| `sig_scale` (`double`) | `sig_scale_s2 = 1.027786` | same | No | Active calibration selector; produces the published 2*sigma^2 convention. |

## Catalog indices (`LensingConfig`)
Compiled in `config/LensingConfig.hpp`. These derived schema positions must remain consistent with writers and readers. Do not edit them to change raw catalog column mapping; use `[extcat]` instead.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `iid` (`int`) | `1 - 1 = 0` | same | No | PSF polynomial chi-square; zero-based index, not a selection threshold. |
| `ipixx` (`int`) | `2 - 1 = 1` | same | No | Source center x; zero-based index, not a selection threshold. |
| `ipixy` (`int`) | `3 - 1 = 2` | same | No | Source center y; zero-based index, not a selection threshold. |
| `isig` (`int`) | `4 - 1 = 3` | same | No | Local noise sigma; zero-based index, not a selection threshold. |
| `istar` (`int`) | `5 - 1 = 4` | same | No | Available PSF-star count; zero-based index, not a selection threshold. |
| `ipeak` (`int`) | `5 - 1 = 4` | same | No | Historical peak alias; zero-based index, not a selection threshold. |
| `i_imax` (`int`) | `6 - 1 = 5` | same | No | Peak x; zero-based index, not a selection threshold. |
| `i_jmax` (`int`) | `7 - 1 = 6` | same | No | Peak y; zero-based index, not a selection threshold. |
| `ih_flux` (`int`) | `8 - 1 = 7` | same | No | Half-light flux; zero-based index, not a selection threshold. |
| `ih_area` (`int`) | `9 - 1 = 8` | same | No | Source area; zero-based index, not a selection threshold. |
| `iflag` (`int`) | `10 - 1 = 9` | same | No | Quality flag; zero-based index, not a selection threshold. |
| `iPSF` (`int`) | `11 - 1 = 10` | same | No | Local PSF size; zero-based index, not a selection threshold. |
| `iSNR_F` (`int`) | `12 - 1 = 11` | same | No | Fourier S/N; zero-based index, not a selection threshold. |
| `ira` (`int`) | `13 - 1 = 12` | same | No | Right Ascension; zero-based index, not a selection threshold. |
| `idec` (`int`) | `14 - 1 = 13` | same | No | Declination; zero-based index, not a selection threshold. |
| `igf1` (`int`) | `15 - 1 = 14` | same | No | Field-distortion g1; zero-based index, not a selection threshold. |
| `igf2` (`int`) | `16 - 1 = 15` | same | No | Field-distortion g2; zero-based index, not a selection threshold. |
| `ig1` (`int`) | `17 - 1 = 16` | same | No | Fourier_Quad g1 estimator; zero-based index, not a selection threshold. |
| `ig2` (`int`) | `18 - 1 = 17` | same | No | Fourier_Quad g2 estimator; zero-based index, not a selection threshold. |
| `ide` (`int`) | `19 - 1 = 18` | same | No | Shear response estimator; zero-based index, not a selection threshold. |
| `ih1` (`int`) | `20 - 1 = 19` | same | No | Higher-order h1 estimator; zero-based index, not a selection threshold. |
| `ih2` (`int`) | `21 - 1 = 20` | same | No | Higher-order h2 estimator; zero-based index, not a selection threshold. |
| `icos2` (`int`) | `22 - 1 = 21` | same | No | Spin-2 cosine; zero-based index, not a selection threshold. |
| `isin2` (`int`) | `23 - 1 = 22` | same | No | Spin-2 sine; zero-based index, not a selection threshold. |
| `iparity` (`int`) | `24 - 1 = 23` | same | No | WCS parity; zero-based index, not a selection threshold. |
| `igalsizeT` (`int`) | `25 - 1 = 24` | same | No | Curvature galaxy size; zero-based index, not a selection threshold. |
| `ipsfsizeT` (`int`) | `26 - 1 = 25` | same | No | Curvature PSF size; zero-based index, not a selection threshold. |
| `idelta_chi2` (`int`) | `27 - 1 = 26` | same | No | Point-template delta chi-square; zero-based index, not a selection threshold. |
| `iorth_ext` (`int`) | `28 - 1 = 27` | same | No | Orthogonal point-template extension; zero-based index, not a selection threshold. |
| `shear_cat_ncols` (`int`) | `iorth_ext + 1 = 28` | same | No | Derived Stage-7 row width. |
| `expo_cat_ncols` (`int`) | `shear_cat_ncols + 1 = 29` | same | No | Derived per-exposure suffix width including Chi2. |
| `ichi2` (`int`) | `shear_cat_ncols = 28` | same | No | Zero-based exposure chi-square index within the pipeline suffix. |

## Calibration and PCA (`LensingConfig`)
Compiled in `config/LensingConfig.hpp`; rebuild after any deliberate change. Lite has no PCA/multi-scale branch.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `g1_c` (`double`) | `0.0` | same | No | Additive field-distortion g1 correction; recalibrate for another dataset/band. |
| `g2_c` (`double`) | `0.0` | same | No | Additive field-distortion g2 correction; recalibrate for another dataset/band. |
| `chi2_thresh` (`double`) | `0.1` | same | No | Nonnegative maximum exposure PSF chi-square; distinct from FDConfig::chi2_thresh. |
| `rescale_size` (`double`) | `1.2` | Not present | No | Positive target PSF residual rescaling size for PSF_Ms=1. |
| `procs_pn` (`int`) | `40` | Not present | No | Positive MPI ranks per PCA scheduling group; coordinate with resource allocation. |
| `work_pn` (`int`) | `10` | Not present | No | Positive concurrent PCA workers per group. |
| `nblocks` (`int`) | `2` | Not present | No | Positive PCA spatial block count per CCD axis. |
| `n_pcs` (`int`) | `100` | Not present | No | Positive maximum PCA principal components. |
| `npp6th` (`int`) | `28` | Not present | No | Sixth-degree two-dimensional polynomial term count; derived from the basis. |
| `pca_negative_eigenvalue_threshold` (`double`) | `-1.0e-5` | Not present | No | Nonpositive cutoff for invalid PCA eigenvalues. |
| `pi` (`double`) | `3.14159265358979323846` | same | No | Mathematical pi; do not tune. |
| `arc_convert` (`double`) | `pi / 180.0` | same | No | Derived degrees-to-radians conversion; do not tune. |

## Spatial rearrangement (`ProcessRearrConfig`)
Compiled in `config/ProcessRearrConfig.hpp`; fixed names are physically defined in `config/pathconfig.hpp`. Every entry requires rebuilding after an edit.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `SKY_GRID_DEGREES` (`double`) | `0.1` | same | No | Positive sky tile width in degrees; coordinate grid count changes. |
| `RA_BIN_COUNT` (`int`) | `3600` | same | No | Full-sky RA bin count consistent with SKY_GRID_DEGREES. |
| `DEC_BIN_COUNT` (`int`) | `1800` | same | No | Full-sky Dec bin count consistent with SKY_GRID_DEGREES. |
| `SKY_TILE_COUNT` (`std::size_t`) | `RA_BIN_COUNT * DEC_BIN_COUNT = 6480000` | same | No | Derived total full-sky tile count. |
| `TARGET_SUBCAT_ROWS` (`std::uint64_t`) | `500000` | same | No | Positive target partition row count; controls file size and memory use. |
| `SKIP_DIRECTORY_NAME` (`std::string_view`) | `"Large_Field"` | same | No | Directory name excluded from scans. |
| `SUBCAT_PREFIX` (`std::string_view`) | `"subcat_"` | same | No | Published partition filename prefix. |
| `SUBCAT_EXTENSION` (`std::string_view`) | `".cat"` | same | No | Published partition filename extension. |
| `SUBCAT_ID_WIDTH` (`int`) | `6` | same | No | Positive minimum zero-padded partition-ID width. |
| `SUMMARY_FILENAME` (`std::string_view`) | `"catalog_summary.txt"` | same | No | Published summary filename. |
| `OUTPUT_PRECISION` (`int`) | `10` | same | No | Positive significant digits in catalog output. |
| `SUMMARY_PRECISION` (`int`) | `4` | same | No | Nonnegative decimal places in summary bounds. |
| `SKIP_MISSING_CATALOGS` (`bool`) | `true` | same | No | Continue past missing catalogs; false enforces strict completeness. |
| `SKIP_MALFORMED_ROWS` (`bool`) | `true` | same | No | Continue past malformed rows; false enforces strict row validity. |

Catalog widths and offsets are resolved from the effective runtime schema, not
configured in this header. Complete rows contain the external prefix,
`EXPO_NUM`, `ccD_NUM`, 28 Stage-7 fields and exposure `Chi2`. The default is
**18 + 2 + 29 = 49 columns** (`EXPO_NUM=18`, `ccD_NUM=19`, source base 20,
`Chi2=48`, all zero-based). Standard with `ext_cat=0` has 31 columns. Regenerate
older incompatible rows before rearrangement or FD processing.
## Output layout (`OutputLayout`)
Both arrays are identical across variants and defined in `config/pathconfig.hpp`. They contain relative names below the dataset root, not independent input/output roots. One dataset-root mount normally covers them all. Changing names requires a rebuild and coordinated product regeneration.
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `NON_CHIP_BASE_DIRECTORIES` (`std::array<const char*, 16>`) | `16 paths below` | same | No | Complete fixed dataset directories without a chip suffix. |
| `CHIP_PRODUCT_DIRECTORIES` (`std::array<const char*, 16>`) | `16 paths below` | same | No | Complete fixed per-chip product-directory prefixes. |

| NON_CHIP_BASE_DIRECTORIES | CHIP_PRODUCT_DIRECTORIES |
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
| Parameter (type) | Standard default | Lite default | INI override | Meaning and constraints |
| --- | --- | --- | --- | --- |
| `FD_STATIC_MODE` (`StaticMode`) | `StaticMode::PDF_SIGMA` | same | No | PDF_SIGMA: PDF chi-square estimate with analytic sigma; PDF_JACK: PDF estimate with jackknife sigma; SWSE_JACK: ratio estimate with jackknife sigma. |
| `FD_USE_PDF_STATIS` (`bool`) | `mode is PDF_SIGMA or PDF_JACK` | same | No | Derived PDF sign-test switch; change FD_STATIC_MODE, not this expression. |
| `FD_USE_JACKKNIFE` (`bool`) | `mode is PDF_JACK or SWSE_JACK` | same | No | Derived jackknife-uncertainty switch. |
| `FD_USE_SWSE_DATA` (`bool`) | `mode is SWSE_JACK` | same | No | Derived SWSE-data-model switch. |
| `FD_PER_EXPOSURE_STAR_BAR` (`bool`) | `false` | same | No | false: single global N-to-1 star bar; true: per-exposure N-to-N star bars. |
| `nmax_per_core` (`int`) | `20000000` | same | No | Maximum sources reserved per MPI rank; review per-rank memory use. |
| `fd_num` (`int`) | `21` | same | No | Positive field-distortion spatial bin count. |
| `PDF_BINS` (`int`) | `4` | same | No | Positive equal-probability inner-bin count. |
| `gf_lim` (`float`) | `0.0015` | same | No | Positive spatial field-distortion half-range, from -gf_lim to +gf_lim. |
| `NMAX` (`int`) | `200` | same | No | Positive fine-grid sampling count. |
| `MAX_DUP` (`int`) | `5` | same | No | Positive maximum duplicate measurements. |
| `N_jack` (`int`) | `50` | same | No | Positive jackknife region count. |
| `nmax_total` (`int`) | `1000000` | same | No | Positive maximum sources used in k-means. |
| `Km_iter` (`int`) | `100` | same | No | Positive k-means iteration count. |
| `snrfcut` (`float`) | `4.0` | same | No | Minimum Fourier signal-to-noise ratio; nonnegative. |
| `snrlow` (`float`) | `0.0` | same | No | Optional lower source-S/N bound; zero disables it. |
| `snrhigh` (`float`) | `0.0` | same | No | Optional upper source-S/N bound; zero disables it; keep above the lower bound when both active. |
| `starcut` (`float`) | `20.0` | same | No | Point-source size threshold; recalibrate with the data. |
| `chi2_thresh` (`float`) | `0.01` | same | No | Maximum exposure chi-square; separate from the LensingConfig threshold. |
| `flagcut` (`float`) | `0.0` | same | No | Maximum accepted source quality flag. |
| `imaxcut` (`float`) | `64.0` | same | No | Upper source-peak x coordinate cut; coordinate with stamp geometry. |
| `jmaxcut` (`float`) | `64.0` | same | No | Upper source-peak y coordinate cut; coordinate with stamp geometry. |
| `zplow` (`float`) | `0.0` | same | No | Minimum photometric redshift. |
| `zphigh` (`float`) | `3.0` | same | No | Maximum photometric redshift; greater than zplow. |
| `r_half_thresh` (`float`) | `0.0` | same | No | Optional half-light-radius threshold; zero disables it. |
| `star_bar_mltp` (`float`) | `3.0` | same | No | Positive stellar-locus sigma multiplier. |
| `psf_chi2_mltp` (`float`) | `3.0` | same | No | Positive PSF chi-square sigma multiplier. |
| `n_size_bins` (`int`) | `100` | same | No | Positive stellar-size histogram bin count. |
| `n_mag_bins` (`int`) | `20` | same | No | Positive magnitude histogram bin count. |
| `size_min` (`float`) | `-2.0` | same | No | Lower size-histogram bound. |
| `size_max` (`float`) | `2.0` | same | No | Upper size-histogram bound, greater than size_min. |
| `mag_min_val` (`float`) | `10.0` | same | No | Lower magnitude-histogram bound. |
| `mag_max_val` (`float`) | `30.0` | same | No | Upper magnitude-histogram bound, greater than mag_min_val. |
| `min_bin_count` (`int`) | `100` | same | No | Positive minimum sample count for a usable bin. |
| `peak_match_tol` (`float`) | `0.05` | same | No | Positive stellar-peak matching tolerance. |
| `min_concentration` (`float`) | `0.6` | same | No | Minimum stellar-locus concentration fraction in [0,1]. |
| `star_phy_min` (`float`) | `-0.5` | same | No | Lower physical stellar-size bound. |
| `star_phy_max` (`float`) | `0.2` | same | No | Upper physical stellar-size bound, above star_phy_min. |
| `stage1_snr` (`float`) | `40.0` | same | No | Nonnegative S/N threshold for histogram accumulation. |
| `stage2_snr` (`float`) | `0.0` | same | No | Optional second-pass S/N threshold; zero disables it. |
| `init_win_active` (`float`) | `0.1` | same | No | Positive initial size window for an active exposure. |
| `init_win_fallback` (`float`) | `0.15` | same | No | Positive fallback initial size window. |
| `default_s_init` (`float`) | `0.5` | same | No | Default initial stellar-size center. |
| `clip_nsigma` (`float`) | `3.0` | same | No | Positive iterative clipping sigma multiplier. |
| `min_clip_limit` (`float`) | `0.015` | same | No | Positive minimum clipping half-width. |
| `default_s_std` (`float`) | `0.05` | same | No | Positive default stellar-size scatter. |
| `fallback_scut_default` (`float`) | `0.6` | same | No | Fallback stellar-size cut. |
| `bad_ccds` (`int[]`) | `{-1}` | same | No | Backing CCD-ID array. With n_bad_ccds=0 no CCD is excluded; the -1 entry is inactive storage, not an active detector ID. |
| `n_bad_ccds` (`int`) | `0` | same | No | Active count of excluded CCD IDs. Keep this explicit count within the actual bad_ccds array length; it is not automatically inferred. |
| `chip_mask_edge` (`int`) | `50` | same | No | Nonnegative symmetric edge mask in pixels; inclusive bounds derive as [edge, chipnx-edge] and [edge, chipny-edge]. |

Per-exposure statistics use serialized `EXPO_NUM` identities, not file order or
terminal `Chi2`. The runtime allocation follows the MPI-global maximum valid
exposure ID; there is no separate configured exposure cap. Scientific FD
measurement remains part of the C++ pipeline and is unrelated to removed
container smoke-test utilities.
## Lite `Initialize.hpp`

This is Lite's convenient compiled-default entry point. Editing any entry below
requires rebuilding. The corresponding namespace headers alias these values;
INI overrides change runtime copies only. The table lists **every declaration**
in `Initialize.hpp`; `DatasetSpec` is the structure containing `target` and
`prefix`, not a third configuration mechanism.

| Initialize entry | Compiled default | Effective setting / INI |
| --- | --- | --- |
| `N_CCD` | `62` | LensingConfig::DEFAULT_CHIP_COUNT; `[lensing].nmax_chip` |
| `pixel_size` | `0.2628` | `[lensing].pixel_size` |
| `chipnx` | `2046` | `[lensing].chipnx` |
| `chipny` | `4094` | `[lensing].chipny` |
| `saturation_thresh` | `25000.0` | LensingConfig::saturation_thresh; compile only |
| `CCD_split` | `2` | `[lensing].ccd_split` |
| `PROCESS_stage` | `223092870` | `[lensing].process_stage` |
| `RUN_PROCESS_ASTROCAT` | `false` | `[process].run_process_astrocat` |
| `RUN_PROCESS_EXTCAT` | `false` | `[process].run_process_extcat` |
| `RUN_PROCESS_INIT` | `true` | `[process].run_process_init` |
| `RUN_PROCESS_MAIN` | `true` | `[process].run_process_main` |
| `RUN_PROCESS_REARR` | `false` | `[process].run_process_rearr` |
| `RUN_PROCESS_FD` | `false` | `[process].run_process_fd` |
| `DatasetSpec::target` | `"gband" in DATASETS` | `[init].datasets` |
| `DatasetSpec::prefix` | `"c4d_" in DATASETS` | `[init].datasets` |
| `DATASETS` | `{{"gband", "c4d_"}}` | `[init].datasets` |
| `CONTAINS` | `{"v1"}` | `[init].contains` |
| `EXPO_LIST` | `""` | `[process].expo_list` |
| `ASTROMETRY_CAT` | `"/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted"` | `[lensing].astrometry_cat` |
| `SOURCE_CAT_DEFAULT` | `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"` | `[lensing].source_cat` |
| `SCIENCE_ROOT` | `"/lustre/home/acct-phyzj/share/DES/g"` | `[init].science_root` |
| `DQ_ROOT` | `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"` | `[init].dq_root` |
| `OUTPUT_ROOT` | `"/lustre/home/acct-phyzj/share/DES/g_band_v1"` | `[init].output_root` |
| `ASTROMETRY_TILE_PREFIX` | `"astra_"` | AstroCatConfig::ASTROMETRY_TILE_PREFIX; compile only |
| `ASTROCAT_INPUT_DIRECTORY` | `""` | `[astrocat].input_directory` |
| `ASTROCAT_ADD_HEADER` | `true` | `[astrocat].add_header` |
| `SOURCE_CAT_TILE_PREFIX` | `"extern_"` | ExtCatConfig::SOURCE_CAT_TILE_PREFIX; compile only |
| `EXTCAT_INPUT_DIRECTORY` | `""` | `[extcat].input_directory` |
| `EXTCAT_TOTAL_COLUMNS` | `18` | `[extcat].total_columns` |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `5` | `[extcat].ra_column` |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `6` | `[extcat].dec_column` |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `17` | `[extcat].zp_column` |
| `ARCHIVE_SUFFIX` | `".fits.fz"` | InitConfig::ARCHIVE_SUFFIX; compile only |
| `CCDNUM_KEYWORD` | `"CCDNUM"` | InitConfig::CCDNUM_KEYWORD; compile only |
| `DQ_STEM_REPLACE_FROM` | `"ood"` | InitConfig::DQ_STEM_REPLACE_FROM; compile only |
| `DQ_STEM_REPLACE_TO` | `"ooi"` | InitConfig::DQ_STEM_REPLACE_TO; compile only |

## Runtime selection and coupled values

The five INI sections are `[process]`, `[astrocat]`, `[extcat]`, `[init]` and
`[lensing]`. All supported keys appear in the tables above. A sample with all
keys is provided in each variant's `pipeline.example.ini`; its entries are
commented until you opt in.

| Relationship | User action |
| --- | --- |
| Main stages | Select prime factors in `process_stage`; include 19 whenever including 23. |
| Gaia producer and consumer | Set `[astrocat].output_directory` and `[lensing].astrometry_cat` together when reusing newly produced tiles. Select Type 2 in Standard; Lite is fixed Type 2. |
| Standard source tiles | `[extcat].output_directory` and `[lensing].source_cat` address the same field. Set only one, or make both identical. |
| Lite source tiles | Producer output and main input are independent. Set both to the same directory when appropriate. |
| Projection | Set `use_explicit_columns = true` together with `input_columns`; set `use_explicit_coordinate_columns = true` for explicit RA/Dec positions. |
| Downstream-only lists | An explicit `expo_list` is an existing file. Without it, the program derives one list per dataset from `output_root`. |
| Main output | Dataset roots come from the listed image layout; changing only `output_root` does not relocate an already prepared dataset. |
| Rebuild | Any header edit changes compiled settings only after a rebuild. An explicit INI still overrides mapped values. |

Select scientific thresholds and branches using representative observations.
A syntactically valid configuration does not certify data completeness or
scientific suitability.

