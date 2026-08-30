# C++ Pipeline Parameter Reference

This reference follows the live Standard and Lite configuration namespaces.
Each namespace has one section and one complete table. Defaults are compiled into
the executable; only entries named in **INI / CLI override** can change at run time.

Runtime precedence is `config header default < --config INI < CLI`. An INI or CLI
override applies only to fields represented by `RuntimeConfig`; all other edits require
rebuilding the selected variant. `N/A — removed in Lite` means the alternate branch is
physically absent from Lite and cannot be restored by adding a constant.

Derived parameters remain listed for source coverage. Do not edit them directly; change
their source parameter and preserve the associated assertions and consumers.

## Centralized path configuration

Each variant has its own `config/pathconfig.hpp`. It is the sole physical source
for fixed input/output paths, workflow list/output names, rearrangement filenames,
and fixed relative output-directory layouts. The established namespaces remain
unchanged, so existing call sites and INI mapping continue to use the same
namespace-qualified symbols.

| Namespace | Definitions physically owned by `pathconfig.hpp` | Coupling and RuntimeConfig rule |
|---|---|---|
| `LensingConfig` | `ASTROMETRY_CAT`, `SOURCE_CAT`; Standard only: `FLAT_PATH`, `PSF_PATH` | These seed RuntimeConfig. `[lensing].astrometry_cat`, `source_cat`, `flat_path`, and `psf_path` can replace the surviving values at run time; Lite rejects the removed flat/PSF keys. |
| `AstroCatConfig` | `ASTROCAT_INPUT_DIRECTORY`, `ASTROCAT_OUTPUT_DIRECTORY` | The compiled output default is deliberately initialized as `LensingConfig::ASTROMETRY_CAT`. Once RuntimeConfig is created, producer output and Stage-1 consumer input are separate fields. |
| `ExtCatConfig` | `EXTCAT_INPUT_DIRECTORY`, `EXTCAT_OUTPUT_DIRECTORY` | `EXTCAT_OUTPUT_DIRECTORY` deliberately remains a `const std::string&` to `LensingConfig::SOURCE_CAT`. `[lensing].source_cat` and `[extcat].output_directory` address the same effective RuntimeConfig field; CLI has final precedence. |
| `InitConfig` | `SCIENCE_ROOT`, `DQ_ROOT`, `OUTPUT_ROOT` | These seed RuntimeConfig and have INI/CLI overrides. |
| `ProcessConfig` | `EXPO_LIST`, `REARR_OUTPUT_DIRECTORY`, `REARR_OUTPUT_BASE_DIRECTORY`, `REARRANGED_EXPO_LIST_FILENAME`, `REARRANGED_EXPO_LIST_DIRECTORY`, `FD_EXPO_LIST`, `FD_OUTPUT_DIRECTORY`, `FD_OUTPUT_BASE_DIRECTORY` | These seed RuntimeConfig and have INI/CLI overrides. |
| `ProcessRearrConfig` | `SKIP_DIRECTORY_NAME`, `SUBCAT_PREFIX`, `SUBCAT_EXTENSION`, `SUMMARY_FILENAME` | No RuntimeConfig field; edit the selected variant and rebuild. |
| `OutputLayout` | `NON_CHIP_BASE_DIRECTORIES`, `CHIP_PRODUCT_DIRECTORIES` | No RuntimeConfig field; these are fixed relative directory contracts used by initialization and processing. |

Editing `pathconfig.hpp` changes compiled defaults and therefore requires a clean
rebuild. INI and CLI overrides change only RuntimeConfig copies according to
`compiled default < INI < CLI`; they do not mutate the header constants.

## `ProcessConfig` (`config/ProcessConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `RUN_PROCESS_ASTROCAT` | `bool` | `false` | same | `[process].run_process_astrocat`; `--run-astrocat` | Boolean | Run Gaia-catalog tiling by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `RUN_PROCESS_EXTCAT` | `bool` | `false` | same | `[process].run_process_extcat`; `--run-extcat` | Boolean | Run external-catalog tiling by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `RUN_PROCESS_INIT` | `bool` | `true` | same | `[process].run_process_init`; `--run-init` | Boolean | Run archive initialization by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `RUN_PROCESS_MAIN` | `bool` | `true` | same | `[process].run_process_main`; `--run-main` | Boolean | Run the nine-stage pipeline by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `RUN_PROCESS_REARR` | `bool` | `true` | `false` | `[process].run_process_rearr`; `--run-rearr` | Boolean | Run catalog rearrangement by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `RUN_PROCESS_FD` | `bool` | `true` | `false` | `[process].run_process_fd`; `--run-fd` | Boolean | Run the field-distortion test by default. | Select phases per run with CLI. | No at run time; yes if editing the header |
| `EXPO_LIST` | `const char*` | empty | same | `[process].expo_list`; `--expo-list` or one positional path | Exposure-list path | Default top-level exposure-list path. | Prefer CLI for each run. | No at run time; yes if editing the header |
| `REARR_OUTPUT_DIRECTORY` | `const char*` | `"baked"` | same | `[process].rearr_output_directory`; `--rearr-output-dir` | Directory name/path | Rearranged catalog directory. | Change when output layout changes. | No at run time; yes if editing the header |
| `REARR_OUTPUT_BASE_DIRECTORY` | `const char*` | empty | same | `[process].rearr_output_base_directory`; `--rearr-output-base` | Empty = dataset root | Empty uses the dataset root. | Change when outputs live outside the dataset root. | No at run time; yes if editing the header |
| `REARRANGED_EXPO_LIST_FILENAME` | `const char*` | `"cat_gband_ori.list"` | same | `[process].rearranged_expo_list_filename`; `--rearr-list-name` | Filename | FD input list name. | Change for another naming convention. | No at run time; yes if editing the header |
| `REARRANGED_EXPO_LIST_DIRECTORY` | `const char*` | empty | same | `[process].rearranged_expo_list_directory`; `--rearr-list-dir` | Empty = input-list parent | Empty uses the input-list directory. | Change for another list location. | No at run time; yes if editing the header |
| `FD_EXPO_LIST` | `const char*` | empty | same | `[process].fd_expo_list`; `--fd-expo-list` | Exposure-list path | Optional explicit FD exposure-list path. | Override when FD must use a different list. | No at run time; yes if editing the header |
| `FD_OUTPUT_DIRECTORY` | `const char*` | `"fdout"` | same | `[process].fd_output_directory`; `--fd-output-dir` | Directory name/path | FD result directory name. | Change when output layout changes. | No at run time; yes if editing the header |
| `FD_OUTPUT_BASE_DIRECTORY` | `const char*` | empty | same | `[process].fd_output_base_directory`; `--fd-output-base` | Empty = dataset root | Empty uses the dataset root. | Change when outputs live outside the dataset root. | No at run time; yes if editing the header |

The runtime option structs are mutable copies of these defaults, not a second
source of defaults. At least one top-level phase must be enabled.

## `InitConfig` (`config/InitConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `DatasetSpec::target` | `std::string` | `"gband" in DATASETS` | same | `[init].datasets`; `--dataset` or `--target` | Non-empty output-directory name | Names the dataset output directory. | Change for each dataset. | No at run time; yes if editing the header |
| `DatasetSpec::prefix` | `std::string` | `"c4d_" in DATASETS` | same | `[init].datasets`; `--dataset` or `--prefix` | Non-empty archive-basename prefix | Selects matching archive basenames. | Change for each dataset. | No at run time; yes if editing the header |
| `SCIENCE_ROOT` | `const char*` | `"/lustre/home/acct-phyzj/share/DES/g"` | same | `[init].science_root`; `--science-root` | Readable directory path | Science archive root. | Change for another site or archive. | No at run time; yes if editing the header |
| `DQ_ROOT` | `const char*` | `"/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"` | same | `[init].dq_root`; `--dq-root` | Readable directory path | DQ archive root. | Change for another site; irrelevant only when DQ access is disabled. | No at run time; yes if editing the header |
| `OUTPUT_ROOT` | `const char*` | `"/lustre/home/acct-phyzj/share/DES/g_band_v1"` | same | `[init].output_root`; `--output-root` | Writable directory path | Pipeline output root. | Change for every deployment. | No at run time; yes if editing the header |
| `DATASETS` | `std::vector<DatasetSpec>` | `{ {"gband", "c4d_"} }` | same | `[init].datasets`; repeatable `--dataset` | Unique targets with non-empty prefixes | Datasets processed sequentially. | Change for another dataset set. | No at run time; yes if editing the header |
| `CONTAINS` | `std::vector<std::string>` | `{"v1"}` | same | `[init].contains`; repeatable `--contains` | Non-empty case-sensitive basename tokens; OR matching | OR-matched archive basename tokens. | Change when archive naming changes. | No at run time; yes if editing the header |
| `EXISTING` | `const char*` | `"fail"` | same | `[init].existing`; `--existing` | fail, resume, or overwrite | Existing-output policy. | Select intentionally per run. | No at run time; yes if editing the header |
| `F77_MAX_PATH` | `int` | `0` | same | `[init].f77_max_path`; `--f77-max-path` | Non-negative; 0 disables the guard | Generated-path compatibility limit; zero disables it. | Change only for path-policy compatibility. | No at run time; yes if editing the header |

## `AstroCatConfig` (`config/AstroCatConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `ASTROCAT_INPUT_DIRECTORY` | `const char*` | empty | same | `[astrocat].input_directory`; `--astrocat-input` | Readable flat directory | Raw Gaia files; each data row begins with RA and Dec. | Set when running `process_astrocat`. | No at run time; yes if editing the header |
| `ASTROCAT_OUTPUT_DIRECTORY` | `std::string` | `LensingConfig::ASTROMETRY_CAT` | same | `[astrocat].output_directory`; `--astrocat-output` | Writable directory that does not equal, contain, or sit below the input directory | Compiled default destination for one-degree Type-2 Gaia tiles. | Override independently for each publication. | No at run time; yes if editing the header |
| `ASTROCAT_ADD_HEADER` | `bool` | `true` | same | `[astrocat].add_header`; `--astrocat-add-header` | `true` starts at the first line; `false` skips exactly one line per input file | Controls raw-input header handling; output tiles always contain `RA    DEC`. | Change to match the raw files. | No at run time; yes if editing the header |
| `ASTROCAT_EXISTING_POLICY` | `const char*` | `"fail"` | same | `[astrocat].existing_policy`; `--astrocat-existing` | fail or overwrite | Existing generated-tile policy. | Select intentionally for reruns. | No at run time; yes if editing the header |

In `pathconfig.hpp`, `ASTROCAT_OUTPUT_DIRECTORY` intentionally derives from
`LensingConfig::ASTROMETRY_CAT`; keep both symbols and that expression rather
than merging them. After RuntimeConfig is constructed,
`[astrocat].output_directory` and `--astrocat-output` change only the
`process_astrocat` destination. They are not compared with or propagated back
to `[lensing].astrometry_cat`.

The phase discovers only direct regular children of the input directory; it
does not recurse. It reads each complete file through dynamic MPI scheduling,
optionally skips exactly one first line, replaces commas with spaces, consumes
the first two parseable doubles, and silently skips rows without two doubles.
The input contract is finite sky coordinates with `0 <= RA <= 360` and
`-90 <= Dec <= 90`; exactly `RA=360` is stored as zero and exactly `Dec=90`
belongs to the last Dec tile. Exact and one-ULP duplicates in both coordinates
are removed, including duplicates across tile boundaries. Output files use
`des_y6_RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`, always begin with `RA    DEC`,
and contain round-trip-precision doubles. `overwrite` removes only files that
match this generated basename contract and preserves unrelated directory
content.

## `ExtCatConfig` (`config/ExtCatConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `EXTCAT_INPUT_DIRECTORY` | `const char*` | empty | same | `[extcat].input_directory`; `--extcat-input` | Readable directory path | Root containing raw catalog files. | Set when running `process_extcat`. | No at run time; yes if editing the header |
| `EXTCAT_OUTPUT_DIRECTORY` | `std::string&` | `LensingConfig::SOURCE_CAT` | same | `[extcat].output_directory`; `--extcat-output` | Writable tile directory | Tile output and pipeline input directory. | Change for each catalog deployment. | No at run time; yes if editing the header |
| `EXTCAT_FILENAME_TOKENS` | `std::vector<std::string>` | `{}` | same | `[extcat].filename_tokens`; `--extcat-contains` | Non-empty basename tokens; OR matching | OR-matched basename filters. | Change when filenames need filtering. | No at run time; yes if editing the header |
| `EXTCAT_RECURSIVE` | `bool` | `true` | same | `[extcat].recursive`; `--extcat-recursive` | Boolean | Recurse below the input directory. | Disable for a flat directory only. | No at run time; yes if editing the header |
| `EXTCAT_DELIMITER` | `const char*` | `"auto"` | same | `[extcat].delimiter`; `--extcat-delimiter` | auto, whitespace, comma, or tab | Input delimiter detection mode. | Change when auto-detection is unsuitable. | No at run time; yes if editing the header |
| `EXTCAT_HEADER_MODE` | `const char*` | `"auto"` | same | `[extcat].header_mode`; `--extcat-header` | auto, present, or absent | Input header handling mode. | Change when auto-detection is unsuitable. | No at run time; yes if editing the header |
| `EXTCAT_MALFORMED_POLICY` | `const char*` | `"fail"` | same | `[extcat].malformed_policy`; `--extcat-malformed` | fail or skip | Malformed-row handling policy. | Use `skip` only after accepting data loss. | No at run time; yes if editing the header |
| `EXTCAT_EXISTING_POLICY` | `const char*` | `"fail"` | same | `[extcat].existing_policy`; `--extcat-existing` | fail or overwrite | Existing-tile handling policy. | Select intentionally for reruns. | No at run time; yes if editing the header |
| `EXTCAT_CHUNK_MIB` | `std::uint64_t` | `64` | same | `[extcat].chunk_mib`; `--extcat-chunk-mib` | Positive MiB value | MPI byte-range task size in MiB. | Tune for storage and rank count. | No at run time; yes if editing the header |
| `EXTCAT_TOTAL_COLUMNS` | `std::size_t` | `18` | same | `[extcat].total_columns` | Positive canonical width | Pass-through external catalog width. | Change only with coordinated schema consumers. | No at run time; yes if editing the header |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `bool` | `false` | same | `[extcat].use_explicit_columns`; enabled by `--extcat-columns` | Boolean | Enable ordered column projection. | Use when raw tables contain extra/reordered fields. | No at run time; yes if editing the header |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `std::vector<std::size_t>` | `{1, ..., 18}` | same | `[extcat].input_columns`; `--extcat-columns` | Non-empty positive one-based indices | Raw one-based columns emitted in output order. | Change for another raw schema. | No at run time; yes if editing the header |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` | `bool` | `false` | same | `[extcat].use_explicit_coordinate_columns`; enabled by coordinate-column CLI | Boolean | Override named-field discovery. | Use for known positional schemas. | No at run time; yes if editing the header |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `std::size_t` | `5` | same | `[extcat].ra_column`; `--extcat-ra-column` | Positive one-based raw column | Raw one-based RA column. | Change for another catalog schema. | No at run time; yes if editing the header |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `std::size_t` | `6` | same | `[extcat].dec_column`; `--extcat-dec-column` | Positive one-based raw column | Raw one-based Dec column. | Change for another catalog schema. | No at run time; yes if editing the header |
| `EXTCAT_MAG_G_COLUMN_ONE_BASED` | `std::size_t` | `7` | same | `[extcat].mag_g_column` | 0 = band absent; otherwise a positive one-based raw column | Raw one-based g-magnitude column. | Change when the raw external-catalog schema changes. | No at run time; yes if editing the header |
| `EXTCAT_MAG_R_COLUMN_ONE_BASED` | `std::size_t` | `9` | same | `[extcat].mag_r_column` | 0 = band absent; otherwise a positive one-based raw column | Raw one-based r-magnitude column. | Change when the raw external-catalog schema changes. | No at run time; yes if editing the header |
| `EXTCAT_MAG_I_COLUMN_ONE_BASED` | `std::size_t` | `11` | same | `[extcat].mag_i_column` | 0 = band absent; otherwise a positive one-based raw column | Raw one-based i-magnitude column. | Change when the raw external-catalog schema changes. | No at run time; yes if editing the header |
| `EXTCAT_MAG_Z_COLUMN_ONE_BASED` | `std::size_t` | `13` | same | `[extcat].mag_z_column` | 0 = band absent; otherwise a positive one-based raw column | Raw one-based z-magnitude column. | Change when the raw external-catalog schema changes. | No at run time; yes if editing the header |
| `EXTCAT_MAG_Y_COLUMN_ONE_BASED` | `std::size_t` | `15` | same | `[extcat].mag_y_column` | 0 = band absent; otherwise a positive one-based raw column | Raw one-based y-magnitude column. | Change when the raw external-catalog schema changes. | No at run time; yes if editing the header |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `std::size_t` | `17` | same | `[extcat].zp_column`; `--extcat-zp-column` | Positive one-based raw column | Raw one-based photo-z column. | Change for another catalog schema. | No at run time; yes if editing the header |

`[extcat].output_directory`, `[lensing].source_cat`, and `--extcat-output`
address the same effective external-tile directory. Explicit projection must
retain RA, Dec, photo-z, and every nonzero magnitude field. FD needs at least one
magnitude and selects the first available band in i, z, r, g, y order.

## `LensingConfig` (`config/LensingConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `pi` | `double` | `3.14159265358979323846` | same | No | Mathematical constant | Mathematical pi. | Derived parameter — do not edit directly. | Yes |
| `arc_convert` | `double` | `pi / 180.0` | same | No | Radians per degree | Degrees-to-radians conversion factor. | Derived parameter — do not edit directly. | Yes |
| `ASTROMETRY_trivial` | `int` | `0` | N/A — removed in Lite | `[lensing].astrometry_trivial` | 0 = Gaia; 1 = identity mapping | Use Gaia astrometry; one selects identity mapping. | Debug or deliberately bypass Gaia only. | No at run time; yes if editing the header |
| `AstroCatType` | `int` | `1` | same | `[lensing].astrometry_cat_type` | 1 = legacy large Gaia tiles; 2 = 1-degree Gaia tiles | Select the Stage-1 Gaia catalog filename layout without changing its two-column row format. | Set to 2 when `[lensing].astrometry_cat` points to `process_astrocat` output. | No at run time; yes if editing the header |
| `PROCESS_stage` | `int` | `223092870` | same | `[lensing].process_stage` | Product of selected stage primes 2 through 23; Stage 9 requires Stage 8 | Selects numerical stages by prime-product divisibility. | Change for staged/restart runs. | No at run time; yes if editing the header |
| `include_FLAT` | `int` | `0` | N/A — removed in Lite | `[lensing].include_flat` | 0 = off; 1 = on | Apply super-flat correction when one. | Enable only with valid flat files. | No at run time; yes if editing the header |
| `include_Mask` | `int` | `2` | N/A — removed in Lite | `[lensing].include_mask` | 0 to 3; selects the Standard mask-input mode | Select the DQ-mask input mode. | Change when DQ masks are unavailable or mask mode changes. | No at run time; yes if editing the header |
| `include_BGsub` | `int` | `1` | same | No | `0` off, `1` on | Subtract the fitted science-image background. | Change only for controlled preprocessing experiments. | Yes |
| `ASTROMETRY_CAT` | `std::string` | `"/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted"` | same | `[lensing].astrometry_cat` | Readable Gaia tile directory | Gaia tile directory. | Change for every site/catalog deployment. | No at run time; yes if editing the header |
| `SOURCE_CAT` | `std::string` | `"/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat"` | same | `[lensing].source_cat` or `[extcat].output_directory`; `--extcat-output` | String/path value; see the source comment and validation | Seeds the effective external-tile input/output directory. | Change for another deployment, data set, or path layout. | No at run time; yes if editing the header |
| `FLAT_PATH` | `std::string` | `"/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014"` | N/A — removed in Lite | `[lensing].flat_path` | Readable flat FITS directory | Super-flat FITS directory. | Change when `include_FLAT=1`. | No at run time; yes if editing the header |
| `PSF_PATH` | `std::string` | `"hahahaha"` | N/A — removed in Lite | `[lensing].psf_path` | Readable PSF image directory | External PSF image directory. | Replace before `ext_PSF=1`. | No at run time; yes if editing the header |
| `ext_cat` | `int` | `1` | N/A — removed in Lite | `[lensing].ext_cat` | 0 = off; 1 = use external catalog | Use the external source catalog when one. | Change only for a Standard no-catalog run. | No at run time; yes if editing the header |
| `ext_PSF` | `int` | `0` | N/A — removed in Lite | `[lensing].ext_psf` | 0 = frame stars; 1 = external PSF | Use externally supplied PSF images when one. | Change only with valid external PSFs. | No at run time; yes if editing the header |
| `CCD_split` | `int` | `2` | same | `[lensing].ccd_split` | 1 = whole CCD; 2 = amplifier split | Split each CCD into one or two amplifier regions. | Change for another detector/readout model. | No at run time; yes if editing the header |
| `nct` | `int` | `12` | same | No | Positive rectangle count | Number of background rectangles. | Tune only with background-model validation. | Yes |
| `ncx` | `int` | `3` | same | No | Positive x count | Number of background rectangles along x. | Keep consistent with `nct` geometry. | Yes |
| `npo` | `int` | `64` | same | No | Positive sample count | Exposure PSF sample count. | Adjust PSF sampling only. | Yes |
| `nstar_min` | `int` | `npo * 3 / 2 = 96` | same | No | Derived positive count | Minimum stars for exposure PSF fitting. | Derived parameter — do not edit directly. | Yes |
| `npl` | `int` | `10` | same | No | Non-negative coefficient-count offset | Local PSF polynomial coefficient count minus one. | Adjust local PSF model only. | Yes |
| `nstar_min_local` | `int` | `16` | same | No | Positive count | Minimum retained stars for a local fit. | Tune only with PSF fit validation. | Yes |
| `psf_exposure_min_candidates` | `int` | `60` | same | No | Positive count | Minimum exposure-wide PSF candidates. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_hist_bins` | `int` | `128` | same | No | Integer ≥ 3 | FWHM histogram bin count. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_locus_sigma` | `double` | `4.0` | same | No | Positive sigma multiplier | Exposure FWHM-locus sigma window. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_locus_min_samples` | `int` | `30` | same | No | Positive count | Minimum FWHM-locus samples. | Tune sparse-exposure handling. | Yes |
| `PsfGroupingType` | `int` | `2` | same | No | 1 = threshold graph; 2 = mutual KNN | One selects threshold graph; two selects mutual KNN. | Change for controlled algorithm comparison. | Yes |
| `psf_minchi_reference_fraction` | `double` | `1.0 / 3.0` | same | No | `(0, 1]` | Exposure top-size reference fraction. | Tune Stage 5 threshold estimation. | Yes |
| `psf_minchi_reference_max_per_chip` | `int` | `5` | same | No | Positive count | Reference-star cap per chip. | Tune Stage 5 threshold estimation. | Yes |
| `psf_minchi_sigma_cut` | `double` | `4.0` | same | No | Positive sigma multiplier | Minimum-chi rejection sigma. | Tune Stage 5 selection. | Yes |
| `psf_knn_k` | `int` | `8` | same | No | Positive neighbor count | Neighbors retained by the PSF KNN graph. | Change with grouping validation. | Yes |
| `psf_group_merge_ratio` | `double` | `0.30` | same | No | Non-negative group-size ratio | Secondary-group relative-size threshold. | Tune Stage 5 grouping. | Yes |
| `psf_group_merge_min_gaia` | `int` | `2` | same | No | Non-negative match count | Minimum Gaia matches in a merged group. | Tune Stage 5 grouping. | Yes |
| `psf_gaia_match_radius_pix` | `double` | `2.5` | same | No | Positive pixels | Gaia match radius in pixels. | Change for astrometric precision/pixel scale. | Yes |
| `psf_gaia_locus_min_matches` | `int` | `10` | same | No | Positive match count | Minimum Gaia matches for locus support. | Tune sparse fields. | Yes |
| `psf_press_rejection_enabled` | `bool` | `true` | same | No | Boolean | Enable optional post-fit PRESS cleanup. | Disable for controlled fallback testing. | Yes |
| `psf_press_sigma_cut` | `double` | `4.0` | same | No | Positive sigma multiplier | Standardized PRESS rejection sigma. | Tune only with PSF residual validation. | Yes |
| `psf_press_max_removals` | `int` | `5` | same | No | Non-negative count | Maximum PRESS removals permitted per chip. | Tune only with PSF residual validation. | Yes |
| `psf_loo_min_denom` | `double` | `1.0e-6` | same | No | `(0, 1)` | Minimum leave-one-out denominator. | Numerical guard; normally unchanged. | Yes |
| `step_psf` | `int` | `100` | N/A — removed in Lite | No | Positive pixels | PSF star spatial sampling step. | Change only for the corresponding Standard branch. | Yes |
| `deblending` | `int` | `1` | N/A — removed in Lite | `[lensing].deblending` | 0 = off; 1 = on | Enable source deblending when one. | Change only for controlled source tests. | No at run time; yes if editing the header |
| `n_neighbor` | `int` | `5` | N/A — removed in Lite | No | Positive neighbor count | Neighbor count used by deblending. | Change with deblending validation. | Yes |
| `PSF_type` | `int` | `1` | N/A — removed in Lite | `[lensing].psf_type` | 1 = local polynomial; 2 = hybrid | One selects local polynomial; two selects hybrid PSF. | Change only for Standard hybrid PSF. | No at run time; yes if editing the header |
| `PSF_Ms` | `int` | `0` | N/A — removed in Lite | `[lensing].psf_ms` | 0 = off; 1 = PCA/multi-scale | Enable PCA/multi-scale PSF reconstruction when one. | Change only with PCA inputs/resources. | No at run time; yes if editing the header |
| `ns` | `int` | `64` | same | No | Positive even pixels | Science stamp and Fourier-grid side length. | Change only with coordinated stamp/FFT validation. | Yes |
| `nsns` | `int` | `ns * ns = 4096` | same | No | Derived pixels | Pixels in one science stamp. | Derived parameter — do not edit directly. | Yes |
| `chip_margin` | `int` | `8` | same | No | Non-negative pixels | Extra chip-edge extraction margin. | Change with stamp geometry. | Yes |
| `ns_2` | `int` | `ns / 2 = 32` | same | No | Derived pixels | Half science-stamp side length. | Derived parameter — do not edit directly. | Yes |
| `nl_2` | `int` | `ns_2 + chip_margin = 40` | same | No | Derived pixels | Half expanded extraction side. | Derived parameter — do not edit directly. | Yes |
| `nl` | `int` | `nl_2 * 2 = 80` | same | No | Derived pixels | Full expanded extraction side. | Derived parameter — do not edit directly. | Yes |
| `flag_thresh` | `int` | `3` | same | No | Non-negative extraction flag | Maximum accepted source extraction flag. | Tune source quality selection. | Yes |
| `chip_edge_margin` | `int` | `chip_margin = 8` | same | No | Derived pixels | Alias used by chip-edge checks. | Derived parameter — do not edit directly. | Yes |
| `dz_thresh` | `double` | `0.1` | same | No | Non-negative redshift difference | Redshift tolerance for catalog matching. | Tune matching for another catalog/error model. | Yes |
| `n_user_max` | `int` | `200` | same | No | Positive count | Bright detections used for astrometric matching. | Tune astrometric pattern matching. | Yes |
| `ngal_max` | `int` | `4000` | same | No | Positive reservation hint | Initial galaxy-vector reservation hint. | Change only for allocation tuning. | Yes |
| `nstar_max` | `int` | `2000` | same | No | Positive reservation hint | Initial star-vector reservation hint. | Change only for allocation tuning. | Yes |
| `src_npara` | `int` | `12` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Sets the Stage-6 source-marker field count. | Change only for a validated configuration change. | Yes |
| `npd` | `int` | `33` | same | No | Positive coefficient count | PU astrometric distortion coefficient count. | Change only with astrometric model code. | Yes |
| `blocksize` | `int` | `200` | same | No | Positive pixels | Target background block side length. | Tune background modeling. | Yes |
| `bg_rough_grid_x` | `int` | `32` | same | No | Positive grid count | Rough background grid columns. | Tune background modeling. | Yes |
| `bg_rough_grid_y` | `int` | `32` | same | No | Positive grid count | Rough background grid rows. | Tune background modeling. | Yes |
| `bg_min_block_pixels` | `int` | `1000` | same | No | Positive pixel count | Minimum pixels in a background block. | Tune masked/small images. | Yes |
| `bg_min_clipped_pixels` | `int` | `200` | same | No | Positive pixel count | Minimum pixels after block clipping. | Tune masked/small images. | Yes |
| `bg_min_valid_frac` | `double` | `0.25` | same | No | Fraction in `(0, 1]` | Minimum valid fraction per background block. | Tune masking tolerance. | Yes |
| `bg_clip_low` | `double` | `4.0` | same | No | Positive sigma multiplier | Lower background clipping sigma. | Tune background robustness. | Yes |
| `bg_clip_high` | `double` | `2.5` | same | No | Positive sigma multiplier | Upper background clipping sigma. | Tune background robustness. | Yes |
| `bg_fit_clip_sigma` | `double` | `3.0` | same | No | Positive sigma multiplier | Background-plane fit clipping sigma. | Tune background robustness. | Yes |
| `bg_fit_max_iter` | `int` | `4` | same | No | Non-negative iterations | Maximum background-plane clipping iterations. | Tune convergence/runtime. | Yes |
| `bg_min_fit_factor` | `int` | `3` | same | No | Positive samples-per-coefficient factor | Minimum samples per fitted coefficient factor. | Tune fit stability. | Yes |
| `source_thresh` | `double` | `2.0` | same | No | Positive S/N threshold | Source-detection SNR threshold. | Adjust scientific source selection. | Yes |
| `core_thresh` | `double` | `4.0` | same | No | Positive threshold | Source-core detection threshold. | Adjust scientific source selection. | Yes |
| `NstampType` | `int` | `1` | same | No | 1 = physical blank stamp; 2 = local covariance power | One uses blank stamps; two uses covariance power. | Change for controlled noise-method runs. | Yes |
| `noise_sigma_ratio_min` | `double` | `0.80` | same | No | Positive lower ratio | Minimum blank-to-source sigma ratio. | Tune blank-stamp quality. | Yes |
| `noise_sigma_ratio_max` | `double` | `1.25` | same | No | Above minimum | Maximum blank-to-source sigma ratio. | Tune blank-stamp quality. | Yes |
| `noise_mad_ratio_min` | `double` | `0.70` | same | No | Positive lower ratio | Minimum blank-to-source MAD ratio. | Tune blank-stamp quality. | Yes |
| `noise_mad_ratio_max` | `double` | `1.30` | same | No | Above minimum | Maximum blank-to-source MAD ratio. | Tune blank-stamp quality. | Yes |
| `noise_tail_sigma` | `double` | `2.5` | same | No | Positive sigma | Tail-count sigma threshold. | Tune blank-stamp quality. | Yes |
| `noise_max_tail_fraction` | `double` | `0.05` | same | No | Fraction `[0, 1]` | Maximum blank-stamp tail fraction. | Tune blank-stamp quality. | Yes |
| `noise_max_mask_fraction` | `double` | `0.02` | same | No | Fraction `[0, 1]` | Maximum blank-stamp masked fraction. | Tune blank-stamp quality. | Yes |
| `noise_region_size` | `int` | `192` | same | No | Positive even pixels; greater than inner size | Outer local-noise square side length. | Change with covariance geometry. | Yes |
| `noise_inner_size` | `int` | `96` | same | No | Even pixels; at least `nl` | Central exclusion square side length. | Change with source/stamp geometry. | Yes |
| `noise_plane_min_valid_fraction` | `double` | `0.30` | same | No | Fraction `(0, 1]` | Minimum plane-fit shell fraction. | Tune masking tolerance. | Yes |
| `noise_cov_padding_factor` | `double` | `2.0` | same | No | Positive; padded side must support linear autocorrelation | Covariance FFT padding multiplier. | Change only with FFT validation. | Yes |
| `noise_cov_fft_size` | `int` | `ceil(noise_region_size * noise_cov_padding_factor) = 384` | same | No | Derived padded side | Padded covariance FFT side. | Derived parameter — do not edit directly. | Yes |
| `noise_cov_max_lag` | `int` | `8` | same | No | `0 <= lag < noise_region_size` | Maximum retained signed covariance lag. | Tune covariance model. | Yes |
| `noise_cov_min_valid_pixels` | `int` | `4096` | same | No | Positive count | Minimum covariance-mask pixels. | Tune masking tolerance. | Yes |
| `noise_cov_min_pair_fraction` | `double` | `0.50` | same | No | Fraction `(0, 1]` | Minimum lag pair-count fraction. | Tune covariance reliability. | Yes |
| `noise_cov_sigma_ratio_min` | `double` | `0.80` | same | No | Positive lower ratio | Minimum covariance sigma ratio. | Tune covariance quality. | Yes |
| `noise_cov_sigma_ratio_max` | `double` | `1.25` | same | No | Above minimum | Maximum covariance sigma ratio. | Tune covariance quality. | Yes |
| `noise_cov_max_negative_fraction` | `double` | `0.25` | same | No | Fraction `[0, 1]` | Maximum negative power fraction. | Tune covariance quality. | Yes |
| `noise_cov_imag_tolerance` | `double` | `1.0e-10` | same | No | Non-negative numerical tolerance | Imaginary FFT residual tolerance. | Numerical guard; normally unchanged. | Yes |
| `sig_blocksize` | `int` | `200` | same | No | Positive pixels | Noise-estimator block side length. | Tune only with noise-estimator calibration. | Yes |
| `sig_block_max` | `int` | `sig_blocksize * sig_blocksize = 40000` | same | No | Derived pixel count | Maximum pixels per block. | Derived parameter — do not edit directly. | Yes |
| `sig_max_blocks` | `int` | `2048` | same | No | Positive count | Maximum sampled noise blocks. | Tune memory/runtime only. | Yes |
| `sig_min_block_pixels` | `int` | `1000` | same | No | Positive count | Minimum pixels in one block. | Tune sparse/masked data. | Yes |
| `sig_min_block_triples` | `int` | `1000` | same | No | Positive count | Minimum valid triples per block. | Tune sparse/masked data. | Yes |
| `sig_min_blocks` | `int` | `4` | same | No | Positive count | Minimum blocks for a plane fit. | Tune sparse data only. | Yes |
| `sig_hist_nbin` | `int` | `256` | same | No | Positive bin count | Mode-finding histogram bins. | Tune estimator resolution. | Yes |
| `sig_hist_range` | `double` | `6.0` | same | No | Positive sigma range | Histogram range in sigma units. | Tune estimator robustness. | Yes |
| `sig_min_mode_count` | `int` | `500` | same | No | Positive count | Minimum samples defining the mode. | Tune sparse data only. | Yes |
| `sig_min_lower_count` | `int` | `1000` | same | No | Positive count | Minimum lower-side samples. | Tune sparse data only. | Yes |
| `sig_lower_quantile` | `double` | `0.3173105` | same | No | Quantile in `(0, 1)` | Lower-side width quantile. | Calibration constant; normally unchanged. | Yes |
| `sig_clip_k` | `double` | `3.0` | same | No | Positive sigma multiplier | Symmetric clipping sigma. | Tune estimator robustness. | Yes |
| `sig_rdil` | `int` | `2` | same | No | Positive pixel stride | Pixel stride used by the estimator. | Tune sampling/runtime only. | Yes |
| `sig_clip_niter` | `int` | `2` | same | No | Non-negative iterations | Number of clipping iterations. | Tune convergence/runtime. | Yes |
| `sig_min_fit_triples` | `int` | `1000` | same | No | Positive count | Minimum triples for the final fit. | Tune sparse data only. | Yes |
| `sig_min_fit_frac` | `double` | `0.20` | same | No | Fraction `(0, 1]` | Minimum retained fit fraction. | Tune robustness only. | Yes |
| `sig_median_ratio` | `double` | `1.2678405` | same | No | Positive calibration factor | Median-to-sigma conversion. | Calibration constant; normally unchanged. | Yes |
| `sig_plane_min` | `double` | `1.0e-8` | same | No | Positive floor | Minimum positive sigma-plane value. | Numerical guard; normally unchanged. | Yes |
| `sig_max_plane_ratio` | `double` | `4.0` | same | No | Ratio ≥ 1 | Maximum plane variation ratio. | Tune rejection only with validation. | Yes |
| `sig_pivot_min` | `double` | `1.0e-8` | same | No | Positive floor | Minimum linear-solve pivot. | Numerical guard; normally unchanged. | Yes |
| `sig_scale_s1` | `double` | `0.673475` | same | No | Positive calibration candidate | Stage-1 noise calibration candidate. | Calibration experiments only. | Yes |
| `sig_scale_s2` | `double` | `1.027786` | same | No | Positive calibration value | Stage-2 noise calibration. | Calibration experiments only. | Yes |
| `sig_scale` | `double` | `sig_scale_s2 = 1.027786` | same | No | Derived active selector | Active noise calibration selector. | Derived parameter — do not edit directly. | Yes |
| `area_max` | `int` | `ns * ns = 4096` | same | No | Derived pixels | Maximum connected source area. | Derived parameter — do not edit directly. | Yes |
| `area_thresh` | `int` | `6` | same | No | Positive pixels | Minimum connected source area. | Adjust scientific source selection. | Yes |
| `gal_smooth` | `int` | `0` | same | `[lensing].gal_smooth` | Supported smoothing selector | Galaxy-stamp smoothing radius. | Change for controlled processing tests. | No at run time; yes if editing the header |
| `star_smooth` | `int` | `2` | same | `[lensing].star_smooth` | Supported smoothing selector | Star-stamp smoothing radius. | Change for controlled PSF tests. | No at run time; yes if editing the header |
| `size_fit_rmax` | `int` | `4` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Sets the maximum radius for curvature-size fitting. | Change only with scientific or numerical validation. | Yes |
| `point_stat_beta` | `double` | `0.10` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Sets the survey-wide point-source template shape. | Change only for a validated configuration change. | Yes |
| `point_stat_k_frac` | `double` | `0.90` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Sets the Fourier-radius fraction used by point statistics. | Change only for a validated configuration change. | Yes |
| `point_stat_eps` | `double` | `1.0e-20` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Provides the point-statistic denominator floor. | Change only for a validated configuration change. | Yes |
| `point_stat_min_corr` | `double` | `1.0e-6` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Sets the minimum accepted point-template correlation. | Change only with scientific or numerical validation. | Yes |
| `SNR_PSF` | `double` | `100.0` | same | No | Positive S/N threshold | Minimum PSF-star signal-to-noise ratio. | Adjust PSF star selection. | Yes |
| `saturation_thresh` | `double` | `25000.0` | same | No | Detector-count threshold | Saturated pixel threshold. | Change for detector/gain regime. | Yes |
| `pixel_size` | `double` | `0.2628` | same | `[lensing].pixel_size` | Positive arcsec/pixel | arcsec. | Change only for another detector. | No at run time; yes if editing the header |
| `iid` | `int` | `1 - 1` | same | No | Derived zero-based index | PSF polynomial chi-square field index. | Derived parameter — do not edit directly. | Yes |
| `ipixx` | `int` | `2 - 1` | same | No | Derived zero-based index | Source-center x field index. | Derived parameter — do not edit directly. | Yes |
| `ipixy` | `int` | `3 - 1` | same | No | Derived zero-based index | Source-center y field index. | Derived parameter — do not edit directly. | Yes |
| `isig` | `int` | `4 - 1` | same | No | Derived zero-based index | Local noise sigma field index. | Derived parameter — do not edit directly. | Yes |
| `istar` | `int` | `5 - 1` | same | No | Derived zero-based index | Available PSF-star count field index. | Derived parameter — do not edit directly. | Yes |
| `ipeak` | `int` | `5 - 1` | same | No | Derived legacy alias | Historical peak alias field index. | Derived parameter — do not edit directly. | Yes |
| `i_imax` | `int` | `6 - 1` | same | No | Derived zero-based index | Peak x field index. | Derived parameter — do not edit directly. | Yes |
| `i_jmax` | `int` | `7 - 1` | same | No | Derived zero-based index | Peak y field index. | Derived parameter — do not edit directly. | Yes |
| `ih_flux` | `int` | `8 - 1` | same | No | Derived zero-based index | Half-light flux field index. | Derived parameter — do not edit directly. | Yes |
| `ih_area` | `int` | `9 - 1` | same | No | Derived zero-based index | Source area field index. | Derived parameter — do not edit directly. | Yes |
| `iflag` | `int` | `10 - 1` | same | No | Derived zero-based index | Quality flag field index. | Derived parameter — do not edit directly. | Yes |
| `iPSF` | `int` | `11 - 1` | same | No | Derived zero-based index | Local PSF size field index. | Derived parameter — do not edit directly. | Yes |
| `iSNR_F` | `int` | `12 - 1` | same | No | Derived zero-based index | Fourier SNR field index. | Derived parameter — do not edit directly. | Yes |
| `ira` | `int` | `13 - 1` | same | No | Derived zero-based index | Right-ascension field index. | Derived parameter — do not edit directly. | Yes |
| `idec` | `int` | `14 - 1` | same | No | Derived zero-based index | Declination field index. | Derived parameter — do not edit directly. | Yes |
| `igf1` | `int` | `15 - 1` | same | No | Derived zero-based index | Field-distortion g1 index. | Derived parameter — do not edit directly. | Yes |
| `igf2` | `int` | `16 - 1` | same | No | Derived zero-based index | Field-distortion g2 index. | Derived parameter — do not edit directly. | Yes |
| `ig1` | `int` | `17 - 1` | same | No | Derived zero-based index | Fourier_Quad g1 estimator index. | Derived parameter — do not edit directly. | Yes |
| `ig2` | `int` | `18 - 1` | same | No | Derived zero-based index | Fourier_Quad g2 estimator index. | Derived parameter — do not edit directly. | Yes |
| `ide` | `int` | `19 - 1` | same | No | Derived zero-based index | Shear response estimator index. | Derived parameter — do not edit directly. | Yes |
| `ih1` | `int` | `20 - 1` | same | No | Derived zero-based index | Higher-order h1 estimator index. | Derived parameter — do not edit directly. | Yes |
| `ih2` | `int` | `21 - 1` | same | No | Derived zero-based index | Higher-order h2 estimator index. | Derived parameter — do not edit directly. | Yes |
| `icos2` | `int` | `22 - 1` | same | No | Derived zero-based index | Spin-2 cosine field index. | Derived parameter — do not edit directly. | Yes |
| `isin2` | `int` | `23 - 1` | same | No | Derived zero-based index | Spin-2 sine field index. | Derived parameter — do not edit directly. | Yes |
| `iparity` | `int` | `24 - 1` | same | No | Derived zero-based index | WCS parity field index. | Derived parameter — do not edit directly. | Yes |
| `igalsizeT` | `int` | `25 - 1` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Curvature galaxy-size field index. | Derived parameter — do not edit directly. | Yes |
| `ipsfsizeT` | `int` | `26 - 1` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Curvature PSF-size field index. | Derived parameter — do not edit directly. | Yes |
| `idelta_chi2` | `int` | `27 - 1` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Point-template delta-chi-square field index. | Derived parameter — do not edit directly. | Yes |
| `iorth_ext` | `int` | `28 - 1` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Orthogonal point-template extension field index. | Derived parameter — do not edit directly. | Yes |
| `shear_cat_ncols` | `int` | `iorth_ext + 1 = 28` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Derives the 28-field Stage-7 catalog width. | Derived parameter — do not edit directly. | Yes |
| `expo_cat_ncols` | `int` | `shear_cat_ncols + 1 = 29` | same | No | Numeric or derived value; keep source assertions and consumers consistent | Derives the 29-field exposure-catalog width. | Derived parameter — do not edit directly. | Yes |
| `ichi2` | `int` | `shear_cat_ncols = 28` | same | No | Derived zero-based index | Exposure chi-square field index. | Derived parameter — do not edit directly. | Yes |
| `DEFAULT_CHIP_COUNT` | `int` | `62` | same | `[lensing].nmax_chip` | Numeric or derived value; keep source assertions and consumers consistent | Seeds RuntimeConfig's CCD count limit. | Change only for a validated configuration change. | No at run time; yes if editing the header |
| `g1_c` | `double` | `-0.001` | same | No | Additive calibration | Additive field-distortion g1 correction. | Recalibrate for another dataset/band. | Yes |
| `g2_c` | `double` | `-0.0003` | same | No | Additive calibration | Additive field-distortion g2 correction. | Recalibrate for another dataset/band. | Yes |
| `chi2_thresh` | `double` | `0.01` | same | No | Non-negative threshold | Maximum exposure PSF chi-square. | Adjust scientific quality selection. | Yes |
| `chipnx` | `int` | `2046` | same | `[lensing].chipnx` | Positive pixels | Science CCD width used for PSF coordinates. | Change only for another detector. | No at run time; yes if editing the header |
| `chipny` | `int` | `4094` | same | `[lensing].chipny` | Positive pixels | Science CCD height used for PSF coordinates. | Change only for another detector. | No at run time; yes if editing the header |
| `rescale_size` | `double` | `1.2` | N/A — removed in Lite | No | Positive target size | Target PSF residual rescaling size. | Change only for `PSF_Ms=1`. | Yes |
| `procs_pn` | `int` | `40` | N/A — removed in Lite | No | Positive rank count | MPI ranks per PCA scheduling group. | Tune Standard PCA scheduling. | Yes |
| `work_pn` | `int` | `10` | N/A — removed in Lite | No | Positive worker count | Concurrent PCA workers per group. | Tune Standard PCA scheduling. | Yes |
| `nblocks` | `int` | `2` | N/A — removed in Lite | No | Positive blocks per CCD axis | PCA spatial blocks per CCD axis. | Tune Standard PCA modeling. | Yes |
| `n_pcs` | `int` | `100` | N/A — removed in Lite | No | Positive component count | Maximum PCA principal components. | Tune Standard PCA modeling. | Yes |
| `npp6th` | `int` | `28` | N/A — removed in Lite | No | Polynomial term count | Sixth-degree 2D polynomial term count. | Derived by polynomial basis; normally unchanged. | Yes |
| `pca_negative_eigenvalue_threshold` | `double` | `-1.0e-5` | N/A — removed in Lite | No | Non-positive tolerance | Invalid PCA eigenvalue cutoff. | Numerical guard; normally unchanged. | Yes |

Lite accepts only the `[lensing]` keys whose symbols survive in its table. The
implementation-local `PSFr_ratio = 0.75` in `ShearMeasurement.cpp` is not a
configuration-header parameter and is intentionally excluded.

## `ProcessRearrConfig` (`config/ProcessRearrConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `SKY_GRID_DEGREES` | `double` | `0.1` | same | No | Positive degrees | Full-sky tile width in degrees. | Change for another spatial partition resolution. | Yes |
| `RA_BIN_COUNT` | `int` | `3600` | same | No | Positive full-sky bin count | Number of right-ascension bins. | Keep consistent with grid width. | Yes |
| `DEC_BIN_COUNT` | `int` | `1800` | same | No | Positive full-sky bin count | Number of declination bins. | Keep consistent with grid width. | Yes |
| `SKY_TILE_COUNT` | `std::size_t` | `RA_BIN_COUNT * DEC_BIN_COUNT = 6480000` | same | No | `RA_BIN_COUNT * DEC_BIN_COUNT` | Total full-sky tile count. | Derived parameter — do not edit directly. | Yes |
| `TARGET_SUBCAT_ROWS` | `std::uint64_t` | `500000` | same | No | Positive row count | Target rows per partition. | Tune output file size and memory. | Yes |
| `SKIP_DIRECTORY_NAME` | `std::string_view` | `"Large_Field"` | same | No | Directory name | Directory excluded from scans. | Change only when layout conventions change. | Yes |
| `SUBCAT_PREFIX` | `std::string_view` | `"subcat_"` | same | No | Filename prefix | Partition filename prefix. | Change for another naming convention. | Yes |
| `SUBCAT_EXTENSION` | `std::string_view` | `".cat"` | same | No | Filename extension | Partition filename extension. | Change for another naming convention. | Yes |
| `SUBCAT_ID_WIDTH` | `int` | `6` | same | No | Positive digit count | Minimum zero-padded partition ID width. | Change for another naming convention. | Yes |
| `SUMMARY_FILENAME` | `std::string_view` | `"catalog_summary.txt"` | same | No | Filename | Summary report name. | Change for another naming convention. | Yes |
| `OUTPUT_PRECISION` | `int` | `10` | same | No | Positive significant digits | Significant digits in catalog rows. | Change only after precision/size review. | Yes |
| `SUMMARY_PRECISION` | `int` | `4` | same | No | Non-negative decimal places | Decimal places in summary bounds. | Change for reporting needs. | Yes |
| `SKIP_MISSING_CATALOGS` | `bool` | `true` | same | No | Boolean | Continue past absent input catalogs. | Set false for strict completeness checks. | Yes |
| `SKIP_MALFORMED_ROWS` | `bool` | `true` | same | No | Boolean | Continue past malformed catalog rows. | Set false for strict schema checks. | Yes |

Catalog widths and column offsets are resolved by runtime `CatalogLayout`; they
are not parameters in this header. The canonical complete-row order is the
effective external prefix, `EXPO_NUM`, `ccD_NUM`, 28 Stage-7 fields, then the
exposure `Chi2`. With the default 18-field external prefix this is 49 fields:
zero-based `EXPO_NUM=18`, `ccD_NUM=19`, source base `20`, and `Chi2=48`.
Standard's `ext_cat=0` rearrangement schema has 31 fields using the same two
identity fields. Catalogs produced by the former 48/30-field schemas must be
regenerated before rearrangement or FD processing.

## `OutputLayout` (`config/pathconfig.hpp`)

These arrays are identical in Standard and Lite and have no RuntimeConfig, INI,
or CLI override. They contain relative directory names, not deployment roots;
the runtime dataset root is prepended by the existing path helpers.

| Parameter | Type | Standard / Lite compiled value | Function | Rebuild after change |
|---|---|---|---|---|
| `NON_CHIP_BASE_DIRECTORIES` | `std::array<const char*, 14>` | `science`, `dqmask`, `stamps`, `result`, `stamps/dat_StarInfo`, `stamps/fits_StarP`, `stamps/fits_PsfSrc`, `stamps/dat_ExpoInfo`, `stamps/dat_StarComp`, `stamps/dat_Rescale`, `stamps/dat_Pcs`, `stamps/dat_StarCompV2`, `astrometry/Head`, `astrometry/dat_Chk` | Complete fixed base-directory contract created without a chip suffix. | Yes |
| `CHIP_PRODUCT_DIRECTORIES` | `std::array<const char*, 16>` | `stamps/Norm`, `stamps/cat_Orig`, `stamps/dat_StarCanInfo`, `stamps/fits_StarCan`, `stamps/fits_StarCanN`, `stamps/fits_StarCanP`, `stamps/dat_SrcInfo`, `stamps/fits_Src`, `stamps/fits_Noise`, `stamps/fits_SrcP`, `stamps/dat_PsfFit`, `stamps/fits_PsfLocal`, `stamps/dat_Shear`, `stamps/dat_StarXY`, `stamps/fits_PsfResi`, `astrometry/dat_Astro` | Complete fixed per-chip product-directory contract. | Yes |

`include/general/OutputLayout.hpp` now contains only the functions that derive
exposure and chip paths from these centralized arrays.

## `config/FDConfig.hpp`

| Parameter | Type | Standard default | Lite default | INI / CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `FD_STATIC_MODE` | `StaticMode` | `StaticMode::PDF_SIGMA` | same | No | PDF_SIGMA, PDF_JACK, or SWSE_JACK | Active FD estimator mode. | Change for another validated statistical method. | Yes |
| `FD_USE_PDF_STATIS` | `bool` | `(FD_STATIC_MODE == StaticMode::PDF_SIGMA \|\| FD_STATIC_MODE == StaticMode::PDF_JACK)` | same | No | Derived from `FD_STATIC_MODE` | Enables PDF sign-test path. | Derived parameter — do not edit directly. | Yes |
| `FD_USE_JACKKNIFE` | `bool` | `(FD_STATIC_MODE == StaticMode::PDF_JACK \|\| FD_STATIC_MODE == StaticMode::SWSE_JACK)` | same | No | Derived from `FD_STATIC_MODE` | Enables jackknife uncertainty. | Derived parameter — do not edit directly. | Yes |
| `FD_USE_SWSE_DATA` | `bool` | `(FD_STATIC_MODE == StaticMode::SWSE_JACK)` | same | No | Derived from `FD_STATIC_MODE` | Enables SWSE data model. | Derived parameter — do not edit directly. | Yes |
| `FD_PER_EXPOSURE_STAR_BAR` | `bool` | `false` | same | No | `false` global N-to-1, `true` per-exposure N-to-N | Fit one star bar per exposure when true. | Change for another validated FD mode. | Yes |
| `nmax_per_core` | `int` | `20000000` | same | No | Positive source capacity | Maximum sources reserved per MPI rank. | Tune memory/capacity for workload. | Yes |
| `fd_num` | `int` | `21` | same | No | Positive bin count | spatial bins by field distortion. | Tune FD spatial resolution. | Yes |
| `PDF_BINS` | `int` | `4` | same | No | Positive inner-bin count | equal-probability inner bins. | Tune PDF estimator resolution. | Yes |
| `gf_lim` | `float` | `0.0015` | same | No | Positive distortion half-range | spatial bin range ±gf_lim. | Change for another distortion range. | Yes |
| `NMAX` | `int` | `200` | same | No | Positive sampling count | fine grid sampling points. | Tune fit resolution/runtime. | Yes |
| `MAX_DUP` | `int` | `5` | same | No | Positive duplicate limit | max duplicate measurements. | Change for another catalog duplication policy. | Yes |
| `N_jack` | `int` | `50` | same | No | Positive region count | jackknife regions. | Tune covariance resolution. | Yes |
| `nmax_total` | `int` | `1000000` | same | No | Positive source count | max total sources for k-means. | Tune memory/runtime. | Yes |
| `Km_iter` | `int` | `100` | same | No | Positive iterations | k-means iterations. | Tune convergence/runtime. | Yes |
| `snrfcut` | `float` | `4.0` | same | No | Non-negative S/N | Minimum Fourier signal-to-noise ratio. | Adjust scientific selection. | Yes |
| `snrlow` | `float` | `0.0` | same | No | `0` disables or lower S/N bound | Optional lower source-SNR bound; zero disables it. | Enable for a bounded S/N sample. | Yes |
| `snrhigh` | `float` | `0.0` | same | No | `0` disables or upper S/N bound | Optional upper source-SNR bound; zero disables it. | Enable for a bounded S/N sample. | Yes |
| `starcut` | `float` | `20.0` | same | No | Size threshold | Point-source size cut. | Recalibrate stellar selection. | Yes |
| `chi2_thresh` | `float` | `0.01` | same | No | Non-negative threshold | Maximum exposure chi-square. | Adjust scientific quality selection. | Yes |
| `flagcut` | `float` | `0.0` | same | No | Maximum accepted flag | Maximum accepted source quality flag. | Adjust scientific quality selection. | Yes |
| `imaxcut` | `float` | `64.0` | same | No | Pixel-coordinate upper cut | Maximum source peak x coordinate. | Change with stamp geometry. | Yes |
| `jmaxcut` | `float` | `64.0` | same | No | Pixel-coordinate upper cut | Maximum source peak y coordinate. | Change with stamp geometry. | Yes |
| `zplow` | `float` | `0.0` | same | No | Lower `zp` bound | Minimum photometric redshift. | Adjust redshift/sample selection. | Yes |
| `zphigh` | `float` | `3.0` | same | No | Upper `zp` bound above `zplow` | Maximum photometric redshift. | Adjust redshift/sample selection. | Yes |
| `r_half_thresh` | `float` | `0.0` | same | No | `0` disables or size threshold | Optional half-light-radius threshold. | Enable for a size-selected sample. | Yes |
| `star_bar_mltp` | `float` | `3.0` | same | No | Positive sigma multiplier | Stellar-locus sigma multiplier. | Recalibrate stellar selection. | Yes |
| `psf_chi2_mltp` | `float` | `3.0` | same | No | Positive sigma multiplier | PSF chi-square sigma multiplier. | Recalibrate PSF quality selection. | Yes |
| `n_size_bins` | `int` | `100` | same | No | Positive bin count | Stellar-size histogram bins. | Tune stellar-locus resolution. | Yes |
| `n_mag_bins` | `int` | `20` | same | No | Positive bin count | Magnitude histogram bins. | Tune stellar-locus resolution. | Yes |
| `size_min` | `float` | `-2.0` | same | No | Lower size bound | Minimum histogram size coordinate. | Recalibrate stellar locus. | Yes |
| `size_max` | `float` | `2.0` | same | No | Above `size_min` | Maximum histogram size coordinate. | Recalibrate stellar locus. | Yes |
| `mag_min_val` | `float` | `10.0` | same | No | Lower magnitude bound | Minimum histogram magnitude. | Change for another band/depth. | Yes |
| `mag_max_val` | `float` | `30.0` | same | No | Above `mag_min_val` | Maximum histogram magnitude. | Change for another band/depth. | Yes |
| `min_bin_count` | `int` | `100` | same | No | Positive sample count | Minimum samples in a usable bin. | Tune sparse samples. | Yes |
| `peak_match_tol` | `float` | `0.05` | same | No | Positive size tolerance | Stellar-peak matching tolerance. | Recalibrate stellar locus. | Yes |
| `min_concentration` | `float` | `0.6` | same | No | Fraction `[0, 1]` | Minimum stellar-locus concentration. | Recalibrate stellar locus. | Yes |
| `star_phy_min` | `float` | `-0.5` | same | No | Lower physical-size bound | Minimum physical stellar size. | Recalibrate stellar locus. | Yes |
| `star_phy_max` | `float` | `0.2` | same | No | Above `star_phy_min` | Maximum physical stellar size. | Recalibrate stellar locus. | Yes |
| `stage1_snr` | `float` | `40.0` | same | No | Non-negative S/N | SNR for histogram accumulation. | Tune per-exposure star cuts. | Yes |
| `stage2_snr` | `float` | `0.0` | same | No | `0` disables or S/N threshold | Second-pass SNR threshold; zero disables it. | Enable only for a validated second pass. | Yes |
| `init_win_active` | `float` | `0.1` | same | No | Positive window | Active exposure initial size window. | Tune per-exposure star cuts. | Yes |
| `init_win_fallback` | `float` | `0.15` | same | No | Positive window | Fallback exposure initial size window. | Tune fallback behavior. | Yes |
| `default_s_init` | `float` | `0.5` | same | No | Size coordinate | Default initial stellar-size center. | Recalibrate stellar locus. | Yes |
| `clip_nsigma` | `float` | `3.0` | same | No | Positive sigma multiplier | Iterative clipping sigma. | Tune per-exposure robustness. | Yes |
| `min_clip_limit` | `float` | `0.015` | same | No | Positive half-width | Minimum clipping half-width. | Tune per-exposure robustness. | Yes |
| `default_s_std` | `float` | `0.05` | same | No | Positive scatter | Default stellar-size scatter. | Recalibrate stellar locus. | Yes |
| `fallback_scut_default` | `float` | `0.6` | same | No | Size threshold | Fallback stellar-size cut. | Recalibrate fallback behavior. | Yes |
| `bad_ccds` | `int[]` | `{2, 31, 53, 61}` | same | No | Detector-specific CCD IDs | DES CCD numbers excluded from analysis. | Change for another detector/quality list. | Yes |
| `n_bad_ccds` | `int` | `4` | same | No | Derived array length | Number of excluded CCDs. | Derived parameter — do not edit directly. | Yes |
| `chip_xmin` | `int` | `50` | same | No | Pixel lower bound | Minimum accepted chip x coordinate. | Change for detector/edge-mask policy. | Yes |
| `chip_xmax` | `int` | `1990` | same | No | Pixel upper bound | Maximum accepted chip x coordinate. | Change for detector/edge-mask policy. | Yes |
| `chip_ymin` | `int` | `100` | same | No | Pixel lower bound | Minimum accepted chip y coordinate. | Change for detector/edge-mask policy. | Yes |
| `chip_ymax` | `int` | `3990` | same | No | Pixel upper bound | Maximum accepted chip y coordinate. | Change for detector/edge-mask policy. | Yes |

FD consumes the runtime catalog layout. External-prefix and pipeline-column
offsets are therefore not duplicated in this header; per-exposure mode reads
the serialized `EXPO_NUM` identity rather than inferring it from file order or
the terminal `Chi2` value.
