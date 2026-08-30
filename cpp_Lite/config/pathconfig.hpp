#ifndef PATHCONFIG_HPP
#define PATHCONFIG_HPP

// ==========================================
// Configuration: Compile-time input and output path defaults
// Method: Keep every user-selectable path and fixed output-layout name in one
//         header while preserving the established configuration namespaces and
//         derived relationships copied into RuntimeConfig at startup.
// ==========================================

#include <array>
#include <string>
#include <string_view>

namespace LensingConfig {

// Catalog inputs originally supplied by para.inc; optional flat/PSF branches are absent in Lite.
inline const std::string ASTROMETRY_CAT =
    "/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted";  // Gaia tiles.
inline const std::string SOURCE_CAT =
    "/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat";  // External source tiles.

}  // namespace LensingConfig

namespace AstroCatConfig {

inline constexpr const char* ASTROCAT_INPUT_DIRECTORY = "";  // Raw two-column Gaia files.
inline const std::string ASTROCAT_OUTPUT_DIRECTORY =
    LensingConfig::ASTROMETRY_CAT;  // Generated one-degree tiles.

}  // namespace AstroCatConfig

namespace ExtCatConfig {

inline constexpr const char* EXTCAT_INPUT_DIRECTORY = "";  // Raw catalog root.
inline const std::string& EXTCAT_OUTPUT_DIRECTORY =
    LensingConfig::SOURCE_CAT;  // Tile output and process_main input.

}  // namespace ExtCatConfig

namespace InitConfig {

inline constexpr const char* SCIENCE_ROOT =
    "/lustre/home/acct-phyzj/share/DES/g";  // Science archive root.
inline constexpr const char* DQ_ROOT =
    "/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask";  // DQ archive root.
inline constexpr const char* OUTPUT_ROOT =
    "/lustre/home/acct-phyzj/share/DES/g_band_v1";  // Pipeline output root.

}  // namespace InitConfig

namespace ProcessConfig {

// Workflow path defaults copied into RuntimeConfig before INI/CLI overrides.
inline constexpr const char* EXPO_LIST = "";  // Top-level exposure-list path.
inline constexpr const char* REARR_OUTPUT_DIRECTORY = "baked";  // Rearranged catalogs.
inline constexpr const char* REARR_OUTPUT_BASE_DIRECTORY = "";  // Empty uses dataset root.
inline constexpr const char* REARRANGED_EXPO_LIST_FILENAME =
    "cat_gband_ori.list";  // Published rearranged list name.
inline constexpr const char* REARRANGED_EXPO_LIST_DIRECTORY = "";  // Input-list parent.
inline constexpr const char* FD_EXPO_LIST = "";  // Optional FD-specific list.
inline constexpr const char* FD_OUTPUT_DIRECTORY = "fdout";  // FD results.
inline constexpr const char* FD_OUTPUT_BASE_DIRECTORY = "";  // Empty uses dataset root.

}  // namespace ProcessConfig

namespace ProcessRearrConfig {

inline constexpr std::string_view SKIP_DIRECTORY_NAME = "Large_Field";  // Excluded scan root.
inline constexpr std::string_view SUBCAT_PREFIX = "subcat_";  // Partition prefix.
inline constexpr std::string_view SUBCAT_EXTENSION = ".cat";  // Partition extension.
inline constexpr std::string_view SUMMARY_FILENAME = "catalog_summary.txt";  // Summary output.

}  // namespace ProcessRearrConfig

namespace OutputLayout {

// Complete process_init base-directory contract without per-chip products.
inline constexpr std::array<const char*, 14> NON_CHIP_BASE_DIRECTORIES = {
    "science",
    "dqmask",
    "stamps",
    "result",
    "stamps/dat_StarInfo",
    "stamps/fits_StarP",
    "stamps/fits_PsfSrc",
    "stamps/dat_ExpoInfo",
    "stamps/dat_StarComp",
    "stamps/dat_Rescale",
    "stamps/dat_Pcs",
    "stamps/dat_StarCompV2",
    "astrometry/Head",
    "astrometry/dat_Chk",
};

// Complete process_main per-chip product-directory contract.
inline constexpr std::array<const char*, 16> CHIP_PRODUCT_DIRECTORIES = {
    "stamps/Norm",
    "stamps/cat_Orig",
    "stamps/dat_StarCanInfo",
    "stamps/fits_StarCan",
    "stamps/fits_StarCanN",
    "stamps/fits_StarCanP",
    "stamps/dat_SrcInfo",
    "stamps/fits_Src",
    "stamps/fits_Noise",
    "stamps/fits_SrcP",
    "stamps/dat_PsfFit",
    "stamps/fits_PsfLocal",
    "stamps/dat_Shear",
    "stamps/dat_StarXY",
    "stamps/fits_PsfResi",
    "astrometry/dat_Astro",
};

}  // namespace OutputLayout

#endif  // PATHCONFIG_HPP
