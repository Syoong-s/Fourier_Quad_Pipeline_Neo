#ifndef PATHCONFIG_HPP
#define PATHCONFIG_HPP

// ==========================================
// Configuration: Compile-time input and output path defaults
// Method: Keep every user-selectable path and fixed output-layout name in one
//         header while preserving the established configuration namespaces and
//         derived relationships copied into RuntimeConfig at startup.
// ==========================================

#include "Initialize.hpp"

#include <array>
#include <string>
#include <string_view>


namespace ProcessConfig {

// Workflow path defaults copied into RuntimeConfig before INI/CLI overrides.
inline constexpr const char* EXPO_LIST = Initialize::EXPO_LIST;  // Top-level exposure-list path.
inline constexpr const char* REARR_OUTPUT_DIRECTORY = "baked";  // Rearranged catalogs.
inline constexpr const char* REARR_OUTPUT_BASE_DIRECTORY = "";  // Empty uses dataset root.
inline constexpr const char* REARRANGED_EXPO_LIST_FILENAME =
    "cat_gband_ori.list";  // Published rearranged list name.
inline constexpr const char* REARRANGED_EXPO_LIST_DIRECTORY = "";  // Input-list parent.
inline constexpr const char* FD_EXPO_LIST = "";  // Optional FD-specific list.
inline constexpr const char* FD_OUTPUT_DIRECTORY = "fdout";  // FD results.
inline constexpr const char* FD_OUTPUT_BASE_DIRECTORY = "";  // Empty uses dataset root.

}  // namespace ProcessConfig

namespace LensingConfig {

// Catalog inputs originally supplied by para.inc; optional flat/PSF branches are absent in Lite.
inline const std::string ASTROMETRY_CAT = Initialize::ASTROMETRY_CAT;  // Gaia tiles.
inline constexpr const char* SOURCE_CAT_DEFAULT =
    Initialize::SOURCE_CAT_DEFAULT;  // External source tiles.

}  // namespace LensingConfig

namespace InitConfig {

inline constexpr const char* SCIENCE_ROOT = Initialize::SCIENCE_ROOT;
inline constexpr const char* DQ_ROOT = Initialize::DQ_ROOT;
inline constexpr const char* OUTPUT_ROOT = Initialize::OUTPUT_ROOT;

}  // namespace InitConfig

namespace AstroCatConfig {

inline constexpr std::string_view ASTROMETRY_TILE_PREFIX = Initialize::ASTROMETRY_TILE_PREFIX;
inline constexpr const char* ASTROCAT_INPUT_DIRECTORY = Initialize::ASTROCAT_INPUT_DIRECTORY;
inline constexpr const char* ASTROCAT_OUTPUT_DIRECTORY = "";  // Advanced producer output.

}  // namespace AstroCatConfig

namespace ExtCatConfig {

inline constexpr std::string_view SOURCE_CAT_TILE_PREFIX = Initialize::SOURCE_CAT_TILE_PREFIX;
inline constexpr const char* EXTCAT_INPUT_DIRECTORY = Initialize::EXTCAT_INPUT_DIRECTORY;
inline constexpr const char* EXTCAT_OUTPUT_DIRECTORY = "";  // Advanced producer output.

}  // namespace ExtCatConfig

namespace ProcessRearrConfig {

inline constexpr std::string_view SKIP_DIRECTORY_NAME = "Large_Field";  // Excluded scan root.
inline constexpr std::string_view SUBCAT_PREFIX = "subcat_";  // Partition prefix.
inline constexpr std::string_view SUBCAT_EXTENSION = ".cat";  // Partition extension.
inline constexpr std::string_view SUMMARY_FILENAME = "catalog_summary.txt";  // Summary output.

}  // namespace ProcessRearrConfig

namespace OutputLayout {

// Complete process_init base-directory contract without per-chip products.
inline constexpr std::array<const char*, 15> NON_CHIP_BASE_DIRECTORIES = {
    "science",
    "dqmask",
    "stamps",
    "result",
    "stamps/dat_StarInfo",
    "stamps/svg_StarLocus",
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
