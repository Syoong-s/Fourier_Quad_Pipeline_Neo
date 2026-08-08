#ifndef PROCESS_REARR_CONFIG_HPP
#define PROCESS_REARR_CONFIG_HPP

#include "ProcessConfig.hpp"
#include "ExtCatConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace ProcessRearrConfig {

// ==========================================
// Configuration: Derived _all.cat column layout
// Method: Treat ichi2 as the 25 fields appended by process_main after CCD_NUM,
//         then apply external columns + 1 CCD column + ichi2 exactly once.
// ==========================================
inline constexpr std::size_t ichi2 =
    static_cast<std::size_t>(LensingConfig::ichi2) + 1;
inline constexpr std::size_t CCD_COLUMN_COUNT = 1;
inline constexpr std::size_t ALL_CAT_TOTAL_COLUMNS =
    ExtCatConfig::EXTCAT_TOTAL_COLUMNS + CCD_COLUMN_COUNT + ichi2;

// ==========================================
// Configuration: Spatial partitioning and output defaults
// Method: Preserve the F77 0.1-degree full-sky grid and approximately
//         500,000 rows per deterministic weighted k-d partition.
// ==========================================
inline constexpr double SKY_GRID_DEGREES = 0.1;
inline constexpr int RA_BIN_COUNT = 3600;
inline constexpr int DEC_BIN_COUNT = 1800;
inline constexpr std::size_t SKY_TILE_COUNT =
    static_cast<std::size_t>(RA_BIN_COUNT) * DEC_BIN_COUNT;
inline constexpr std::uint64_t TARGET_SUBCAT_ROWS = 500000;
inline constexpr std::string_view SKIP_DIRECTORY_NAME = "Large_Field";
inline constexpr std::string_view SUBCAT_PREFIX = "subcat_";
inline constexpr std::string_view SUBCAT_EXTENSION = ".cat";
inline constexpr int SUBCAT_ID_WIDTH = 6;
inline constexpr std::string_view SUMMARY_FILENAME = "catalog_summary.txt";
inline constexpr int OUTPUT_PRECISION = 10;
inline constexpr int SUMMARY_PRECISION = 4;
inline constexpr bool SKIP_MISSING_CATALOGS = true;
inline constexpr bool SKIP_MALFORMED_ROWS = true;

// ==========================================
// Function: Determine the external-field width emitted by process_extcat
// Method: Use the configured pass-through width, or the explicit projection
//         length because projected output contains exactly those fields.
// ==========================================
inline std::size_t externalCatalogColumns(
    const ProcessConfig::RuntimeOptions& options) {
    return options.extcat_use_explicit_columns
               ? options.extcat_input_columns_one_based.size()
               : ExtCatConfig::EXTCAT_TOTAL_COLUMNS;
}

// ==========================================
// Function: Compute the complete _all.cat row width
// Method: Apply external catalog total columns + 1 + ichi2 in the
//         process_rearr-specific parameter header as required.
// ==========================================
inline std::size_t allCatalogColumns(
    const ProcessConfig::RuntimeOptions& options) {
    return externalCatalogColumns(options) + CCD_COLUMN_COUNT + ichi2;
}

static_assert(ichi2 == 25, "process_main must append 25 fields through exposure Chi2");
static_assert(ExtCatConfig::EXTCAT_TOTAL_COLUMNS > 0,
              "the external catalog column count must be positive");
static_assert(SKY_TILE_COUNT
                  <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "the full-sky MPI reduction count must fit int");

}  // namespace ProcessRearrConfig

#endif  // PROCESS_REARR_CONFIG_HPP
