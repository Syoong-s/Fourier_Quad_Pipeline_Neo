#ifndef PROCESS_REARR_CONFIG_HPP
#define PROCESS_REARR_CONFIG_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace ProcessRearrConfig {

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

static_assert(SKY_TILE_COUNT
                  <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "the full-sky MPI reduction count must fit int");

}  // namespace ProcessRearrConfig

#endif  // PROCESS_REARR_CONFIG_HPP
