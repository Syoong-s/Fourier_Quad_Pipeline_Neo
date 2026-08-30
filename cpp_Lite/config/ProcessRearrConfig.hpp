#ifndef PROCESS_REARR_CONFIG_HPP
#define PROCESS_REARR_CONFIG_HPP

#include "pathconfig.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace ProcessRearrConfig {

// ==========================================
// Configuration: Spatial partitioning and output defaults
// Method: Preserve the F77 0.1-degree full-sky grid and approximately
//         500,000 rows per deterministic weighted k-d partition.
// ==========================================
inline constexpr double SKY_GRID_DEGREES = 0.1;  // Full-sky tile width in degrees.
inline constexpr int RA_BIN_COUNT = 3600;  // Number of right-ascension bins.
inline constexpr int DEC_BIN_COUNT = 1800;  // Number of declination bins.
inline constexpr std::size_t SKY_TILE_COUNT =
    static_cast<std::size_t>(RA_BIN_COUNT) * DEC_BIN_COUNT;  // Total full-sky tile count.
inline constexpr std::uint64_t TARGET_SUBCAT_ROWS = 500000;  // Target rows per partition.
inline constexpr int SUBCAT_ID_WIDTH = 6;  // Minimum zero-padded partition ID width.
inline constexpr int OUTPUT_PRECISION = 10;  // Significant digits in catalog rows.
inline constexpr int SUMMARY_PRECISION = 4;  // Decimal places in summary bounds.
inline constexpr bool SKIP_MISSING_CATALOGS = true;  // Continue past absent input catalogs.
inline constexpr bool SKIP_MALFORMED_ROWS = true;  // Continue past malformed catalog rows.

static_assert(SKY_TILE_COUNT
                  <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "the full-sky MPI reduction count must fit int");

}  // namespace ProcessRearrConfig

#endif  // PROCESS_REARR_CONFIG_HPP
