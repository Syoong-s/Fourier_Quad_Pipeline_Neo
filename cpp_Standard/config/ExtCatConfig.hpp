#ifndef EXTCAT_CONFIG_HPP
#define EXTCAT_CONFIG_HPP

// ==========================================
// ExtCatConfig - External source-catalog repartitioning defaults
// Method: Derive the tile output from the main pipeline SOURCE_CAT while keeping raw-input,
//         parsing, and optional projection controls together for the first workflow phase.
//
// Note: CatalogLayout requires RA, Dec, and ZP. Each magnitude is optional:
//       set its raw one-based column to zero when that band is unavailable.
//       Non-required fields may remain in the physical row but are not
//       represented in the layout or consumed by downstream algorithms.
// ==========================================

#include "LensingConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ExtCatConfig {

inline constexpr const char* EXTCAT_INPUT_DIRECTORY = "";
inline const std::string& EXTCAT_OUTPUT_DIRECTORY = LensingConfig::SOURCE_CAT;
inline const std::vector<std::string> EXTCAT_FILENAME_TOKENS = {};
inline constexpr bool EXTCAT_RECURSIVE = true;
inline constexpr const char* EXTCAT_DELIMITER = "auto";
inline constexpr const char* EXTCAT_HEADER_MODE = "auto";
inline constexpr const char* EXTCAT_MALFORMED_POLICY = "fail";
inline constexpr const char* EXTCAT_EXISTING_POLICY = "fail";
inline constexpr std::uint64_t EXTCAT_CHUNK_MIB = 64;
// Pass-through external-catalog width only. CatalogLayout uses this value only
// when explicit projection is disabled; downstream phases must consume the
// resolved runtime layout instead of deriving offsets from this constant.
inline constexpr std::size_t EXTCAT_TOTAL_COLUMNS = 18;
inline constexpr bool EXTCAT_USE_EXPLICIT_COLUMNS = false;
inline const std::vector<std::size_t> EXTCAT_INPUT_COLUMNS_ONE_BASED = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
};
inline constexpr bool EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS = false;
// Named field identities in the raw input schema (one-based). RA, Dec, and ZP
// must be positive. A zero magnitude position means that band is unavailable.
// Edit positive identities together with EXTCAT_TOTAL_COLUMNS when adopting a
// different pass-through extcat schema.
inline constexpr std::size_t EXTCAT_RA_COLUMN_ONE_BASED = 5;
inline constexpr std::size_t EXTCAT_DEC_COLUMN_ONE_BASED = 6;
inline constexpr std::size_t EXTCAT_MAG_G_COLUMN_ONE_BASED = 7;
inline constexpr std::size_t EXTCAT_MAG_R_COLUMN_ONE_BASED = 9;
inline constexpr std::size_t EXTCAT_MAG_I_COLUMN_ONE_BASED = 11;
inline constexpr std::size_t EXTCAT_MAG_Z_COLUMN_ONE_BASED = 13;
inline constexpr std::size_t EXTCAT_MAG_Y_COLUMN_ONE_BASED = 15;
inline constexpr std::size_t EXTCAT_ZP_COLUMN_ONE_BASED = 17;

}  // namespace ExtCatConfig

#endif  // EXTCAT_CONFIG_HPP
