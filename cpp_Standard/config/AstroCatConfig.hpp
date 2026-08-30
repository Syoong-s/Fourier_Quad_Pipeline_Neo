#ifndef ASTROCAT_CONFIG_HPP
#define ASTROCAT_CONFIG_HPP

#include "LensingConfig.hpp"

#include <string>

namespace AstroCatConfig {

inline constexpr const char* ASTROCAT_INPUT_DIRECTORY = "";

inline const std::string ASTROCAT_OUTPUT_DIRECTORY =
    LensingConfig::ASTROMETRY_CAT;

// true means raw files begin with data; false skips exactly one input header.
inline constexpr bool ASTROCAT_ADD_HEADER = true;

inline constexpr const char* ASTROCAT_EXISTING_POLICY = "fail";

}  // namespace AstroCatConfig

#endif  // ASTROCAT_CONFIG_HPP
