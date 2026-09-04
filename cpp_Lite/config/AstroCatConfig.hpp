#ifndef ASTROCAT_CONFIG_HPP
#define ASTROCAT_CONFIG_HPP

#include "Initialize.hpp"
#include "LensingConfig.hpp"
#include "pathconfig.hpp"

namespace AstroCatConfig {

// true means raw files begin with data; false skips exactly one input header.
inline constexpr bool ASTROCAT_ADD_HEADER = Initialize::ASTROCAT_ADD_HEADER;

inline constexpr const char* ASTROCAT_EXISTING_POLICY = "fail";

}  // namespace AstroCatConfig

#endif  // ASTROCAT_CONFIG_HPP
