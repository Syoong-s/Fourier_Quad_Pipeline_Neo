#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include "pathconfig.hpp"

namespace ProcessConfig {

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_ASTROCAT = false;  // Run Gaia-catalog tiling by default.
inline constexpr bool RUN_PROCESS_EXTCAT = false;  // Run external-catalog tiling by default.
inline constexpr bool RUN_PROCESS_INIT = true;  // Run archive initialization by default.
inline constexpr bool RUN_PROCESS_MAIN = true;  // Run the nine-stage pipeline by default.
inline constexpr bool RUN_PROCESS_REARR = false;  // Run catalog rearrangement by default.
inline constexpr bool RUN_PROCESS_FD = false;  // Run the field-distortion test by default.

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
