#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include "Initialize.hpp"
#include "pathconfig.hpp"

namespace ProcessConfig {

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_ASTROCAT = Initialize::RUN_PROCESS_ASTROCAT;
inline constexpr bool RUN_PROCESS_EXTCAT = Initialize::RUN_PROCESS_EXTCAT;
inline constexpr bool RUN_PROCESS_INIT = Initialize::RUN_PROCESS_INIT;
inline constexpr bool RUN_PROCESS_MAIN = Initialize::RUN_PROCESS_MAIN;
inline constexpr bool RUN_PROCESS_REARR = Initialize::RUN_PROCESS_REARR;
inline constexpr bool RUN_PROCESS_FD = Initialize::RUN_PROCESS_FD;

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
