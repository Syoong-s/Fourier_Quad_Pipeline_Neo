#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

namespace ProcessConfig {

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_EXTCAT = false;
inline constexpr bool RUN_PROCESS_INIT = true;
inline constexpr bool RUN_PROCESS_MAIN = true;
inline constexpr bool RUN_PROCESS_REARR = true;
inline constexpr bool RUN_PROCESS_FD = true;

// ==========================================
// Configuration: Path interface defaults for process_rearr and process_fd
// Method: These I/O path constants are the compile-time defaults seeded into
//         RuntimeConfig. CLI and INI options override them without rebuilding.
// ==========================================
inline constexpr const char* EXPO_LIST = "";
inline constexpr const char* REARR_OUTPUT_DIRECTORY = "baked";
inline constexpr const char* REARR_OUTPUT_BASE_DIRECTORY = "";
inline constexpr const char* REARRANGED_EXPO_LIST_FILENAME = "cat_gband_ori.list";
inline constexpr const char* REARRANGED_EXPO_LIST_DIRECTORY = "";
inline constexpr const char* FD_EXPO_LIST = "";
inline constexpr const char* FD_OUTPUT_DIRECTORY = "fdout";
inline constexpr const char* FD_OUTPUT_BASE_DIRECTORY = "";

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
