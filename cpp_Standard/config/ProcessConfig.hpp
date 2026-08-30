#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

namespace ProcessConfig {

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_ASTROCAT = false;  // Run Gaia-catalog tiling by default.
inline constexpr bool RUN_PROCESS_EXTCAT = false;  // Run external-catalog tiling by default.
inline constexpr bool RUN_PROCESS_INIT = true;  // Run archive initialization by default.
inline constexpr bool RUN_PROCESS_MAIN = true;  // Run the nine-stage pipeline by default.
inline constexpr bool RUN_PROCESS_REARR = true;  // Run catalog rearrangement by default.
inline constexpr bool RUN_PROCESS_FD = true;  // Run the field-distortion test by default.

// ==========================================
// Configuration: Path interface defaults for process_rearr and process_fd
// Method: These I/O path constants are the compile-time defaults seeded into
//         RuntimeConfig. CLI and INI options override them without rebuilding.
// ==========================================
inline constexpr const char* EXPO_LIST = "";  // Default top-level exposure-list path.
inline constexpr const char* REARR_OUTPUT_DIRECTORY = "baked";  // Rearranged catalog directory.
inline constexpr const char* REARR_OUTPUT_BASE_DIRECTORY = "";  // Empty uses the dataset root.
inline constexpr const char* REARRANGED_EXPO_LIST_FILENAME = "cat_gband_ori.list";  // FD input list name.
inline constexpr const char* REARRANGED_EXPO_LIST_DIRECTORY = "";  // Empty uses the input-list directory.
inline constexpr const char* FD_EXPO_LIST = "";  // Optional explicit FD exposure-list path.
inline constexpr const char* FD_OUTPUT_DIRECTORY = "fdout";  // FD result directory name.
inline constexpr const char* FD_OUTPUT_BASE_DIRECTORY = "";  // Empty uses the dataset root.

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
