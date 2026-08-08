#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include "ExtCatConfig.hpp"
#include "InitConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ProcessConfig {

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_EXTCAT = false;
inline constexpr bool RUN_PROCESS_INIT = true;
inline constexpr bool RUN_PROCESS_MAIN = true;
inline constexpr bool RUN_PROCESS_REARR = false;
inline constexpr bool RUN_PROCESS_FD = false;

// ==========================================
// Configuration: Path interface defaults for process_rearr and process_fd
// Method: These I/O path constants are the compile-time defaults seeded into
//         RuntimeOptions below.  CLI options override them at runtime without rebuild.
// ==========================================
inline constexpr const char* EXPO_LIST = "";
inline constexpr const char* REARR_OUTPUT_DIRECTORY = "baked";
inline constexpr const char* REARR_OUTPUT_BASE_DIRECTORY = "";
inline constexpr const char* REARRANGED_EXPO_LIST_FILENAME = "cat_gband_ori.list";
inline constexpr const char* REARRANGED_EXPO_LIST_DIRECTORY = "";
inline constexpr const char* FD_EXPO_LIST = "";
inline constexpr const char* FD_OUTPUT_DIRECTORY = "fdout";
inline constexpr const char* FD_OUTPUT_BASE_DIRECTORY = "";

// ==========================================
// Configuration: Runtime workflow options
// Method: Seed every optional command-line value from the per-process config defaults, then
//         let the unified parser override only explicitly supplied options.
// ==========================================
struct RuntimeOptions {
    bool run_process_extcat = RUN_PROCESS_EXTCAT;
    bool run_process_init = RUN_PROCESS_INIT;
    bool run_process_main = RUN_PROCESS_MAIN;
    bool run_process_rearr = RUN_PROCESS_REARR;
    bool run_process_fd = RUN_PROCESS_FD;
    std::string extcat_input_directory = ExtCatConfig::EXTCAT_INPUT_DIRECTORY;
    std::string extcat_output_directory = ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY;
    std::vector<std::string> extcat_filename_tokens = ExtCatConfig::EXTCAT_FILENAME_TOKENS;
    bool extcat_recursive = ExtCatConfig::EXTCAT_RECURSIVE;
    std::string extcat_delimiter = ExtCatConfig::EXTCAT_DELIMITER;
    std::string extcat_header_mode = ExtCatConfig::EXTCAT_HEADER_MODE;
    std::string extcat_malformed_policy = ExtCatConfig::EXTCAT_MALFORMED_POLICY;
    std::string extcat_existing_policy = ExtCatConfig::EXTCAT_EXISTING_POLICY;
    std::uint64_t extcat_chunk_mib = ExtCatConfig::EXTCAT_CHUNK_MIB;
    bool extcat_use_explicit_columns = ExtCatConfig::EXTCAT_USE_EXPLICIT_COLUMNS;
    std::vector<std::size_t> extcat_input_columns_one_based =
        ExtCatConfig::EXTCAT_INPUT_COLUMNS_ONE_BASED;
    bool extcat_use_explicit_coordinate_columns =
        ExtCatConfig::EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS;
    std::size_t extcat_ra_column_one_based = ExtCatConfig::EXTCAT_RA_COLUMN_ONE_BASED;
    std::size_t extcat_dec_column_one_based = ExtCatConfig::EXTCAT_DEC_COLUMN_ONE_BASED;
    std::size_t extcat_zp_column_one_based = ExtCatConfig::EXTCAT_ZP_COLUMN_ONE_BASED;
    std::string science_root = InitConfig::SCIENCE_ROOT;
    std::string dq_root = InitConfig::DQ_ROOT;
    std::string output_root = InitConfig::OUTPUT_ROOT;
    std::vector<InitConfig::DatasetSpec> datasets = InitConfig::DATASETS;
    std::vector<std::string> contains = InitConfig::CONTAINS;
    std::string existing = InitConfig::EXISTING;
    int f77_max_path = InitConfig::F77_MAX_PATH;
    std::string expo_list = EXPO_LIST;
    std::string rearr_output_directory = REARR_OUTPUT_DIRECTORY;
    std::string rearr_output_base_directory = REARR_OUTPUT_BASE_DIRECTORY;
    std::string rearranged_expo_list_filename = REARRANGED_EXPO_LIST_FILENAME;
    std::string rearranged_expo_list_directory = REARRANGED_EXPO_LIST_DIRECTORY;
    std::string fd_expo_list = FD_EXPO_LIST;
    std::string fd_output_directory = FD_OUTPUT_DIRECTORY;
    std::string fd_output_base_directory = FD_OUTPUT_BASE_DIRECTORY;
    bool external_expo_list_supplied = false;
    bool help_requested = false;
};

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
