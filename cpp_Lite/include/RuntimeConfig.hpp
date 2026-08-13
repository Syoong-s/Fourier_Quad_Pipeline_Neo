#ifndef RUNTIME_CONFIG_HPP
#define RUNTIME_CONFIG_HPP

#include "InitConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// ==========================================
// Structure: ProcessRuntimeConfig
// Method: Hold user-facing workflow switches and downstream path controls while
//         leaving parser-only state outside the persisted process section.
// ==========================================
struct ProcessRuntimeConfig {
    bool run_process_extcat = false;
    bool run_process_init = true;
    bool run_process_main = true;
    bool run_process_rearr = false;
    bool run_process_fd = false;
    std::string expo_list;
    std::string rearr_output_directory;
    std::string rearr_output_base_directory;
    std::string rearranged_expo_list_filename;
    std::string rearranged_expo_list_directory;
    std::string fd_expo_list;
    std::string fd_output_directory;
    std::string fd_output_base_directory;
};

// ==========================================
// Structure: ExtCatRuntimeConfig
// Method: Represent the complete external-catalog parser, projection, schema,
//         and publication configuration with one authoritative tile directory.
// ==========================================
struct ExtCatRuntimeConfig {
    std::string input_directory;
    std::string output_directory;
    std::vector<std::string> filename_tokens;
    bool recursive = true;
    std::string delimiter;
    std::string header_mode;
    std::string malformed_policy;
    std::string existing_policy;
    std::uint64_t chunk_mib = 64;
    std::size_t total_columns = 18;
    bool use_explicit_columns = false;
    std::vector<std::size_t> input_columns_one_based;
    bool use_explicit_coordinate_columns = false;
    std::size_t ra_column_one_based = 5;
    std::size_t dec_column_one_based = 6;
    std::size_t mag_g_column_one_based = 7;
    std::size_t mag_r_column_one_based = 9;
    std::size_t mag_i_column_one_based = 11;
    std::size_t mag_z_column_one_based = 13;
    std::size_t mag_y_column_one_based = 15;
    std::size_t zp_column_one_based = 17;
};

// ==========================================
// Structure: InitRuntimeConfig
// Method: Keep archive roots, dataset pairs, filename filters, publication
//         policy, and the optional legacy path guard together.
// ==========================================
struct InitRuntimeConfig {
    std::string science_root;
    std::string dq_root;
    std::string output_root;
    std::vector<InitConfig::DatasetSpec> datasets;
    std::vector<std::string> contains;
    std::string existing;
    int f77_max_path = 0;
};

// ==========================================
// Structure: LensingRuntimeConfig
// Method: Store only parameters that still have live consumers in Lite; the
//         physically deleted Standard branches remain absent from this API.
// ==========================================
struct LensingRuntimeConfig {
    std::int64_t process_stage = 1;
    std::string astrometry_cat;
    int ccd_split = 2;
    int gal_smooth = 0;
    int star_smooth = 2;
    double pixel_size = 0.2628;
    int nmax_chip = 62;
    int chipnx = 2046;
    int chipny = 4094;
};

// ==========================================
// Structure: RuntimeConfig
// Method: Aggregate every run-selectable domain plus non-config-file parser
//         state into one value that becomes immutable before pipeline dispatch.
// ==========================================
struct RuntimeConfig {
    ProcessRuntimeConfig process;
    ExtCatRuntimeConfig extcat;
    InitRuntimeConfig init;
    LensingRuntimeConfig lensing;
    bool external_expo_list_supplied = false;
    bool help_requested = false;
};

// ==========================================
// Function: makeDefaultRuntimeConfig
// Method: Copy the four compiled configuration-header domains into one mutable
//         startup value before file and CLI overrides are applied.
// ==========================================
RuntimeConfig makeDefaultRuntimeConfig();

// ==========================================
// Function: findRuntimeConfigPath
// Method: Perform the phase-A CLI scan for --config and --help without applying
//         any other option or touching the filesystem.
// ==========================================
bool findRuntimeConfigPath(int argc, char* argv[], std::string& config_path,
                           bool& help_requested, std::string& error);

// ==========================================
// Function: applyRuntimeConfigText
// Method: Parse minimal INI text transactionally and commit recognized typed
//         overrides only when every section/key/value is valid.
// ==========================================
bool applyRuntimeConfigText(const std::string& text,
                            const std::string& source_name,
                            RuntimeConfig& config,
                            std::string& error);

// ==========================================
// Function: parseRuntimeCommandLine
// Method: Apply all existing CLI controls after file values, accepting both
//         --name value and --name=value while skipping the phase-A --config.
// ==========================================
bool parseRuntimeCommandLine(int argc, char* argv[], RuntimeConfig& config,
                             std::string& error);

// ==========================================
// Function: validateRuntimeConfig
// Method: Enforce executable branch domains, stage dependencies, positive
//         geometry, and workflow input invariants without instrument hard limits.
// ==========================================
bool validateRuntimeConfig(const RuntimeConfig& config, std::string& error);

namespace RuntimeConfigStore {

// ==========================================
// Function: initialize
// Method: Copy the finalized per-rank configuration exactly once and reject
//         accidental attempts to change it after pipeline startup.
// ==========================================
bool initialize(RuntimeConfig config, std::string& error);

// ==========================================
// Function: isInitialized
// Method: Let compatibility entry points determine whether compiled defaults
//         still need to be installed before accessing deep-stage configuration.
// ==========================================
bool isInitialized();

// ==========================================
// Function: get
// Method: Return the immutable finalized configuration or throw when startup
//         failed to initialize the store.
// ==========================================
const RuntimeConfig& get();

}  // namespace RuntimeConfigStore

#endif  // RUNTIME_CONFIG_HPP
