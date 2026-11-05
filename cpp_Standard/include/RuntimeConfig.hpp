#ifndef RUNTIME_CONFIG_HPP
#define RUNTIME_CONFIG_HPP
#include "InitConfig.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct ProcessRuntimeConfig {
    bool run_process_astrocat = false;
    bool run_process_extcat = false;
    bool run_process_init = true;
    bool run_process_main = true;
    bool run_process_rearr = true;
    bool run_process_fd = true;
    std::string expo_list;
    std::string rearr_output_directory;
    std::string rearr_output_base_directory;
    std::string rearranged_expo_list_filename;
    std::string rearranged_expo_list_directory;
    std::string fd_expo_list;
    std::string fd_output_directory;
    std::string fd_output_base_directory;
};

struct AstroCatRuntimeConfig {
    std::string input_directory;
    std::string output_directory;
    bool add_header = true;
    std::string existing_policy;
};

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

struct InitRuntimeConfig {
    std::string science_root;
    std::string dq_root;
    std::string output_root;
    std::vector<InitConfig::DatasetSpec> datasets;
    std::vector<std::string> contains;
    std::string existing;
    int f77_max_path = 0;
};

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
    int astrometry_trivial = 0;
    int astrometry_cat_type = 1;
    int include_flat = 0;
    int include_mask = 2;
    std::string flat_path;
    std::string psf_path;
    int ext_cat = 1;
    int ext_psf = 0;
    int deblending = 1;
    int psf_type = 1;
    int psf_ms = 0;
};

struct RuntimeConfig {
    ProcessRuntimeConfig process;
    AstroCatRuntimeConfig astrocat;
    ExtCatRuntimeConfig extcat;
    InitRuntimeConfig init;
    LensingRuntimeConfig lensing;
    bool external_expo_list_supplied = false;
    bool help_requested = false;  // Internal compatibility field; no --help option.
};

RuntimeConfig makeDefaultRuntimeConfig();
// Validate only --config PATH / --config=PATH. No implicit INI lookup.
bool findRuntimeConfigPath(int argc, char* argv[], std::string& path, bool& help_requested, std::string& error);
bool applyRuntimeConfigText(const std::string& text, const std::string& source, RuntimeConfig& config, std::string& error);
// Internal compatibility API: validates arguments without applying overrides.
bool parseRuntimeCommandLine(int argc, char* argv[], RuntimeConfig& config, std::string& error);
bool validateRuntimeConfig(const RuntimeConfig& config, std::string& error);
namespace RuntimeConfigStore {
bool initialize(RuntimeConfig config, std::string& error);
bool isInitialized();
const RuntimeConfig& get();
}
#endif
