#ifndef RUNTIME_CONFIG_HPP
#define RUNTIME_CONFIG_HPP
#include "InitConfig.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Value objects do not supply scientific defaults. Initialize through
// makeDefaultRuntimeConfig() before applying optional INI overrides.
struct ProcessRuntimeConfig {
    bool run_process_astrocat{};
    bool run_process_extcat{};
    bool run_process_init{};
    bool run_process_main{};
    bool run_process_rearr{};
    bool run_process_fd{};
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
    bool add_header{};
    std::string existing_policy;
};

struct ExtCatRuntimeConfig {
    std::string input_directory;
    std::string output_directory;
    std::vector<std::string> filename_tokens;
    bool recursive{};
    std::string delimiter;
    std::string header_mode;
    std::string malformed_policy;
    std::string existing_policy;
    std::uint64_t chunk_mib{};
    std::size_t total_columns{};
    bool use_explicit_columns{};
    std::vector<std::size_t> input_columns_one_based;
    bool use_explicit_coordinate_columns{};
    std::size_t ra_column_one_based{};
    std::size_t dec_column_one_based{};
    std::size_t mag_g_column_one_based{};
    std::size_t mag_r_column_one_based{};
    std::size_t mag_i_column_one_based{};
    std::size_t mag_z_column_one_based{};
    std::size_t mag_y_column_one_based{};
    std::size_t zp_column_one_based{};
};

struct InitRuntimeConfig {
    std::string science_root;
    std::string dq_root;
    std::string output_root;
    std::vector<InitConfig::DatasetSpec> datasets;
    std::vector<std::string> contains;
    std::string existing;
    int f77_max_path{};
};

struct LensingRuntimeConfig {
    std::int64_t process_stage{};
    std::string astrometry_cat;
    int ccd_split{};
    int gal_smooth{};
    int star_smooth{};
    double pixel_size{};
    int nmax_chip{};
    int chipnx{};
    int chipny{};
    int astrometry_trivial{};
    int astrometry_cat_type{};
    int include_flat{};
    int include_mask{};
    std::string flat_path;
    std::string psf_path;
    int ext_cat{};
    int ext_psf{};
    int deblending{};
    int psf_type{};
    int psf_ms{};
};

struct RuntimeConfig {
    ProcessRuntimeConfig process;
    AstroCatRuntimeConfig astrocat;
    ExtCatRuntimeConfig extcat;
    InitRuntimeConfig init;
    LensingRuntimeConfig lensing;
};

RuntimeConfig makeDefaultRuntimeConfig();
// Validate only --config PATH / --config=PATH. No implicit INI lookup.
bool findRuntimeConfigPath(int argc, char* argv[], std::string& path, std::string& error);
bool applyRuntimeConfigText(const std::string& text, const std::string& source, RuntimeConfig& config, std::string& error);
bool validateRuntimeConfig(const RuntimeConfig& config, std::string& error);
namespace RuntimeConfigStore {
bool initialize(RuntimeConfig config, std::string& error);
bool isInitialized();
const RuntimeConfig& get();
}
#endif
