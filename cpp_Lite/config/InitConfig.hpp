#ifndef INIT_CONFIG_HPP
#define INIT_CONFIG_HPP

// ==========================================
// InitConfig - Initializer and exposure-list defaults
// Method: Edit this file for one site's usual dataset.
// ==========================================

#include <string>
#include <vector>

namespace InitConfig {

struct DatasetSpec {
    std::string target;
    std::string prefix;
};

inline constexpr const char* SCIENCE_ROOT = "/lustre/home/acct-phyzj/share/DES/g";
inline constexpr const char* DQ_ROOT = "/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask";
inline constexpr const char* OUTPUT_ROOT = "/lustre/home/acct-phyzj/share/DES/g_band_v1";
// {"Target1", "Prefix1"}, {"Target2", "Prefix2"} ...
inline const std::vector<DatasetSpec> DATASETS = {
    {"gband", "c4d_"}
};
// "Contains1", "Contains2" ...
inline const std::vector<std::string> CONTAINS = {"v1"};
inline constexpr const char* EXISTING = "fail";
inline constexpr int F77_MAX_PATH = 150;

}  // namespace InitConfig

#endif  // INIT_CONFIG_HPP
