#ifndef INIT_CONFIG_HPP
#define INIT_CONFIG_HPP

#include "pathconfig.hpp"

// ==========================================
// InitConfig - Initializer and exposure-list defaults
// Method: Edit this file for one site's usual dataset.
// ==========================================

#include <string>
#include <vector>

namespace InitConfig {

struct DatasetSpec {
    std::string target;  // Dataset output directory name.
    std::string prefix;  // Archive basename prefix.
};

// {"Target1", "Prefix1"}, {"Target2", "Prefix2"} ...
inline const std::vector<DatasetSpec> DATASETS = {
    {"gband", "c4d_"}
};  // Datasets processed sequentially.
// "Contains1", "Contains2" ...
inline const std::vector<std::string> CONTAINS = {"v1"};  // OR-matched archive basename tokens.
inline constexpr const char* EXISTING = "fail";  // Existing-output policy.
inline constexpr int F77_MAX_PATH = 0;  // Generated-path compatibility limit; zero disables it.

}  // namespace InitConfig

#endif  // INIT_CONFIG_HPP
