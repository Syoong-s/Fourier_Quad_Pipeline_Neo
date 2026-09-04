#ifndef INIT_CONFIG_HPP
#define INIT_CONFIG_HPP

#include "Initialize.hpp"
#include "pathconfig.hpp"

// ==========================================
// InitConfig - Initializer and exposure-list defaults
// Method: Edit this file for one site's usual dataset.
// ==========================================

#include <string>
#include <vector>

namespace InitConfig {

using DatasetSpec = Initialize::DatasetSpec;
inline const auto& DATASETS = Initialize::DATASETS;  // Datasets processed sequentially.
inline const auto& CONTAINS = Initialize::CONTAINS;  // OR-matched archive basename tokens.
inline constexpr const char* EXISTING = "fail";  // Existing-output policy.
inline constexpr int F77_MAX_PATH = 0;  // Generated-path compatibility limit; zero disables it.

}  // namespace InitConfig

#endif  // INIT_CONFIG_HPP
