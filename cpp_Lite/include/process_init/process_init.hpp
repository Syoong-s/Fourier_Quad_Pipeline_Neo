#ifndef PROCESS_INIT_PROCESS_INIT_HPP
#define PROCESS_INIT_PROCESS_INIT_HPP

#include "RuntimeConfig.hpp"

#include <string>

// ==========================================
// Function: Run the integrated MPI pipeline initializer
// Method: Convert unified runtime options into the preserved initializer modules
//         and return the generated absolute exposure-list path on success.
// ==========================================
int process_init(const RuntimeConfig& runtime_config,
                 const InitConfig::DatasetSpec& dataset,
                 std::string& generated_expo_list);

#endif  // PROCESS_INIT_PROCESS_INIT_HPP
