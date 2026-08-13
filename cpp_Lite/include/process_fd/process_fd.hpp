#ifndef PROCESS_FD_PROCESS_FD_HPP
#define PROCESS_FD_PROCESS_FD_HPP

#include "CatalogLayout.hpp"

#include <string>

// ==========================================
// Function: Run the FD (field-distortion) shear test on one exposure list
// Method: Consume the startup catalog layout while reading per-exposure rows,
//         applying star removal, recovering shear, and writing FD_test_comb.dat.
// ==========================================
int process_fd(const std::string& exposure_list,
               const RuntimeConfig& runtime_config,
               const std::string& dataset_root,
               const PipelineCatalog::CatalogLayout& layout);

#endif  // PROCESS_FD_PROCESS_FD_HPP
