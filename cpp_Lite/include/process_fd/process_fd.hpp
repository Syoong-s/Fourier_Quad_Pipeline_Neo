#ifndef PROCESS_FD_PROCESS_FD_HPP
#define PROCESS_FD_PROCESS_FD_HPP

#include "ProcessConfig.hpp"

#include <string>

// ==========================================
// Function: Run the FD (field-distortion) shear test on one exposure list
// Method: Read per-exposure shear catalogs, compute star-bar point-source
//         removal, perform spatial binning and shear-recovery (PDF or
//         jackknife), and write FD_test_comb.dat.  Acts as the fifth
//         pipeline stage alongside process_extcat / process_init /
//         process_main / process_rearr.
// ==========================================
int process_fd(const std::string& exposure_list,
               const ProcessConfig::RuntimeOptions& options,
               const std::string& dataset_root);

#endif  // PROCESS_FD_PROCESS_FD_HPP
