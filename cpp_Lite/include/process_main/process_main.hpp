#ifndef PROCESS_MAIN_PROCESS_MAIN_HPP
#define PROCESS_MAIN_PROCESS_MAIN_HPP

#include "ProcessConfig.hpp"

#include <string>

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline on one exposure list
// Method: Load and broadcast the list, then execute the configured MPI stages
//         without owning MPI initialization or finalization.
// ==========================================
int process_main(const std::string& exposure_list);

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline with unified runtime options
// Method: Resolve external-catalog projection and RA/Dec/ZP columns before loading and
//         broadcasting the exposure list, preserving the one-argument compatibility entry point.
// ==========================================
int process_main(const std::string& exposure_list,
                 const ProcessConfig::RuntimeOptions& options);

#endif  // PROCESS_MAIN_PROCESS_MAIN_HPP
