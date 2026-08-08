#ifndef PROCESS_REARR_PROCESS_REARR_HPP
#define PROCESS_REARR_PROCESS_REARR_HPP

#include "ProcessConfig.hpp"

#include <mpi.h>

#include <string>

// ==========================================
// Function: Rearrange exposure _all.cat files into spatial subcatalogs
// Method: Read catalogs across MPI ranks, build a global weighted k-d sky
//         partition, redistribute complete rows, and write sorted outputs.
// ==========================================
int process_rearr(const std::string& exposure_list,
                  const ProcessConfig::RuntimeOptions& options,
                  MPI_Comm communicator = MPI_COMM_WORLD);

#endif  // PROCESS_REARR_PROCESS_REARR_HPP
