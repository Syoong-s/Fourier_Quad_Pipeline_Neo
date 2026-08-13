#ifndef PROCESS_REARR_PROCESS_REARR_HPP
#define PROCESS_REARR_PROCESS_REARR_HPP

#include "CatalogLayout.hpp"

#include <mpi.h>

#include <string>

// ==========================================
// Function: Rearrange exposure _all.cat files into spatial subcatalogs
// Method: Consume the startup layout while reading catalogs across MPI ranks,
//         partitioning the sky, redistributing complete rows, and writing output.
// ==========================================
int process_rearr(const std::string& exposure_list,
                  const RuntimeConfig& runtime_config,
                  const PipelineCatalog::RearrCatalogSchema& schema,
                  MPI_Comm communicator = MPI_COMM_WORLD);

#endif  // PROCESS_REARR_PROCESS_REARR_HPP
