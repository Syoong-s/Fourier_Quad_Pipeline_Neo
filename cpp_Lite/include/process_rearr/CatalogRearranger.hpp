#ifndef PROCESS_REARR_CATALOG_REARRANGER_HPP
#define PROCESS_REARR_CATALOG_REARRANGER_HPP

#include "ProcessConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ProcessRearr {

// ==========================================
// Structure: Describe the effective numeric _all.cat layout
// Method: Store the external and complete widths plus zero-based coordinate
//         positions after applying process_main's explicit projection rules.
// ==========================================
struct CatalogLayout {
    std::size_t external_columns = 0;
    std::size_t all_columns = 0;
    std::size_t ra_column = 0;
    std::size_t dec_column = 0;
};

// ==========================================
// Function: Resolve and validate the _all.cat layout
// Method: Reuse process_main's raw-to-projected coordinate mapping, then
//         validate both coordinates against the effective external width.
// ==========================================
bool resolveCatalogLayout(const ProcessConfig::RuntimeOptions& options,
                          CatalogLayout& layout,
                          std::string& error);

// ==========================================
// Function: Parse one complete numeric _all.cat row
// Method: Require exactly the configured number of finite numeric tokens and
//         retain double precision for MPI redistribution and output.
// ==========================================
bool parseCatalogRow(const std::string& line,
                     std::size_t column_count,
                     std::vector<double>& values,
                     std::string& error);

// ==========================================
// Function: Map one celestial coordinate to the legacy full-sky grid
// Method: Apply the F77 single-step RA wrap, pole clamping, floor-based
//         0.1-degree bins, and final integer bounds checks.
// ==========================================
bool skyTileIndex(double ra,
                  double dec,
                  std::size_t& tile_index,
                  std::string& error);

// ==========================================
// Function: Choose the requested number of target subcatalogs
// Method: Compute ceil(total rows / target rows), with one empty partition for
//         an otherwise valid header-only input set.
// ==========================================
std::size_t partitionCount(std::uint64_t total_rows);

// ==========================================
// Function: Assign populated sky tiles to deterministic k-d partitions
// Method: Alternate Dec/RA sorting and split cumulative tile weight nearest
//         the proportional target used by the F77 recursive bisection.
// ==========================================
bool buildTilePartitions(const std::vector<std::uint64_t>& tile_counts,
                         std::size_t partition_count,
                         std::vector<int>& tile_partitions,
                         std::string& error);

// ==========================================
// Function: Sort selected flat-buffer rows by celestial position
// Method: Order by Dec then RA, using exposure and source-row keys only to make
//         equal-coordinate output deterministic across MPI process counts.
// ==========================================
void sortRowIndices(const std::vector<double>& row_values,
                    std::size_t column_count,
                    std::size_t ra_column,
                    std::size_t dec_column,
                    const std::vector<std::uint64_t>& source_exposures,
                    const std::vector<std::uint64_t>& source_rows,
                    std::vector<std::size_t>& row_indices);

}  // namespace ProcessRearr

#endif  // PROCESS_REARR_CATALOG_REARRANGER_HPP
