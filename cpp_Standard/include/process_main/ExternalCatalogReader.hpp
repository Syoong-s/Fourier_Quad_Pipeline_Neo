#ifndef EXTERNAL_CATALOG_READER_HPP
#define EXTERNAL_CATALOG_READER_HPP

#include "ProcessConfig.hpp"

#include <cstddef>
#include <string>

namespace ExternalCatalogReader {

// ==========================================
// Structure: Locate the three numerical fields consumed from each generated catalog row
// Method: Store one-based output-tile positions so configuration matches user-facing columns.
// ==========================================
struct ColumnSelection {
    std::size_t ra_column_one_based = 0;
    std::size_t dec_column_one_based = 0;
    std::size_t zp_column_one_based = 0;
};

// ==========================================
// Structure: Hold the external-catalog values used by source extraction
// Method: Parse only configured RA, Dec, and photometric-redshift fields as finite doubles.
// ==========================================
struct Record {
    double ra = 0.0;
    double dec = 0.0;
    double zp = 0.0;
};

// ==========================================
// Function: Resolve raw configured columns to generated-tile columns
// Method: Preserve positions in pass-through mode or locate each raw index in the ordered
//         explicit projection, rejecting missing, zero, or overlapping field selections.
// ==========================================
bool resolveColumnSelection(const ProcessConfig::RuntimeOptions& options,
                            ColumnSelection& selection,
                            std::string& error);

// ==========================================
// Function: Configure the numerical external-catalog reader
// Method: Resolve and store the effective generated-tile positions before MPI stage scheduling.
// ==========================================
bool configure(const ProcessConfig::RuntimeOptions& options, std::string& error);

// ==========================================
// Function: Read one external-catalog row at the configured positions
// Method: Scan whitespace-delimited tokens through the highest requested column and convert
//         only RA, Dec, and ZP, allowing arbitrary text in every unselected field.
// ==========================================
bool parseRecord(const std::string& line, Record& record);

}  // namespace ExternalCatalogReader

#endif  // EXTERNAL_CATALOG_READER_HPP
