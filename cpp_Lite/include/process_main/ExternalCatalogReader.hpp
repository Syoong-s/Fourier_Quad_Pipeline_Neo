#ifndef EXTERNAL_CATALOG_READER_HPP
#define EXTERNAL_CATALOG_READER_HPP

#include "CatalogLayout.hpp"

#include <cstddef>
#include <string>

namespace ExternalCatalogReader {

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
// Function: Configure the numerical external-catalog reader
// Method: Consume and store effective RA, Dec, and ZP positions from the shared layout.
// ==========================================
bool configure(const PipelineCatalog::CatalogLayout& layout,
               std::string& error);

// ==========================================
// Function: Read one external-catalog row at the configured positions
// Method: Scan whitespace-delimited tokens through the highest requested column and convert
//         only RA, Dec, and ZP, allowing arbitrary text in every unselected field.
// ==========================================
bool parseRecord(const std::string& line, Record& record);

}  // namespace ExternalCatalogReader

#endif  // EXTERNAL_CATALOG_READER_HPP
