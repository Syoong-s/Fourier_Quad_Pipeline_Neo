#ifndef CATALOG_ROW_COUNT_HPP
#define CATALOG_ROW_COUNT_HPP

#include "MPIFailure.hpp"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>

namespace CatalogCombiner {
namespace Internal {

// ==========================================
// Function: Count physical catalog data rows after the required header
// Method: Open an independent input stream, skip exactly one header line, and
//         abort the MPI world on open, header-read, or stream I/O failure.
// ==========================================
inline std::size_t countCatalogDataRows(const std::string& filename) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        MPIFailure::abortWorld("count catalog rows", filename);
    }

    std::string line;
    if (!std::getline(input, line)) {
        MPIFailure::abortWorld("read catalog header", filename);
    }

    std::size_t row_count = 0;
    while (std::getline(input, line)) {
        ++row_count;
    }
    if (input.bad()) {
        MPIFailure::abortWorld("count catalog rows", filename);
    }
    return row_count;
}

// ==========================================
// Function: Require matching shear and original catalog data-row counts
// Method: Count each file independently and abort the MPI world with both
//         paths and counts when the physical row totals differ.
// ==========================================
inline void requireMatchingCatalogDataRows(const std::string& shear_filename,
                                           const std::string& orig_filename) {
    const std::size_t shear_rows = countCatalogDataRows(shear_filename);
    const std::size_t orig_rows = countCatalogDataRows(orig_filename);
    if (shear_rows == orig_rows) {
        return;
    }

    std::ostringstream detail;
    detail << "shear_rows=" << shear_rows
           << " orig_rows=" << orig_rows
           << " shear=" << shear_filename
           << " orig=" << orig_filename;
    MPIFailure::abortWorld(
        "combine catalog row-count mismatch", detail.str());
}

}  // namespace Internal
}  // namespace CatalogCombiner

#endif  // CATALOG_ROW_COUNT_HPP
