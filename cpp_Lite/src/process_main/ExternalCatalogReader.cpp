#include "ExternalCatalogReader.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace ExternalCatalogReader {
namespace {

// ==========================================
// Structure: Locate the numerical fields consumed from generated catalog rows
// Method: Store one-based positions derived only from the shared runtime layout.
// ==========================================
struct ColumnSelection {
    std::size_t ra_column_one_based = 0;
    std::size_t dec_column_one_based = 0;
    std::size_t zp_column_one_based = 0;
};

ColumnSelection active_columns;

// ==========================================
// Function: Convert one selected catalog token to a finite double
// Method: Use strtod with full-token and range validation so malformed values are rejected.
// ==========================================
bool parseFiniteDouble(const std::string& token, double& value) {
    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(token.c_str(), &end);
    if (end != token.c_str() + token.size() || errno == ERANGE
        || !std::isfinite(parsed)) {
        return false;
    }
    value = parsed;
    return true;
}

}  // namespace

// ==========================================
// Function: Configure the numerical external-catalog reader
// Method: Consume and store effective RA, Dec, and ZP positions from the shared layout.
// ==========================================
bool configure(const PipelineCatalog::CatalogLayout& layout,
               std::string& error) {
    if (layout.external_columns == 0
        || layout.external.ra >= layout.external_columns
        || layout.external.dec >= layout.external_columns
        || layout.external.zp >= layout.external_columns) {
        error = "external-catalog reader positions lie outside the shared layout";
        return false;
    }
    if (layout.external.ra == layout.external.dec
        || layout.external.ra == layout.external.zp
        || layout.external.dec == layout.external.zp) {
        error = "external-catalog reader RA, Dec, and ZP positions must be distinct";
        return false;
    }
    active_columns.ra_column_one_based = layout.external.ra + 1;
    active_columns.dec_column_one_based = layout.external.dec + 1;
    active_columns.zp_column_one_based = layout.external.zp + 1;
    error.clear();
    return true;
}

// ==========================================
// Function: Read one external-catalog row at the configured positions
// Method: Scan whitespace-delimited tokens through the highest requested column and convert
//         only RA, Dec, and ZP, allowing arbitrary text in every unselected field.
// ==========================================
bool parseRecord(const std::string& line, Record& record) {
    const std::size_t final_column = std::max(
        active_columns.ra_column_one_based,
        std::max(active_columns.dec_column_one_based,
                 active_columns.zp_column_one_based));
    if (final_column == 0) {
        return false;
    }

    std::istringstream input(line);
    Record parsed;
    for (std::size_t column = 1; column <= final_column; ++column) {
        std::string token;
        if (!(input >> token)) {
            return false;
        }
        if (column == active_columns.ra_column_one_based
            && !parseFiniteDouble(token, parsed.ra)) {
            return false;
        }
        if (column == active_columns.dec_column_one_based
            && !parseFiniteDouble(token, parsed.dec)) {
            return false;
        }
        if (column == active_columns.zp_column_one_based
            && !parseFiniteDouble(token, parsed.zp)) {
            return false;
        }
    }
    record = parsed;
    return true;
}

}  // namespace ExternalCatalogReader
