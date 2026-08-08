#include "ExternalCatalogReader.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace ExternalCatalogReader {
namespace {

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

// ==========================================
// Function: Map one raw input column to its generated-tile position
// Method: Keep its index for pass-through output or use the first matching projection entry.
// ==========================================
bool resolveProjectedColumn(std::size_t raw_column_one_based,
                            const std::string& field_name,
                            const ProcessConfig::RuntimeOptions& options,
                            std::size_t& output_column_one_based,
                            std::string& error) {
    if (raw_column_one_based == 0) {
        error = "external-catalog " + field_name + " column must be a positive one-based index";
        return false;
    }
    if (!options.extcat_use_explicit_columns) {
        output_column_one_based = raw_column_one_based;
        return true;
    }

    const std::vector<std::size_t>& projection =
        options.extcat_input_columns_one_based;
    const auto match = std::find(projection.begin(), projection.end(), raw_column_one_based);
    if (match == projection.end()) {
        error = "external-catalog explicit projection omits the configured "
                + field_name + " raw column " + std::to_string(raw_column_one_based);
        return false;
    }
    output_column_one_based =
        static_cast<std::size_t>(std::distance(projection.begin(), match)) + 1;
    return true;
}

}  // namespace

// ==========================================
// Function: Resolve raw configured columns to generated-tile columns
// Method: Preserve positions in pass-through mode or locate each raw index in the ordered
//         explicit projection, rejecting missing, zero, or overlapping field selections.
// ==========================================
bool resolveColumnSelection(const ProcessConfig::RuntimeOptions& options,
                            ColumnSelection& selection,
                            std::string& error) {
    ColumnSelection resolved;
    if (!resolveProjectedColumn(options.extcat_ra_column_one_based, "RA", options,
                                resolved.ra_column_one_based, error)
        || !resolveProjectedColumn(options.extcat_dec_column_one_based, "Dec", options,
                                   resolved.dec_column_one_based, error)
        || !resolveProjectedColumn(options.extcat_zp_column_one_based, "ZP", options,
                                   resolved.zp_column_one_based, error)) {
        return false;
    }
    if (resolved.ra_column_one_based == resolved.dec_column_one_based
        || resolved.ra_column_one_based == resolved.zp_column_one_based
        || resolved.dec_column_one_based == resolved.zp_column_one_based) {
        error = "external-catalog RA, Dec, and ZP columns must be distinct";
        return false;
    }
    selection = resolved;
    error.clear();
    return true;
}

// ==========================================
// Function: Configure the numerical external-catalog reader
// Method: Resolve and store the effective generated-tile positions before MPI stage scheduling.
// ==========================================
bool configure(const ProcessConfig::RuntimeOptions& options, std::string& error) {
    ColumnSelection resolved;
    if (!resolveColumnSelection(options, resolved, error)) {
        return false;
    }
    active_columns = resolved;
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
