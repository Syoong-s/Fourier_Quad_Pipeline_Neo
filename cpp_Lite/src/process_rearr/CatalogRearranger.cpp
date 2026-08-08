#include "process_rearr/CatalogRearranger.hpp"

#include "ProcessRearrConfig.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace ProcessRearr {
namespace {

// ==========================================
// Structure: Hold one populated 0.1-degree sky tile
// Method: Retain its center, row weight, and original flat-grid position for
//         deterministic recursive sorting and lookup reconstruction.
// ==========================================
struct ActiveTile {
    double dec = 0.0;
    double ra = 0.0;
    std::uint64_t weight = 0;
    std::size_t flat_index = 0;
};

// ==========================================
// Function: Map one configured raw coordinate to its external-output position
// Method: Preserve the one-based raw index in pass-through mode or locate it
//         in the ordered explicit projection exactly as process_main does.
// ==========================================
bool resolveProjectedCoordinate(
    std::size_t raw_column_one_based,
    const std::string& field_name,
    const ProcessConfig::RuntimeOptions& options,
    std::size_t& output_column_one_based,
    std::string& error) {
    if (raw_column_one_based == 0) {
        error = "process_rearr " + field_name
                + " column must be a positive one-based index";
        return false;
    }
    if (!options.extcat_use_explicit_columns) {
        output_column_one_based = raw_column_one_based;
        return true;
    }

    const auto match = std::find(options.extcat_input_columns_one_based.begin(),
                                 options.extcat_input_columns_one_based.end(),
                                 raw_column_one_based);
    if (match == options.extcat_input_columns_one_based.end()) {
        error = "process_rearr explicit projection omits the configured "
                + field_name + " raw column "
                + std::to_string(raw_column_one_based);
        return false;
    }
    output_column_one_based =
        static_cast<std::size_t>(std::distance(
            options.extcat_input_columns_one_based.begin(), match))
        + 1;
    return true;
}

// ==========================================
// Function: Convert one catalog token to a finite double
// Method: Use strtod with full-token and range checks so NaN, infinity,
//         overflow, suffixes, and empty values are rejected uniformly.
// ==========================================
bool parseFiniteDouble(const std::string& token, double& value) {
    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(token.c_str(), &end);
    if (token.empty() || end != token.c_str() + token.size()
        || errno == ERANGE || !std::isfinite(parsed)) {
        return false;
    }
    value = parsed;
    return true;
}

// ==========================================
// Function: Compare two active tiles along one k-d split dimension
// Method: Match the F77 primary/secondary coordinate ordering and use the
//         full-sky flat index only as a deterministic final tie breaker.
// ==========================================
bool tileLess(const ActiveTile& left, const ActiveTile& right, int dimension) {
    if (dimension == 0) {
        if (left.dec != right.dec) {
            return left.dec < right.dec;
        }
        if (left.ra != right.ra) {
            return left.ra < right.ra;
        }
    } else {
        if (left.ra != right.ra) {
            return left.ra < right.ra;
        }
        if (left.dec != right.dec) {
            return left.dec < right.dec;
        }
    }
    return left.flat_index < right.flat_index;
}

// ==========================================
// Function: Recursively assign active tiles to weighted k-d partitions
// Method: Sort the current slice, split k into floor(k/2) and the remainder,
//         then choose the closest cumulative-weight boundary before recursing.
// ==========================================
void partitionRecursive(const std::vector<ActiveTile>& active_tiles,
                        std::vector<std::size_t>& order,
                        std::size_t left,
                        std::size_t right,
                        std::size_t partition_count,
                        int dimension,
                        int first_partition,
                        std::vector<int>& active_partitions) {
    if (left > right) {
        return;
    }
    if (partition_count <= 1 || left >= right) {
        for (std::size_t index = left; index <= right; ++index) {
            active_partitions[order[index]] = first_partition;
        }
        return;
    }

    std::sort(order.begin() + static_cast<std::ptrdiff_t>(left),
              order.begin() + static_cast<std::ptrdiff_t>(right + 1),
              [&](std::size_t first, std::size_t second) {
                  return tileLess(active_tiles[first], active_tiles[second], dimension);
              });

    long double total_weight = 0.0L;
    for (std::size_t index = left; index <= right; ++index) {
        total_weight += static_cast<long double>(active_tiles[order[index]].weight);
    }

    const std::size_t left_partition_count = partition_count / 2;
    const std::size_t right_partition_count =
        partition_count - left_partition_count;
    const long double target_weight =
        total_weight * static_cast<long double>(left_partition_count)
        / static_cast<long double>(partition_count);

    long double running_weight = 0.0L;
    std::size_t middle = left;
    for (std::size_t index = left; index < right; ++index) {
        const long double previous_weight = running_weight;
        running_weight +=
            static_cast<long double>(active_tiles[order[index]].weight);
        if (running_weight >= target_weight) {
            if (index > left
                && std::fabs(previous_weight - target_weight)
                       < std::fabs(running_weight - target_weight)) {
                middle = index - 1;
            } else {
                middle = index;
            }
            break;
        }
        middle = index;
    }

    partitionRecursive(active_tiles, order, left, middle,
                       left_partition_count, 1 - dimension, first_partition,
                       active_partitions);
    partitionRecursive(active_tiles, order, middle + 1, right,
                       right_partition_count, 1 - dimension,
                       first_partition + static_cast<int>(left_partition_count),
                       active_partitions);
}

}  // namespace

// ==========================================
// Function: Resolve and validate the _all.cat layout
// Method: Reuse process_main's raw-to-projected coordinate mapping, then
//         validate both coordinates against the effective external width.
// ==========================================
bool resolveCatalogLayout(const ProcessConfig::RuntimeOptions& options,
                          CatalogLayout& layout,
                          std::string& error) {
    std::size_t ra_column_one_based = 0;
    std::size_t dec_column_one_based = 0;
    if (!resolveProjectedCoordinate(options.extcat_ra_column_one_based, "RA",
                                    options, ra_column_one_based, error)
        || !resolveProjectedCoordinate(options.extcat_dec_column_one_based, "Dec",
                                       options, dec_column_one_based, error)) {
        return false;
    }
    if (ra_column_one_based == dec_column_one_based) {
        error = "process_rearr RA and Dec columns must be distinct";
        return false;
    }

    const std::size_t external_columns =
        ProcessRearrConfig::externalCatalogColumns(options);
    if (external_columns == 0) {
        error = "process_rearr external catalog column count must be positive";
        return false;
    }
    if (ra_column_one_based > external_columns
        || dec_column_one_based > external_columns) {
        error = "process_rearr RA/Dec columns must lie within the effective external "
                "catalog width "
                + std::to_string(external_columns);
        return false;
    }
    if (external_columns
        > std::numeric_limits<std::size_t>::max()
              - ProcessRearrConfig::CCD_COLUMN_COUNT
              - ProcessRearrConfig::ichi2) {
        error = "process_rearr _all.cat column count overflows size_t";
        return false;
    }

    layout.external_columns = external_columns;
    layout.all_columns = ProcessRearrConfig::allCatalogColumns(options);
    layout.ra_column = ra_column_one_based - 1;
    layout.dec_column = dec_column_one_based - 1;
    error.clear();
    return true;
}

// ==========================================
// Function: Parse one complete numeric _all.cat row
// Method: Require exactly the configured number of finite numeric tokens and
//         retain double precision for MPI redistribution and output.
// ==========================================
bool parseCatalogRow(const std::string& line,
                     std::size_t column_count,
                     std::vector<double>& values,
                     std::string& error) {
    if (column_count == 0) {
        error = "catalog row width must be positive";
        return false;
    }

    std::istringstream input(line);
    std::vector<double> parsed_values;
    parsed_values.reserve(column_count);
    for (std::size_t column = 0; column < column_count; ++column) {
        std::string token;
        if (!(input >> token)) {
            error = "catalog row has fewer than " + std::to_string(column_count)
                    + " columns";
            return false;
        }
        double value = 0.0;
        if (!parseFiniteDouble(token, value)) {
            error = "catalog row column " + std::to_string(column + 1)
                    + " is not a finite number";
            return false;
        }
        parsed_values.push_back(value);
    }

    std::string extra;
    if (input >> extra) {
        error = "catalog row has more than " + std::to_string(column_count)
                + " columns";
        return false;
    }

    values = std::move(parsed_values);
    error.clear();
    return true;
}

// ==========================================
// Function: Map one celestial coordinate to the legacy full-sky grid
// Method: Apply the F77 single-step RA wrap, pole clamping, floor-based
//         0.1-degree bins, and final integer bounds checks.
// ==========================================
bool skyTileIndex(double ra,
                  double dec,
                  std::size_t& tile_index,
                  std::string& error) {
    if (!std::isfinite(ra) || !std::isfinite(dec)) {
        error = "RA and Dec must be finite";
        return false;
    }

    if (ra >= 360.0) {
        ra -= 360.0;
    }
    if (ra < 0.0) {
        ra += 360.0;
    }
    if (dec >= 90.0) {
        dec = 89.999;
    }
    if (dec < -90.0) {
        dec = -90.0;
    }

    const double bins_per_degree = 1.0 / ProcessRearrConfig::SKY_GRID_DEGREES;
    int dec_bin = static_cast<int>(std::floor(dec * bins_per_degree));
    int ra_bin = static_cast<int>(std::floor(ra * bins_per_degree));
    dec_bin = std::clamp(dec_bin, -900, 899);
    ra_bin = std::clamp(ra_bin, 0, ProcessRearrConfig::RA_BIN_COUNT - 1);
    tile_index = static_cast<std::size_t>(dec_bin + 900)
                     * ProcessRearrConfig::RA_BIN_COUNT
                 + static_cast<std::size_t>(ra_bin);
    error.clear();
    return true;
}

// ==========================================
// Function: Choose the requested number of target subcatalogs
// Method: Compute ceil(total rows / target rows), with one empty partition for
//         an otherwise valid header-only input set.
// ==========================================
std::size_t partitionCount(std::uint64_t total_rows) {
    static_assert(ProcessRearrConfig::TARGET_SUBCAT_ROWS > 0,
                  "TARGET_SUBCAT_ROWS must be positive");
    if (total_rows == 0) {
        return 1;
    }
    return static_cast<std::size_t>(
        1 + (total_rows - 1) / ProcessRearrConfig::TARGET_SUBCAT_ROWS);
}

// ==========================================
// Function: Assign populated sky tiles to deterministic k-d partitions
// Method: Alternate Dec/RA sorting and split cumulative tile weight nearest
//         the proportional target used by the F77 recursive bisection.
// ==========================================
bool buildTilePartitions(const std::vector<std::uint64_t>& tile_counts,
                         std::size_t partition_count,
                         std::vector<int>& tile_partitions,
                         std::string& error) {
    if (tile_counts.size() != ProcessRearrConfig::SKY_TILE_COUNT) {
        error = "tile count array must contain exactly "
                + std::to_string(ProcessRearrConfig::SKY_TILE_COUNT)
                + " full-sky bins";
        return false;
    }
    if (partition_count == 0
        || partition_count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        error = "partition count must fit a positive int";
        return false;
    }

    tile_partitions.assign(ProcessRearrConfig::SKY_TILE_COUNT, 0);
    std::vector<ActiveTile> active_tiles;
    for (std::size_t flat_index = 0; flat_index < tile_counts.size(); ++flat_index) {
        if (tile_counts[flat_index] == 0) {
            continue;
        }
        const std::size_t dec_index =
            flat_index / static_cast<std::size_t>(ProcessRearrConfig::RA_BIN_COUNT);
        const std::size_t ra_index =
            flat_index % static_cast<std::size_t>(ProcessRearrConfig::RA_BIN_COUNT);
        ActiveTile tile;
        tile.dec = (static_cast<double>(dec_index) - 900.0)
                       * ProcessRearrConfig::SKY_GRID_DEGREES
                   + ProcessRearrConfig::SKY_GRID_DEGREES / 2.0;
        tile.ra = static_cast<double>(ra_index)
                      * ProcessRearrConfig::SKY_GRID_DEGREES
                  + ProcessRearrConfig::SKY_GRID_DEGREES / 2.0;
        tile.weight = tile_counts[flat_index];
        tile.flat_index = flat_index;
        active_tiles.push_back(tile);
    }

    if (active_tiles.empty()) {
        error.clear();
        return true;
    }

    std::vector<std::size_t> order(active_tiles.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::vector<int> active_partitions(active_tiles.size(), 0);
    partitionRecursive(active_tiles, order, 0, active_tiles.size() - 1,
                       partition_count, 0, 1, active_partitions);

    for (std::size_t index = 0; index < active_tiles.size(); ++index) {
        if (active_partitions[index] <= 0) {
            error = "weighted k-d partition left an active sky tile unassigned";
            return false;
        }
        tile_partitions[active_tiles[index].flat_index] = active_partitions[index];
    }
    error.clear();
    return true;
}

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
                    std::vector<std::size_t>& row_indices) {
    if (column_count == 0 || ra_column >= column_count || dec_column >= column_count) {
        return;
    }
    const std::size_t row_count = row_values.size() / column_count;
    if (source_exposures.size() < row_count || source_rows.size() < row_count) {
        return;
    }

    std::sort(row_indices.begin(), row_indices.end(),
              [&](std::size_t first, std::size_t second) {
                  const double first_dec = row_values[first * column_count + dec_column];
                  const double second_dec = row_values[second * column_count + dec_column];
                  if (first_dec != second_dec) {
                      return first_dec < second_dec;
                  }
                  const double first_ra = row_values[first * column_count + ra_column];
                  const double second_ra = row_values[second * column_count + ra_column];
                  if (first_ra != second_ra) {
                      return first_ra < second_ra;
                  }
                  if (source_exposures[first] != source_exposures[second]) {
                      return source_exposures[first] < source_exposures[second];
                  }
                  if (source_rows[first] != source_rows[second]) {
                      return source_rows[first] < source_rows[second];
                  }
                  return first < second;
              });
}

}  // namespace ProcessRearr
