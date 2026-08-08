#include "process_rearr/process_rearr.hpp"

#include "process_rearr/CatalogRearranger.hpp"
#include "ProcessRearrConfig.hpp"

#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

// ==========================================
// Structure: Hold root-prepared input and output paths
// Method: Resolve every exposure catalog once, retain one shared header, and
//         broadcast only normalized strings needed by distributed readers.
// ==========================================
struct PreparedInputs {
    std::vector<std::string> catalog_paths;
    std::string dataset_root;
    std::string output_directory;
    std::string header;
};

// ==========================================
// Structure: Store locally read rows in contiguous column-major-independent form
// Method: Keep row-major values plus deterministic source keys while counting
//         skipped files and malformed records for the collective report.
// ==========================================
struct LocalRows {
    std::vector<double> values;
    std::vector<std::uint64_t> source_exposures;
    std::vector<std::uint64_t> source_rows;
    std::uint64_t missing_catalogs = 0;
    std::uint64_t malformed_rows = 0;
};

// ==========================================
// Structure: Describe row and flattened-element MPI Alltoallv layouts
// Method: Build row counts first, then validate int-safe value and metadata
//         multipliers before entering the variable-count collectives.
// ==========================================
struct TransferPlan {
    std::vector<int> send_rows;
    std::vector<int> receive_rows;
    std::vector<int> send_row_displacements;
    std::vector<int> receive_row_displacements;
    std::vector<int> send_values;
    std::vector<int> receive_values;
    std::vector<int> send_value_displacements;
    std::vector<int> receive_value_displacements;
    std::vector<int> send_metadata;
    std::vector<int> receive_metadata;
    std::vector<int> send_metadata_displacements;
    std::vector<int> receive_metadata_displacements;
    int total_send_rows = 0;
    int total_receive_rows = 0;
};

// ==========================================
// Structure: Store rows owned by this rank after MPI redistribution
// Method: Keep the numeric table flat while unpacking partition and stable
//         source keys into parallel arrays used by sorting and output.
// ==========================================
struct ReceivedRows {
    std::vector<double> values;
    std::vector<int> partitions;
    std::vector<std::uint64_t> source_exposures;
    std::vector<std::uint64_t> source_rows;
};

// ==========================================
// Function: Remove leading and trailing ASCII whitespace
// Method: Use unsigned-char ctype calls so non-ASCII path bytes are not passed
//         to std::isspace with undefined signed-char values.
// ==========================================
std::string trimWhitespace(const std::string& value) {
    std::size_t first = 0;
    while (first < value.size()
           && std::isspace(static_cast<unsigned char>(value[first])) != 0) {
        ++first;
    }
    std::size_t last = value.size();
    while (last > first
           && std::isspace(static_cast<unsigned char>(value[last - 1])) != 0) {
        --last;
    }
    return value.substr(first, last - first);
}

// ==========================================
// Function: Remove one matching pair of path quotes
// Method: Accept the quoted path convention already supported by process_main
//         while leaving unquoted and mismatched values untouched.
// ==========================================
std::string stripMatchingQuotes(const std::string& value) {
    if (value.size() >= 2
        && ((value.front() == '"' && value.back() == '"')
            || (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

// ==========================================
// Function: Load the top-level exposure list
// Method: Match process_main's path-and-chip-count records while validating a
//         nonempty, fully parsed input before any catalog discovery begins.
// ==========================================
bool loadExposureList(const std::string& exposure_list,
                      std::vector<std::string>& exposure_paths,
                      std::string& error) {
    std::ifstream input(exposure_list);
    if (!input.is_open()) {
        error = "process_rearr cannot open exposure list: " + exposure_list;
        return false;
    }

    std::string path;
    int chip_count = 0;
    while (input >> path >> chip_count) {
        exposure_paths.push_back(stripMatchingQuotes(path));
    }
    if (!input.eof()) {
        error = "process_rearr exposure list contains an invalid record: "
                + exposure_list;
        return false;
    }
    if (exposure_paths.empty()) {
        error = "process_rearr exposure list contains no exposures: "
                + exposure_list;
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Resolve one exposure _all.cat path via getDir(image, 3)
// Method: Open the per-exposure list, read the first non-empty image line,
//         compute the great-grandparent directory (getDir level 3) as the dataset
//         root for this exposure, and construct the _all.cat path. Called
//         once per exposure, matching process_main's per-exposure derivation.
// ==========================================
bool resolveCatalogPathFromImage(const std::string& exposure_list_path,
                                 fs::path& dataset_root,
                                 fs::path& catalog_path,
                                 std::string& error) {
    std::ifstream input(exposure_list_path);
    if (!input.is_open()) {
        error = "process_rearr cannot open per-exposure list: "
                + exposure_list_path;
        return false;
    }
    std::string line;
    while (std::getline(input, line)) {
        line = stripMatchingQuotes(trimWhitespace(line));
        if (line.empty()) {
            continue;
        }
        const fs::path image_path(line);
        const fs::path parent = image_path.parent_path();
        const fs::path grandparent = parent.parent_path();
        const fs::path great_grandparent = grandparent.parent_path();
        if (parent.empty() || grandparent.empty() || great_grandparent.empty()) {
            error = "process_rearr image path has fewer than three parent "
                    "levels: " + line;
            return false;
        }
        dataset_root = fs::absolute(great_grandparent).lexically_normal();
        const std::string basename = image_path.filename().string();
        const std::size_t underscore = basename.find_last_of('_');
        if (underscore == std::string::npos || underscore == 0) {
            error = "process_rearr image basename lacks an exposure suffix: "
                    + basename;
            return false;
        }
        const std::string prefix = basename.substr(0, underscore);
        catalog_path = dataset_root / "result" / (prefix + "_all.cat");
        error.clear();
        return true;
    }
    error = "process_rearr per-exposure list is empty: " + exposure_list_path;
    return false;
}

// ==========================================
// Function: Prepare every catalog path, output path, and shared header on rank zero
// Method: Resolve a single dataset root, select the first readable catalog
//         header, and preserve missing-catalog skip behavior for distributed reads.
// ==========================================
bool prepareInputs(const std::string& exposure_list,
                   const ProcessConfig::RuntimeOptions& options,
                   PreparedInputs& prepared,
                   std::string& error) {
    std::vector<std::string> exposure_paths;
    if (!loadExposureList(exposure_list, exposure_paths, error)) {
        return false;
    }

    fs::path output_base_root;
    for (const std::string& exposure_path : exposure_paths) {
        fs::path dataset_root;
        fs::path catalog_path;
        if (!resolveCatalogPathFromImage(exposure_path, dataset_root,
                                         catalog_path, error)) {
            return false;
        }
        if (output_base_root.empty()) {
            output_base_root = dataset_root;
        }
        prepared.catalog_paths.push_back(catalog_path.string());
    }

    const std::string base_dir_str(options.rearr_output_base_directory);
    const fs::path base_dir =
        base_dir_str.empty() ? output_base_root : fs::path(base_dir_str);
    fs::path configured_output(options.rearr_output_directory);
    if (configured_output.empty()) {
        configured_output = base_dir;
    } else if (configured_output.is_relative()) {
        configured_output = base_dir / configured_output;
    }
    prepared.dataset_root = output_base_root.string();
    prepared.output_directory = fs::absolute(configured_output).lexically_normal().string();

    for (const std::string& catalog_path : prepared.catalog_paths) {
        std::ifstream input(catalog_path);
        if (!input.is_open()) {
            if (ProcessRearrConfig::SKIP_MISSING_CATALOGS) {
                continue;
            }
            error = "process_rearr cannot open catalog: " + catalog_path;
            return false;
        }
        if (!std::getline(input, prepared.header)) {
            error = "process_rearr catalog has no header: " + catalog_path;
            return false;
        }
        prepared.header = trimWhitespace(prepared.header);
        break;
    }
    if (prepared.header.empty()) {
        error = "process_rearr found no readable catalog header";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Broadcast one mutable string
// Method: Send an int-safe byte length followed by the exact character payload.
// ==========================================
bool broadcastString(std::string& value,
                     int root,
                     int rank,
                     MPI_Comm communicator,
                     std::string& error) {
    int length = 0;
    if (rank == root) {
        if (value.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            length = -1;
        } else {
            length = static_cast<int>(value.size());
        }
    }
    if (MPI_Bcast(&length, 1, MPI_INT, root, communicator) != MPI_SUCCESS) {
        error = "process_rearr failed to broadcast a string length";
        return false;
    }
    if (length < 0) {
        error = "process_rearr broadcast string exceeds MPI int length";
        return false;
    }
    if (rank != root) {
        value.resize(static_cast<std::size_t>(length));
    }
    if (length > 0
        && MPI_Bcast(value.data(), length, MPI_CHAR, root, communicator)
               != MPI_SUCCESS) {
        error = "process_rearr failed to broadcast a string payload";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Broadcast a root-prepared string vector
// Method: Send the vector length followed by length-prefixed strings so every
//         rank receives identical catalog ordering and stable exposure keys.
// ==========================================
bool broadcastStringVector(std::vector<std::string>& values,
                           int root,
                           int rank,
                           MPI_Comm communicator,
                           std::string& error) {
    int count = 0;
    if (rank == root) {
        if (values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            count = -1;
        } else {
            count = static_cast<int>(values.size());
        }
    }
    if (MPI_Bcast(&count, 1, MPI_INT, root, communicator) != MPI_SUCCESS) {
        error = "process_rearr failed to broadcast catalog path count";
        return false;
    }
    if (count < 0) {
        error = "process_rearr catalog path count exceeds MPI int range";
        return false;
    }
    if (rank != root) {
        values.resize(static_cast<std::size_t>(count));
    }
    bool success = true;
    for (std::string& value : values) {
        std::string item_error;
        const bool item_success =
            broadcastString(value, root, rank, communicator, item_error);
        if (!item_success && success) {
            error = item_error;
        }
        success = success && item_success;
    }
    if (!success) {
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Report and combine one rank-local validation result
// Method: Reduce with MPI_MIN and serialize only failing rank diagnostics to
//         keep collective error exits readable and deadlock-free.
// ==========================================
bool collectiveSuccess(bool local_success,
                       const std::string& local_error,
                       const std::string& stage,
                       int rank,
                       int world_size,
                       MPI_Comm communicator) {
    const int local_value = local_success ? 1 : 0;
    int global_value = 0;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_MIN, communicator);
    if (global_value != 0) {
        return true;
    }
    for (int reporting_rank = 0; reporting_rank < world_size; ++reporting_rank) {
        if (rank == reporting_rank && !local_success) {
            std::cerr << "process_rearr [" << stage << "] rank " << rank << ": "
                      << (local_error.empty() ? "unspecified error" : local_error)
                      << std::endl;
        }
        MPI_Barrier(communicator);
    }
    return false;
}

// ==========================================
// Function: Read this rank's static exposure subset
// Method: Use exposure-index striding, require the shared header and exact row
//         width, and skip only policies explicitly enabled in the phase config.
// ==========================================
bool readLocalCatalogs(const PreparedInputs& prepared,
                       const ProcessRearr::CatalogLayout& layout,
                       int rank,
                       int world_size,
                       LocalRows& rows,
                       std::string& error) {
    for (std::size_t exposure = static_cast<std::size_t>(rank);
         exposure < prepared.catalog_paths.size();
         exposure += static_cast<std::size_t>(world_size)) {
        const std::string& catalog_path = prepared.catalog_paths[exposure];
        std::ifstream input(catalog_path);
        if (!input.is_open()) {
            ++rows.missing_catalogs;
            if (ProcessRearrConfig::SKIP_MISSING_CATALOGS) {
                continue;
            }
            error = "cannot open catalog " + catalog_path;
            return false;
        }

        std::string header;
        if (!std::getline(input, header)) {
            error = "catalog has no header: " + catalog_path;
            return false;
        }
        if (trimWhitespace(header) != prepared.header) {
            error = "catalog header differs from the shared schema: " + catalog_path;
            return false;
        }

        std::string line;
        std::uint64_t line_number = 1;
        while (std::getline(input, line)) {
            ++line_number;
            if (trimWhitespace(line).empty()) {
                continue;
            }
            std::vector<double> parsed;
            std::string row_error;
            if (!ProcessRearr::parseCatalogRow(line, layout.all_columns,
                                               parsed, row_error)) {
                ++rows.malformed_rows;
                if (ProcessRearrConfig::SKIP_MALFORMED_ROWS) {
                    continue;
                }
                error = catalog_path + ":" + std::to_string(line_number)
                        + ": " + row_error;
                return false;
            }
            rows.values.insert(rows.values.end(), parsed.begin(), parsed.end());
            rows.source_exposures.push_back(static_cast<std::uint64_t>(exposure));
            rows.source_rows.push_back(line_number);
        }
        if (input.bad()) {
            error = "I/O error while reading catalog: " + catalog_path;
            return false;
        }
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Bin all local rows on the full-sky grid
// Method: Reuse the validated layout and preserve each row's tile index for
//         later partition lookup without recomputing coordinate boundaries.
// ==========================================
bool binLocalRows(const LocalRows& rows,
                  const ProcessRearr::CatalogLayout& layout,
                  std::vector<std::uint64_t>& tile_counts,
                  std::vector<std::size_t>& row_tiles,
                  std::string& error) {
    tile_counts.assign(ProcessRearrConfig::SKY_TILE_COUNT, 0);
    const std::size_t row_count = rows.source_rows.size();
    row_tiles.resize(row_count);
    if (rows.values.size() != row_count * layout.all_columns) {
        error = "local row buffer has inconsistent dimensions";
        return false;
    }

    for (std::size_t row = 0; row < row_count; ++row) {
        const double ra = rows.values[row * layout.all_columns + layout.ra_column];
        const double dec = rows.values[row * layout.all_columns + layout.dec_column];
        std::size_t tile = 0;
        if (!ProcessRearr::skyTileIndex(ra, dec, tile, error)) {
            error = "local row " + std::to_string(row) + ": " + error;
            return false;
        }
        if (tile_counts[tile] == std::numeric_limits<std::uint64_t>::max()) {
            error = "one sky tile exceeds uint64 row capacity";
            return false;
        }
        ++tile_counts[tile];
        row_tiles[row] = tile;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Sum global tile weights without unsigned wraparound
// Method: Check each addition explicitly before deriving the target partition count.
// ==========================================
bool sumTileCounts(const std::vector<std::uint64_t>& counts,
                   std::uint64_t& total,
                   std::string& error) {
    total = 0;
    for (const std::uint64_t count : counts) {
        if (count > std::numeric_limits<std::uint64_t>::max() - total) {
            error = "global catalog row count exceeds uint64 capacity";
            return false;
        }
        total += count;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Build row counts for destination ranks
// Method: Map partition p to rank (p-1) mod world_size and reject any invalid
//         assignment before the first all-to-all collective.
// ==========================================
bool buildSendRowCounts(const std::vector<int>& row_partitions,
                        int world_size,
                        TransferPlan& plan,
                        std::string& error) {
    plan.send_rows.assign(static_cast<std::size_t>(world_size), 0);
    for (const int partition : row_partitions) {
        if (partition <= 0) {
            error = "a local row has no positive partition assignment";
            return false;
        }
        const int destination = (partition - 1) % world_size;
        if (plan.send_rows[static_cast<std::size_t>(destination)]
            == std::numeric_limits<int>::max()) {
            error = "one MPI destination exceeds int row-count capacity";
            return false;
        }
        ++plan.send_rows[static_cast<std::size_t>(destination)];
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Build int-safe displacements for one MPI count vector
// Method: Accumulate in long long and reject any total beyond Alltoallv's
//         standard int count/displacement interface.
// ==========================================
bool buildDisplacements(const std::vector<int>& counts,
                        std::vector<int>& displacements,
                        int& total,
                        std::string& error) {
    displacements.resize(counts.size());
    long long running = 0;
    for (std::size_t index = 0; index < counts.size(); ++index) {
        if (counts[index] < 0
            || running > static_cast<long long>(std::numeric_limits<int>::max())) {
            error = "MPI row counts or displacements exceed int range";
            return false;
        }
        displacements[index] = static_cast<int>(running);
        running += counts[index];
    }
    if (running > static_cast<long long>(std::numeric_limits<int>::max())) {
        error = "MPI total row count exceeds int range";
        return false;
    }
    total = static_cast<int>(running);
    error.clear();
    return true;
}

// ==========================================
// Function: Multiply MPI row counts and displacements by a record width
// Method: Validate every product against int limits before generating the
//         flattened value or metadata Alltoallv layout.
// ==========================================
bool scaleTransferLayout(const std::vector<int>& row_counts,
                         const std::vector<int>& row_displacements,
                         std::size_t width,
                         std::vector<int>& element_counts,
                         std::vector<int>& element_displacements,
                         std::string& error) {
    if (width == 0
        || width > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        error = "MPI record width must fit a positive int";
        return false;
    }
    element_counts.resize(row_counts.size());
    element_displacements.resize(row_displacements.size());
    for (std::size_t index = 0; index < row_counts.size(); ++index) {
        const long long count = static_cast<long long>(row_counts[index])
                                * static_cast<long long>(width);
        const long long displacement =
            static_cast<long long>(row_displacements[index])
            * static_cast<long long>(width);
        if (count > std::numeric_limits<int>::max()
            || displacement > std::numeric_limits<int>::max()) {
            error = "flattened MPI transfer exceeds int count/displacement range";
            return false;
        }
        element_counts[index] = static_cast<int>(count);
        element_displacements[index] = static_cast<int>(displacement);
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Complete receive counts and flattened transfer layouts
// Method: Exchange row counts, derive displacements, then build value-width
//         and three-field metadata layouts with collective-safe validation.
// ==========================================
bool completeTransferPlan(std::size_t column_count,
                          TransferPlan& plan,
                          MPI_Comm communicator,
                          std::string& error) {
    const std::size_t world_size = plan.send_rows.size();
    plan.receive_rows.assign(world_size, 0);
    if (MPI_Alltoall(plan.send_rows.data(), 1, MPI_INT,
                     plan.receive_rows.data(), 1, MPI_INT,
                     communicator)
        != MPI_SUCCESS) {
        error = "MPI_Alltoall failed for process_rearr row counts";
        return false;
    }
    if (!buildDisplacements(plan.send_rows, plan.send_row_displacements,
                            plan.total_send_rows, error)
        || !buildDisplacements(plan.receive_rows, plan.receive_row_displacements,
                               plan.total_receive_rows, error)
        || !scaleTransferLayout(plan.send_rows, plan.send_row_displacements,
                                column_count, plan.send_values,
                                plan.send_value_displacements, error)
        || !scaleTransferLayout(plan.receive_rows, plan.receive_row_displacements,
                                column_count, plan.receive_values,
                                plan.receive_value_displacements, error)
        || !scaleTransferLayout(plan.send_rows, plan.send_row_displacements, 3,
                                plan.send_metadata,
                                plan.send_metadata_displacements, error)
        || !scaleTransferLayout(plan.receive_rows, plan.receive_row_displacements, 3,
                                plan.receive_metadata,
                                plan.receive_metadata_displacements, error)) {
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Pack and exchange complete rows plus partition/source metadata
// Method: Group row-major values by destination, send doubles and three uint64
//         metadata fields with matching Alltoallv layouts, then unpack keys.
// ==========================================
bool exchangeRows(const LocalRows& local_rows,
                  const std::vector<int>& row_partitions,
                  const ProcessRearr::CatalogLayout& layout,
                  const TransferPlan& plan,
                  int world_size,
                  MPI_Comm communicator,
                  ReceivedRows& received,
                  std::string& error) {
    const std::size_t local_count = local_rows.source_rows.size();
    if (local_rows.values.size() != local_count * layout.all_columns
        || local_rows.source_exposures.size() != local_count
        || row_partitions.size() != local_count
        || plan.total_send_rows != static_cast<int>(local_count)) {
        error = "local buffers do not match the MPI transfer plan";
        return false;
    }

    std::vector<double> send_values(local_rows.values.size());
    std::vector<std::uint64_t> send_metadata(local_count * 3);
    std::vector<int> cursor = plan.send_row_displacements;
    for (std::size_t row = 0; row < local_count; ++row) {
        const int destination = (row_partitions[row] - 1) % world_size;
        const std::size_t packed_row = static_cast<std::size_t>(
            cursor[static_cast<std::size_t>(destination)]++);
        std::copy_n(local_rows.values.begin()
                        + static_cast<std::ptrdiff_t>(row * layout.all_columns),
                    layout.all_columns,
                    send_values.begin()
                        + static_cast<std::ptrdiff_t>(packed_row * layout.all_columns));
        send_metadata[packed_row * 3] =
            static_cast<std::uint64_t>(row_partitions[row]);
        send_metadata[packed_row * 3 + 1] = local_rows.source_exposures[row];
        send_metadata[packed_row * 3 + 2] = local_rows.source_rows[row];
    }

    received.values.resize(static_cast<std::size_t>(plan.total_receive_rows)
                           * layout.all_columns);
    std::vector<std::uint64_t> receive_metadata(
        static_cast<std::size_t>(plan.total_receive_rows) * 3);
    const int value_status =
        MPI_Alltoallv(send_values.data(), plan.send_values.data(),
                      plan.send_value_displacements.data(), MPI_DOUBLE,
                      received.values.data(), plan.receive_values.data(),
                      plan.receive_value_displacements.data(), MPI_DOUBLE,
                      communicator);
    const int metadata_status =
        MPI_Alltoallv(send_metadata.data(), plan.send_metadata.data(),
                      plan.send_metadata_displacements.data(), MPI_UINT64_T,
                      receive_metadata.data(), plan.receive_metadata.data(),
                      plan.receive_metadata_displacements.data(), MPI_UINT64_T,
                      communicator);
    if (value_status != MPI_SUCCESS || metadata_status != MPI_SUCCESS) {
        error = "MPI_Alltoallv failed while redistributing process_rearr rows";
        return false;
    }

    const std::size_t receive_count =
        static_cast<std::size_t>(plan.total_receive_rows);
    received.partitions.resize(receive_count);
    received.source_exposures.resize(receive_count);
    received.source_rows.resize(receive_count);
    for (std::size_t row = 0; row < receive_count; ++row) {
        if (receive_metadata[row * 3] == 0
            || receive_metadata[row * 3]
                   > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            error = "received partition identifier is outside positive int range";
            return false;
        }
        received.partitions[row] =
            static_cast<int>(receive_metadata[row * 3]);
        received.source_exposures[row] = receive_metadata[row * 3 + 1];
        received.source_rows[row] = receive_metadata[row * 3 + 2];
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Format one six-digit-compatible subcatalog filename
// Method: Preserve the legacy zero-padded partition identifier and allow
//         larger identifiers to expand rather than truncate.
// ==========================================
std::string subcatalogFilename(std::size_t partition) {
    std::ostringstream name;
    name << ProcessRearrConfig::SUBCAT_PREFIX
         << std::setw(ProcessRearrConfig::SUBCAT_ID_WIDTH)
         << std::setfill('0') << partition
         << ProcessRearrConfig::SUBCAT_EXTENSION;
    return name.str();
}

// ==========================================
// Function: Write this rank's owned, coordinate-sorted partitions
// Method: Group received rows by partition, sort Dec then RA, emit the shared
//         header and complete rows, and fill global-summary reduction arrays.
// ==========================================
bool writeLocalPartitions(const ReceivedRows& received,
                          const ProcessRearr::CatalogLayout& layout,
                          const PreparedInputs& prepared,
                          std::size_t partition_count,
                          int rank,
                          int world_size,
                          std::vector<std::uint64_t>& counts,
                          std::vector<double>& dec_min,
                          std::vector<double>& dec_max,
                          std::vector<double>& ra_min,
                          std::vector<double>& ra_max,
                          std::string& error) {
    const std::size_t row_count = received.partitions.size();
    if (received.values.size() != row_count * layout.all_columns
        || received.source_exposures.size() != row_count
        || received.source_rows.size() != row_count) {
        error = "received row buffers have inconsistent dimensions";
        return false;
    }

    std::vector<std::vector<std::size_t>> rows_by_partition(partition_count + 1);
    for (std::size_t row = 0; row < row_count; ++row) {
        const int partition = received.partitions[row];
        if (partition <= 0
            || static_cast<std::size_t>(partition) > partition_count
            || (partition - 1) % world_size != rank) {
            error = "received row belongs to an invalid or foreign partition";
            return false;
        }
        rows_by_partition[static_cast<std::size_t>(partition)].push_back(row);
    }

    counts.assign(partition_count, 0);
    dec_min.assign(partition_count, std::numeric_limits<double>::infinity());
    dec_max.assign(partition_count, -std::numeric_limits<double>::infinity());
    ra_min.assign(partition_count, std::numeric_limits<double>::infinity());
    ra_max.assign(partition_count, -std::numeric_limits<double>::infinity());

    for (std::size_t partition = static_cast<std::size_t>(rank) + 1;
         partition <= partition_count;
         partition += static_cast<std::size_t>(world_size)) {
        std::vector<std::size_t>& indices = rows_by_partition[partition];
        if (indices.empty()) {
            continue;
        }
        ProcessRearr::sortRowIndices(received.values, layout.all_columns,
                                     layout.ra_column, layout.dec_column,
                                     received.source_exposures,
                                     received.source_rows, indices);

        const fs::path output_path =
            fs::path(prepared.output_directory) / subcatalogFilename(partition);
        std::ofstream output(output_path, std::ios::trunc);
        if (!output.is_open()) {
            error = "cannot write subcatalog: " + output_path.string();
            return false;
        }
        output << prepared.header << '\n';
        output << std::setprecision(ProcessRearrConfig::OUTPUT_PRECISION);

        const std::size_t summary_index = partition - 1;
        counts[summary_index] = static_cast<std::uint64_t>(indices.size());
        for (const std::size_t row : indices) {
            for (std::size_t column = 0; column < layout.all_columns; ++column) {
                if (column != 0) {
                    output << ' ';
                }
                output << received.values[row * layout.all_columns + column];
            }
            output << '\n';

            const double dec =
                received.values[row * layout.all_columns + layout.dec_column];
            const double ra =
                received.values[row * layout.all_columns + layout.ra_column];
            dec_min[summary_index] = std::min(dec_min[summary_index], dec);
            dec_max[summary_index] = std::max(dec_max[summary_index], dec);
            ra_min[summary_index] = std::min(ra_min[summary_index], ra);
            ra_max[summary_index] = std::max(ra_max[summary_index], ra);
        }
        output.close();
        if (!output) {
            error = "I/O failure while writing subcatalog: " + output_path.string();
            return false;
        }
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Reduce per-partition summaries and write the ordered root report
// Method: Use SUM for unique-owner counts and MIN/MAX for coordinate bounds,
//         then emit only nonempty partitions in ascending identifier order.
// ==========================================
bool reduceAndWriteSummary(const PreparedInputs& prepared,
                           const std::vector<std::uint64_t>& local_counts,
                           const std::vector<double>& local_dec_min,
                           const std::vector<double>& local_dec_max,
                           const std::vector<double>& local_ra_min,
                           const std::vector<double>& local_ra_max,
                           int rank,
                           MPI_Comm communicator,
                           std::string& error) {
    if (local_counts.size()
        > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        error = "summary partition count exceeds MPI int range";
        return false;
    }
    const int count = static_cast<int>(local_counts.size());
    std::vector<std::uint64_t> global_counts(rank == 0 ? local_counts.size() : 0);
    std::vector<double> global_dec_min(rank == 0 ? local_counts.size() : 0);
    std::vector<double> global_dec_max(rank == 0 ? local_counts.size() : 0);
    std::vector<double> global_ra_min(rank == 0 ? local_counts.size() : 0);
    std::vector<double> global_ra_max(rank == 0 ? local_counts.size() : 0);

    const int count_status =
        MPI_Reduce(local_counts.data(), global_counts.data(), count, MPI_UINT64_T,
                   MPI_SUM, 0, communicator);
    const int dec_min_status =
        MPI_Reduce(local_dec_min.data(), global_dec_min.data(), count, MPI_DOUBLE,
                   MPI_MIN, 0, communicator);
    const int dec_max_status =
        MPI_Reduce(local_dec_max.data(), global_dec_max.data(), count, MPI_DOUBLE,
                   MPI_MAX, 0, communicator);
    const int ra_min_status =
        MPI_Reduce(local_ra_min.data(), global_ra_min.data(), count, MPI_DOUBLE,
                   MPI_MIN, 0, communicator);
    const int ra_max_status =
        MPI_Reduce(local_ra_max.data(), global_ra_max.data(), count, MPI_DOUBLE,
                   MPI_MAX, 0, communicator);
    if (count_status != MPI_SUCCESS || dec_min_status != MPI_SUCCESS
        || dec_max_status != MPI_SUCCESS || ra_min_status != MPI_SUCCESS
        || ra_max_status != MPI_SUCCESS) {
        error = "MPI_Reduce failed for process_rearr summary";
        return false;
    }

    if (rank == 0) {
        const fs::path summary_path =
            fs::path(prepared.output_directory)
            / std::string(ProcessRearrConfig::SUMMARY_FILENAME);
        std::ofstream output(summary_path, std::ios::trunc);
        if (!output.is_open()) {
            error = "cannot write process_rearr summary: " + summary_path.string();
            return false;
        }
        output << "part_id  count  dec_min  dec_max  ra_min   ra_max\n";
        output << std::fixed << std::setprecision(ProcessRearrConfig::SUMMARY_PRECISION);
        for (std::size_t index = 0; index < global_counts.size(); ++index) {
            if (global_counts[index] == 0) {
                continue;
            }
            output << std::setw(6) << index + 1 << "  "
                   << std::setw(8) << global_counts[index] << "  "
                   << std::setw(9) << global_dec_min[index] << "  "
                   << std::setw(9) << global_dec_max[index] << "  "
                   << std::setw(9) << global_ra_min[index] << "  "
                   << std::setw(9) << global_ra_max[index] << '\n';
        }
        output.close();
        if (!output) {
            error = "I/O failure while writing process_rearr summary: "
                    + summary_path.string();
            return false;
        }
    }
    error.clear();
    return true;
}


// ==========================================
// Function: Generate a new expo list from the rearranged output directory
// Method: Walk the output directory recursively, skip the configured skip
//         directory (Large_Field), collect every .cat file, and write their
//         absolute paths as quoted lines to the caller-supplied list file path.
// ==========================================
bool generateRearrangedExpoList(const std::string& output_directory,
                                const std::string& list_file_path,
                                std::string& error) {
    const fs::path output_dir(output_directory);
    const fs::path list_path(list_file_path);

    std::vector<std::string> catalog_paths;
    const std::string skip_dir(ProcessRearrConfig::SKIP_DIRECTORY_NAME);

    std::error_code walk_error;
    for (auto iterator = fs::recursive_directory_iterator(
             output_dir, fs::directory_options::skip_permission_denied,
             walk_error);
         iterator != fs::recursive_directory_iterator(); ++iterator) {
        const fs::directory_entry& entry = *iterator;
        if (entry.is_directory()) {
            if (entry.path().filename() == skip_dir) {
                iterator.disable_recursion_pending();
            }
            continue;
        }
        if (!entry.is_regular_file()) {
            continue;
        }
        if (entry.path().extension()
            == ProcessRearrConfig::SUBCAT_EXTENSION) {
            catalog_paths.push_back(
                fs::absolute(entry.path()).lexically_normal().string());
        }
    }
    if (walk_error) {
        error = "process_rearr failed to walk output directory: "
                + output_directory + ": " + walk_error.message();
        return false;
    }

    std::sort(catalog_paths.begin(), catalog_paths.end());

    std::ofstream output(list_path, std::ios::trunc);
    if (!output.is_open()) {
        error = "process_rearr cannot write rearranged expo list: "
                + list_path.string();
        return false;
    }
    for (const std::string& path : catalog_paths) {
        output << '"' << path << '"' << '\n';
    }
    output.close();
    if (!output) {
        error = "process_rearr I/O failure while writing rearranged expo list: "
                + list_path.string();
        return false;
    }
    error.clear();
    return true;
}

}  // namespace

// ==========================================
// Function: Rearrange exposure _all.cat files into spatial subcatalogs
// Method: Read catalogs across MPI ranks, build a global weighted k-d sky
//         partition, redistribute complete rows, and write sorted outputs.
// ==========================================
int process_rearr(const std::string& exposure_list,
                  const ProcessConfig::RuntimeOptions& options,
                  MPI_Comm communicator) {
    int rank = 0;
    int world_size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &world_size);

    ProcessRearr::CatalogLayout layout;
    std::string local_error;
    bool local_success =
        ProcessRearr::resolveCatalogLayout(options, layout, local_error);
    if (!collectiveSuccess(local_success, local_error, "layout", rank,
                           world_size, communicator)) {
        return 1;
    }

    PreparedInputs prepared;
    int preparation_ok = 1;
    if (rank == 0 && !prepareInputs(exposure_list, options, prepared, local_error)) {
        preparation_ok = 0;
    }
    MPI_Bcast(&preparation_ok, 1, MPI_INT, 0, communicator);
    if (preparation_ok == 0) {
        if (rank == 0) {
            std::cerr << "process_rearr [prepare]: " << local_error << std::endl;
        }
        return 1;
    }

    std::string catalog_paths_error;
    std::string dataset_root_error;
    std::string output_directory_error;
    std::string header_error;
    const bool catalog_paths_ok =
        broadcastStringVector(prepared.catalog_paths, 0, rank, communicator,
                              catalog_paths_error);
    const bool dataset_root_ok =
        broadcastString(prepared.dataset_root, 0, rank, communicator,
                        dataset_root_error);
    const bool output_directory_ok =
        broadcastString(prepared.output_directory, 0, rank, communicator,
                        output_directory_error);
    const bool header_ok =
        broadcastString(prepared.header, 0, rank, communicator, header_error);
    local_success = catalog_paths_ok && dataset_root_ok
                    && output_directory_ok && header_ok;
    if (!catalog_paths_ok) {
        local_error = catalog_paths_error;
    } else if (!dataset_root_ok) {
        local_error = dataset_root_error;
    } else if (!output_directory_ok) {
        local_error = output_directory_error;
    } else if (!header_ok) {
        local_error = header_error;
    }
    if (!collectiveSuccess(local_success, local_error, "broadcast", rank,
                           world_size, communicator)) {
        return 1;
    }

    LocalRows local_rows;
    local_success = readLocalCatalogs(prepared, layout, rank, world_size,
                                      local_rows, local_error);
    if (!collectiveSuccess(local_success, local_error, "read", rank,
                           world_size, communicator)) {
        return 1;
    }

    std::vector<std::uint64_t> local_tile_counts;
    std::vector<std::size_t> row_tiles;
    local_success = binLocalRows(local_rows, layout, local_tile_counts,
                                 row_tiles, local_error);
    if (!collectiveSuccess(local_success, local_error, "bin", rank,
                           world_size, communicator)) {
        return 1;
    }

    std::vector<std::uint64_t> global_tile_counts(
        ProcessRearrConfig::SKY_TILE_COUNT, 0);
    if (MPI_Allreduce(local_tile_counts.data(), global_tile_counts.data(),
                      static_cast<int>(ProcessRearrConfig::SKY_TILE_COUNT),
                      MPI_UINT64_T, MPI_SUM, communicator)
        != MPI_SUCCESS) {
        if (rank == 0) {
            std::cerr << "process_rearr: MPI_Allreduce failed for sky tiles"
                      << std::endl;
        }
        return 1;
    }
    std::vector<std::uint64_t>().swap(local_tile_counts);

    std::uint64_t total_rows = 0;
    local_success = sumTileCounts(global_tile_counts, total_rows, local_error);
    const std::size_t partition_count =
        local_success ? ProcessRearr::partitionCount(total_rows) : 0;
    std::vector<int> tile_partitions;
    if (local_success) {
        local_success = ProcessRearr::buildTilePartitions(
            global_tile_counts, partition_count, tile_partitions, local_error);
    }
    if (!collectiveSuccess(local_success, local_error, "partition", rank,
                           world_size, communicator)) {
        return 1;
    }
    std::vector<std::uint64_t>().swap(global_tile_counts);

    std::vector<int> row_partitions(row_tiles.size(), 0);
    for (std::size_t row = 0; row < row_tiles.size(); ++row) {
        row_partitions[row] = tile_partitions[row_tiles[row]];
        if (row_partitions[row] <= 0
            || static_cast<std::size_t>(row_partitions[row]) > partition_count) {
            local_success = false;
            local_error = "row maps to an invalid weighted k-d partition";
            break;
        }
    }
    std::vector<std::size_t>().swap(row_tiles);
    std::vector<int>().swap(tile_partitions);
    if (!collectiveSuccess(local_success, local_error, "map", rank,
                           world_size, communicator)) {
        return 1;
    }

    TransferPlan transfer_plan;
    local_success = buildSendRowCounts(row_partitions, world_size,
                                       transfer_plan, local_error);
    if (!collectiveSuccess(local_success, local_error, "send-counts", rank,
                           world_size, communicator)) {
        return 1;
    }
    local_success = completeTransferPlan(layout.all_columns, transfer_plan,
                                         communicator, local_error);
    if (!collectiveSuccess(local_success, local_error, "transfer-plan", rank,
                           world_size, communicator)) {
        return 1;
    }

    ReceivedRows received;
    local_success = exchangeRows(local_rows, row_partitions, layout,
                                 transfer_plan, world_size, communicator,
                                 received, local_error);
    if (!collectiveSuccess(local_success, local_error, "redistribute", rank,
                           world_size, communicator)) {
        return 1;
    }
    const std::uint64_t local_missing = local_rows.missing_catalogs;
    const std::uint64_t local_malformed = local_rows.malformed_rows;
    local_rows = LocalRows{};
    std::vector<int>().swap(row_partitions);

    local_success = true;
    if (rank == 0) {
        std::error_code filesystem_error;
        fs::create_directories(prepared.output_directory, filesystem_error);
        if (filesystem_error
            || !fs::is_directory(prepared.output_directory, filesystem_error)) {
            local_success = false;
            local_error = "cannot create output directory "
                          + prepared.output_directory + ": "
                          + filesystem_error.message();
        }
    }
    if (!collectiveSuccess(local_success, local_error, "output-directory", rank,
                           world_size, communicator)) {
        return 1;
    }
    MPI_Barrier(communicator);

    std::vector<std::uint64_t> local_summary_counts;
    std::vector<double> local_dec_min;
    std::vector<double> local_dec_max;
    std::vector<double> local_ra_min;
    std::vector<double> local_ra_max;
    local_success = writeLocalPartitions(
        received, layout, prepared, partition_count, rank, world_size,
        local_summary_counts, local_dec_min, local_dec_max,
        local_ra_min, local_ra_max, local_error);
    if (!collectiveSuccess(local_success, local_error, "write", rank,
                           world_size, communicator)) {
        return 1;
    }

    local_success = reduceAndWriteSummary(
        prepared, local_summary_counts, local_dec_min, local_dec_max,
        local_ra_min, local_ra_max, rank, communicator, local_error);
    if (!collectiveSuccess(local_success, local_error, "summary", rank,
                           world_size, communicator)) {
        return 1;
    }

    std::uint64_t global_missing = 0;
    std::uint64_t global_malformed = 0;
    const int missing_status =
        MPI_Reduce(&local_missing, &global_missing, 1, MPI_UINT64_T,
                   MPI_SUM, 0, communicator);
    const int malformed_status =
        MPI_Reduce(&local_malformed, &global_malformed, 1, MPI_UINT64_T,
                   MPI_SUM, 0, communicator);
    local_success = missing_status == MPI_SUCCESS
                    && malformed_status == MPI_SUCCESS;
    local_error = local_success
                      ? std::string{}
                      : "MPI_Reduce failed for process_rearr skip counters";
    if (!collectiveSuccess(local_success, local_error, "report", rank,
                           world_size, communicator)) {
        return 1;
    }
    MPI_Barrier(communicator);

    // Generate the rearranged expo list for downstream processing.
    local_success = true;
    std::string rearranged_list_path;
    if (rank == 0) {
        const std::string& configured_list_dir(
            options.rearranged_expo_list_directory);
        const fs::path list_dir = configured_list_dir.empty()
            ? fs::path(exposure_list).parent_path()
            : fs::path(configured_list_dir);
        rearranged_list_path =
            fs::absolute(list_dir
                         / options.rearranged_expo_list_filename)
                .lexically_normal().string();
        if (!generateRearrangedExpoList(prepared.output_directory,
                                        rearranged_list_path,
                                        local_error)) {
            local_success = false;
        }
    }
    if (!collectiveSuccess(local_success, local_error, "expo-list", rank,
                           world_size, communicator)) {
        return 1;
    }

    if (rank == 0) {
        std::cout << "process_rearr completed: rows=" << total_rows
                  << " partitions=" << partition_count
                  << " missing_catalogs=" << global_missing
                  << " malformed_rows=" << global_malformed
                  << " output=" << prepared.output_directory
                  << " expo_list=" << rearranged_list_path << std::endl;
    }
    return 0;
}
