#include "process_astrocat/process_astrocat.hpp"

#include "general/MPIUtils.hpp"
#include "general/MPIScheduler.hpp"
#include "general/PathUtils.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

constexpr int kRaTileCount = 360;
constexpr int kDecTileCount = 180;
constexpr std::uint32_t kTileCount =
    static_cast<std::uint32_t>(kRaTileCount * kDecTileCount);

// ==========================================
// Structure: Store one parsed Gaia coordinate and its output tile
// Method: Retain round-trip double coordinates in one standard-layout record
//         that can be transferred with a resized MPI derived datatype.
// ==========================================
struct AstroRecord {
    double ra = 0.0;
    double dec = 0.0;
    std::uint32_t tile_id = 0;
};

// ==========================================
// Structure: Route one edge coordinate to a neighboring tile owner
// Method: Preserve the source tile for canonical selection while naming the
//         target tile whose local records must be compared.
// ==========================================
struct BoundaryRecord {
    AstroRecord record;
    std::uint32_t target_tile_id = 0;
};

static_assert(std::is_standard_layout<AstroRecord>::value,
              "AstroRecord must remain standard-layout for MPI");
static_assert(std::is_standard_layout<BoundaryRecord>::value,
              "BoundaryRecord must remain standard-layout for MPI");

// ==========================================
// Structure: Hold normalized root-prepared workflow inputs
// Method: Discover and sort every regular raw file once on rank zero, then
//         broadcast only path strings needed by dynamic workers.
// ==========================================
struct PreparedInputs {
    std::vector<std::string> catalog_paths;
};

// ==========================================
// Structure: Describe one standard-int MPI Alltoallv layout
// Method: Keep per-rank counts and displacements plus validated aggregate
//         element counts for contiguous record buffers.
// ==========================================
struct TransferPlan {
    std::vector<int> send_counts;
    std::vector<int> receive_counts;
    std::vector<int> send_displacements;
    std::vector<int> receive_displacements;
    int total_send = 0;
    int total_receive = 0;
};

// ==========================================
// Function: Report and combine one rank-local phase result
// Method: Reduce success with MPI_MIN and serialize failing-rank messages so
//         every process follows the same collective exit path.
// ==========================================
bool collectiveSuccess(bool local_success,
                       const std::string& local_error,
                       const std::string& stage) {
    bool global_success = false;
    std::string reduction_error;
    if (!MPIUtils::allRanksSucceeded(
            local_success, global_success, reduction_error)) {
        if (MPIScheduler::state.rank == 0) {
            std::cerr << "process_astrocat [" << stage
                      << "]: " << reduction_error << std::endl;
        }
        return false;
    }
    if (global_success) {
        return true;
    }
    for (int reporting_rank = 0;
         reporting_rank < MPIScheduler::state.size; ++reporting_rank) {
        if (MPIScheduler::state.rank == reporting_rank && !local_success) {
            std::cerr << "process_astrocat [" << stage << "] rank "
                      << MPIScheduler::state.rank << ": "
                      << (local_error.empty() ? "unspecified error"
                                              : local_error)
                      << std::endl;
        }
        MPIScheduler::barrier();
    }
    return false;
}

// ==========================================
// Function: Recognize one generated Type-2 Gaia tile basename
// Method: Validate the fixed prefix, digit/sign fields, separators, and suffix
//         without matching unrelated files in the output directory.
// ==========================================
bool isGeneratedTileFilename(const std::string& filename) {
    constexpr std::size_t kExpectedSize = 33;
    if (filename.size() != kExpectedSize
        || filename.compare(0, 10, "des_y6_RA_") != 0
        || filename.compare(13, 1, "_") != 0
        || filename.compare(17, 5, "_Dec_") != 0
        || filename.compare(25, 1, "_") != 0
        || filename.compare(29, 4, ".dat") != 0) {
        return false;
    }
    const std::array<std::size_t, 10> digit_positions = {
        10, 11, 12, 14, 15, 16, 23, 24, 27, 28};
    for (const std::size_t position : digit_positions) {
        if (filename[position] < '0' || filename[position] > '9') {
            return false;
        }
    }
    return (filename[22] == 'p' || filename[22] == 'm')
           && (filename[26] == 'p' || filename[26] == 'm');
}

// ==========================================
// Function: List generated tile files in one output directory
// Method: Inspect only direct regular children and retain exact normalized
//         Type-2 tile basenames for policy checks and cleanup.
// ==========================================
bool listGeneratedTiles(const fs::path& output_directory,
                        std::vector<fs::path>& files,
                        std::string& error) {
    files.clear();
    std::error_code filesystem_error;
    if (!fs::exists(output_directory, filesystem_error)) {
        if (filesystem_error) {
            error = "cannot inspect output directory "
                    + output_directory.string() + ": "
                    + filesystem_error.message();
            return false;
        }
        return true;
    }
    if (!fs::is_directory(output_directory, filesystem_error)
        || filesystem_error) {
        error = "process_astrocat output is not a directory: "
                + output_directory.string();
        return false;
    }
    for (fs::directory_iterator iterator(output_directory, filesystem_error);
         !filesystem_error && iterator != fs::directory_iterator();
         iterator.increment(filesystem_error)) {
        const fs::directory_entry& entry = *iterator;
        std::error_code type_error;
        if (entry.is_regular_file(type_error)
            && !type_error
            && isGeneratedTileFilename(entry.path().filename().string())) {
            files.push_back(entry.path());
        }
    }
    if (filesystem_error) {
        error = "cannot scan output directory " + output_directory.string()
                + ": " + filesystem_error.message();
        return false;
    }
    std::sort(files.begin(), files.end());
    error.clear();
    return true;
}

// ==========================================
// Function: Remove every generated Type-2 tile
// Method: Re-list exact generated basenames and delete only those files while
//         preserving all unrelated output-directory content.
// ==========================================
bool removeGeneratedTiles(const fs::path& output_directory,
                          std::string& error) {
    std::vector<fs::path> files;
    if (!listGeneratedTiles(output_directory, files, error)) {
        return false;
    }
    for (const fs::path& file : files) {
        std::error_code remove_error;
        if (!fs::remove(file, remove_error) || remove_error) {
            error = "cannot remove generated tile " + file.string()
                    + ": " + remove_error.message();
            return false;
        }
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Discover raw Gaia catalogs on rank zero
// Method: Select direct regular children only, normalize their absolute paths,
//         sort lexically, and reject an empty or int-oversized job list.
// ==========================================
bool discoverInputFiles(const ProcessAstrocat::Config& config,
                        std::vector<std::string>& catalog_paths,
                        std::string& error) {
    std::error_code filesystem_error;
    if (!fs::is_directory(config.input_directory, filesystem_error)
        || filesystem_error) {
        error = "process_astrocat input is not a readable directory: "
                + config.input_directory.string();
        return false;
    }
    catalog_paths.clear();
    for (fs::directory_iterator iterator(config.input_directory,
                                         filesystem_error);
         !filesystem_error && iterator != fs::directory_iterator();
         iterator.increment(filesystem_error)) {
        const fs::directory_entry& entry = *iterator;
        std::error_code type_error;
        if (entry.is_regular_file(type_error) && !type_error) {
            catalog_paths.push_back(
                PathUtils::normalizedAbsolute(entry.path()).string());
        }
    }
    if (filesystem_error) {
        error = "cannot scan process_astrocat input directory: "
                + filesystem_error.message();
        return false;
    }
    std::sort(catalog_paths.begin(), catalog_paths.end());
    if (catalog_paths.empty()) {
        error = "process_astrocat found no regular input files";
        return false;
    }
    if (catalog_paths.size()
        > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        error = "process_astrocat input file count exceeds MPI scheduler capacity";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Prepare rank-zero input and output state
// Method: Discover sorted jobs, enforce the existing-tile policy before costly
//         work, and create the independent output directory if needed.
// ==========================================
bool prepareInputs(const ProcessAstrocat::Config& config,
                   PreparedInputs& prepared,
                   std::string& error) {
    if (!discoverInputFiles(config, prepared.catalog_paths, error)) {
        return false;
    }
    std::error_code create_error;
    fs::create_directories(config.output_directory, create_error);
    if (create_error
        || !fs::is_directory(config.output_directory, create_error)
        || create_error) {
        error = "cannot create process_astrocat output directory "
                + config.output_directory.string() + ": "
                + create_error.message();
        return false;
    }
    std::vector<fs::path> existing_tiles;
    if (!listGeneratedTiles(config.output_directory, existing_tiles, error)) {
        return false;
    }
    if (config.existing_policy == ProcessAstrocat::ExistingPolicy::Fail
        && !existing_tiles.empty()) {
        error = "generated astrometry tiles already exist in "
                + config.output_directory.string();
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Convert one valid sky coordinate to its one-degree tile id
// Method: Wrap exactly RA=360 to zero, place exactly Dec=90 below its upper
//         boundary, then encode Dec-major and RA-minor integer bins.
// ==========================================
std::uint32_t tileId(double ra, double dec) {
    const double normalized_ra = ra == 360.0 ? 0.0 : ra;
    const double bounded_dec =
        dec == 90.0
            ? std::nextafter(
                  90.0, -std::numeric_limits<double>::infinity())
            : dec;
    const int ra_lower = static_cast<int>(std::floor(normalized_ra));
    const int dec_lower = static_cast<int>(std::floor(bounded_dec));
    return static_cast<std::uint32_t>(
        (dec_lower + 90) * kRaTileCount + ra_lower);
}

// ==========================================
// Function: Read one complete raw Gaia catalog
// Method: Optionally skip exactly one header, accept comma or whitespace
//         separators, silently ignore non-two-double rows, and assign tiles.
// ==========================================
void readOneCatalog(const std::string& path,
                    bool add_header,
                    std::vector<AstroRecord>& rows,
                    std::uint64_t& parsed_rows) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("cannot open raw Gaia catalog: " + path);
    }
    if (!add_header) {
        std::string header;
        std::getline(input, header);
    }
    std::string line;
    while (std::getline(input, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream stream(line);
        double ra = 0.0;
        double dec = 0.0;
        if (!(stream >> ra >> dec)) {
            continue;
        }
        if (parsed_rows == std::numeric_limits<std::uint64_t>::max()) {
            throw std::runtime_error(
                "local parsed-row count exceeds uint64 capacity");
        }
        const double stored_ra = ra == 360.0 ? 0.0 : ra;
        rows.push_back(AstroRecord{stored_ra, dec, tileId(ra, dec)});
        ++parsed_rows;
    }
    if (input.bad()) {
        throw std::runtime_error("I/O failure while reading raw Gaia catalog: "
                                 + path);
    }
}

// ==========================================
// Function: Build int-safe displacements for MPI records
// Method: Accumulate signed counts in long long and reject any count,
//         displacement, or total outside the standard Alltoallv int range.
// ==========================================
bool buildDisplacements(const std::vector<int>& counts,
                        std::vector<int>& displacements,
                        int& total,
                        std::string& error) {
    displacements.resize(counts.size());
    long long running = 0;
    for (std::size_t index = 0; index < counts.size(); ++index) {
        if (counts[index] < 0
            || running > std::numeric_limits<int>::max()) {
            error = "MPI record counts or displacements exceed int range";
            return false;
        }
        displacements[index] = static_cast<int>(running);
        running += counts[index];
    }
    if (running > std::numeric_limits<int>::max()) {
        error = "MPI total record count exceeds int range";
        return false;
    }
    total = static_cast<int>(running);
    error.clear();
    return true;
}

// ==========================================
// Function: Exchange per-destination record counts
// Method: Use MPI_Alltoall, then derive int-safe send and receive
//         displacements before any variable-count payload collective.
// ==========================================
bool completeTransferPlan(TransferPlan& plan, std::string& error) {
    plan.receive_counts.assign(plan.send_counts.size(), 0);
    if (MPI_Alltoall(plan.send_counts.data(), 1, MPI_INT,
                     plan.receive_counts.data(), 1, MPI_INT,
                     MPIScheduler::state.communicator)
        != MPI_SUCCESS) {
        error = "MPI_Alltoall failed for process_astrocat record counts";
        return false;
    }
    return buildDisplacements(plan.send_counts, plan.send_displacements,
                              plan.total_send, error)
           && buildDisplacements(plan.receive_counts,
                                 plan.receive_displacements,
                                 plan.total_receive, error);
}

// ==========================================
// Function: Create the MPI datatype for one AstroRecord
// Method: Describe both doubles and the uint32 tile id, then resize the extent
//         to sizeof(AstroRecord) so arrays preserve compiler padding.
// ==========================================
bool createAstroRecordMpiType(MPI_Datatype& record_type,
                              std::string& error) {
    const int block_lengths[3] = {1, 1, 1};
    const MPI_Aint displacements[3] = {
        static_cast<MPI_Aint>(offsetof(AstroRecord, ra)),
        static_cast<MPI_Aint>(offsetof(AstroRecord, dec)),
        static_cast<MPI_Aint>(offsetof(AstroRecord, tile_id))};
    MPI_Datatype member_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_UINT32_T};
    MPI_Datatype raw_type = MPI_DATATYPE_NULL;
    if (MPI_Type_create_struct(3, block_lengths, displacements, member_types,
                               &raw_type)
        != MPI_SUCCESS) {
        error = "cannot create AstroRecord MPI datatype";
        return false;
    }
    const int resize_status = MPI_Type_create_resized(
        raw_type, 0, static_cast<MPI_Aint>(sizeof(AstroRecord)), &record_type);
    MPI_Type_free(&raw_type);
    if (resize_status != MPI_SUCCESS
        || MPI_Type_commit(&record_type) != MPI_SUCCESS) {
        if (record_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&record_type);
        }
        error = "cannot commit AstroRecord MPI datatype";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Redistribute parsed coordinates to unique tile owners
// Method: Count tile_id modulo world size, exchange layouts, pack one AoS send
//         buffer, release local input memory, and call one MPI_Alltoallv.
// ==========================================
bool redistributeRows(std::vector<AstroRecord>& local_rows,
                      std::vector<AstroRecord>& received_rows,
                      std::string& error) {
    const int world_size = MPIScheduler::state.size;
    TransferPlan plan;
    plan.send_counts.assign(static_cast<std::size_t>(world_size), 0);
    bool local_plan_success = true;
    for (const AstroRecord& record : local_rows) {
        const int destination = static_cast<int>(
            record.tile_id % static_cast<std::uint32_t>(world_size));
        int& count = plan.send_counts[static_cast<std::size_t>(destination)];
        if (count == std::numeric_limits<int>::max()) {
            error = "one tile-owner send count exceeds MPI int capacity";
            local_plan_success = false;
            break;
        }
        ++count;
    }
    if (!collectiveSuccess(local_plan_success, error, "send-counts")) {
        return false;
    }
    local_plan_success = completeTransferPlan(plan, error);
    if (!collectiveSuccess(
            local_plan_success, error, "transfer-plan")) {
        return false;
    }

    std::vector<AstroRecord> send_rows(
        static_cast<std::size_t>(plan.total_send));
    std::vector<int> cursors = plan.send_displacements;
    for (const AstroRecord& record : local_rows) {
        const int destination = static_cast<int>(
            record.tile_id % static_cast<std::uint32_t>(world_size));
        send_rows[static_cast<std::size_t>(
            cursors[static_cast<std::size_t>(destination)]++)] = record;
    }
    std::vector<AstroRecord>().swap(local_rows);
    received_rows.resize(static_cast<std::size_t>(plan.total_receive));

    MPI_Datatype record_type = MPI_DATATYPE_NULL;
    const bool datatype_success =
        createAstroRecordMpiType(record_type, error);
    if (!collectiveSuccess(datatype_success, error, "record-datatype")) {
        if (record_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&record_type);
        }
        return false;
    }
    const int transfer_status = MPI_Alltoallv(
        send_rows.data(), plan.send_counts.data(),
        plan.send_displacements.data(), record_type,
        received_rows.data(), plan.receive_counts.data(),
        plan.receive_displacements.data(), record_type,
        MPIScheduler::state.communicator);
    MPI_Type_free(&record_type);
    if (transfer_status != MPI_SUCCESS) {
        error = "MPI_Alltoallv failed while redistributing Gaia rows";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Convert one double to a monotonic IEEE-754 key
// Method: Complement negative bit patterns and set the positive sign offset so
//         unsigned key distance equals representable-step distance.
// ==========================================
std::uint64_t orderedDoubleKey(double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    constexpr std::uint64_t kSignBit = UINT64_C(1) << 63;
    return (bits & kSignBit) != 0 ? ~bits : bits | kSignBit;
}

// ==========================================
// Function: Measure representable double-key separation
// Method: Subtract monotonic unsigned keys in sorted order without signed
//         overflow or a decimal tolerance.
// ==========================================
std::uint64_t ulpDistance(double first, double second) {
    const std::uint64_t first_key = orderedDoubleKey(first);
    const std::uint64_t second_key = orderedDoubleKey(second);
    return first_key > second_key ? first_key - second_key
                                  : second_key - first_key;
}

// ==========================================
// Function: Test the frozen coordinate duplicate predicate
// Method: Require both RA and Dec to differ by at most one representable
//         double step, including exact equality.
// ==========================================
bool isDuplicate(const AstroRecord& first, const AstroRecord& second) {
    return ulpDistance(first.ra, second.ra) <= 1
           && ulpDistance(first.dec, second.dec) <= 1;
}

// ==========================================
// Function: Order rows deterministically for grouping and canonical output
// Method: Compare tile id, monotonic RA key, then monotonic Dec key so signed
//         zero and representable neighbors have stable positions.
// ==========================================
bool recordLess(const AstroRecord& first, const AstroRecord& second) {
    if (first.tile_id != second.tile_id) {
        return first.tile_id < second.tile_id;
    }
    const std::uint64_t first_ra = orderedDoubleKey(first.ra);
    const std::uint64_t second_ra = orderedDoubleKey(second.ra);
    if (first_ra != second_ra) {
        return first_ra < second_ra;
    }
    return orderedDoubleKey(first.dec) < orderedDoubleKey(second.dec);
}

// ==========================================
// Function: Test whether a sorted key group has a one-step neighbor
// Method: Lower-bound the candidate and inspect only that position and its
//         predecessor, which are sufficient for a sorted one-dimensional set.
// ==========================================
bool containsNearKey(const std::vector<std::uint64_t>& keys,
                     std::uint64_t candidate) {
    const auto position = std::lower_bound(keys.begin(), keys.end(), candidate);
    if (position != keys.end()) {
        const std::uint64_t distance =
            *position > candidate ? *position - candidate
                                  : candidate - *position;
        if (distance <= 1) {
            return true;
        }
    }
    if (position != keys.begin()) {
        const std::uint64_t previous = *(position - 1);
        const std::uint64_t distance =
            previous > candidate ? previous - candidate
                                 : candidate - previous;
        return distance <= 1;
    }
    return false;
}

// ==========================================
// Function: De-duplicate sorted records within each owned tile
// Method: Retain sorted Dec keys for the current RA key and its one-step
//         predecessor so non-adjacent two-dimensional duplicates are found.
// ==========================================
std::uint64_t deduplicateOwnedTiles(std::vector<AstroRecord>& rows) {
    std::sort(rows.begin(), rows.end(), recordLess);
    std::vector<AstroRecord> unique_rows;
    unique_rows.reserve(rows.size());
    std::uint64_t duplicate_count = 0;

    std::size_t tile_begin = 0;
    while (tile_begin < rows.size()) {
        std::size_t tile_end = tile_begin + 1;
        while (tile_end < rows.size()
               && rows[tile_end].tile_id == rows[tile_begin].tile_id) {
            ++tile_end;
        }

        std::uint64_t current_ra_key = orderedDoubleKey(rows[tile_begin].ra);
        std::uint64_t previous_ra_key = 0;
        bool previous_group_valid = false;
        std::vector<std::uint64_t> current_dec_keys;
        std::vector<std::uint64_t> previous_dec_keys;

        for (std::size_t index = tile_begin; index < tile_end; ++index) {
            const std::uint64_t ra_key = orderedDoubleKey(rows[index].ra);
            const std::uint64_t dec_key = orderedDoubleKey(rows[index].dec);
            if (ra_key != current_ra_key) {
                previous_ra_key = current_ra_key;
                previous_dec_keys = std::move(current_dec_keys);
                current_dec_keys.clear();
                current_ra_key = ra_key;
                previous_group_valid =
                    current_ra_key > previous_ra_key
                    && current_ra_key - previous_ra_key <= 1;
                if (!previous_group_valid) {
                    previous_dec_keys.clear();
                }
            }
            const bool duplicate_in_current =
                containsNearKey(current_dec_keys, dec_key);
            const bool duplicate_in_previous =
                previous_group_valid
                && containsNearKey(previous_dec_keys, dec_key);
            if (duplicate_in_current || duplicate_in_previous) {
                ++duplicate_count;
                continue;
            }
            current_dec_keys.push_back(dec_key);
            unique_rows.push_back(rows[index]);
        }
        tile_begin = tile_end;
    }
    rows.swap(unique_rows);
    return duplicate_count;
}

// ==========================================
// Function: Encode one validated pair of lower tile bounds
// Method: Convert RA [0,359] and Dec [-90,89] to the shared Dec-major tile id.
// ==========================================
std::uint32_t tileIdFromLower(int ra_lower, int dec_lower) {
    return static_cast<std::uint32_t>(
        (dec_lower + 90) * kRaTileCount + ra_lower);
}

// ==========================================
// Function: Wrap one RA lower bound into the 360-tile circle
// Method: Apply integer modulo with a positive correction for the 359/0 edge.
// ==========================================
int wrapRaLower(int lower) {
    int wrapped = lower % kRaTileCount;
    if (wrapped < 0) {
        wrapped += kRaTileCount;
    }
    return wrapped;
}

// ==========================================
// Function: Find neighboring tile ids within one coordinate ULP
// Method: Tile nextafter values in both directions, combine changed RA and Dec
//         lower bounds, wrap RA, and return up to three distinct neighbors.
// ==========================================
std::vector<std::uint32_t> neighboringTileIds(const AstroRecord& record) {
    const int primary_ra = static_cast<int>(record.tile_id % kRaTileCount);
    const int primary_dec =
        static_cast<int>(record.tile_id / kRaTileCount) - 90;
    std::vector<int> ra_lowers{primary_ra};
    std::vector<int> dec_lowers{primary_dec};

    for (const double direction : {
             -std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()}) {
        const double adjacent_ra = std::nextafter(record.ra, direction);
        const int ra_lower = wrapRaLower(
            static_cast<int>(std::floor(adjacent_ra)));
        ra_lowers.push_back(ra_lower);

        double adjacent_dec = std::nextafter(record.dec, direction);
        if (adjacent_dec == 90.0) {
            adjacent_dec = std::nextafter(
                90.0, -std::numeric_limits<double>::infinity());
        }
        const int dec_lower = static_cast<int>(std::floor(adjacent_dec));
        if (dec_lower >= -90 && dec_lower < 90) {
            dec_lowers.push_back(dec_lower);
        }
    }
    std::sort(ra_lowers.begin(), ra_lowers.end());
    ra_lowers.erase(std::unique(ra_lowers.begin(), ra_lowers.end()),
                    ra_lowers.end());
    std::sort(dec_lowers.begin(), dec_lowers.end());
    dec_lowers.erase(std::unique(dec_lowers.begin(), dec_lowers.end()),
                     dec_lowers.end());

    std::vector<std::uint32_t> neighbors;
    for (const int dec_lower : dec_lowers) {
        for (const int ra_lower : ra_lowers) {
            const std::uint32_t candidate =
                tileIdFromLower(ra_lower, dec_lower);
            if (candidate != record.tile_id && candidate < kTileCount) {
                neighbors.push_back(candidate);
            }
        }
    }
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                    neighbors.end());
    return neighbors;
}

// ==========================================
// Function: Create the MPI datatype for one BoundaryRecord
// Method: Describe nested coordinate fields and both tile identifiers, then
//         resize to the compiler's complete structure extent.
// ==========================================
bool createBoundaryRecordMpiType(MPI_Datatype& record_type,
                                 std::string& error) {
    const int block_lengths[4] = {1, 1, 1, 1};
    const MPI_Aint record_offset =
        static_cast<MPI_Aint>(offsetof(BoundaryRecord, record));
    const MPI_Aint displacements[4] = {
        record_offset + static_cast<MPI_Aint>(offsetof(AstroRecord, ra)),
        record_offset + static_cast<MPI_Aint>(offsetof(AstroRecord, dec)),
        record_offset + static_cast<MPI_Aint>(offsetof(AstroRecord, tile_id)),
        static_cast<MPI_Aint>(offsetof(BoundaryRecord, target_tile_id))};
    MPI_Datatype member_types[4] = {
        MPI_DOUBLE, MPI_DOUBLE, MPI_UINT32_T, MPI_UINT32_T};
    MPI_Datatype raw_type = MPI_DATATYPE_NULL;
    if (MPI_Type_create_struct(4, block_lengths, displacements, member_types,
                               &raw_type)
        != MPI_SUCCESS) {
        error = "cannot create BoundaryRecord MPI datatype";
        return false;
    }
    const int resize_status = MPI_Type_create_resized(
        raw_type, 0, static_cast<MPI_Aint>(sizeof(BoundaryRecord)),
        &record_type);
    MPI_Type_free(&raw_type);
    if (resize_status != MPI_SUCCESS
        || MPI_Type_commit(&record_type) != MPI_SUCCESS) {
        if (record_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&record_type);
        }
        error = "cannot commit BoundaryRecord MPI datatype";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Exchange only integer-boundary coordinate candidates
// Method: Expand each record to distinct neighboring tiles, group by target
//         owner, and perform one compact MPI_Alltoallv with source tile keys.
// ==========================================
bool exchangeBoundaryCandidates(const std::vector<AstroRecord>& rows,
                                std::vector<BoundaryRecord>& received,
                                std::string& error) {
    const int world_size = MPIScheduler::state.size;
    std::vector<BoundaryRecord> local;
    for (const AstroRecord& record : rows) {
        for (const std::uint32_t neighbor : neighboringTileIds(record)) {
            local.push_back(BoundaryRecord{record, neighbor});
        }
    }

    TransferPlan plan;
    plan.send_counts.assign(static_cast<std::size_t>(world_size), 0);
    bool local_plan_success = true;
    for (const BoundaryRecord& candidate : local) {
        const int destination = static_cast<int>(
            candidate.target_tile_id
            % static_cast<std::uint32_t>(world_size));
        int& count = plan.send_counts[static_cast<std::size_t>(destination)];
        if (count == std::numeric_limits<int>::max()) {
            error = "one boundary-owner send count exceeds MPI int capacity";
            local_plan_success = false;
            break;
        }
        ++count;
    }
    if (!collectiveSuccess(
            local_plan_success, error, "boundary-send-counts")) {
        return false;
    }
    local_plan_success = completeTransferPlan(plan, error);
    if (!collectiveSuccess(
            local_plan_success, error, "boundary-transfer-plan")) {
        return false;
    }

    std::vector<BoundaryRecord> send(
        static_cast<std::size_t>(plan.total_send));
    std::vector<int> cursors = plan.send_displacements;
    for (const BoundaryRecord& candidate : local) {
        const int destination = static_cast<int>(
            candidate.target_tile_id
            % static_cast<std::uint32_t>(world_size));
        send[static_cast<std::size_t>(
            cursors[static_cast<std::size_t>(destination)]++)] = candidate;
    }
    std::vector<BoundaryRecord>().swap(local);
    received.resize(static_cast<std::size_t>(plan.total_receive));

    MPI_Datatype record_type = MPI_DATATYPE_NULL;
    const bool datatype_success =
        createBoundaryRecordMpiType(record_type, error);
    if (!collectiveSuccess(
            datatype_success, error, "boundary-datatype")) {
        if (record_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&record_type);
        }
        return false;
    }
    const int transfer_status = MPI_Alltoallv(
        send.data(), plan.send_counts.data(), plan.send_displacements.data(),
        record_type, received.data(), plan.receive_counts.data(),
        plan.receive_displacements.data(), record_type,
        MPIScheduler::state.communicator);
    MPI_Type_free(&record_type);
    if (transfer_status != MPI_SUCCESS) {
        error = "MPI_Alltoallv failed for boundary de-duplication records";
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Compare canonical duplicate ownership keys
// Method: Order by monotonic RA, monotonic Dec, then source tile id so adjacent
//         owners independently reach the same keep/remove decision.
// ==========================================
bool canonicalLess(const AstroRecord& first, const AstroRecord& second) {
    const std::uint64_t first_ra = orderedDoubleKey(first.ra);
    const std::uint64_t second_ra = orderedDoubleKey(second.ra);
    if (first_ra != second_ra) {
        return first_ra < second_ra;
    }
    const std::uint64_t first_dec = orderedDoubleKey(first.dec);
    const std::uint64_t second_dec = orderedDoubleKey(second.dec);
    if (first_dec != second_dec) {
        return first_dec < second_dec;
    }
    return first.tile_id < second.tile_id;
}

// ==========================================
// Function: Remove local rows superseded across tile boundaries
// Method: Index contiguous local tile ranges, compare received edge candidates
//         with the same one-ULP predicate, and retain the canonical smaller key.
// ==========================================
std::uint64_t applyBoundaryDedup(
    std::vector<AstroRecord>& rows,
    const std::vector<BoundaryRecord>& received) {
    std::unordered_map<std::uint32_t, std::pair<std::size_t, std::size_t>> ranges;
    for (std::size_t begin = 0; begin < rows.size();) {
        std::size_t end = begin + 1;
        while (end < rows.size() && rows[end].tile_id == rows[begin].tile_id) {
            ++end;
        }
        ranges.emplace(rows[begin].tile_id, std::make_pair(begin, end));
        begin = end;
    }

    std::vector<unsigned char> remove(rows.size(), 0);
    for (const BoundaryRecord& candidate : received) {
        const auto range = ranges.find(candidate.target_tile_id);
        if (range == ranges.end()) {
            continue;
        }
        for (std::size_t index = range->second.first;
             index < range->second.second; ++index) {
            if (isDuplicate(rows[index], candidate.record)
                && canonicalLess(candidate.record, rows[index])) {
                remove[index] = 1;
            }
        }
    }

    std::vector<AstroRecord> retained;
    retained.reserve(rows.size());
    std::uint64_t removed = 0;
    for (std::size_t index = 0; index < rows.size(); ++index) {
        if (remove[index] != 0) {
            ++removed;
        } else {
            retained.push_back(rows[index]);
        }
    }
    rows.swap(retained);
    return removed;
}

// ==========================================
// Function: Format one signed integer declination boundary
// Method: Emit the Type-2 p/m prefix plus at least two zero-padded digits.
// ==========================================
std::string formatDeclination(int value) {
    std::ostringstream stream;
    stream << (value >= 0 ? 'p' : 'm')
           << std::setw(2) << std::setfill('0') << std::abs(value);
    return stream.str();
}

// ==========================================
// Function: Format one Type-2 Gaia tile basename
// Method: Decode the shared tile id and reproduce generateGalCatFileNames'
//         exact RA/Dec boundary convention.
// ==========================================
std::string tileFilename(std::uint32_t tile_id) {
    const int ra_lower = static_cast<int>(tile_id % kRaTileCount);
    const int dec_lower = static_cast<int>(tile_id / kRaTileCount) - 90;
    std::ostringstream stream;
    stream << "des_y6_RA_"
           << std::setw(3) << std::setfill('0') << ra_lower << '_'
           << std::setw(3) << std::setfill('0') << ra_lower + 1
           << "_Dec_" << formatDeclination(dec_lower) << '_'
           << formatDeclination(dec_lower + 1) << ".dat";
    return stream.str();
}

// ==========================================
// Function: Write this rank's nonempty owned Gaia tiles
// Method: Stream sorted segments with the mandatory RA/DEC header and
//         max_digits10 precision while never sharing a file between ranks.
// ==========================================
bool writeOwnedTiles(const std::vector<AstroRecord>& rows,
                     const fs::path& output_directory,
                     std::uint64_t& generated_tiles,
                     std::string& error) {
    generated_tiles = 0;
    const int rank = MPIScheduler::state.rank;
    const int world_size = MPIScheduler::state.size;
    for (std::size_t begin = 0; begin < rows.size();) {
        std::size_t end = begin + 1;
        while (end < rows.size() && rows[end].tile_id == rows[begin].tile_id) {
            ++end;
        }
        const std::uint32_t tile_id = rows[begin].tile_id;
        if (static_cast<int>(
                tile_id % static_cast<std::uint32_t>(world_size)) != rank) {
            error = "received rows include a tile owned by another rank";
            return false;
        }
        const fs::path output_path =
            output_directory / tileFilename(tile_id);
        std::ofstream output(output_path, std::ios::trunc);
        if (!output.is_open()) {
            error = "cannot write generated Gaia tile: "
                    + output_path.string();
            return false;
        }
        output << "RA    DEC\n"
               << std::setprecision(
                      std::numeric_limits<double>::max_digits10);
        for (std::size_t index = begin; index < end; ++index) {
            output << rows[index].ra << "    " << rows[index].dec << '\n';
        }
        output.close();
        if (!output) {
            error = "I/O failure while writing generated Gaia tile: "
                    + output_path.string();
            return false;
        }
        ++generated_tiles;
        begin = end;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Prepare the final output lifecycle on rank zero
// Method: Delay overwrite deletion until all read, transfer, sort, and
//         de-duplication work succeeds, preserving unrelated files.
// ==========================================
bool prepareFinalWrite(const ProcessAstrocat::Config& config,
                       std::string& error) {
    if (config.existing_policy == ProcessAstrocat::ExistingPolicy::Overwrite
        && !removeGeneratedTiles(config.output_directory, error)) {
        return false;
    }
    error.clear();
    return true;
}

// ==========================================
// Function: Remove partial generated output after a write failure
// Method: Delete every exact generated tile basename while leaving unrelated
//         files and the output directory itself intact.
// ==========================================
bool cleanupFailedWrite(const ProcessAstrocat::Config& config,
                        std::string& error) {
    return removeGeneratedTiles(config.output_directory, error);
}

}  // namespace

namespace ProcessAstrocat {

void normalizeAndValidateConfig(Config& config) {
    if (config.input_directory.empty()) {
        throw std::invalid_argument(
            "process_astrocat input directory must not be empty");
    }
    if (config.output_directory.empty()) {
        throw std::invalid_argument(
            "process_astrocat output directory must not be empty");
    }
    config.input_directory =
        PathUtils::normalizedAbsolute(config.input_directory);
    config.output_directory =
        PathUtils::normalizedAbsolute(config.output_directory);
    if (config.input_directory == config.output_directory
        || PathUtils::isPathWithin(config.output_directory,
                                   config.input_directory)
        || PathUtils::isPathWithin(config.input_directory,
                                   config.output_directory)) {
        throw std::invalid_argument(
            "process_astrocat input and output directories must not overlap");
    }
}

}  // namespace ProcessAstrocat

int process_astrocat(ProcessAstrocat::Config config) {
    const int rank = MPIScheduler::state.rank;
    const MPI_Comm communicator = MPIScheduler::state.communicator;
    std::string local_error;
    bool local_success = true;
    try {
        ProcessAstrocat::normalizeAndValidateConfig(config);
    } catch (const std::exception& exception) {
        local_success = false;
        local_error = exception.what();
    }
    if (!collectiveSuccess(local_success, local_error, "configuration")) {
        return 1;
    }

    PreparedInputs prepared;
    local_success = true;
    if (rank == 0
        && !prepareInputs(config, prepared, local_error)) {
        local_success = false;
    }
    if (!collectiveSuccess(local_success, local_error, "prepare")) {
        return 1;
    }
    local_success = MPIUtils::broadcastStrings(
        prepared.catalog_paths, 0, local_error);
    if (!collectiveSuccess(local_success, local_error, "input-broadcast")) {
        return 1;
    }

    std::vector<AstroRecord> local_rows;
    std::uint64_t local_input_rows = 0;
    bool local_read_ok = true;
    std::string local_read_error;
    MPIScheduler::distribute(
        static_cast<int>(prepared.catalog_paths.size()),
        [&](int one_based_job) {
            if (!local_read_ok) {
                return;
            }
            try {
                readOneCatalog(
                    prepared.catalog_paths[
                        static_cast<std::size_t>(one_based_job - 1)],
                    config.add_header, local_rows, local_input_rows);
            } catch (const std::exception& exception) {
                local_read_ok = false;
                local_read_error = exception.what();
            }
        },
        "process_astrocat read");
    if (!collectiveSuccess(local_read_ok, local_read_error, "read")) {
        return 1;
    }

    std::vector<AstroRecord> owned_rows;
    local_success = redistributeRows(local_rows, owned_rows, local_error);
    if (!collectiveSuccess(local_success, local_error, "redistribute")) {
        return 1;
    }

    std::uint64_t local_duplicates = deduplicateOwnedTiles(owned_rows);
    std::vector<BoundaryRecord> boundary_records;
    local_success = exchangeBoundaryCandidates(
        owned_rows, boundary_records, local_error);
    if (!collectiveSuccess(local_success, local_error, "boundary-exchange")) {
        return 1;
    }
    local_duplicates += applyBoundaryDedup(owned_rows, boundary_records);
    std::vector<BoundaryRecord>().swap(boundary_records);

    local_success = true;
    if (rank == 0 && !prepareFinalWrite(config, local_error)) {
        local_success = false;
    }
    if (!collectiveSuccess(local_success, local_error, "output-policy")) {
        return 1;
    }
    MPI_Barrier(communicator);

    std::uint64_t local_generated_tiles = 0;
    local_success = writeOwnedTiles(
        owned_rows, config.output_directory,
        local_generated_tiles, local_error);
    if (!collectiveSuccess(local_success, local_error, "write")) {
        bool cleanup_success = true;
        std::string cleanup_error;
        if (rank == 0) {
            cleanup_success = cleanupFailedWrite(config, cleanup_error);
        }
        collectiveSuccess(cleanup_success, cleanup_error, "write-cleanup");
        MPI_Barrier(communicator);
        return 1;
    }

    const std::uint64_t local_unique_rows =
        static_cast<std::uint64_t>(owned_rows.size());
    std::uint64_t global_input_rows = 0;
    std::uint64_t global_unique_rows = 0;
    std::uint64_t global_duplicates = 0;
    std::uint64_t global_generated_tiles = 0;
    const int input_status = MPI_Reduce(
        &local_input_rows, &global_input_rows, 1, MPI_UINT64_T,
        MPI_SUM, 0, communicator);
    const int unique_status = MPI_Reduce(
        &local_unique_rows, &global_unique_rows, 1, MPI_UINT64_T,
        MPI_SUM, 0, communicator);
    const int duplicate_status = MPI_Reduce(
        &local_duplicates, &global_duplicates, 1, MPI_UINT64_T,
        MPI_SUM, 0, communicator);
    const int tile_status = MPI_Reduce(
        &local_generated_tiles, &global_generated_tiles, 1, MPI_UINT64_T,
        MPI_SUM, 0, communicator);
    local_success = input_status == MPI_SUCCESS
                    && unique_status == MPI_SUCCESS
                    && duplicate_status == MPI_SUCCESS
                    && tile_status == MPI_SUCCESS;
    local_error = local_success
                      ? std::string{}
                      : "MPI_Reduce failed for process_astrocat summary";
    if (!collectiveSuccess(local_success, local_error, "summary")) {
        return 1;
    }
    if (rank == 0) {
        std::cout << "process_astrocat completed: input_files="
                  << prepared.catalog_paths.size()
                  << " input_rows=" << global_input_rows
                  << " unique_rows=" << global_unique_rows
                  << " duplicate_rows=" << global_duplicates
                  << " generated_tiles=" << global_generated_tiles
                  << " output=" << config.output_directory.string()
                  << std::endl;
    }
    return 0;
}

int process_astrocat(const RuntimeConfig& runtime_config) {
    ProcessAstrocat::Config config;
    config.input_directory = runtime_config.astrocat.input_directory;
    config.output_directory = runtime_config.astrocat.output_directory;
    config.add_header = runtime_config.astrocat.add_header;
    if (runtime_config.astrocat.existing_policy == "overwrite") {
        config.existing_policy = ProcessAstrocat::ExistingPolicy::Overwrite;
    } else if (runtime_config.astrocat.existing_policy == "fail") {
        config.existing_policy = ProcessAstrocat::ExistingPolicy::Fail;
    } else {
        if (MPIScheduler::state.rank == 0) {
            std::cerr << "process_astrocat: existing policy must be fail or overwrite"
                      << std::endl;
        }
        return 1;
    }
    return process_astrocat(std::move(config));
}
