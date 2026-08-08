#ifndef PROCESS_EXTCAT_PROCESS_EXTCAT_HPP
#define PROCESS_EXTCAT_PROCESS_EXTCAT_HPP

#include "ProcessConfig.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace ProcessExtcat {

constexpr std::size_t kCanonicalColumnCount = 18;

enum class Delimiter {
    Auto,
    Whitespace,
    Comma,
    Tab,
};

enum class HeaderMode {
    Auto,
    Present,
    Absent,
};

enum class MalformedPolicy {
    Fail,
    Skip,
};

enum class ExistingPolicy {
    Fail,
    Overwrite,
};

// ==========================================
// Structure: Configure standalone or integrated external-catalog tiling
// Method: Keep filesystem, schema, parsing, task-size, and lifecycle controls in
//         a reusable form below the pipeline RuntimeOptions adapter.
// ==========================================
struct Config {
    std::filesystem::path input_directory;
    std::filesystem::path output_directory;
    std::vector<std::string> filename_tokens;
    bool recursive = true;
    Delimiter delimiter = Delimiter::Auto;
    HeaderMode header_mode = HeaderMode::Auto;
    MalformedPolicy malformed_policy = MalformedPolicy::Fail;
    ExistingPolicy existing_policy = ExistingPolicy::Fail;
    std::uint64_t chunk_bytes = 64ULL * 1024ULL * 1024ULL;
    bool use_explicit_columns = false;
    std::vector<std::size_t> input_columns = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    };
    bool use_explicit_coordinate_columns = false;
    std::size_t ra_column = 4;
    std::size_t dec_column = 5;
};

// ==========================================
// Function: Return the legacy DES Y6 GOLD reference schema
// Method: Retain the immutable 18-field API reference without constraining generated width.
// ==========================================
const std::array<std::string, kCanonicalColumnCount>& canonicalColumnNames();

// ==========================================
// Function: Normalize and validate external-catalog tiler configuration
// Method: Canonicalize paths and reject unsafe overlap, empty filters, invalid mappings,
//         and unsupported task sizes before collective processing begins.
// ==========================================
void normalizeAndValidateConfig(Config& config);

}  // namespace ProcessExtcat

// ==========================================
// Function: Split raw external catalogs into one-degree sky tiles
// Method: Discover and inspect inputs on rank zero, optionally project any ordered column
//         subset, process newline-aligned byte ranges, and publish deterministic tiles.
// ==========================================
int process_extcat(ProcessExtcat::Config config, MPI_Comm communicator = MPI_COMM_WORLD);

// ==========================================
// Function: Run external-catalog tiling from unified pipeline options
// Method: Translate ProcessConfig values into the reusable tiler configuration and
//         participate collectively on the pipeline-owned MPI communicator.
// ==========================================
int process_extcat(const ProcessConfig::RuntimeOptions& options,
                   MPI_Comm communicator = MPI_COMM_WORLD);

#endif  // PROCESS_EXTCAT_PROCESS_EXTCAT_HPP
