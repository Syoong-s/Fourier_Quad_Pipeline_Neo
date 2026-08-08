#include "process_extcat/process_extcat.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

constexpr const char* kTilePrefix = "des_y6_RA_";
constexpr const char* kTileSuffix = ".dat";

// ==========================================
// Structure: Describe one inspected input catalog
// Method: Retain its path, byte range, delimiter, output projection, raw coordinate columns,
//         and root-only output header for deterministic merging.
// ==========================================
struct FileMetadata {
    fs::path path;
    std::uint64_t data_offset = 0;
    std::uint64_t file_size = 0;
    ProcessExtcat::Delimiter delimiter = ProcessExtcat::Delimiter::Whitespace;
    std::size_t input_column_count = 0;
    std::size_t ra_column = 4;
    std::size_t dec_column = 5;
    std::vector<std::size_t> columns;
    std::vector<std::string> output_header;
};

// ==========================================
// Structure: Describe one deterministic input byte-range task
// Method: Identify the source file and half-open byte interval independently of MPI rank.
// ==========================================
struct Task {
    std::uint64_t id = 0;
    std::size_t file_index = 0;
    std::uint64_t begin = 0;
    std::uint64_t end = 0;
};

// ==========================================
// Structure: Accumulate task-level row counters
// Method: Separate accepted and policy-skipped malformed records for MPI reduction.
// ==========================================
struct TaskStats {
    std::uint64_t accepted_rows = 0;
    std::uint64_t malformed_rows = 0;
};

// ==========================================
// Structure: Identify one integral-degree sky tile
// Method: Store lower RA and Dec boundaries and order them deterministically.
// ==========================================
struct TileKey {
    int ra_lower = 0;
    int dec_lower = 0;

    // ==========================================
    // Function: Order sky tile keys
    // Method: Compare RA first and declination second for stable map traversal.
    // ==========================================
    bool operator<(const TileKey& other) const {
        if (ra_lower != other.ra_lower) {
            return ra_lower < other.ra_lower;
        }
        return dec_lower < other.dec_lower;
    }
};

// ==========================================
// Function: Trim ASCII whitespace from both ends of a string
// Method: Find the first and last non-space bytes without locale-dependent mutation.
// ==========================================
std::string trim(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return std::string();
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

// ==========================================
// Function: Convert an ASCII identifier to lowercase
// Method: Lowercase each unsigned byte for case-insensitive schema matching.
// ==========================================
std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        if (character >= 'A' && character <= 'Z') {
            return static_cast<char>(character - 'A' + 'a');
        }
        return static_cast<char>(character);
    });
    return value;
}

// ==========================================
// Function: Remove an optional UTF-8 byte-order mark from the first field
// Method: Erase the three BOM bytes after delimiter parsing so headers and data share behavior.
// ==========================================
void stripUtf8Bom(std::vector<std::string>& tokens) {
    if (!tokens.empty() && tokens.front().size() >= 3
        && static_cast<unsigned char>(tokens.front()[0]) == 0xef
        && static_cast<unsigned char>(tokens.front()[1]) == 0xbb
        && static_cast<unsigned char>(tokens.front()[2]) == 0xbf) {
        tokens.front().erase(0, 3);
    }
}

// ==========================================
// Function: Normalize a path for deterministic discovery and safety checks
// Method: Resolve existing components while permitting a not-yet-created output tail.
// ==========================================
fs::path normalizedAbsolute(const fs::path& path) {
    return fs::weakly_canonical(fs::absolute(path));
}

// ==========================================
// Function: Test whether one normalized path is equal to or below another
// Method: Compare complete path components instead of raw string prefixes.
// ==========================================
bool pathIsWithin(const fs::path& candidate, const fs::path& parent) {
    auto candidate_iterator = candidate.begin();
    auto parent_iterator = parent.begin();
    for (; parent_iterator != parent.end(); ++parent_iterator, ++candidate_iterator) {
        if (candidate_iterator == candidate.end() || *candidate_iterator != *parent_iterator) {
            return false;
        }
    }
    return true;
}

// ==========================================
// Function: Test a basename against repeatable substring filters
// Method: Select every file when no filters exist, otherwise apply case-sensitive OR semantics.
// ==========================================
bool matchesFilenameTokens(const std::string& basename,
                           const std::vector<std::string>& tokens) {
    if (tokens.empty()) {
        return true;
    }
    return std::any_of(tokens.begin(), tokens.end(), [&](const std::string& token) {
        return basename.find(token) != std::string::npos;
    });
}

// ==========================================
// Function: Recognize one pipeline-generated tile basename
// Method: Match the exact des_y6 one-degree naming grammar used by catalog lookup.
// ==========================================
bool isGeneratedTileName(const std::string& basename) {
    static const std::regex pattern(
        R"(^des_y6_RA_[0-9]{3}_[0-9]{3}_Dec_[pm][0-9]{2}_[pm][0-9]{2}\.dat$)");
    return std::regex_match(basename, pattern);
}

// ==========================================
// Function: Broadcast one dynamically sized string
// Method: Send an integer byte count followed by the contiguous character payload.
// ==========================================
void broadcastString(std::string& value, int root_rank, MPI_Comm communicator) {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int length = 0;
    if (rank == root_rank) {
        length = value.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
                     ? -1
                     : static_cast<int>(value.size());
    }
    MPI_Bcast(&length, 1, MPI_INT, root_rank, communicator);
    if (length < 0) {
        throw std::runtime_error("MPI string exceeds int length range");
    }
    if (rank != root_rank) {
        value.resize(static_cast<std::size_t>(length));
    }
    if (length > 0) {
        MPI_Bcast(value.data(), length, MPI_CHAR, root_rank, communicator);
    }
}

// ==========================================
// Function: Report the first processing error from every failed MPI rank
// Method: Broadcast each rank's optional message in rank order and print only on rank zero.
// ==========================================
void reportRankErrors(const std::string& local_error, MPI_Comm communicator) {
    int rank = 0;
    int world_size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &world_size);
    for (int source_rank = 0; source_rank < world_size; ++source_rank) {
        std::string message = rank == source_rank ? local_error : std::string();
        broadcastString(message, source_rank, communicator);
        if (rank == 0 && !message.empty()) {
            std::cerr << "process_extcat rank " << source_rank << ": " << message << '\n';
        }
    }
}

// ==========================================
// Function: Resolve the effective delimiter for one representative line
// Method: Honor explicit configuration or infer comma, tab, then generic whitespace.
// ==========================================
ProcessExtcat::Delimiter resolveDelimiter(ProcessExtcat::Delimiter configured,
                                          const std::string& line) {
    if (configured != ProcessExtcat::Delimiter::Auto) {
        return configured;
    }
    if (line.find(',') != std::string::npos) {
        return ProcessExtcat::Delimiter::Comma;
    }
    if (line.find('\t') != std::string::npos) {
        return ProcessExtcat::Delimiter::Tab;
    }
    return ProcessExtcat::Delimiter::Whitespace;
}

// ==========================================
// Function: Split one delimited catalog record
// Method: Use stream tokenization for whitespace and a quote-aware finite-state parser
//         for comma or tab input.
// ==========================================
std::vector<std::string> splitLine(const std::string& line,
                                   ProcessExtcat::Delimiter delimiter) {
    std::vector<std::string> tokens;
    if (delimiter == ProcessExtcat::Delimiter::Whitespace
        || delimiter == ProcessExtcat::Delimiter::Auto) {
        std::istringstream input(line);
        std::string token;
        while (input >> token) {
            tokens.push_back(token);
        }
        stripUtf8Bom(tokens);
        return tokens;
    }

    const char separator = delimiter == ProcessExtcat::Delimiter::Comma ? ',' : '\t';
    std::string token;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (character == '"') {
            if (quoted && index + 1 < line.size() && line[index + 1] == '"') {
                token.push_back('"');
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (character == separator && !quoted) {
            tokens.push_back(trim(token));
            token.clear();
        } else {
            token.push_back(character);
        }
    }
    if (quoted) {
        throw std::runtime_error("unterminated quoted field");
    }
    tokens.push_back(trim(token));
    stripUtf8Bom(tokens);
    return tokens;
}

// ==========================================
// Function: Parse one finite floating-point catalog token
// Method: Use strtod with complete-consumption, range, and finiteness checks.
// ==========================================
bool parseFiniteDouble(const std::string& token, double& value) {
    if (token.empty()) {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const char* begin = token.c_str();
    value = std::strtod(begin, &end);
    return end != begin && end != nullptr && end == begin + token.size() && errno != ERANGE
           && std::isfinite(value);
}

// ==========================================
// Function: Locate named celestial-coordinate columns
// Method: Match unique ra and dec field names case-insensitively without constraining any
//         other input or output column.
// ==========================================
bool findNamedCoordinateColumns(const std::vector<std::string>& tokens,
                                std::size_t& ra_column,
                                std::size_t& dec_column) {
    std::map<std::string, std::size_t> positions;
    std::set<std::string> duplicate_names;
    for (std::size_t index = 0; index < tokens.size(); ++index) {
        const std::string name = lowercase(trim(tokens[index]));
        if (!name.empty() && !positions.emplace(name, index).second) {
            duplicate_names.insert(name);
        }
    }
    const auto ra_position = positions.find("ra");
    const auto dec_position = positions.find("dec");
    if (ra_position == positions.end() || dec_position == positions.end()
        || duplicate_names.count("ra") != 0 || duplicate_names.count("dec") != 0) {
        return false;
    }
    ra_column = ra_position->second;
    dec_column = dec_position->second;
    return true;
}

// ==========================================
// Function: Build one output projection and header
// Method: Preserve every input column when projection is disabled, otherwise retain the
//         configured order; use input names when available and stable generic names otherwise.
// ==========================================
void configureMetadataSchema(const std::vector<std::string>& tokens,
                             bool tokens_are_header,
                             const ProcessExtcat::Config& config,
                             FileMetadata& metadata) {
    if (tokens.empty()) {
        throw std::runtime_error("catalog schema contains no columns");
    }
    metadata.input_column_count = tokens.size();
    metadata.ra_column = config.ra_column;
    metadata.dec_column = config.dec_column;
    if (tokens_are_header && !config.use_explicit_coordinate_columns) {
        std::size_t named_ra = 0;
        std::size_t named_dec = 0;
        if (findNamedCoordinateColumns(tokens, named_ra, named_dec)) {
            metadata.ra_column = named_ra;
            metadata.dec_column = named_dec;
        }
    }
    if (metadata.ra_column >= tokens.size() || metadata.dec_column >= tokens.size()) {
        throw std::runtime_error(
            "catalog has " + std::to_string(tokens.size())
            + " fields but RA/Dec columns are "
            + std::to_string(metadata.ra_column + 1) + "/"
            + std::to_string(metadata.dec_column + 1));
    }

    metadata.columns.clear();
    if (config.use_explicit_columns) {
        metadata.columns = config.input_columns;
    } else {
        metadata.columns.reserve(tokens.size());
        for (std::size_t index = 0; index < tokens.size(); ++index) {
            metadata.columns.push_back(index);
        }
    }

    metadata.output_header.clear();
    metadata.output_header.reserve(metadata.columns.size());
    for (const std::size_t column : metadata.columns) {
        if (column >= tokens.size()) {
            throw std::runtime_error(
                "catalog has " + std::to_string(tokens.size())
                + " fields but the projection requests field "
                + std::to_string(column + 1));
        }
        if (tokens_are_header) {
            std::string name = trim(tokens[column]);
            std::replace_if(name.begin(), name.end(), [](unsigned char character) {
                return character == ' ' || character == '\t';
            }, '_');
            metadata.output_header.push_back(
                name.empty() ? "column_" + std::to_string(column + 1) : name);
        } else {
            metadata.output_header.push_back("column_" + std::to_string(column + 1));
        }
    }
}

// ==========================================
// Function: Project and validate one raw record
// Method: Read finite RA/Dec directly from raw columns, enforce stable width in pass-through
//         mode, and copy the configured output fields without requiring numeric payloads.
// ==========================================
bool projectRow(const std::vector<std::string>& tokens,
                const FileMetadata& metadata,
                bool require_exact_width,
                std::vector<std::string>& projected,
                double& ra,
                double& dec,
                std::string& error) {
    if (require_exact_width && tokens.size() != metadata.input_column_count) {
        error = "record has " + std::to_string(tokens.size())
                + " fields but the pass-through schema requires "
                + std::to_string(metadata.input_column_count);
        return false;
    }
    if (metadata.ra_column >= tokens.size() || metadata.dec_column >= tokens.size()) {
        error = "record has " + std::to_string(tokens.size())
                + " fields but RA/Dec columns are "
                + std::to_string(metadata.ra_column + 1) + "/"
                + std::to_string(metadata.dec_column + 1);
        return false;
    }
    if (!parseFiniteDouble(trim(tokens[metadata.ra_column]), ra)) {
        error = "RA field " + std::to_string(metadata.ra_column + 1)
                + " is not a finite number: " + trim(tokens[metadata.ra_column]);
        return false;
    }
    if (!parseFiniteDouble(trim(tokens[metadata.dec_column]), dec)) {
        error = "Dec field " + std::to_string(metadata.dec_column + 1)
                + " is not a finite number: " + trim(tokens[metadata.dec_column]);
        return false;
    }

    projected.clear();
    projected.reserve(metadata.columns.size());
    for (const std::size_t column : metadata.columns) {
        if (column >= tokens.size()) {
            error = "record has " + std::to_string(tokens.size())
                    + " fields but the projection requests field "
                    + std::to_string(column + 1);
            return false;
        }
        projected.push_back(trim(tokens[column]));
    }
    return true;
}

// ==========================================
// Function: Join one projected output record
// Method: Emit the selected tokens in exact configured order with whitespace separation.
// ==========================================
std::string joinProjectedRow(const std::vector<std::string>& projected) {
    std::ostringstream output;
    for (std::size_t index = 0; index < projected.size(); ++index) {
        if (index > 0) {
            output << ' ';
        }
        output << projected[index];
    }
    return output.str();
}

// ==========================================
// Function: Map celestial coordinates to one integral-degree tile
// Method: Accept RA in [0,360], wrap the exact 360 boundary to zero, and include
//         the exact north pole in the final [89,90] declination tile.
// ==========================================
bool coordinateTile(double ra, double dec, TileKey& tile, std::string& error) {
    if (ra < 0.0 || ra > 360.0) {
        error = "right ascension is outside [0, 360] degrees: " + std::to_string(ra);
        return false;
    }
    if (dec < -90.0 || dec > 90.0) {
        error = "declination is outside [-90, 90] degrees: " + std::to_string(dec);
        return false;
    }
    const double normalized_ra = ra == 360.0 ? 0.0 : ra;
    const double bounded_dec = dec == 90.0
                                   ? std::nextafter(90.0, -std::numeric_limits<double>::infinity())
                                   : dec;
    tile.ra_lower = static_cast<int>(std::floor(normalized_ra));
    tile.dec_lower = static_cast<int>(std::floor(bounded_dec));
    return true;
}

// ==========================================
// Function: Format one signed declination boundary
// Method: Prefix non-negative integers with p, negative integers with m, and zero-pad two digits.
// ==========================================
std::string formatDeclination(int value) {
    std::ostringstream output;
    output << (value >= 0 ? 'p' : 'm') << std::setw(2) << std::setfill('0')
           << std::abs(value);
    return output.str();
}

// ==========================================
// Function: Build the pipeline filename for one one-degree tile
// Method: Reproduce the des_y6 RA padding and pXX/mXX declination boundary convention.
// ==========================================
std::string tileFilename(const TileKey& tile) {
    std::ostringstream output;
    output << kTilePrefix << std::setw(3) << std::setfill('0') << tile.ra_lower << '_'
           << std::setw(3) << std::setfill('0') << tile.ra_lower + 1 << "_Dec_"
           << formatDeclination(tile.dec_lower) << '_'
           << formatDeclination(tile.dec_lower + 1) << kTileSuffix;
    return output.str();
}

// ==========================================
// Function: Construct one projected commented output header
// Method: Join the validated effective field names with one leading hash marker.
// ==========================================
std::string outputHeader(const std::vector<std::string>& columns) {
    std::ostringstream output;
    output << "# ";
    for (std::size_t index = 0; index < columns.size(); ++index) {
        if (index > 0) {
            output << ' ';
        }
        output << columns[index];
    }
    return output.str();
}

// ==========================================
// Function: Discover all selected raw catalog files
// Method: Traverse one immutable input root, select regular files by basename substring OR,
//         canonicalize and deduplicate paths, then return deterministic lexical order.
// ==========================================
std::vector<fs::path> discoverInputFiles(const ProcessExtcat::Config& config) {
    std::set<fs::path> unique_paths;
    std::error_code iterator_error;

    if (config.recursive) {
        fs::recursive_directory_iterator iterator(
            config.input_directory, fs::directory_options::none, iterator_error);
        const fs::recursive_directory_iterator end;
        if (iterator_error) {
            throw std::runtime_error("cannot scan input directory: " + iterator_error.message());
        }
        while (iterator != end) {
            const fs::directory_entry entry = *iterator;
            std::error_code type_error;
            const bool regular = entry.is_regular_file(type_error);
            if (type_error) {
                throw std::runtime_error("cannot inspect input entry " + entry.path().string()
                                         + ": " + type_error.message());
            }
            if (regular
                && matchesFilenameTokens(entry.path().filename().string(), config.filename_tokens)) {
                unique_paths.insert(normalizedAbsolute(entry.path()));
            }
            iterator.increment(iterator_error);
            if (iterator_error) {
                throw std::runtime_error("input traversal failed: " + iterator_error.message());
            }
        }
    } else {
        fs::directory_iterator iterator(config.input_directory, iterator_error);
        const fs::directory_iterator end;
        if (iterator_error) {
            throw std::runtime_error("cannot scan input directory: " + iterator_error.message());
        }
        while (iterator != end) {
            const fs::directory_entry entry = *iterator;
            std::error_code type_error;
            const bool regular = entry.is_regular_file(type_error);
            if (type_error) {
                throw std::runtime_error("cannot inspect input entry " + entry.path().string()
                                         + ": " + type_error.message());
            }
            if (regular
                && matchesFilenameTokens(entry.path().filename().string(), config.filename_tokens)) {
                unique_paths.insert(normalizedAbsolute(entry.path()));
            }
            iterator.increment(iterator_error);
            if (iterator_error) {
                throw std::runtime_error("input traversal failed: " + iterator_error.message());
            }
        }
    }

    if (unique_paths.empty()) {
        throw std::runtime_error("no regular input files matched the configured basename filters");
    }
    return std::vector<fs::path>(unique_paths.begin(), unique_paths.end());
}

// ==========================================
// Function: Inspect one raw catalog schema and locate its first data record
// Method: Detect delimiters and general headers, resolve coordinate columns separately from
//         output projection, and validate the first celestial-coordinate row.
// ==========================================
FileMetadata inspectFile(const fs::path& path, const ProcessExtcat::Config& config) {
    std::ifstream input(path, std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("cannot open input catalog: " + path.string());
    }

    FileMetadata metadata;
    metadata.path = path;
    metadata.file_size = static_cast<std::uint64_t>(fs::file_size(path));
    bool header_seen = false;
    bool schema_configured = false;
    ProcessExtcat::Delimiter resolved_delimiter = config.delimiter;
    std::string line;
    while (true) {
        const std::streampos line_position = input.tellg();
        if (!std::getline(input, line)) {
            break;
        }
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const std::string stripped = trim(line);
        if (stripped.empty()) {
            continue;
        }

        const bool commented = stripped.front() == '#';
        const std::string candidate = commented ? trim(stripped.substr(1)) : stripped;
        if (candidate.empty()) {
            continue;
        }

        const ProcessExtcat::Delimiter candidate_delimiter =
            resolveDelimiter(config.delimiter, candidate);
        std::vector<std::string> tokens;
        try {
            tokens = splitLine(candidate, candidate_delimiter);
        } catch (const std::exception& exception) {
            throw std::runtime_error("cannot parse leading record in " + path.string()
                                     + ": " + exception.what());
        }

        std::size_t named_ra = 0;
        std::size_t named_dec = 0;
        const bool named_coordinates = findNamedCoordinateColumns(
            tokens, named_ra, named_dec);
        const bool header_allowed = config.header_mode != ProcessExtcat::HeaderMode::Absent;
        const bool header_candidate = header_allowed && !header_seen
                                      && (named_coordinates
                                          || (commented
                                              && config.header_mode
                                                     == ProcessExtcat::HeaderMode::Present));
        if (header_candidate) {
            configureMetadataSchema(tokens, true, config, metadata);
            header_seen = true;
            schema_configured = true;
            resolved_delimiter = candidate_delimiter;
            continue;
        }
        if (commented) {
            continue;
        }

        if (!schema_configured) {
            configureMetadataSchema(tokens, false, config, metadata);
            schema_configured = true;
        }
        std::vector<std::string> projected;
        double ra = 0.0;
        double dec = 0.0;
        std::string row_error;
        const bool coordinate_row = projectRow(
            tokens, metadata, !config.use_explicit_columns, projected, ra, dec, row_error);

        if (!coordinate_row && header_allowed && !header_seen) {
            configureMetadataSchema(tokens, true, config, metadata);
            header_seen = true;
            schema_configured = true;
            resolved_delimiter = candidate_delimiter;
            continue;
        }
        if (!coordinate_row) {
            throw std::runtime_error("cannot identify the first data row in " + path.string()
                                     + ": " + row_error
                                     + "; set coordinate columns for a nonstandard schema");
        }
        if (config.header_mode == ProcessExtcat::HeaderMode::Present && !header_seen) {
            throw std::runtime_error("--header present was requested but no header was found in "
                                     + path.string());
        }
        if (line_position < 0) {
            throw std::runtime_error("cannot determine data offset in " + path.string());
        }
        metadata.data_offset = static_cast<std::uint64_t>(line_position);
        metadata.delimiter = resolveDelimiter(resolved_delimiter, candidate);
        return metadata;
    }

    throw std::runtime_error(
        "input catalog contains no row with finite RA/Dec: " + path.string());
}

// ==========================================
// Function: Validate merged output schema compatibility
// Method: Require every input file to produce the same projected header before any MPI task
//         can concatenate their rows into shared tile files.
// ==========================================
void validateCompatibleOutputSchemas(const std::vector<FileMetadata>& metadata) {
    if (metadata.empty()) {
        throw std::runtime_error("no inspected input metadata is available");
    }
    const std::vector<std::string>& expected = metadata.front().output_header;
    for (std::size_t index = 1; index < metadata.size(); ++index) {
        if (metadata[index].output_header != expected) {
            throw std::runtime_error(
                "projected output schema differs between "
                + metadata.front().path.string() + " and "
                + metadata[index].path.string());
        }
    }
}

// ==========================================
// Function: Broadcast inspected file metadata to every MPI rank
// Method: Length-prefix paths and variable projections, then send offsets, widths,
//         delimiter enums, and raw coordinate indices.
// ==========================================
void broadcastMetadata(std::vector<FileMetadata>& metadata,
                       int root_rank,
                       MPI_Comm communicator) {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int count = 0;
    if (rank == root_rank) {
        count = metadata.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
                    ? -1
                    : static_cast<int>(metadata.size());
    }
    MPI_Bcast(&count, 1, MPI_INT, root_rank, communicator);
    if (count < 0) {
        throw std::runtime_error("input file count exceeds MPI int range");
    }
    if (rank != root_rank) {
        metadata.resize(static_cast<std::size_t>(count));
    }

    for (int item = 0; item < count; ++item) {
        FileMetadata& file = metadata[static_cast<std::size_t>(item)];
        std::string path_text = rank == root_rank ? file.path.string() : std::string();
        broadcastString(path_text, root_rank, communicator);
        if (rank != root_rank) {
            file.path = path_text;
        }

        std::uint64_t sizes[5] = {
            file.data_offset,
            file.file_size,
            static_cast<std::uint64_t>(file.input_column_count),
            static_cast<std::uint64_t>(file.ra_column),
            static_cast<std::uint64_t>(file.dec_column),
        };
        MPI_Bcast(sizes, 5, MPI_UINT64_T, root_rank, communicator);
        if (rank != root_rank) {
            file.data_offset = sizes[0];
            file.file_size = sizes[1];
            for (int index = 2; index < 5; ++index) {
                if (sizes[index]
                    > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
                    throw std::runtime_error(
                        "received catalog metadata outside size_t range");
                }
            }
            file.input_column_count = static_cast<std::size_t>(sizes[2]);
            file.ra_column = static_cast<std::size_t>(sizes[3]);
            file.dec_column = static_cast<std::size_t>(sizes[4]);
        }

        int delimiter = rank == root_rank ? static_cast<int>(file.delimiter) : 0;
        MPI_Bcast(&delimiter, 1, MPI_INT, root_rank, communicator);
        if (rank != root_rank) {
            file.delimiter = static_cast<ProcessExtcat::Delimiter>(delimiter);
        }

        int column_count = 0;
        if (rank == root_rank) {
            column_count = file.columns.size()
                               > static_cast<std::size_t>(std::numeric_limits<int>::max())
                               ? -1
                               : static_cast<int>(file.columns.size());
        }
        MPI_Bcast(&column_count, 1, MPI_INT, root_rank, communicator);
        if (column_count <= 0) {
            throw std::runtime_error("projected column count is outside MPI int range");
        }
        std::vector<std::uint64_t> columns(static_cast<std::size_t>(column_count), 0);
        if (rank == root_rank) {
            for (std::size_t index = 0; index < file.columns.size(); ++index) {
                columns[index] = static_cast<std::uint64_t>(file.columns[index]);
            }
        }
        MPI_Bcast(columns.data(), static_cast<int>(columns.size()), MPI_UINT64_T,
                  root_rank, communicator);
        if (rank != root_rank) {
            file.columns.resize(columns.size());
            for (std::size_t index = 0; index < columns.size(); ++index) {
                if (columns[index]
                    > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
                    throw std::runtime_error("received a column index outside size_t range");
                }
                file.columns[index] = static_cast<std::size_t>(columns[index]);
            }
        }
    }
}

// ==========================================
// Function: Build deterministic byte-range tasks for all input catalogs
// Method: Cap chunks by configuration and shrink the effective size enough to expose
//         at least one approximately balanced task per rank when data volume permits.
// ==========================================
std::vector<Task> buildTasks(const std::vector<FileMetadata>& metadata,
                             std::uint64_t configured_chunk_bytes,
                             int world_size) {
    std::uint64_t total_bytes = 0;
    for (const FileMetadata& file : metadata) {
        if (file.file_size <= file.data_offset) {
            throw std::runtime_error("catalog has no bytes after its first data row: "
                                     + file.path.string());
        }
        const std::uint64_t data_bytes = file.file_size - file.data_offset;
        if (data_bytes > std::numeric_limits<std::uint64_t>::max() - total_bytes) {
            throw std::runtime_error("total catalog byte count exceeds uint64 range");
        }
        total_bytes += data_bytes;
    }
    if (total_bytes == 0 || world_size <= 0) {
        throw std::runtime_error("cannot create MPI tasks for empty data or zero ranks");
    }

    const std::uint64_t balanced_bytes =
        (total_bytes + static_cast<std::uint64_t>(world_size) - 1)
        / static_cast<std::uint64_t>(world_size);
    const std::uint64_t chunk_bytes = std::max<std::uint64_t>(
        1, std::min(configured_chunk_bytes, balanced_bytes));

    std::vector<Task> tasks;
    for (std::size_t file_index = 0; file_index < metadata.size(); ++file_index) {
        const FileMetadata& file = metadata[file_index];
        std::uint64_t begin = file.data_offset;
        while (begin < file.file_size) {
            const std::uint64_t remaining = file.file_size - begin;
            const std::uint64_t length = std::min(chunk_bytes, remaining);
            Task task;
            task.id = static_cast<std::uint64_t>(tasks.size());
            task.file_index = file_index;
            task.begin = begin;
            task.end = begin + length;
            tasks.push_back(task);
            begin += length;
        }
    }
    return tasks;
}

// ==========================================
// Function: Build one zero-padded per-task staging directory
// Method: Encode the deterministic task identifier so lexical and processing order agree.
// ==========================================
fs::path taskDirectory(const fs::path& staging_directory, std::uint64_t task_id) {
    std::ostringstream name;
    name << "task_" << std::setw(12) << std::setfill('0') << task_id;
    return staging_directory / "shards" / name.str();
}

// ==========================================
// Function: Process one newline-aligned catalog byte range
// Method: Seek to the nominal range, discard only an initial partial record, project rows in
//         configured order, group by raw RA/Dec, and write collision-free shard files.
// ==========================================
TaskStats processTask(const Task& task,
                      const FileMetadata& metadata,
                      const ProcessExtcat::Config& config,
                      const fs::path& staging_directory) {
    std::ifstream input(metadata.path, std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("cannot open input catalog: " + metadata.path.string());
    }

    if (task.begin > metadata.data_offset) {
        input.seekg(static_cast<std::streamoff>(task.begin - 1));
        char previous = '\0';
        input.get(previous);
        if (!input.good()) {
            throw std::runtime_error("cannot align byte range in " + metadata.path.string());
        }
        if (previous != '\n') {
            std::string partial;
            std::getline(input, partial);
        }
    } else {
        input.seekg(static_cast<std::streamoff>(task.begin));
    }
    if (!input.good()) {
        throw std::runtime_error("cannot seek input catalog: " + metadata.path.string());
    }

    TaskStats stats;
    std::map<TileKey, std::vector<std::string>> grouped_rows;
    std::string line;
    while (true) {
        const std::streampos line_position = input.tellg();
        if (line_position < 0
            || static_cast<std::uint64_t>(line_position) >= task.end
            || !std::getline(input, line)) {
            break;
        }
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped.front() == '#') {
            continue;
        }

        std::string row_error;
        try {
            const std::vector<std::string> tokens = splitLine(stripped, metadata.delimiter);
            std::vector<std::string> projected;
            double ra = 0.0;
            double dec = 0.0;
            if (!projectRow(tokens, metadata, !config.use_explicit_columns,
                            projected, ra, dec, row_error)) {
                throw std::runtime_error(row_error);
            }
            TileKey tile;
            if (!coordinateTile(ra, dec, tile, row_error)) {
                throw std::runtime_error(row_error);
            }
            grouped_rows[tile].push_back(joinProjectedRow(projected));
            ++stats.accepted_rows;
        } catch (const std::exception& exception) {
            ++stats.malformed_rows;
            if (config.malformed_policy == ProcessExtcat::MalformedPolicy::Fail) {
                throw std::runtime_error(
                    metadata.path.string() + " at byte "
                    + std::to_string(static_cast<std::uint64_t>(line_position))
                    + ": " + exception.what());
            }
        }
    }

    if (!grouped_rows.empty()) {
        const fs::path directory = taskDirectory(staging_directory, task.id);
        fs::create_directories(directory);
        for (const auto& [tile, rows] : grouped_rows) {
            const fs::path shard_path = directory / (tileFilename(tile) + ".part");
            std::ofstream output(shard_path, std::ios::binary | std::ios::trunc);
            if (!output.is_open()) {
                throw std::runtime_error("cannot create shard file: " + shard_path.string());
            }
            for (const std::string& row : rows) {
                output << row << '\n';
            }
            if (!output.good()) {
                throw std::runtime_error("cannot write shard file: " + shard_path.string());
            }
            output.close();
            if (output.fail()) {
                throw std::runtime_error("cannot finalize shard file: " + shard_path.string());
            }
        }
    }
    return stats;
}

// ==========================================
// Function: List existing generated tiles in deterministic order
// Method: Inspect only immediate regular files whose basenames match the fixed tile grammar.
// ==========================================
std::vector<fs::path> listExistingTiles(const fs::path& output_directory) {
    std::vector<fs::path> paths;
    if (!fs::exists(output_directory)) {
        return paths;
    }
    for (const fs::directory_entry& entry : fs::directory_iterator(output_directory)) {
        if (entry.is_regular_file()
            && isGeneratedTileName(entry.path().filename().string())) {
            paths.push_back(entry.path());
        }
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

// ==========================================
// Function: Create one unique staging directory below the output root
// Method: Combine epoch microseconds with the rank-zero process identifier.
// ==========================================
fs::path createStagingDirectory(const fs::path& output_directory) {
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    const auto microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(now).count();
    const fs::path staging = output_directory
                             / (".process_extcat_staging_"
                                + std::to_string(microseconds) + "_"
                                + std::to_string(getpid()));
    if (!fs::create_directories(staging / "shards")) {
        throw std::runtime_error("cannot create unique staging directory: " + staging.string());
    }
    return staging;
}

// ==========================================
// Function: Publish a complete staged tile set with rollback
// Method: Move old generated tiles into a private backup, atomically rename new files on the
//         same filesystem, and restore the old set if any publication step fails.
// ==========================================
void publishTiles(const fs::path& final_directory,
                  const fs::path& output_directory,
                  const fs::path& staging_directory,
                  ProcessExtcat::ExistingPolicy policy) {
    const std::vector<fs::path> existing = listExistingTiles(output_directory);
    if (policy == ProcessExtcat::ExistingPolicy::Fail && !existing.empty()) {
        throw std::runtime_error("output directory already contains generated tiles; "
                                 "use --existing overwrite: " + output_directory.string());
    }

    const fs::path backup_directory = staging_directory / "backup";
    std::vector<fs::path> backed_up;
    std::vector<fs::path> published;
    try {
        if (policy == ProcessExtcat::ExistingPolicy::Overwrite && !existing.empty()) {
            fs::create_directories(backup_directory);
            for (const fs::path& old_path : existing) {
                const fs::path backup_path = backup_directory / old_path.filename();
                fs::rename(old_path, backup_path);
                backed_up.push_back(backup_path);
            }
        }

        std::vector<fs::path> new_tiles;
        for (const fs::directory_entry& entry : fs::directory_iterator(final_directory)) {
            if (entry.is_regular_file()) {
                new_tiles.push_back(entry.path());
            }
        }
        std::sort(new_tiles.begin(), new_tiles.end());
        for (const fs::path& new_path : new_tiles) {
            const fs::path destination = output_directory / new_path.filename();
            if (policy == ProcessExtcat::ExistingPolicy::Fail && fs::exists(destination)) {
                throw std::runtime_error("output tile appeared during processing: "
                                         + destination.string());
            }
            fs::rename(new_path, destination);
            published.push_back(destination);
        }
    } catch (const std::exception& publication_exception) {
        std::vector<std::string> rollback_errors;
        for (const fs::path& path : published) {
            std::error_code remove_error;
            fs::remove(path, remove_error);
            if (remove_error) {
                rollback_errors.push_back(
                    "cannot remove partially published tile " + path.string()
                    + ": " + remove_error.message());
            }
        }
        for (const fs::path& backup_path : backed_up) {
            std::error_code restore_error;
            fs::rename(backup_path, output_directory / backup_path.filename(), restore_error);
            if (restore_error) {
                rollback_errors.push_back(
                    "cannot restore backup tile " + backup_path.string()
                    + ": " + restore_error.message());
            }
        }
        if (!rollback_errors.empty()) {
            std::ostringstream message;
            message << "publication failed: " << publication_exception.what()
                    << "; rollback was incomplete";
            for (const std::string& rollback_error : rollback_errors) {
                message << "; " << rollback_error;
            }
            throw std::runtime_error(message.str());
        }
        throw;
    }
}

// ==========================================
// Function: Merge all deterministic task shards into final projected tiles
// Method: Traverse task directories in task order, concatenate shards by original input byte
//         order, prepend the validated shared header, and publish transactionally.
// ==========================================
std::size_t mergeShards(const std::vector<Task>& tasks,
                        const fs::path& staging_directory,
                        const std::vector<std::string>& projected_header,
                        const ProcessExtcat::Config& config) {
    std::map<std::string, std::vector<fs::path>> shards_by_tile;
    for (const Task& task : tasks) {
        const fs::path directory = taskDirectory(staging_directory, task.id);
        if (!fs::exists(directory)) {
            continue;
        }
        std::vector<fs::path> task_shards;
        for (const fs::directory_entry& entry : fs::directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                task_shards.push_back(entry.path());
            }
        }
        std::sort(task_shards.begin(), task_shards.end());
        for (const fs::path& shard : task_shards) {
            const std::string basename = shard.filename().string();
            constexpr const char* part_suffix = ".part";
            const std::size_t suffix_length = std::char_traits<char>::length(part_suffix);
            if (basename.size() <= suffix_length
                || basename.compare(basename.size() - suffix_length, suffix_length,
                                    part_suffix) != 0) {
                throw std::runtime_error("unexpected shard basename: " + basename);
            }
            const std::string tile_name = basename.substr(0, basename.size() - suffix_length);
            if (!isGeneratedTileName(tile_name)) {
                throw std::runtime_error("invalid staged tile basename: " + tile_name);
            }
            shards_by_tile[tile_name].push_back(shard);
        }
    }
    if (shards_by_tile.empty()) {
        throw std::runtime_error("no valid catalog rows were available to publish");
    }

    const fs::path final_directory = staging_directory / "final";
    fs::create_directories(final_directory);
    for (const auto& [tile_name, shards] : shards_by_tile) {
        const fs::path final_path = final_directory / tile_name;
        std::ofstream output(final_path, std::ios::binary | std::ios::trunc);
        if (!output.is_open()) {
            throw std::runtime_error("cannot create staged tile: " + final_path.string());
        }
        output << outputHeader(projected_header) << '\n';
        for (const fs::path& shard : shards) {
            std::ifstream input(shard, std::ios::binary);
            if (!input.is_open()) {
                throw std::runtime_error("cannot read shard during merge: " + shard.string());
            }
            output << input.rdbuf();
            if (input.bad()) {
                throw std::runtime_error("cannot read complete shard during merge: "
                                         + shard.string());
            }
            if (!output.good()) {
                throw std::runtime_error("cannot append shard to tile: " + final_path.string());
            }
        }
        output.close();
        if (output.fail()) {
            throw std::runtime_error("cannot finalize staged tile: " + final_path.string());
        }
    }

    publishTiles(final_directory, config.output_directory, staging_directory,
                 config.existing_policy);
    return shards_by_tile.size();
}

}  // namespace

namespace ProcessExtcat {

// ==========================================
// Function: Return the legacy DES Y6 GOLD reference schema
// Method: Retain the immutable 18-field API reference without constraining generated width.
// ==========================================
const std::array<std::string, kCanonicalColumnCount>& canonicalColumnNames() {
    static const std::array<std::string, kCanonicalColumnCount> names = {
        "flags_footprint",
        "flags_foreground",
        "flags_gold",
        "ext_mash",
        "ra",
        "dec",
        "bdf_mag_g",
        "bdf_mag_err_g",
        "bdf_mag_r",
        "bdf_mag_err_r",
        "bdf_mag_i",
        "bdf_mag_err_i",
        "bdf_mag_z",
        "bdf_mag_err_z",
        "bdf_mag_y",
        "bdf_mag_err_y",
        "dnf_z",
        "dnf_zsigma",
    };
    return names;
}

// ==========================================
// Function: Normalize and validate external-catalog tiler configuration
// Method: Canonicalize paths and reject unsafe overlap, empty filters, invalid mappings,
//         and unsupported task sizes before collective processing begins.
// ==========================================
void normalizeAndValidateConfig(Config& config) {
    if (config.input_directory.empty()) {
        throw std::runtime_error("input directory must be provided");
    }
    if (config.output_directory.empty()) {
        throw std::runtime_error("output directory must be provided");
    }
    config.input_directory = normalizedAbsolute(config.input_directory);
    config.output_directory = normalizedAbsolute(config.output_directory);

    if (!fs::exists(config.input_directory) || !fs::is_directory(config.input_directory)) {
        throw std::runtime_error("input path is not a directory: "
                                 + config.input_directory.string());
    }
    if (config.input_directory == config.output_directory
        || pathIsWithin(config.output_directory, config.input_directory)) {
        throw std::runtime_error("output directory must not equal or be below the input directory");
    }
    if (fs::exists(config.output_directory) && !fs::is_directory(config.output_directory)) {
        throw std::runtime_error("output path exists but is not a directory: "
                                 + config.output_directory.string());
    }
    if (config.chunk_bytes == 0) {
        throw std::runtime_error("chunk byte size must be positive");
    }
    for (const std::string& token : config.filename_tokens) {
        if (token.empty()) {
            throw std::runtime_error("filename filter tokens must not be empty");
        }
    }
    if (config.use_explicit_columns && config.input_columns.empty()) {
        throw std::runtime_error("explicit column list must not be empty");
    }
    if (config.ra_column == config.dec_column) {
        throw std::runtime_error("RA and Dec input columns must be distinct");
    }
}

}  // namespace ProcessExtcat

namespace {

// ==========================================
// Function: Convert configured delimiter text
// Method: Map the documented ProcessConfig strings to the reusable parser enum.
// ==========================================
ProcessExtcat::Delimiter configuredDelimiter(const std::string& value) {
    if (value == "auto") {
        return ProcessExtcat::Delimiter::Auto;
    }
    if (value == "whitespace") {
        return ProcessExtcat::Delimiter::Whitespace;
    }
    if (value == "comma") {
        return ProcessExtcat::Delimiter::Comma;
    }
    if (value == "tab") {
        return ProcessExtcat::Delimiter::Tab;
    }
    throw std::runtime_error(
        "extcat delimiter must be auto, whitespace, comma, or tab");
}

// ==========================================
// Function: Convert configured header-mode text
// Method: Map automatic, required, and absent header modes to the reusable enum.
// ==========================================
ProcessExtcat::HeaderMode configuredHeaderMode(const std::string& value) {
    if (value == "auto") {
        return ProcessExtcat::HeaderMode::Auto;
    }
    if (value == "present") {
        return ProcessExtcat::HeaderMode::Present;
    }
    if (value == "absent") {
        return ProcessExtcat::HeaderMode::Absent;
    }
    throw std::runtime_error("extcat header mode must be auto, present, or absent");
}

// ==========================================
// Function: Convert configured malformed-row policy
// Method: Restrict the integrated phase to the reusable fail and skip behaviors.
// ==========================================
ProcessExtcat::MalformedPolicy configuredMalformedPolicy(const std::string& value) {
    if (value == "fail") {
        return ProcessExtcat::MalformedPolicy::Fail;
    }
    if (value == "skip") {
        return ProcessExtcat::MalformedPolicy::Skip;
    }
    throw std::runtime_error("extcat malformed policy must be fail or skip");
}

// ==========================================
// Function: Convert configured existing-output policy
// Method: Restrict publication to fail-safe refusal or transactional overwrite.
// ==========================================
ProcessExtcat::ExistingPolicy configuredExtcatExistingPolicy(const std::string& value) {
    if (value == "fail") {
        return ProcessExtcat::ExistingPolicy::Fail;
    }
    if (value == "overwrite") {
        return ProcessExtcat::ExistingPolicy::Overwrite;
    }
    throw std::runtime_error("extcat existing policy must be fail or overwrite");
}

// ==========================================
// Function: Build reusable tiler configuration from pipeline options
// Method: Convert variable one-based output and coordinate columns plus MiB task sizes while
//         preserving every validated parser and publication setting.
// ==========================================
ProcessExtcat::Config buildIntegratedConfig(const ProcessConfig::RuntimeOptions& options) {
    constexpr std::uint64_t bytes_per_mib = 1024ULL * 1024ULL;
    if (options.extcat_chunk_mib == 0
        || options.extcat_chunk_mib
               > std::numeric_limits<std::uint64_t>::max() / bytes_per_mib) {
        throw std::runtime_error("extcat chunk MiB must be a positive uint64 value");
    }

    ProcessExtcat::Config config;
    config.input_directory = options.extcat_input_directory;
    config.output_directory = options.extcat_output_directory;
    config.filename_tokens = options.extcat_filename_tokens;
    config.recursive = options.extcat_recursive;
    config.delimiter = configuredDelimiter(options.extcat_delimiter);
    config.header_mode = configuredHeaderMode(options.extcat_header_mode);
    config.malformed_policy = configuredMalformedPolicy(options.extcat_malformed_policy);
    config.existing_policy = configuredExtcatExistingPolicy(
        options.extcat_existing_policy);
    config.chunk_bytes = options.extcat_chunk_mib * bytes_per_mib;
    config.use_explicit_columns = options.extcat_use_explicit_columns;
    if (config.use_explicit_columns) {
        config.input_columns.clear();
        config.input_columns.reserve(options.extcat_input_columns_one_based.size());
        for (const std::size_t one_based : options.extcat_input_columns_one_based) {
            if (one_based == 0) {
                throw std::runtime_error(
                    "extcat explicit columns must use positive one-based indices");
            }
            config.input_columns.push_back(one_based - 1);
        }
    }
    config.use_explicit_coordinate_columns =
        options.extcat_use_explicit_coordinate_columns;
    if (options.extcat_ra_column_one_based == 0
        || options.extcat_dec_column_one_based == 0) {
        throw std::runtime_error(
            "extcat coordinate columns must use positive one-based indices");
    }
    config.ra_column = options.extcat_ra_column_one_based - 1;
    config.dec_column = options.extcat_dec_column_one_based - 1;
    return config;
}

}  // namespace

// ==========================================
// Function: Split raw external catalogs into one-degree sky tiles
// Method: Discover and inspect inputs on rank zero, preserve or project arbitrary columns,
//         process byte ranges across MPI ranks, and publish deterministic tiles.
// ==========================================
int process_extcat(ProcessExtcat::Config config, MPI_Comm communicator) {
    int rank = 0;
    int world_size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &world_size);

    int local_config_ok = 1;
    std::string local_error;
    try {
        ProcessExtcat::normalizeAndValidateConfig(config);
    } catch (const std::exception& exception) {
        local_config_ok = 0;
        local_error = exception.what();
    }
    int global_config_ok = 0;
    MPI_Allreduce(&local_config_ok, &global_config_ok, 1, MPI_INT, MPI_MIN, communicator);
    if (global_config_ok == 0) {
        reportRankErrors(local_error, communicator);
        return 2;
    }

    std::vector<FileMetadata> metadata;
    fs::path staging_directory;
    int root_preparation_ok = 1;
    std::string root_error;
    if (rank == 0) {
        try {
            const std::vector<fs::path> inputs = discoverInputFiles(config);
            const std::vector<fs::path> existing = listExistingTiles(config.output_directory);
            if (config.existing_policy == ProcessExtcat::ExistingPolicy::Fail
                && !existing.empty()) {
                throw std::runtime_error("output directory already contains generated tiles; "
                                         "use --existing overwrite: "
                                         + config.output_directory.string());
            }
            fs::create_directories(config.output_directory);
            staging_directory = createStagingDirectory(config.output_directory);
            metadata.reserve(inputs.size());
            for (const fs::path& path : inputs) {
                metadata.push_back(inspectFile(path, config));
            }
            validateCompatibleOutputSchemas(metadata);
        } catch (const std::exception& exception) {
            root_preparation_ok = 0;
            root_error = exception.what();
            if (!staging_directory.empty()) {
                std::error_code cleanup_error;
                fs::remove_all(staging_directory, cleanup_error);
            }
        }
    }

    MPI_Bcast(&root_preparation_ok, 1, MPI_INT, 0, communicator);
    broadcastString(root_error, 0, communicator);
    if (root_preparation_ok == 0) {
        if (rank == 0) {
            std::cerr << "process_extcat preparation error: " << root_error << '\n';
        }
        return 2;
    }

    std::string staging_text = rank == 0 ? staging_directory.string() : std::string();
    broadcastString(staging_text, 0, communicator);
    if (rank != 0) {
        staging_directory = staging_text;
    }
    broadcastMetadata(metadata, 0, communicator);

    std::vector<Task> tasks;
    int local_task_build_ok = 1;
    local_error.clear();
    try {
        tasks = buildTasks(metadata, config.chunk_bytes, world_size);
    } catch (const std::exception& exception) {
        local_task_build_ok = 0;
        local_error = exception.what();
    }
    int global_task_build_ok = 0;
    MPI_Allreduce(&local_task_build_ok, &global_task_build_ok, 1, MPI_INT, MPI_MIN,
                  communicator);
    if (global_task_build_ok == 0) {
        reportRankErrors(local_error, communicator);
        if (rank == 0) {
            std::error_code cleanup_error;
            fs::remove_all(staging_directory, cleanup_error);
        }
        return 2;
    }

    TaskStats local_stats;
    int local_processing_ok = 1;
    local_error.clear();
    for (const Task& task : tasks) {
        if (task.id % static_cast<std::uint64_t>(world_size)
            != static_cast<std::uint64_t>(rank)) {
            continue;
        }
        try {
            const TaskStats task_stats = processTask(
                task, metadata[task.file_index], config, staging_directory);
            local_stats.accepted_rows += task_stats.accepted_rows;
            local_stats.malformed_rows += task_stats.malformed_rows;
        } catch (const std::exception& exception) {
            local_processing_ok = 0;
            local_error = exception.what();
            break;
        }
    }

    int global_processing_ok = 0;
    MPI_Allreduce(&local_processing_ok, &global_processing_ok, 1, MPI_INT, MPI_MIN,
                  communicator);
    if (global_processing_ok == 0) {
        reportRankErrors(local_error, communicator);
        if (rank == 0) {
            std::error_code cleanup_error;
            fs::remove_all(staging_directory, cleanup_error);
        }
        return 1;
    }

    std::uint64_t local_counts[2] = {
        local_stats.accepted_rows,
        local_stats.malformed_rows,
    };
    std::uint64_t global_counts[2] = {0, 0};
    MPI_Reduce(local_counts, global_counts, 2, MPI_UINT64_T, MPI_SUM, 0, communicator);
    MPI_Barrier(communicator);

    int merge_ok = 1;
    std::size_t tile_count = 0;
    std::string merge_error;
    if (rank == 0) {
        try {
            tile_count = mergeShards(
                tasks, staging_directory, metadata.front().output_header, config);
            fs::remove_all(staging_directory);
        } catch (const std::exception& exception) {
            merge_ok = 0;
            merge_error = std::string(exception.what())
                          + "; staging retained at " + staging_directory.string();
        }
    }
    MPI_Bcast(&merge_ok, 1, MPI_INT, 0, communicator);
    broadcastString(merge_error, 0, communicator);
    if (merge_ok == 0) {
        if (rank == 0) {
            std::cerr << "process_extcat merge error: " << merge_error << '\n';
        }
        return 1;
    }

    if (rank == 0) {
        std::cout << "process_extcat completed: files=" << metadata.size()
                  << " tasks=" << tasks.size()
                  << " accepted_rows=" << global_counts[0]
                  << " skipped_malformed_rows=" << global_counts[1]
                  << " tiles=" << tile_count
                  << " output=" << config.output_directory << '\n';
    }
    return 0;
}

// ==========================================
// Function: Run external-catalog tiling from unified pipeline options
// Method: Translate configuration collectively, then call the reusable implementation
//         without taking ownership of MPI initialization or finalization.
// ==========================================
int process_extcat(const ProcessConfig::RuntimeOptions& options, MPI_Comm communicator) {
    ProcessExtcat::Config config;
    int local_adapter_ok = 1;
    std::string local_error;
    try {
        config = buildIntegratedConfig(options);
    } catch (const std::exception& exception) {
        local_adapter_ok = 0;
        local_error = exception.what();
    }

    int global_adapter_ok = 0;
    MPI_Allreduce(&local_adapter_ok, &global_adapter_ok, 1, MPI_INT, MPI_MIN,
                  communicator);
    if (global_adapter_ok == 0) {
        reportRankErrors(local_error, communicator);
        return 2;
    }
    return process_extcat(config, communicator);
}
