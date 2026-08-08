#include "process_init/Initializer.hpp"
#include "OutputLayout.hpp"

#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace fqinit {
namespace {

namespace fs = std::filesystem;

constexpr const char* kArchiveSuffix = ".fits.fz";

struct Task {
    ProductKind kind = ProductKind::Science;
    fs::path source;
};

// ==========================================
// Function: Test whether a string ends with one exact suffix
// Method: Compare the final suffix-sized span without locale transformations.
// ==========================================
bool endsWith(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size()
           && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// ==========================================
// Function: Test whether a filename starts with the configured archive prefix
// Method: Compare the prefix at offset zero to reproduce the original find rule.
// ==========================================
bool startsWith(const std::string& value, const std::string& prefix) {
    return value.size() >= prefix.size() && value.compare(0, prefix.size(), prefix) == 0;
}

// ==========================================
// Function: Read the value following one command-line option
// Method: Advance the shared argument index and reject a missing value.
// ==========================================
std::string optionValue(int& index, int argc, char** argv, const std::string& option) {
    if (index + 1 >= argc) {
        throw std::runtime_error("missing value after " + option);
    }
    ++index;
    return argv[index];
}

// ==========================================
// Function: Parse the existing-output policy name
// Method: Accept only the three explicit safe lifecycle modes.
// ==========================================
ExistingPolicy parseExistingPolicy(const std::string& value) {
    if (value == "fail") {
        return ExistingPolicy::Fail;
    }
    if (value == "resume") {
        return ExistingPolicy::Resume;
    }
    if (value == "overwrite") {
        return ExistingPolicy::Overwrite;
    }
    throw std::runtime_error("--existing must be fail, resume, or overwrite");
}

// ==========================================
// Function: Convert the existing-output policy to manifest text
// Method: Return one stable lowercase schema value for each policy.
// ==========================================
const char* existingPolicyName(ExistingPolicy policy) {
    switch (policy) {
        case ExistingPolicy::Fail:
            return "fail";
        case ExistingPolicy::Resume:
            return "resume";
        case ExistingPolicy::Overwrite:
            return "overwrite";
    }
    return "unknown";
}

// ==========================================
// Function: Normalize a user path for deterministic lists and boundary checks
// Method: Resolve existing symlink components while allowing a new output tail.
// ==========================================
fs::path normalizedAbsolute(const fs::path& path) {
    return fs::weakly_canonical(fs::absolute(path));
}

// ==========================================
// Function: Test whether one normalized path is equal to or below another
// Method: Compare path components so output creation cannot enter a source root.
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
// Function: Match one archive basename against configured tokens
// Method: Disable token filtering for an empty list; otherwise accept when any
//         non-empty token occurs in the basename.
// ==========================================
bool matchesAnyToken(const std::string& filename,
                     const std::vector<std::string>& tokens) {
    if (tokens.empty()) {
        return true;
    }
    return std::any_of(tokens.begin(), tokens.end(), [&filename](const std::string& token) {
        return filename.find(token) != std::string::npos;
    });
}

// ==========================================
// Function: Discover matching archives below one immutable source root
// Method: Recursively collect regular .fits.fz files, apply basename filters,
//         normalize absolute paths, and sort the in-memory result.
// ==========================================
std::vector<fs::path> discoverArchives(const fs::path& root,
                                       const std::string& prefix,
                                       const std::vector<std::string>& tokens) {
    if (!fs::exists(root) || !fs::is_directory(root)) {
        throw std::runtime_error("archive root is not a directory: " + root.string());
    }

    std::vector<fs::path> paths;
    std::error_code iterator_error;
    fs::recursive_directory_iterator iterator(root, fs::directory_options::none, iterator_error);
    const fs::recursive_directory_iterator end;
    if (iterator_error) {
        throw std::runtime_error("cannot scan archive root " + root.string() + ": "
                                 + iterator_error.message());
    }
    while (iterator != end) {
        const fs::directory_entry entry = *iterator;
        std::error_code type_error;
        const bool regular = entry.is_regular_file(type_error);
        if (type_error) {
            throw std::runtime_error("cannot inspect archive entry " + entry.path().string()
                                     + ": " + type_error.message());
        }
        if (regular) {
            const std::string filename = entry.path().filename().string();
            if (endsWith(filename, kArchiveSuffix)
                && startsWith(filename, prefix)
                && matchesAnyToken(filename, tokens)) {
                paths.push_back(normalizedAbsolute(entry.path()));
            }
        }
        iterator.increment(iterator_error);
        if (iterator_error) {
            throw std::runtime_error("archive traversal failed below " + root.string() + ": "
                                     + iterator_error.message());
        }
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

// ==========================================
// Function: Reject source archives that would target the same output exposure
// Method: Compare science stems directly and DQ stems after ood-to-ooi mapping.
// ==========================================
void validateUniqueStems(const std::vector<fs::path>& paths, ProductKind kind) {
    std::map<std::string, fs::path> seen;
    for (const fs::path& path : paths) {
        const std::string stem = kind == ProductKind::Science
                                     ? archiveStem(path)
                                     : dqOutputStem(path);
        const auto [position, inserted] = seen.emplace(stem, path);
        if (!inserted) {
            throw std::runtime_error("duplicate output exposure stem " + stem + " from "
                                     + position->second.string() + " and " + path.string());
        }
    }
}

// ==========================================
// Function: Broadcast one dynamically sized string
// Method: Broadcast an int length followed by contiguous bytes.
// ==========================================
void broadcastString(std::string& value, int root_rank, MPI_Comm communicator) {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int length = rank == root_rank ? static_cast<int>(value.size()) : 0;
    MPI_Bcast(&length, 1, MPI_INT, root_rank, communicator);
    if (length < 0) {
        throw std::runtime_error("received a negative MPI string length");
    }
    if (rank != root_rank) {
        value.resize(static_cast<std::size_t>(length));
    }
    if (length > 0) {
        MPI_Bcast(value.data(), length, MPI_CHAR, root_rank, communicator);
    }
}

// ==========================================
// Function: Broadcast an in-memory vector of filesystem paths
// Method: Send the vector count and reuse the length-prefixed string broadcast.
// ==========================================
void broadcastPaths(std::vector<fs::path>& paths, int root_rank, MPI_Comm communicator) {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int count = rank == root_rank ? static_cast<int>(paths.size()) : 0;
    MPI_Bcast(&count, 1, MPI_INT, root_rank, communicator);
    if (count < 0) {
        throw std::runtime_error("received a negative MPI path count");
    }
    if (rank != root_rank) {
        paths.assign(static_cast<std::size_t>(count), fs::path());
    }
    for (int index = 0; index < count; ++index) {
        std::string text = rank == root_rank ? paths[static_cast<std::size_t>(index)].string()
                                             : std::string();
        broadcastString(text, root_rank, communicator);
        if (rank != root_rank) {
            paths[static_cast<std::size_t>(index)] = fs::path(text);
        }
    }
}

// ==========================================
// Function: Create the complete base-directory contract consumed by all pipeline modes
// Method: Apply the two shared OutputLayout directory sets idempotently.
// ==========================================
void createPipelineDirectories(const fs::path& target_root) {
    for (const char* name : OutputLayout::NON_CHIP_BASE_DIRECTORIES) {
        fs::create_directories(target_root / name);
    }
    for (const char* name : OutputLayout::CHIP_PRODUCT_DIRECTORIES) {
        fs::create_directories(target_root / name);
    }
}

// ==========================================
// Function: Create every chip-product directory for one exposure
// Method: Apply the shared OutputLayout contract idempotently after list publication.
// ==========================================
void createExposureProductDirectories(const fs::path& target_root,
                                      const std::string& exposure) {
    for (const char* name : OutputLayout::CHIP_PRODUCT_DIRECTORIES) {
        fs::create_directories(target_root / name / exposure);
    }
}

// ==========================================
// Function: Construct one collision-resistant staging run identifier
// Method: Combine Unix epoch microseconds with the rank-zero process id.
// ==========================================
std::string makeRunToken() {
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    const auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
    return "run_" + std::to_string(microseconds) + "_" + std::to_string(getpid());
}

// ==========================================
// Function: Escape one string for the initializer JSON manifest
// Method: Encode JSON control characters and preserve all other UTF-8 bytes.
// ==========================================
std::string jsonEscape(const std::string& value) {
    std::ostringstream escaped;
    for (const unsigned char character : value) {
        switch (character) {
            case '"':
                escaped << "\\\"";
                break;
            case '\\':
                escaped << "\\\\";
                break;
            case '\b':
                escaped << "\\b";
                break;
            case '\f':
                escaped << "\\f";
                break;
            case '\n':
                escaped << "\\n";
                break;
            case '\r':
                escaped << "\\r";
                break;
            case '\t':
                escaped << "\\t";
                break;
            default:
                if (character < 0x20U) {
                    escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                            << static_cast<int>(character) << std::dec;
                } else {
                    escaped << static_cast<char>(character);
                }
        }
    }
    return escaped.str();
}

// ==========================================
// Function: Publish one small metadata or list file atomically
// Method: Write a sibling temporary file, flush it, then replace via POSIX rename.
// ==========================================
void writeAtomic(const fs::path& path, const std::string& content) {
    fs::create_directories(path.parent_path());
    const fs::path temporary = path.string() + ".tmp." + std::to_string(getpid());
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output.is_open()) {
            throw std::runtime_error("cannot create temporary metadata file: "
                                     + temporary.string());
        }
        output << content;
        output.flush();
        if (!output.good()) {
            throw std::runtime_error("cannot flush temporary metadata file: "
                                     + temporary.string());
        }
    }
    if (std::rename(temporary.c_str(), path.c_str()) != 0) {
        const std::string error = std::strerror(errno);
        std::error_code cleanup_error;
        fs::remove(temporary, cleanup_error);
        throw std::runtime_error("cannot publish " + path.string() + ": " + error);
    }
}

// ==========================================
// Function: Reject paths the legacy F77 reader cannot consume safely
// Method: Enforce no whitespace and the configured fixed-character limit.
// ==========================================
void validatePipelinePath(const fs::path& path, int f77_max_path) {
    const std::string text = path.string();
    if (std::any_of(text.begin(), text.end(), [](unsigned char character) {
            return std::isspace(character) != 0;
        })) {
        throw std::runtime_error("pipeline path contains whitespace: " + text);
    }
    if (f77_max_path > 0 && text.size() > static_cast<std::size_t>(f77_max_path)) {
        throw std::runtime_error("pipeline path exceeds the configured F77 limit: " + text);
    }
}

// ==========================================
// Function: Publish top-level exposure and compatibility FITS lists
// Method: After every rank has written its per-exposure chip lists into
//         stamps/ during extraction, rank zero scans stamps/, sorts the
//         per-exposure lists, and atomically writes expo_<target>.list
//         ("<list path>" <chip count>) and the flat fits_<target>.list.
//         No re-stat of science images is required.
// ==========================================
void publishPipelineLists(const Config& config, const fs::path& target_root) {
    const fs::path stamps_dir = target_root / "stamps";
    std::vector<fs::path> exposure_lists;
    for (const auto& entry : fs::directory_iterator(stamps_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".list") {
            exposure_lists.push_back(entry.path());
        }
    }
    std::sort(exposure_lists.begin(), exposure_lists.end());

    std::ostringstream top_list;
    std::ostringstream fits_list;
    for (const fs::path& list_path : exposure_lists) {
        validatePipelinePath(list_path, config.f77_max_path);
        std::ifstream list_input(list_path);
        int chip_count = 0;
        std::string line;
        while (std::getline(list_input, line)) {
            while (!line.empty()
                   && (line.back() == '\r' || line.back() == ' '
                       || line.back() == '\t')) {
                line.pop_back();
            }
            if (line.empty()) {
                continue;
            }
            fits_list << line << '\n';
            ++chip_count;
        }
        top_list << '"' << list_path.string() << '"' << "     " << chip_count << '\n';
    }

    const fs::path top_path = config.output_root / ("expo_" + config.target + ".list");
    const fs::path fits_path = config.output_root / ("fits_" + config.target + ".list");
    validatePipelinePath(top_path, config.f77_max_path);
    validatePipelinePath(fits_path, config.f77_max_path);
    writeAtomic(top_path, top_list.str());
    writeAtomic(fits_path, fits_list.str());
}

// ==========================================
// Function: Materialize chip-product exposure directories from the published expo list
// Method: Parse expo_<target>.list after its atomic publication, validate each generated
//         per-exposure list path, derive the exposure from its .list basename, and create
//         every shared chip-product directory idempotently.
// ==========================================
void createExposureDirectoriesFromPublishedList(const Config& config,
                                                const fs::path& target_root) {
    const fs::path top_path = config.output_root / ("expo_" + config.target + ".list");
    std::ifstream input(top_path);
    if (!input.is_open()) {
        throw std::runtime_error("cannot read published exposure list: "
                                 + top_path.string());
    }

    const fs::path expected_parent = (target_root / "stamps").lexically_normal();
    std::set<std::string> exposures;
    std::string list_path_text;
    int chip_count = 0;
    while (input >> std::quoted(list_path_text) >> chip_count) {
        const fs::path list_path = fs::path(list_path_text).lexically_normal();
        if (chip_count < 0 || list_path.extension() != ".list"
            || list_path.parent_path() != expected_parent) {
            throw std::runtime_error("invalid record in published exposure list: "
                                     + list_path_text);
        }
        const std::string exposure = list_path.stem().string();
        if (exposure.empty()) {
            throw std::runtime_error("empty exposure name in published exposure list: "
                                     + list_path_text);
        }
        exposures.insert(exposure);
    }
    if (!input.eof()) {
        throw std::runtime_error("malformed published exposure list: "
                                 + top_path.string());
    }

    for (const std::string& exposure : exposures) {
        createExposureProductDirectories(target_root, exposure);
    }
}

// ==========================================
// Function: Serialize one durable initialization manifest
// Method: Record fixed naming/order decisions, direct-read provenance, counts,
//         existing-output policy, and every failed source path.
// ==========================================
std::string makeManifest(const Config& config,
                         const std::vector<Task>& tasks,
                         const std::vector<int>& statuses,
                         const std::vector<int>& image_counts,
                         const std::vector<int>& resumed_flags,
                         const std::vector<int>& skipped_hdus,
                         bool lists_published,
                         bool exposure_directories_created,
                         const std::string& final_error) {
    int science_sources = 0;
    int dq_sources = 0;
    int science_images = 0;
    int dq_images = 0;
    int resumed_sources = 0;
    int skipped_hdus_total = 0;
    std::vector<std::string> failed_sources;
    std::vector<std::string> partial_sources;
    for (std::size_t index = 0; index < tasks.size(); ++index) {
        if (tasks[index].kind == ProductKind::Science) {
            ++science_sources;
            science_images += std::max(0, image_counts[index]);
        } else {
            ++dq_sources;
            dq_images += std::max(0, image_counts[index]);
        }
        resumed_sources += resumed_flags[index] > 0 ? 1 : 0;
        skipped_hdus_total += skipped_hdus[index];
        if (statuses[index] <= 0) {
            failed_sources.push_back(tasks[index].source.string());
        } else if (skipped_hdus[index] > 0) {
            partial_sources.push_back(tasks[index].source.string());
        }
    }

    const bool has_failures = !failed_sources.empty() || !partial_sources.empty();
    const std::string status_text =
        (!lists_published || !exposure_directories_created)
            ? "failed"
            : has_failures ? "partial" : "success";
    std::ostringstream manifest;
    manifest << "{\n"
             << "  \"schema_version\": 2,\n"
             << "  \"status\": \"" << status_text << "\",\n"
             << "  \"direct_source_read\": true,\n"
             << "  \"copy_staging\": false,\n"
             << "  \"exposure_order\": \"corrected_lexical_no_rotation\",\n"
             << "  \"science_numbering\": \"two_dimensional_hdu_occurrence\",\n"
             << "  \"dq_numbering\": \"CCDNUM\",\n"
             << "  \"science_root\": \"" << jsonEscape(config.science_root.string()) << "\",\n"
             << "  \"dq_root\": \"" << jsonEscape(config.dq_root.string()) << "\",\n"
             << "  \"output_root\": \"" << jsonEscape(config.output_root.string()) << "\",\n"
             << "  \"target\": \"" << jsonEscape(config.target) << "\",\n"
             << "  \"filename_prefix\": \"" << jsonEscape(config.filename_prefix) << "\",\n"
             << "  \"filename_tokens\": [";
    for (std::size_t index = 0; index < config.filename_tokens.size(); ++index) {
        if (index != 0) {
            manifest << ", ";
        }
        manifest << '"' << jsonEscape(config.filename_tokens[index]) << '"';
    }
    manifest << "],\n"
             << "  \"existing_policy\": \"" << existingPolicyName(config.existing_policy) << "\",\n"
             << "  \"f77_max_path\": " << config.f77_max_path << ",\n"
             << "  \"science_sources\": " << science_sources << ",\n"
             << "  \"dq_sources\": " << dq_sources << ",\n"
             << "  \"science_images\": " << science_images << ",\n"
             << "  \"dq_images\": " << dq_images << ",\n"
             << "  \"resumed_sources\": " << resumed_sources << ",\n"
             << "  \"lists_published\": " << (lists_published ? "true" : "false") << ",\n"
             << "  \"exposure_directories_created\": "
             << (exposure_directories_created ? "true" : "false") << ",\n"
             << "  \"error\": \"" << jsonEscape(final_error) << "\",\n"
             << "  \"failed_sources\": [";
    for (std::size_t index = 0; index < failed_sources.size(); ++index) {
        if (index != 0) {
            manifest << ", ";
        }
        manifest << '"' << jsonEscape(failed_sources[index]) << '"';
    }
    manifest << "],\n"
             << "  \"skipped_hdus\": " << skipped_hdus_total << ",\n"
             << "  \"partial_sources\": [";
    for (std::size_t index = 0; index < partial_sources.size(); ++index) {
        if (index != 0) {
            manifest << ", ";
        }
        manifest << '"' << jsonEscape(partial_sources[index]) << '"';
    }
    manifest << "]\n}\n";
    return manifest.str();
}

}  // namespace

// ==========================================
// Function: Parse the standalone initializer command line
// Method: Require explicit archive/output roots and matching keys while
//         retaining safe defaults for existing outputs and F77 path limits.
// ==========================================
Config parseArguments(int argc, char** argv) {
    Config config;
    bool contains_option_seen = false;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--science-root") {
            config.science_root = optionValue(index, argc, argv, option);
        } else if (option == "--dq-root") {
            config.dq_root = optionValue(index, argc, argv, option);
        } else if (option == "--output-root") {
            config.output_root = optionValue(index, argc, argv, option);
        } else if (option == "--target") {
            config.target = optionValue(index, argc, argv, option);
        } else if (option == "--prefix") {
            config.filename_prefix = optionValue(index, argc, argv, option);
        } else if (option == "--contains") {
            const std::string token = optionValue(index, argc, argv, option);
            if (token.empty()) {
                throw std::runtime_error("--contains must not be empty");
            }
            if (!contains_option_seen) {
                config.filename_tokens.clear();
                contains_option_seen = true;
            }
            config.filename_tokens.push_back(token);
        } else if (option == "--existing") {
            config.existing_policy = parseExistingPolicy(optionValue(index, argc, argv, option));
        } else if (option == "--f77-max-path") {
            const std::string value = optionValue(index, argc, argv, option);
            std::size_t parsed_characters = 0;
            try {
                config.f77_max_path = std::stoi(value, &parsed_characters);
            } catch (const std::exception&) {
                throw std::runtime_error("--f77-max-path must be an integer");
            }
            if (parsed_characters != value.size()) {
                throw std::runtime_error("--f77-max-path must be an integer");
            }
            if (config.f77_max_path < 0) {
                throw std::runtime_error("--f77-max-path must be zero or positive");
            }
        } else if (option != "--help") {
            throw std::runtime_error("unknown option: " + option);
        }
    }

    normalizeAndValidateConfig(config);
    return config;
}

// ==========================================
// Function: Normalize and validate one initializer configuration
// Method: Resolve all paths and enforce the same required-field and target-name
//         contract for standalone parsing and integrated workflow execution.
// ==========================================
void normalizeAndValidateConfig(Config& config) {
    if (config.science_root.empty() || config.dq_root.empty() || config.output_root.empty()
        || config.target.empty() || config.filename_prefix.empty()) {
        throw std::runtime_error("science root, DQ root, output root, target, and prefix are required");
    }
    if (config.target == "." || config.target == ".."
        || config.target.find('/') != std::string::npos
        || config.target.find('\\') != std::string::npos) {
        throw std::runtime_error("--target must be one directory name");
    }
    if (config.f77_max_path < 0) {
        throw std::runtime_error("--f77-max-path must be zero or positive");
    }
    for (const std::string& token : config.filename_tokens) {
        if (token.empty()) {
            throw std::runtime_error("filename tokens must be non-empty");
        }
    }

    config.science_root = normalizedAbsolute(config.science_root);
    config.dq_root = normalizedAbsolute(config.dq_root);
    config.output_root = normalizedAbsolute(config.output_root);
    validatePipelinePath(config.output_root / config.target, 0);
}

// ==========================================
// Function: Print the portable initializer command-line contract
// Method: Describe required paths, filters, existing-output policy, and limits.
// ==========================================
void printUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " [options]\n"
        << "  --science-root PATH    Original read-only science FITS/FZ repository\n"
        << "  --dq-root PATH         Original read-only DQ FITS/FZ repository\n"
        << "  --output-root PATH     Parent of target and expo_<target>.list\n"
        << "  --target NAME          Target directory, for example z2015\n"
        << "  --prefix TEXT          Required filename prefix, for example c4d_15\n"
        << "  --contains TEXT        Repeatable filename token; matches any token (default: v1)\n"
        << "  --existing MODE        fail (default), resume, or overwrite\n"
        << "  --f77-max-path N       Maximum generated path length; 0 disables (default: "
        << InitConfig::F77_MAX_PATH << ")\n"
        << "  --help                  Show this help\n";
}

// ==========================================
// Function: Build the pipeline input layout from the original FITS/FZ archives
// Method: Discover on rank zero, broadcast in-memory paths, extract archives
//         directly in parallel, and publish corrected deterministic lists.
// ==========================================
int runInitializer(const Config& input_config, MPI_Comm communicator) {
    int rank = 0;
    int process_count = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &process_count);

    Config config = input_config;
    const fs::path target_root = config.output_root / config.target;
    std::vector<fs::path> science_sources;
    std::vector<fs::path> dq_sources;
    std::string setup_error;
    std::string run_token;
    int setup_ok = 1;

    if (rank == 0) {
        try {
            if (pathIsWithin(target_root, config.science_root)
                || pathIsWithin(target_root, config.dq_root)) {
                throw std::runtime_error(
                    "target output must not be inside either read-only source repository");
            }
            science_sources = discoverArchives(
                config.science_root, config.filename_prefix, config.filename_tokens);
            dq_sources = discoverArchives(
                config.dq_root, config.filename_prefix, config.filename_tokens);
            if (science_sources.empty()) {
                throw std::runtime_error("no matching science .fits.fz archives were found");
            }
            if (dq_sources.empty()) {
                throw std::runtime_error("no matching DQ .fits.fz archives were found");
            }
            validateUniqueStems(science_sources, ProductKind::Science);
            validateUniqueStems(dq_sources, ProductKind::DqMask);
            fs::create_directories(config.output_root);
            createPipelineDirectories(target_root);
            run_token = makeRunToken();
            fs::create_directories(target_root / ".fq_init_tmp" / run_token);
        } catch (const std::exception& exception) {
            setup_ok = 0;
            setup_error = exception.what();
        }
    }

    MPI_Bcast(&setup_ok, 1, MPI_INT, 0, communicator);
    broadcastString(setup_error, 0, communicator);
    if (setup_ok == 0) {
        if (rank == 0) {
            std::cerr << "Initializer setup failed: " << setup_error << std::endl;
        }
        return 1;
    }
    broadcastPaths(science_sources, 0, communicator);
    broadcastPaths(dq_sources, 0, communicator);
    broadcastString(run_token, 0, communicator);

    std::vector<Task> tasks;
    tasks.reserve(science_sources.size() + dq_sources.size());
    for (const fs::path& source : science_sources) {
        tasks.push_back({ProductKind::Science, source});
    }
    for (const fs::path& source : dq_sources) {
        tasks.push_back({ProductKind::DqMask, source});
    }

    std::vector<int> local_status(tasks.size(), 0);
    std::vector<int> local_counts(tasks.size(), 0);
    std::vector<int> local_resumed(tasks.size(), 0);
    std::vector<int> local_skipped(tasks.size(), 0);
    const fs::path staging_root = target_root / ".fq_init_tmp" / run_token;

    for (std::size_t task_index = static_cast<std::size_t>(rank);
         task_index < tasks.size();
         task_index += static_cast<std::size_t>(process_count)) {
        const Task& task = tasks[task_index];
        const std::string exposure_stem = (task.kind == ProductKind::Science)
                                              ? archiveStem(task.source)
                                              : dqOutputStem(task.source);
        const fs::path final_directory = target_root
                                          / (task.kind == ProductKind::Science ? "science" : "dqmask")
                                          / exposure_stem;
        const fs::path task_staging = staging_root / ("rank_" + std::to_string(rank))
                                      / ("task_" + std::to_string(task_index));
        const ExtractionResult result = extractArchive(
            task.source, task.kind, final_directory, task_staging, config.existing_policy);
        local_status[task_index] = result.success ? 1 : -1;
        local_counts[task_index] = result.success
                                       ? static_cast<int>(result.output_paths.size())
                                       : 0;
        local_resumed[task_index] = result.resumed ? 1 : 0;
        local_skipped[task_index] = result.skipped_hdus;
        if (!result.success) {
            std::cerr << "[rank " << rank << "] failed " << task.source << ": "
                      << result.error << std::endl;
        }
        if (task.kind == ProductKind::Science && result.success
            && !result.output_paths.empty()) {
            try {
                std::vector<fs::path> sorted_paths = result.output_paths;
                std::sort(sorted_paths.begin(), sorted_paths.end());
                const std::string exposure = archiveStem(task.source);
                const fs::path list_path = target_root / "stamps" / (exposure + ".list");
                validatePipelinePath(list_path, config.f77_max_path);
                for (const fs::path& image_path : sorted_paths) {
                    validatePipelinePath(image_path, config.f77_max_path);
                }
                std::ostringstream exposure_list;
                for (const fs::path& image_path : sorted_paths) {
                    exposure_list << image_path.string() << '\n';
                }
                writeAtomic(list_path, exposure_list.str());
            } catch (const std::exception& exception) {
                local_status[task_index] = -1;
                std::cerr << "[rank " << rank << "] failed to publish per-exposure list for "
                          << task.source << ": " << exception.what() << std::endl;
            }
        }
    }

    std::vector<int> global_status(tasks.size(), 0);
    std::vector<int> global_counts(tasks.size(), 0);
    std::vector<int> global_resumed(tasks.size(), 0);
    std::vector<int> global_skipped(tasks.size(), 0);
    if (!tasks.empty()) {
        const int mpi_count = static_cast<int>(tasks.size());
        MPI_Allreduce(local_status.data(), global_status.data(), mpi_count,
                      MPI_INT, MPI_SUM, communicator);
        MPI_Allreduce(local_counts.data(), global_counts.data(), mpi_count,
                      MPI_INT, MPI_SUM, communicator);
        MPI_Allreduce(local_resumed.data(), global_resumed.data(), mpi_count,
                      MPI_INT, MPI_SUM, communicator);
        MPI_Allreduce(local_skipped.data(), global_skipped.data(), mpi_count,
                      MPI_INT, MPI_SUM, communicator);
    }
    MPI_Barrier(communicator);

    bool lists_published = false;
    bool exposure_directories_created = false;
    std::string final_error;
    int final_status = 0;
    if (rank == 0) {
        const bool any_failed = std::any_of(
            global_status.begin(), global_status.end(), [](int value) { return value < 0; });
        try {
            publishPipelineLists(config, target_root);
            lists_published = true;
            createExposureDirectoriesFromPublishedList(config, target_root);
            exposure_directories_created = true;
            if (any_failed) {
                final_error = "some source archives failed; pipeline lists published with available data";
            }
        } catch (const std::exception& exception) {
            final_status = 1;
            final_error = exception.what();
        }

        try {
            const fs::path manifest_path = config.output_root
                                           / ("init_" + config.target + "_manifest.json");
            writeAtomic(manifest_path, makeManifest(
                config, tasks, global_status, global_counts, global_resumed,
                global_skipped, lists_published, exposure_directories_created,
                final_error));
        } catch (const std::exception& exception) {
            final_status = 1;
            if (!final_error.empty()) {
                final_error += "; ";
            }
            final_error += exception.what();
        }

        std::error_code cleanup_error;
        fs::remove_all(staging_root, cleanup_error);
        const fs::path staging_parent = target_root / ".fq_init_tmp";
        if (!cleanup_error && fs::is_empty(staging_parent, cleanup_error)) {
            fs::remove(staging_parent, cleanup_error);
        }

        if (final_status == 0) {
            std::cout << "Initialization complete: " << science_sources.size()
                      << " science archives, " << dq_sources.size()
                      << " DQ archives, corrected exposure order, no source copy staging."
                      << std::endl;
            if (!final_error.empty()) {
                std::cerr << "Initialization warning: " << final_error << std::endl;
            }
        } else {
            std::cerr << "Initialization failed: " << final_error << std::endl;
        }
    }

    MPI_Bcast(&final_status, 1, MPI_INT, 0, communicator);
    MPI_Barrier(communicator);
    return final_status;
}

}  // namespace fqinit
