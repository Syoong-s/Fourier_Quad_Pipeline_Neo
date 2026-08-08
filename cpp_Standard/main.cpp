#include "process_main/MPIScheduler.hpp"
#include "process_main/NumericalRecipes.hpp"
#include "process_main/ExternalCatalogReader.hpp"
#include "LensingConfig.hpp"
#include "ProcessConfig.hpp"
#include "ExtCatConfig.hpp"
#include "InitConfig.hpp"
#include "process_extcat/process_extcat.hpp"
#include "process_init/process_init.hpp"
#include "process_main/process_main.hpp"
#include "process_rearr/CatalogRearranger.hpp"
#include "process_rearr/process_rearr.hpp"
#include "process_fd/process_fd.hpp"
#include "FDConfig.hpp"
#include "ProcessRearrConfig.hpp"

#include <mpi.h>

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParserState {
    bool extcat_contains_option_seen = false;
    bool dataset_option_seen = false;
    bool legacy_dataset_option_seen = false;
    bool contains_option_seen = false;
};

// ==========================================
// Function: Parse one command-line boolean
// Method: Accept common textual and numeric forms without locale-dependent conversion.
// ==========================================
bool parseBoolean(const std::string& value, bool& parsed) {
    if (value == "true" || value == "1" || value == "on") {
        parsed = true;
        return true;
    }
    if (value == "false" || value == "0" || value == "off") {
        parsed = false;
        return true;
    }
    return false;
}

// ==========================================
// Function: Parse one non-negative command-line integer
// Method: Require full decimal consumption and guard the C++ int range.
// ==========================================
bool parseNonNegativeInteger(const std::string& value, int& parsed) {
    std::size_t consumed = 0;
    long long number = 0;
    try {
        number = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number < 0
        || number > std::numeric_limits<int>::max()) {
        return false;
    }
    parsed = static_cast<int>(number);
    return true;
}

// ==========================================
// Function: Parse one positive unsigned command-line integer
// Method: Require full decimal consumption and guard the uint64 range.
// ==========================================
bool parsePositiveUnsigned(const std::string& value, std::uint64_t& parsed) {
    if (value.empty() || value.front() == '-') {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number == 0
        || number > std::numeric_limits<std::uint64_t>::max()) {
        return false;
    }
    parsed = static_cast<std::uint64_t>(number);
    return true;
}

// ==========================================
// Function: Parse the external-catalog input column projection
// Method: Read one or more comma-separated positive one-based indices in output order.
// ==========================================
bool parseExtcatColumns(const std::string& value,
                        std::vector<std::size_t>& columns) {
    std::stringstream stream(value);
    std::string token;
    std::vector<std::size_t> parsed;
    while (std::getline(stream, token, ',')) {
        if (token.empty() || token.front() == '-') {
            return false;
        }
        std::size_t consumed = 0;
        unsigned long long number = 0;
        try {
            number = std::stoull(token, &consumed);
        } catch (const std::exception&) {
            return false;
        }
        if (consumed != token.size() || number == 0
            || number > std::numeric_limits<std::size_t>::max()) {
            return false;
        }
        parsed.push_back(static_cast<std::size_t>(number));
    }
    if (parsed.empty()) {
        return false;
    }
    columns = parsed;
    return true;
}

// ==========================================
// Function: Parse one external-catalog field column
// Method: Require one positive one-based RA, Dec, or ZP index within the platform size_t range.
// ==========================================
bool parseExtcatFieldColumn(const std::string& value, std::size_t& column) {
    if (value.empty() || value.front() == '-') {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number == 0
        || number > std::numeric_limits<std::size_t>::max()) {
        return false;
    }
    column = static_cast<std::size_t>(number);
    return true;
}

// ==========================================
// Function: Parse one paired dataset value
// Method: Split exactly one TARGET:PREFIX value and reject missing components.
// ==========================================
bool parseDataset(const std::string& value,
                  InitConfig::DatasetSpec& dataset,
                  std::string& error) {
    const std::size_t separator = value.find(':');
    if (separator == std::string::npos || separator == 0 || separator + 1 == value.size()
        || value.find(':', separator + 1) != std::string::npos) {
        error = "--dataset must use TARGET:PREFIX with both components non-empty";
        return false;
    }
    dataset.target = value.substr(0, separator);
    dataset.prefix = value.substr(separator + 1);
    return true;
}

// ==========================================
// Function: Enter legacy single-dataset override mode
// Method: Preserve the first configured pair, collapse the list to one entry,
//         and forbid mixing legacy target/prefix flags with --dataset.
// ==========================================
bool prepareLegacyDataset(ProcessConfig::RuntimeOptions& options,
                          ParserState& state,
                          std::string& error) {
    if (state.dataset_option_seen) {
        error = "--dataset cannot be combined with --target or --prefix";
        return false;
    }
    if (!state.legacy_dataset_option_seen) {
        InitConfig::DatasetSpec dataset;
        if (!options.datasets.empty()) {
            dataset = options.datasets.front();
        }
        options.datasets.assign(1, dataset);
        state.legacy_dataset_option_seen = true;
    }
    return true;
}

// ==========================================
// Function: Apply one normalized named option
// Method: Override scalar defaults, accumulate batch values, and validate typed
//         values immediately while tracking incompatible dataset syntaxes.
// ==========================================
bool applyNamedOption(const std::string& name,
                      const std::string& value,
                      ProcessConfig::RuntimeOptions& options,
                      ParserState& state,
                      std::string& error) {
    if (name == "--run-extcat") {
        if (!parseBoolean(value, options.run_process_extcat)) {
            error = "--run-extcat must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--extcat-input") {
        options.extcat_input_directory = value;
    } else if (name == "--extcat-output") {
        options.extcat_output_directory = value;
    } else if (name == "--extcat-contains") {
        if (value.empty()) {
            error = "--extcat-contains must not be empty";
            return false;
        }
        if (!state.extcat_contains_option_seen) {
            options.extcat_filename_tokens.clear();
            state.extcat_contains_option_seen = true;
        }
        options.extcat_filename_tokens.push_back(value);
    } else if (name == "--extcat-recursive") {
        if (!parseBoolean(value, options.extcat_recursive)) {
            error = "--extcat-recursive must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--extcat-delimiter") {
        if (value != "auto" && value != "whitespace" && value != "comma"
            && value != "tab") {
            error = "--extcat-delimiter must be auto, whitespace, comma, or tab";
            return false;
        }
        options.extcat_delimiter = value;
    } else if (name == "--extcat-header") {
        if (value != "auto" && value != "present" && value != "absent") {
            error = "--extcat-header must be auto, present, or absent";
            return false;
        }
        options.extcat_header_mode = value;
    } else if (name == "--extcat-columns") {
        if (!parseExtcatColumns(value, options.extcat_input_columns_one_based)) {
            error = "--extcat-columns must contain one or more positive one-based indices";
            return false;
        }
        options.extcat_use_explicit_columns = true;
    } else if (name == "--extcat-ra-column") {
        if (!parseExtcatFieldColumn(value, options.extcat_ra_column_one_based)) {
            error = "--extcat-ra-column must be a positive one-based index";
            return false;
        }
        options.extcat_use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-dec-column") {
        if (!parseExtcatFieldColumn(value, options.extcat_dec_column_one_based)) {
            error = "--extcat-dec-column must be a positive one-based index";
            return false;
        }
        options.extcat_use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-zp-column") {
        if (!parseExtcatFieldColumn(value, options.extcat_zp_column_one_based)) {
            error = "--extcat-zp-column must be a positive one-based index";
            return false;
        }
    } else if (name == "--extcat-chunk-mib") {
        if (!parsePositiveUnsigned(value, options.extcat_chunk_mib)) {
            error = "--extcat-chunk-mib must be a positive integer";
            return false;
        }
    } else if (name == "--extcat-malformed") {
        if (value != "fail" && value != "skip") {
            error = "--extcat-malformed must be fail or skip";
            return false;
        }
        options.extcat_malformed_policy = value;
    } else if (name == "--extcat-existing") {
        if (value != "fail" && value != "overwrite") {
            error = "--extcat-existing must be fail or overwrite";
            return false;
        }
        options.extcat_existing_policy = value;
    } else if (name == "--run-init") {
        if (!parseBoolean(value, options.run_process_init)) {
            error = "--run-init must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-main") {
        if (!parseBoolean(value, options.run_process_main)) {
            error = "--run-main must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-rearr") {
        if (!parseBoolean(value, options.run_process_rearr)) {
            error = "--run-rearr must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-fd") {
        if (!parseBoolean(value, options.run_process_fd)) {
            error = "--run-fd must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--science-root") {
        options.science_root = value;
    } else if (name == "--dq-root") {
        options.dq_root = value;
    } else if (name == "--output-root") {
        options.output_root = value;
    } else if (name == "--dataset") {
        if (state.legacy_dataset_option_seen) {
            error = "--dataset cannot be combined with --target or --prefix";
            return false;
        }
        InitConfig::DatasetSpec dataset;
        if (!parseDataset(value, dataset, error)) {
            return false;
        }
        if (!state.dataset_option_seen) {
            options.datasets.clear();
            state.dataset_option_seen = true;
        }
        options.datasets.push_back(dataset);
    } else if (name == "--target") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.datasets.front().target = value;
    } else if (name == "--prefix") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.datasets.front().prefix = value;
    } else if (name == "--contains") {
        if (value.empty()) {
            error = "--contains must not be empty";
            return false;
        }
        if (!state.contains_option_seen) {
            options.contains.clear();
            state.contains_option_seen = true;
        }
        options.contains.push_back(value);
    } else if (name == "--existing") {
        if (value != "fail" && value != "resume" && value != "overwrite") {
            error = "--existing must be fail, resume, or overwrite";
            return false;
        }
        options.existing = value;
    } else if (name == "--f77-max-path") {
        if (!parseNonNegativeInteger(value, options.f77_max_path)) {
            error = "--f77-max-path must be a non-negative integer";
            return false;
        }
    } else if (name == "--expo-list") {
        options.expo_list = value;
        options.external_expo_list_supplied = true;
    } else if (name == "--rearr-output-dir") {
        options.rearr_output_directory = value;
    } else if (name == "--rearr-output-base") {
        options.rearr_output_base_directory = value;
    } else if (name == "--rearr-list-name") {
        options.rearranged_expo_list_filename = value;
    } else if (name == "--rearr-list-dir") {
        options.rearranged_expo_list_directory = value;
    } else if (name == "--fd-expo-list") {
        options.fd_expo_list = value;
    } else if (name == "--fd-output-dir") {
        options.fd_output_directory = value;
    } else if (name == "--fd-output-base") {
        options.fd_output_base_directory = value;
    } else {
        error = "unknown option: " + name;
        return false;
    }
    return true;
}

// ==========================================
// Function: Validate the effective workflow configuration
// Method: Enforce executable modes, unique safe targets, paired prefixes, valid
//         token lists, and unambiguous exposure-list input for downstream-only runs.
// ==========================================
bool validateOptions(const ProcessConfig::RuntimeOptions& options, std::string& error) {
    if (options.help_requested) {
        return true;
    }
    if (!options.run_process_extcat && !options.run_process_init
        && !options.run_process_fd
        && !options.run_process_main && !options.run_process_rearr) {
        error = "--run-extcat, --run-init, --run-main, --run-rearr, and --run-fd cannot all be false";
        return false;
    }
    if ((options.run_process_extcat || options.run_process_main)
        && options.extcat_output_directory.empty()) {
        error = "external source-catalog output directory must not be empty";
        return false;
    }
    if (options.run_process_extcat && options.extcat_input_directory.empty()) {
        error = "external source-catalog input directory must not be empty";
        return false;
    }
    if ((options.run_process_init || options.run_process_main || options.run_process_fd
         || options.run_process_rearr)
        && options.datasets.empty()) {
        error = "at least one dataset must be configured or supplied with --dataset";
        return false;
    }

    std::set<std::string> targets;
    if (options.run_process_init || options.run_process_main || options.run_process_fd
        || options.run_process_rearr) {
        for (const InitConfig::DatasetSpec& dataset : options.datasets) {
            if (dataset.target.empty() || dataset.target == "." || dataset.target == ".."
                || dataset.target.find('/') != std::string::npos
                || dataset.target.find('\\') != std::string::npos) {
                error = "each dataset target must be one non-empty directory name";
                return false;
            }
            if (dataset.prefix.empty()) {
                error = "each dataset prefix must be non-empty";
                return false;
            }
            if (!targets.insert(dataset.target).second) {
                error = "dataset target is duplicated: " + dataset.target;
                return false;
            }
        }
    }

    if (options.run_process_init) {
        for (const std::string& token : options.contains) {
            if (token.empty()) {
                error = "contains tokens must be non-empty";
                return false;
            }
        }
    }
    if (options.run_process_extcat) {
        for (const std::string& token : options.extcat_filename_tokens) {
            if (token.empty()) {
                error = "external-catalog contains tokens must be non-empty";
                return false;
            }
        }
        if (options.extcat_use_explicit_columns
            && options.extcat_input_columns_one_based.empty()) {
            error = "external-catalog explicit column list must not be empty";
            return false;
        }
        if (options.extcat_use_explicit_coordinate_columns
            && (options.extcat_ra_column_one_based == 0
                || options.extcat_dec_column_one_based == 0
                || options.extcat_ra_column_one_based
                       == options.extcat_dec_column_one_based)) {
            error = "external-catalog RA and Dec columns must be distinct positive indices";
            return false;
        }
    }
    if (options.run_process_main) {
        ExternalCatalogReader::ColumnSelection selection;
        if (!ExternalCatalogReader::resolveColumnSelection(options, selection, error)) {
            return false;
        }
    }
    if (options.run_process_rearr) {
        ProcessRearr::CatalogLayout layout;
        if (!ProcessRearr::resolveCatalogLayout(options, layout, error)) {
            return false;
        }
    }

    if ((options.run_process_main || options.run_process_rearr)
        && !options.run_process_init
        && options.datasets.size() > 1 && !options.expo_list.empty()) {
        error = "one --expo-list cannot serve multiple datasets in downstream-only mode; "
                "omit it to derive expo_<target>.list for each dataset";
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the unified variable-length workflow command line
// Method: Accept --name value and --name=value in any order, accumulate dataset
//         and token options, and retain one positional exposure-list alias.
// ==========================================
bool parseCommandLine(int argc,
                      char** argv,
                      ProcessConfig::RuntimeOptions& options,
                      std::string& error) {
    std::string legacy_exposure_list;
    bool legacy_exposure_list_supplied = false;
    ParserState state;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            options.help_requested = true;
            continue;
        }

        if (argument.rfind("--", 0) != 0) {
            if (legacy_exposure_list_supplied) {
                error = "only one positional exposure-list compatibility argument is allowed";
                return false;
            }
            legacy_exposure_list = argument;
            legacy_exposure_list_supplied = true;
            continue;
        }

        const std::size_t equals = argument.find('=');
        const std::string name = argument.substr(0, equals);
        std::string value;
        if (equals != std::string::npos) {
            value = argument.substr(equals + 1);
        } else {
            if (index + 1 >= argc || std::string(argv[index + 1]).rfind("--", 0) == 0) {
                error = "missing value after " + name;
                return false;
            }
            value = argv[++index];
        }

        if (!applyNamedOption(name, value, options, state, error)) {
            return false;
        }
    }

    if (!options.external_expo_list_supplied && legacy_exposure_list_supplied) {
        options.expo_list = legacy_exposure_list;
        options.external_expo_list_supplied = true;
    }
    return validateOptions(options, error);
}

// ==========================================
// Function: Resolve one dataset exposure list used by downstream-only mode
// Method: Prefer the single external/default list, otherwise derive the current
//         dataset's expo_<target>.list and normalize it to an absolute path.
// ==========================================
std::string resolveExposureList(const ProcessConfig::RuntimeOptions& options,
                                const InitConfig::DatasetSpec& dataset) {
    std::filesystem::path path;
    if (!options.expo_list.empty()) {
        path = options.expo_list;
    } else {
        if (options.output_root.empty()) {
            throw std::runtime_error(
                "--expo-list is absent and output-root cannot derive its default");
        }
        path = std::filesystem::path(options.output_root)
               / ("expo_" + dataset.target + ".list");
    }
    return std::filesystem::weakly_canonical(std::filesystem::absolute(path)).string();
}

// ==========================================
// Function: Format configured dataset defaults for help output
// Method: Join every target/prefix pair without mutating the configured list.
// ==========================================
std::string configuredDatasetsText() {
    if (InitConfig::DATASETS.empty()) {
        return "none";
    }
    std::string text;
    for (const InitConfig::DatasetSpec& dataset : InitConfig::DATASETS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += dataset.target + ":" + dataset.prefix;
    }
    return text;
}

// ==========================================
// Function: Format configured token defaults for help output
// Method: Join every OR-matched basename token, or state that filtering is disabled.
// ==========================================
std::string configuredContainsText() {
    if (InitConfig::CONTAINS.empty()) {
        return "none (no token filter)";
    }
    std::string text;
    for (const std::string& token : InitConfig::CONTAINS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += token;
    }
    return text;
}

// ==========================================
// Function: Format external-catalog filename-token defaults
// Method: Join OR-matched raw-catalog basename tokens or state that all files match.
// ==========================================
std::string configuredExtcatContainsText() {
    if (ExtCatConfig::EXTCAT_FILENAME_TOKENS.empty()) {
        return "none (all files)";
    }
    std::string text;
    for (const std::string& token : ExtCatConfig::EXTCAT_FILENAME_TOKENS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += token;
    }
    return text;
}

// ==========================================
// Function: Print the unified workflow command-line contract
// Method: Describe batch/list accumulation, initializer values, list precedence,
//         compatibility input, and the configured defaults.
// ==========================================
void printUsage(const char* program_name) {
    std::cout
        << "Usage: " << program_name << " [options] [LEGACY_EXPO_LIST]\n"
        << "  --run-extcat BOOL     Repartition raw external catalogs first (default: "
        << (ProcessConfig::RUN_PROCESS_EXTCAT ? "true" : "false") << ")\n"
        << "  --run-init BOOL       Run initializer (default: "
        << (ProcessConfig::RUN_PROCESS_INIT ? "true" : "false") << ")\n"
        << "  --run-main BOOL       Run numerical pipeline (default: "
        << (ProcessConfig::RUN_PROCESS_MAIN ? "true" : "false") << ")\n"
        << "  --run-rearr BOOL      Rearrange _all.cat; follows process_main if both run (default: "
        << (ProcessConfig::RUN_PROCESS_REARR ? "true" : "false") << ")\n"
        << "  --run-fd BOOL         Run FD (field-distortion) shear test; follows process_main (default: "
        << (ProcessConfig::RUN_PROCESS_FD ? "true" : "false") << ")\n"
        << "  --extcat-input PATH   Directory containing raw external catalogs\n"
        << "  --extcat-output PATH  SOURCE_CAT tile directory used by process_main\n"
        << "  --extcat-contains T   Repeatable raw basename token, matched with OR (default: "
        << configuredExtcatContainsText() << ")\n"
        << "  --extcat-recursive B  Recurse below extcat input (default: "
        << (ExtCatConfig::EXTCAT_RECURSIVE ? "true" : "false") << ")\n"
        << "  --extcat-delimiter M  auto, whitespace, comma, or tab\n"
        << "  --extcat-header M     auto, present, or absent\n"
        << "  --extcat-columns LIST Ordered one-based input indices; output width follows LIST\n"
        << "  --extcat-ra-column N  Raw one-based RA column; overrides header discovery\n"
        << "  --extcat-dec-column N Raw one-based Dec column; overrides header discovery\n"
        << "  --extcat-zp-column N  Raw one-based ZP column consumed by process_main\n"
        << "  --extcat-chunk-mib N  MPI byte-range task size in MiB (default: "
        << ExtCatConfig::EXTCAT_CHUNK_MIB << ")\n"
        << "  --extcat-malformed P  fail or skip malformed rows\n"
        << "  --extcat-existing P   fail or overwrite generated tiles\n"
        << "  --science-root PATH   Original Science FITS/FZ repository\n"
        << "  --dq-root PATH        Original DQ FITS/FZ repository\n"
        << "  --output-root PATH    Parent of targets and generated exposure lists\n"
        << "  --dataset T:P         Repeatable TARGET:PREFIX dataset pair (default: "
        << configuredDatasetsText() << ")\n"
        << "  --target NAME         Legacy single-dataset target; do not mix with --dataset\n"
        << "  --prefix TEXT         Legacy single-dataset prefix; do not mix with --dataset\n"
        << "  --contains TEXT       Repeatable basename token, matched with OR (default: "
        << configuredContainsText() << ")\n"
        << "  --existing MODE       fail, resume, or overwrite (default: "
        << InitConfig::EXISTING << ")\n"
        << "  --f77-max-path N      Generated path limit; zero disables (default: "
        << InitConfig::F77_MAX_PATH << ")\n"
        << "  --expo-list PATH      Single exposure list for main/rearr-only mode\n"
        << "  --rearr-output-dir D  process_rearr output sub-directory (default: baked)\n"
        << "  --rearr-output-base P process_rearr output base path (default: dataset root)\n"
        << "  --rearr-list-name F   Rearranged expo-list filename (default: expo_rearranged.list)\n"
        << "  --rearr-list-dir P    Rearranged expo-list directory (default: expo-list parent)\n"
        << "  --fd-expo-list PATH   process_fd exposure list (default: rearranged list)\n"
        << "  --fd-output-dir D    process_fd output sub-directory (default: fdout)\n"
        << "  --fd-output-base P    process_fd output base path (default: dataset root)\n"
        << "  --help                Show this help\n"
        << "Options accept both --name value and --name=value. The first explicit "
           "--dataset, --contains, or --extcat-contains replaces its corresponding "
           "configured list; repeats append.\n"
        << "Other duplicate scalar options use the last value. Downstream-only batches derive "
           "expo_<target>.list per dataset when --expo-list is omitted.\n"
        << "Without --extcat-columns, all raw catalog fields keep their original order; "
           "otherwise output width and order follow LIST exactly.\n"
        << "When process_main runs, LIST must contain the configured RA, Dec, and ZP "
           "raw columns. A rearr-only run requires RA and Dec but not ZP; output "
           "positions follow LIST order.\n"
        << "process_extcat runs once before the dataset loop. When later phases run, each "
           "process_init-generated absolute exposure-list "
           "path overrides external input for its dataset.\n";
}


// ==========================================
// Function: Derive dataset root from the first image path in an expo list
// Method: Read the first per-exposure list entry, open it, read the first
//         non-empty image line, and compute the great-grandparent directory
//         (getDir level 3), matching process_main's dir_output derivation.
// ==========================================
std::string deriveDatasetRootFromExpoList(const std::string& exposure_list) {
    std::ifstream expo_input(exposure_list);
    if (!expo_input.is_open()) {
        return "";
    }
    std::string path;
    int chip_count = 0;
    while (expo_input >> path >> chip_count) {
        if (path.size() >= 2 && path.front() == '"' && path.back() == '"') {
            path = path.substr(1, path.size() - 2);
        }
        std::ifstream list_input(path);
        if (!list_input.is_open()) {
            continue;
        }
        std::string line;
        while (std::getline(list_input, line)) {
            std::size_t first = 0;
            while (first < line.size()
                   && (line[first] == ' ' || line[first] == '\t'
                       || line[first] == '\r' || line[first] == '\n')) {
                ++first;
            }
            std::size_t last = line.size();
            while (last > first
                   && (line[last - 1] == ' ' || line[last - 1] == '\t'
                       || line[last - 1] == '\r' || line[last - 1] == '\n')) {
                --last;
            }
            if (first == last) {
                continue;
            }
            line = line.substr(first, last - first);
            if (line.size() >= 2 && line.front() == '"' && line.back() == '"') {
                line = line.substr(1, line.size() - 2);
            }
            if (line.empty()) {
                continue;
            }
            const std::filesystem::path image_path(line);
            const std::filesystem::path parent = image_path.parent_path();
            const std::filesystem::path grandparent = parent.parent_path();
            const std::filesystem::path great_grandparent = grandparent.parent_path();
            if (parent.empty() || grandparent.empty() || great_grandparent.empty()) {
                return "";
            }
            return std::filesystem::absolute(great_grandparent)
                .lexically_normal().string();
        }
    }
    return "";
}

}  // namespace

// ==========================================
// Function: Dispatch all four external-catalog, initializer, Fourier_Quad, and rearrangement phases
// Method: Own MPI exactly once, run process_extcat before the dataset loop, process datasets
//         sequentially, and invoke process_rearr only after any enabled process_main call.
// ==========================================
int main(int argc, char* argv[]) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::my_id;
    if (rank == 0){
        std::cout << "MPI Init Done..." << std::endl;
    }
    int return_code = 0;

    ProcessConfig::RuntimeOptions options;
    std::string parse_error;
    const int local_parse_ok = parseCommandLine(argc, argv, options, parse_error) ? 1 : 0;
    int global_parse_ok = 0;
    MPI_Allreduce(&local_parse_ok, &global_parse_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (global_parse_ok == 0) {
        if (rank == 0) {
            std::cerr << "Argument error: "
                      << (parse_error.empty() ? "parsing failed on another rank" : parse_error)
                      << std::endl;
            printUsage(argv[0]);
        }
        return_code = 2;
    } else if (options.help_requested) {
        if (rank == 0) {
            printUsage(argv[0]);
        }
    } else {
        LensingConfig::SOURCE_CAT = options.extcat_output_directory;
        if (options.run_process_extcat) {
            if (rank == 0) {
                std::cout << "Running process_extcat before all dataset phases" << std::endl;
            }
            return_code = process_extcat(options, MPI_COMM_WORLD);
            if (return_code == 0) {
                MPIScheduler::barrier();
            }
        }

        bool rng_initialized = false;
        for (std::size_t index = 0;
             index < options.datasets.size() && return_code == 0
                 && (options.run_process_init || options.run_process_main || options.run_process_fd
                     || options.run_process_rearr);
             ++index) {
            const InitConfig::DatasetSpec& dataset = options.datasets[index];
            if (rank == 0) {
                std::cout << "Dataset " << (index + 1) << "/" << options.datasets.size()
                          << ": target=" << dataset.target
                          << " prefix=" << dataset.prefix << std::endl;
            }

            std::string generated_exposure_list;
            if (options.run_process_init) {
                return_code = process_init(options, dataset, generated_exposure_list);
            }

            std::string selected_exposure_list;
            if (return_code == 0
                && (options.run_process_main || options.run_process_rearr || options.run_process_fd)) {
                if (options.run_process_init) {
                    selected_exposure_list = generated_exposure_list;
                    if (rank == 0) {
                        if (!options.expo_list.empty()) {
                            std::cout << "External exposure list overridden by process_init output: ";
                        } else {
                            std::cout << "Downstream phases will use process_init output: ";
                        }
                        std::cout << selected_exposure_list << std::endl;
                    }
                    MPIScheduler::barrier();
                } else {
                    int local_path_ok = 1;
                    std::string path_error;
                    try {
                        selected_exposure_list = resolveExposureList(options, dataset);
                    } catch (const std::exception& exception) {
                        local_path_ok = 0;
                        path_error = exception.what();
                    }
                    int global_path_ok = 0;
                    MPI_Allreduce(&local_path_ok, &global_path_ok, 1, MPI_INT, MPI_MIN,
                                  MPI_COMM_WORLD);
                    if (global_path_ok == 0) {
                        if (rank == 0) {
                            std::cerr << "Exposure-list argument error: "
                                      << (path_error.empty()
                                              ? "path resolution failed on another rank"
                                              : path_error)
                                      << std::endl;
                        }
                        return_code = 2;
                    }
                }

            }

            if (return_code == 0 && options.run_process_main) {
                if (!rng_initialized) {
                    const unsigned int rng_seed = NumericalRecipes::initializeRan1Seed(
                        rank, MPIScheduler::num_procs);
                    std::cout << "RNG_SEED rank seed: " << rank << " " << rng_seed << std::endl;
                    MPIScheduler::barrier();
                    rng_initialized = true;
                }
                if (return_code == 0) {
                    return_code = process_main(selected_exposure_list, options);
                }
            }

            if (return_code == 0 && options.run_process_rearr) {
                MPIScheduler::barrier();
                if (rank == 0) {
                    std::cout << "Running process_rearr"
                              << (options.run_process_main ? " after process_main" : "")
                              << std::endl;
                }
                return_code = process_rearr(selected_exposure_list, options,
                                            MPI_COMM_WORLD);
            }

            if (return_code == 0 && options.run_process_fd) {
                MPIScheduler::barrier();
                if (rank == 0) {
                    std::cout << "Running process_fd" << std::endl;
                }
                if (!rng_initialized) {
                    NumericalRecipes::initializeRan1Seed(rank, MPIScheduler::num_procs);
                    MPIScheduler::barrier();
                    rng_initialized = true;
                }
                std::string fd_expo_list;
                if (options.fd_expo_list.empty()) {
                    const std::filesystem::path list_dir =
                        options.rearranged_expo_list_directory.empty()
                            ? std::filesystem::path(selected_exposure_list)
                                  .parent_path()
                            : std::filesystem::path(
                                  options.rearranged_expo_list_directory);
                    fd_expo_list =
                        std::filesystem::absolute(
                            list_dir / options.rearranged_expo_list_filename)
                            .lexically_normal()
                            .string();
                } else {
                    fd_expo_list = options.fd_expo_list;
                }
                const std::string fd_dataset_root =
                    deriveDatasetRootFromExpoList(selected_exposure_list);
                if (rank == 0) {
                    std::cout << "FD expo list: " << fd_expo_list
                              << "  dataset_root: " << fd_dataset_root << std::endl;
                }
                return_code = process_fd(fd_expo_list, options,
                                         fd_dataset_root, MPI_COMM_WORLD);
            }

            if (return_code == 0 && index + 1 < options.datasets.size()) {
                MPIScheduler::barrier();
            }
        }
    }

    MPIScheduler::finalize();
    return return_code;
}
