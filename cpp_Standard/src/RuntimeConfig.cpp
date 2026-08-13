#include "RuntimeConfig.hpp"

#include "ExtCatConfig.hpp"
#include "LensingConfig.hpp"
#include "ProcessConfig.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

struct CommandLineState {
    bool extcat_contains_seen = false;
    bool dataset_seen = false;
    bool legacy_dataset_seen = false;
    bool contains_seen = false;
};

std::unique_ptr<const RuntimeConfig> stored_config;

// ==========================================
// Function: trim
// Method: Remove leading and trailing ASCII whitespace without locale-specific
//         transformations or changing embedded path characters.
// ==========================================
std::string trim(const std::string& value) {
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
// Function: lowercase
// Method: Normalize INI section/key names while leaving user values untouched.
// ==========================================
std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    return value;
}

// ==========================================
// Function: stripInlineComment
// Method: Remove unquoted # or ; comments and preserve either quoted string form.
// ==========================================
std::string stripInlineComment(const std::string& value) {
    char quote = '\0';
    bool escaped = false;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const char character = value[index];
        if (escaped) {
            escaped = false;
            continue;
        }
        if (character == '\\' && quote != '\0') {
            escaped = true;
            continue;
        }
        if (character == '\'' || character == '"') {
            if (quote == '\0') {
                quote = character;
            } else if (quote == character) {
                quote = '\0';
            }
            continue;
        }
        if (quote == '\0' && (character == '#' || character == ';')) {
            return value.substr(0, index);
        }
    }
    return value;
}

// ==========================================
// Function: decodeString
// Method: Trim one INI value, remove a matching quote pair, and decode only the
//         quote/backslash escapes needed for paths with delimiters or comments.
// ==========================================
bool decodeString(const std::string& input, std::string& output,
                  std::string& reason) {
    const std::string value = trim(stripInlineComment(input));
    if (value.empty()) {
        output.clear();
        return true;
    }
    const bool starts_quoted = value.front() == '\'' || value.front() == '"';
    if (!starts_quoted) {
        output = value;
        return true;
    }
    if (value.size() < 2 || value.back() != value.front()) {
        reason = "unterminated quoted string";
        return false;
    }
    output.clear();
    const char quote = value.front();
    bool escaped = false;
    for (std::size_t index = 1; index + 1 < value.size(); ++index) {
        const char character = value[index];
        if (escaped) {
            if (character != quote && character != '\\') {
                output.push_back('\\');
            }
            output.push_back(character);
            escaped = false;
        } else if (character == '\\') {
            escaped = true;
        } else {
            output.push_back(character);
        }
    }
    if (escaped) {
        reason = "quoted string ends with an incomplete escape";
        return false;
    }
    return true;
}

// ==========================================
// Function: splitList
// Method: Split comma-separated INI lists outside quotes and decode every item.
// ==========================================
bool splitList(const std::string& input, std::vector<std::string>& items,
               std::string& reason) {
    const std::string value = trim(stripInlineComment(input));
    items.clear();
    if (value.empty()) {
        return true;
    }
    char quote = '\0';
    bool escaped = false;
    std::size_t start = 0;
    for (std::size_t index = 0; index <= value.size(); ++index) {
        const bool at_end = index == value.size();
        const char character = at_end ? ',' : value[index];
        if (!at_end && escaped) {
            escaped = false;
            continue;
        }
        if (!at_end && character == '\\' && quote != '\0') {
            escaped = true;
            continue;
        }
        if (!at_end && (character == '\'' || character == '"')) {
            if (quote == '\0') {
                quote = character;
            } else if (quote == character) {
                quote = '\0';
            }
            continue;
        }
        if (character == ',' && quote == '\0') {
            std::string decoded;
            if (!decodeString(value.substr(start, index - start), decoded,
                              reason)) {
                return false;
            }
            if (decoded.empty()) {
                reason = "list entries must not be empty";
                return false;
            }
            items.push_back(std::move(decoded));
            start = index + 1;
        }
    }
    if (quote != '\0') {
        reason = "unterminated quote in list";
        return false;
    }
    return true;
}

// ==========================================
// Function: parseBoolean
// Method: Accept the documented numeric and common textual boolean spellings.
// ==========================================
bool parseBoolean(const std::string& input, bool& parsed) {
    std::string value;
    std::string reason;
    if (!decodeString(input, value, reason)) {
        return false;
    }
    value = lowercase(value);
    if (value == "true" || value == "1" || value == "yes" || value == "on") {
        parsed = true;
        return true;
    }
    if (value == "false" || value == "0" || value == "no" || value == "off") {
        parsed = false;
        return true;
    }
    return false;
}

template <typename Integer>
bool parseSignedInteger(const std::string& input, Integer& parsed,
                        bool positive_only, bool non_negative) {
    std::string value;
    std::string reason;
    if (!decodeString(input, value, reason) || value.empty()) {
        return false;
    }
    std::size_t consumed = 0;
    long long number = 0;
    try {
        number = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || (positive_only && number <= 0)
        || (non_negative && number < 0)
        || number < static_cast<long long>(std::numeric_limits<Integer>::min())
        || number > static_cast<long long>(std::numeric_limits<Integer>::max())) {
        return false;
    }
    parsed = static_cast<Integer>(number);
    return true;
}

template <typename Integer>
bool parseUnsignedInteger(const std::string& input, Integer& parsed,
                          bool positive_only) {
    std::string value;
    std::string reason;
    if (!decodeString(input, value, reason) || value.empty()
        || value.front() == '-') {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || (positive_only && number == 0)
        || number > static_cast<unsigned long long>(
                        std::numeric_limits<Integer>::max())) {
        return false;
    }
    parsed = static_cast<Integer>(number);
    return true;
}

// ==========================================
// Function: parseDouble
// Method: Require full conversion and reject non-finite floating-point values.
// ==========================================
bool parseDouble(const std::string& input, double& parsed) {
    std::string value;
    std::string reason;
    if (!decodeString(input, value, reason) || value.empty()) {
        return false;
    }
    std::size_t consumed = 0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    return consumed == value.size() && std::isfinite(parsed);
}

// ==========================================
// Function: parseDatasetValue
// Method: Split exactly one TARGET:PREFIX pair with two non-empty components.
// ==========================================
bool parseDatasetValue(const std::string& value,
                       InitConfig::DatasetSpec& dataset,
                       std::string& reason) {
    const std::size_t separator = value.find(':');
    if (separator == std::string::npos || separator == 0
        || separator + 1 == value.size()
        || value.find(':', separator + 1) != std::string::npos) {
        reason = "dataset must use TARGET:PREFIX with both components non-empty";
        return false;
    }
    dataset.target = value.substr(0, separator);
    dataset.prefix = value.substr(separator + 1);
    return true;
}

// ==========================================
// Function: parseDatasets
// Method: Parse a comma-separated dataset list into initializer-compatible pairs.
// ==========================================
bool parseDatasets(const std::string& input,
                   std::vector<InitConfig::DatasetSpec>& datasets,
                   std::string& reason) {
    std::vector<std::string> values;
    if (!splitList(input, values, reason)) {
        return false;
    }
    std::vector<InitConfig::DatasetSpec> parsed;
    for (const std::string& value : values) {
        InitConfig::DatasetSpec dataset;
        if (!parseDatasetValue(value, dataset, reason)) {
            return false;
        }
        parsed.push_back(std::move(dataset));
    }
    datasets = std::move(parsed);
    return true;
}

// ==========================================
// Function: parseColumnList
// Method: Convert a non-empty comma-separated list to positive one-based indices.
// ==========================================
bool parseColumnList(const std::string& input,
                     std::vector<std::size_t>& columns,
                     std::string& reason) {
    std::vector<std::string> values;
    if (!splitList(input, values, reason) || values.empty()) {
        reason = "expected one or more positive one-based indices";
        return false;
    }
    std::vector<std::size_t> parsed;
    for (const std::string& value : values) {
        std::size_t column = 0;
        if (!parseUnsignedInteger(value, column, true)) {
            reason = "expected one or more positive one-based indices";
            return false;
        }
        parsed.push_back(column);
    }
    columns = std::move(parsed);
    return true;
}

// ==========================================
// Function: formatConfigError
// Method: Produce a stable source/section/key/value diagnostic for one INI error.
// ==========================================
std::string formatConfigError(const std::string& source_name,
                              const std::string& section,
                              const std::string& key,
                              const std::string& value,
                              const std::string& reason) {
    std::ostringstream output;
    output << "file: " << source_name << "\n"
           << "  section: " << (section.empty() ? "<none>" : section) << "\n"
           << "  key: " << (key.empty() ? "<none>" : key) << "\n"
           << "  value: " << value << "\n"
           << "  reason: " << reason;
    return output.str();
}

// ==========================================
// Function: applyIniValue
// Method: Route one normalized section/key to its typed runtime field and
//         report the expected domain for unknown or malformed values.
// ==========================================
bool applyIniValue(const std::string& section, const std::string& key,
                   const std::string& raw_value, RuntimeConfig& config,
                   std::string& reason) {
    std::string string_value;
    const auto parse_string = [&]() {
        return decodeString(raw_value, string_value, reason);
    };

    if (section == "process") {
        if (key == "run_process_extcat") return parseBoolean(raw_value, config.process.run_process_extcat) || (reason = "expected boolean", false);
        if (key == "run_process_init") return parseBoolean(raw_value, config.process.run_process_init) || (reason = "expected boolean", false);
        if (key == "run_process_main") return parseBoolean(raw_value, config.process.run_process_main) || (reason = "expected boolean", false);
        if (key == "run_process_rearr") return parseBoolean(raw_value, config.process.run_process_rearr) || (reason = "expected boolean", false);
        if (key == "run_process_fd") return parseBoolean(raw_value, config.process.run_process_fd) || (reason = "expected boolean", false);
        if (!parse_string()) return false;
        if (key == "expo_list") config.process.expo_list = string_value;
        else if (key == "rearr_output_directory") config.process.rearr_output_directory = string_value;
        else if (key == "rearr_output_base_directory") config.process.rearr_output_base_directory = string_value;
        else if (key == "rearranged_expo_list_filename") config.process.rearranged_expo_list_filename = string_value;
        else if (key == "rearranged_expo_list_directory") config.process.rearranged_expo_list_directory = string_value;
        else if (key == "fd_expo_list") config.process.fd_expo_list = string_value;
        else if (key == "fd_output_directory") config.process.fd_output_directory = string_value;
        else if (key == "fd_output_base_directory") config.process.fd_output_base_directory = string_value;
        else { reason = "unknown key"; return false; }
        return true;
    }

    if (section == "extcat") {
        if (key == "recursive") return parseBoolean(raw_value, config.extcat.recursive) || (reason = "expected boolean", false);
        if (key == "use_explicit_columns") return parseBoolean(raw_value, config.extcat.use_explicit_columns) || (reason = "expected boolean", false);
        if (key == "use_explicit_coordinate_columns") return parseBoolean(raw_value, config.extcat.use_explicit_coordinate_columns) || (reason = "expected boolean", false);
        if (key == "chunk_mib") return parseUnsignedInteger(raw_value, config.extcat.chunk_mib, true) || (reason = "expected positive integer", false);
        if (key == "total_columns") return parseUnsignedInteger(raw_value, config.extcat.total_columns, true) || (reason = "expected positive integer", false);
        if (key == "ra_column") return parseUnsignedInteger(raw_value, config.extcat.ra_column_one_based, true) || (reason = "expected positive one-based index", false);
        if (key == "dec_column") return parseUnsignedInteger(raw_value, config.extcat.dec_column_one_based, true) || (reason = "expected positive one-based index", false);
        if (key == "zp_column") return parseUnsignedInteger(raw_value, config.extcat.zp_column_one_based, true) || (reason = "expected positive one-based index", false);
        if (key == "mag_g_column") return parseUnsignedInteger(raw_value, config.extcat.mag_g_column_one_based, false) || (reason = "expected non-negative one-based index or zero", false);
        if (key == "mag_r_column") return parseUnsignedInteger(raw_value, config.extcat.mag_r_column_one_based, false) || (reason = "expected non-negative one-based index or zero", false);
        if (key == "mag_i_column") return parseUnsignedInteger(raw_value, config.extcat.mag_i_column_one_based, false) || (reason = "expected non-negative one-based index or zero", false);
        if (key == "mag_z_column") return parseUnsignedInteger(raw_value, config.extcat.mag_z_column_one_based, false) || (reason = "expected non-negative one-based index or zero", false);
        if (key == "mag_y_column") return parseUnsignedInteger(raw_value, config.extcat.mag_y_column_one_based, false) || (reason = "expected non-negative one-based index or zero", false);
        if (key == "input_columns") return parseColumnList(raw_value, config.extcat.input_columns_one_based, reason);
        if (key == "filename_tokens") return splitList(raw_value, config.extcat.filename_tokens, reason);
        if (!parse_string()) return false;
        if (key == "input_directory") config.extcat.input_directory = string_value;
        else if (key == "output_directory") config.extcat.output_directory = string_value;
        else if (key == "delimiter") config.extcat.delimiter = lowercase(string_value);
        else if (key == "header_mode") config.extcat.header_mode = lowercase(string_value);
        else if (key == "malformed_policy") config.extcat.malformed_policy = lowercase(string_value);
        else if (key == "existing_policy") config.extcat.existing_policy = lowercase(string_value);
        else { reason = "unknown key"; return false; }
        return true;
    }

    if (section == "init") {
        if (key == "datasets") return parseDatasets(raw_value, config.init.datasets, reason);
        if (key == "contains") return splitList(raw_value, config.init.contains, reason);
        if (key == "f77_max_path") return parseSignedInteger(raw_value, config.init.f77_max_path, false, true) || (reason = "expected non-negative integer", false);
        if (!parse_string()) return false;
        if (key == "science_root") config.init.science_root = string_value;
        else if (key == "dq_root") config.init.dq_root = string_value;
        else if (key == "output_root") config.init.output_root = string_value;
        else if (key == "existing") config.init.existing = lowercase(string_value);
        else { reason = "unknown key"; return false; }
        return true;
    }

    if (section == "lensing") {
        if (key == "astrometry_trivial") return parseSignedInteger(raw_value, config.lensing.astrometry_trivial, false, false) || (reason = "expected integer", false);
        if (key == "process_stage") return parseSignedInteger(raw_value, config.lensing.process_stage, true, false) || (reason = "expected positive integer", false);
        if (key == "include_flat") return parseSignedInteger(raw_value, config.lensing.include_flat, false, false) || (reason = "expected integer", false);
        if (key == "include_mask") return parseSignedInteger(raw_value, config.lensing.include_mask, false, false) || (reason = "expected integer", false);
        if (key == "ext_cat") return parseSignedInteger(raw_value, config.lensing.ext_cat, false, false) || (reason = "expected integer", false);
        if (key == "ext_psf") return parseSignedInteger(raw_value, config.lensing.ext_psf, false, false) || (reason = "expected integer", false);
        if (key == "ccd_split") return parseSignedInteger(raw_value, config.lensing.ccd_split, false, false) || (reason = "expected integer", false);
        if (key == "deblending") return parseSignedInteger(raw_value, config.lensing.deblending, false, false) || (reason = "expected integer", false);
        if (key == "psf_type") return parseSignedInteger(raw_value, config.lensing.psf_type, false, false) || (reason = "expected integer", false);
        if (key == "psf_ms") return parseSignedInteger(raw_value, config.lensing.psf_ms, false, false) || (reason = "expected integer", false);
        if (key == "gal_smooth") return parseSignedInteger(raw_value, config.lensing.gal_smooth, false, true) || (reason = "expected non-negative integer", false);
        if (key == "star_smooth") return parseSignedInteger(raw_value, config.lensing.star_smooth, false, true) || (reason = "expected non-negative integer", false);
        if (key == "nmax_chip") return parseSignedInteger(raw_value, config.lensing.nmax_chip, true, false) || (reason = "expected positive integer", false);
        if (key == "chipnx") return parseSignedInteger(raw_value, config.lensing.chipnx, true, false) || (reason = "expected positive integer", false);
        if (key == "chipny") return parseSignedInteger(raw_value, config.lensing.chipny, true, false) || (reason = "expected positive integer", false);
        if (key == "pixel_size") return parseDouble(raw_value, config.lensing.pixel_size) || (reason = "expected finite positive number", false);
        if (!parse_string()) return false;
        if (key == "astrometry_cat") config.lensing.astrometry_cat = string_value;
        else if (key == "source_cat") config.extcat.output_directory = string_value;
        else if (key == "flat_path") config.lensing.flat_path = string_value;
        else if (key == "psf_path") config.lensing.psf_path = string_value;
        else { reason = "unknown key"; return false; }
        return true;
    }

    reason = "unknown section";
    return false;
}

// ==========================================
// Function: prepareLegacyDataset
// Method: Collapse configured datasets to one pair for --target/--prefix while
//         forbidding mixed legacy and repeatable dataset syntaxes.
// ==========================================
bool prepareLegacyDataset(RuntimeConfig& config, CommandLineState& state,
                          std::string& error) {
    if (state.dataset_seen) {
        error = "--dataset cannot be combined with --target or --prefix";
        return false;
    }
    if (!state.legacy_dataset_seen) {
        InitConfig::DatasetSpec dataset;
        if (!config.init.datasets.empty()) {
            dataset = config.init.datasets.front();
        }
        config.init.datasets.assign(1, dataset);
        state.legacy_dataset_seen = true;
    }
    return true;
}

// ==========================================
// Function: applyNamedOption
// Method: Apply one normalized CLI option to the nested runtime config and
//         preserve existing list replacement/append semantics.
// ==========================================
bool applyNamedOption(const std::string& name, const std::string& value,
                      RuntimeConfig& config, CommandLineState& state,
                      std::string& error) {
    if (name == "--run-extcat") {
        if (!parseBoolean(value, config.process.run_process_extcat)) error = "--run-extcat must be true, false, 1, 0, yes, no, on, or off";
    } else if (name == "--run-init") {
        if (!parseBoolean(value, config.process.run_process_init)) error = "--run-init must be a boolean";
    } else if (name == "--run-main") {
        if (!parseBoolean(value, config.process.run_process_main)) error = "--run-main must be a boolean";
    } else if (name == "--run-rearr") {
        if (!parseBoolean(value, config.process.run_process_rearr)) error = "--run-rearr must be a boolean";
    } else if (name == "--run-fd") {
        if (!parseBoolean(value, config.process.run_process_fd)) error = "--run-fd must be a boolean";
    } else if (name == "--extcat-input") {
        config.extcat.input_directory = value;
    } else if (name == "--extcat-output") {
        config.extcat.output_directory = value;
    } else if (name == "--extcat-contains") {
        if (value.empty()) error = "--extcat-contains must not be empty";
        else {
            if (!state.extcat_contains_seen) {
                config.extcat.filename_tokens.clear();
                state.extcat_contains_seen = true;
            }
            config.extcat.filename_tokens.push_back(value);
        }
    } else if (name == "--extcat-recursive") {
        if (!parseBoolean(value, config.extcat.recursive)) error = "--extcat-recursive must be a boolean";
    } else if (name == "--extcat-delimiter") {
        config.extcat.delimiter = lowercase(value);
    } else if (name == "--extcat-header") {
        config.extcat.header_mode = lowercase(value);
    } else if (name == "--extcat-columns") {
        std::string reason;
        if (!parseColumnList(value, config.extcat.input_columns_one_based, reason)) error = "--extcat-columns " + reason;
        else config.extcat.use_explicit_columns = true;
    } else if (name == "--extcat-ra-column") {
        if (!parseUnsignedInteger(value, config.extcat.ra_column_one_based, true)) error = "--extcat-ra-column must be a positive one-based index";
        else config.extcat.use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-dec-column") {
        if (!parseUnsignedInteger(value, config.extcat.dec_column_one_based, true)) error = "--extcat-dec-column must be a positive one-based index";
        else config.extcat.use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-zp-column") {
        if (!parseUnsignedInteger(value, config.extcat.zp_column_one_based, true)) error = "--extcat-zp-column must be a positive one-based index";
    } else if (name == "--extcat-chunk-mib") {
        if (!parseUnsignedInteger(value, config.extcat.chunk_mib, true)) error = "--extcat-chunk-mib must be a positive integer";
    } else if (name == "--extcat-malformed") {
        config.extcat.malformed_policy = lowercase(value);
    } else if (name == "--extcat-existing") {
        config.extcat.existing_policy = lowercase(value);
    } else if (name == "--science-root") {
        config.init.science_root = value;
    } else if (name == "--dq-root") {
        config.init.dq_root = value;
    } else if (name == "--output-root") {
        config.init.output_root = value;
    } else if (name == "--dataset") {
        if (state.legacy_dataset_seen) error = "--dataset cannot be combined with --target or --prefix";
        else {
            InitConfig::DatasetSpec dataset;
            std::string reason;
            if (!parseDatasetValue(value, dataset, reason)) error = "--" + reason;
            else {
                if (!state.dataset_seen) {
                    config.init.datasets.clear();
                    state.dataset_seen = true;
                }
                config.init.datasets.push_back(std::move(dataset));
            }
        }
    } else if (name == "--target") {
        if (prepareLegacyDataset(config, state, error)) config.init.datasets.front().target = value;
    } else if (name == "--prefix") {
        if (prepareLegacyDataset(config, state, error)) config.init.datasets.front().prefix = value;
    } else if (name == "--contains") {
        if (value.empty()) error = "--contains must not be empty";
        else {
            if (!state.contains_seen) {
                config.init.contains.clear();
                state.contains_seen = true;
            }
            config.init.contains.push_back(value);
        }
    } else if (name == "--existing") {
        config.init.existing = lowercase(value);
    } else if (name == "--f77-max-path") {
        if (!parseSignedInteger(value, config.init.f77_max_path, false, true)) error = "--f77-max-path must be a non-negative integer";
    } else if (name == "--expo-list") {
        config.process.expo_list = value;
        config.external_expo_list_supplied = true;
    } else if (name == "--rearr-output-dir") {
        config.process.rearr_output_directory = value;
    } else if (name == "--rearr-output-base") {
        config.process.rearr_output_base_directory = value;
    } else if (name == "--rearr-list-name") {
        config.process.rearranged_expo_list_filename = value;
    } else if (name == "--rearr-list-dir") {
        config.process.rearranged_expo_list_directory = value;
    } else if (name == "--fd-expo-list") {
        config.process.fd_expo_list = value;
    } else if (name == "--fd-output-dir") {
        config.process.fd_output_directory = value;
    } else if (name == "--fd-output-base") {
        config.process.fd_output_base_directory = value;
    } else {
        error = "unknown option: " + name;
    }
    return error.empty();
}

// ==========================================
// Function: valueInSet
// Method: Test one small integer branch selector against its implemented domain.
// ==========================================
bool valueInSet(int value, std::initializer_list<int> allowed) {
    return std::find(allowed.begin(), allowed.end(), value) != allowed.end();
}

}  // namespace

RuntimeConfig makeDefaultRuntimeConfig() {
    RuntimeConfig config;
    config.process.run_process_extcat = ProcessConfig::RUN_PROCESS_EXTCAT;
    config.process.run_process_init = ProcessConfig::RUN_PROCESS_INIT;
    config.process.run_process_main = ProcessConfig::RUN_PROCESS_MAIN;
    config.process.run_process_rearr = ProcessConfig::RUN_PROCESS_REARR;
    config.process.run_process_fd = ProcessConfig::RUN_PROCESS_FD;
    config.process.expo_list = ProcessConfig::EXPO_LIST;
    config.process.rearr_output_directory = ProcessConfig::REARR_OUTPUT_DIRECTORY;
    config.process.rearr_output_base_directory = ProcessConfig::REARR_OUTPUT_BASE_DIRECTORY;
    config.process.rearranged_expo_list_filename = ProcessConfig::REARRANGED_EXPO_LIST_FILENAME;
    config.process.rearranged_expo_list_directory = ProcessConfig::REARRANGED_EXPO_LIST_DIRECTORY;
    config.process.fd_expo_list = ProcessConfig::FD_EXPO_LIST;
    config.process.fd_output_directory = ProcessConfig::FD_OUTPUT_DIRECTORY;
    config.process.fd_output_base_directory = ProcessConfig::FD_OUTPUT_BASE_DIRECTORY;

    config.extcat.input_directory = ExtCatConfig::EXTCAT_INPUT_DIRECTORY;
    config.extcat.output_directory = ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY;
    config.extcat.filename_tokens = ExtCatConfig::EXTCAT_FILENAME_TOKENS;
    config.extcat.recursive = ExtCatConfig::EXTCAT_RECURSIVE;
    config.extcat.delimiter = ExtCatConfig::EXTCAT_DELIMITER;
    config.extcat.header_mode = ExtCatConfig::EXTCAT_HEADER_MODE;
    config.extcat.malformed_policy = ExtCatConfig::EXTCAT_MALFORMED_POLICY;
    config.extcat.existing_policy = ExtCatConfig::EXTCAT_EXISTING_POLICY;
    config.extcat.chunk_mib = ExtCatConfig::EXTCAT_CHUNK_MIB;
    config.extcat.total_columns = ExtCatConfig::EXTCAT_TOTAL_COLUMNS;
    config.extcat.use_explicit_columns = ExtCatConfig::EXTCAT_USE_EXPLICIT_COLUMNS;
    config.extcat.input_columns_one_based = ExtCatConfig::EXTCAT_INPUT_COLUMNS_ONE_BASED;
    config.extcat.use_explicit_coordinate_columns = ExtCatConfig::EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS;
    config.extcat.ra_column_one_based = ExtCatConfig::EXTCAT_RA_COLUMN_ONE_BASED;
    config.extcat.dec_column_one_based = ExtCatConfig::EXTCAT_DEC_COLUMN_ONE_BASED;
    config.extcat.mag_g_column_one_based = ExtCatConfig::EXTCAT_MAG_G_COLUMN_ONE_BASED;
    config.extcat.mag_r_column_one_based = ExtCatConfig::EXTCAT_MAG_R_COLUMN_ONE_BASED;
    config.extcat.mag_i_column_one_based = ExtCatConfig::EXTCAT_MAG_I_COLUMN_ONE_BASED;
    config.extcat.mag_z_column_one_based = ExtCatConfig::EXTCAT_MAG_Z_COLUMN_ONE_BASED;
    config.extcat.mag_y_column_one_based = ExtCatConfig::EXTCAT_MAG_Y_COLUMN_ONE_BASED;
    config.extcat.zp_column_one_based = ExtCatConfig::EXTCAT_ZP_COLUMN_ONE_BASED;

    config.init.science_root = InitConfig::SCIENCE_ROOT;
    config.init.dq_root = InitConfig::DQ_ROOT;
    config.init.output_root = InitConfig::OUTPUT_ROOT;
    config.init.datasets = InitConfig::DATASETS;
    config.init.contains = InitConfig::CONTAINS;
    config.init.existing = InitConfig::EXISTING;
    config.init.f77_max_path = InitConfig::F77_MAX_PATH;

    config.lensing.astrometry_trivial = LensingConfig::ASTROMETRY_trivial;
    config.lensing.process_stage = LensingConfig::PROCESS_stage;
    config.lensing.include_flat = LensingConfig::include_FLAT;
    config.lensing.include_mask = LensingConfig::include_Mask;
    config.lensing.astrometry_cat = LensingConfig::ASTROMETRY_CAT;
    config.lensing.flat_path = LensingConfig::FLAT_PATH;
    config.lensing.psf_path = LensingConfig::PSF_PATH;
    config.lensing.ext_cat = LensingConfig::ext_cat;
    config.lensing.ext_psf = LensingConfig::ext_PSF;
    config.lensing.ccd_split = LensingConfig::CCD_split;
    config.lensing.deblending = LensingConfig::deblending;
    config.lensing.psf_type = LensingConfig::PSF_type;
    config.lensing.psf_ms = LensingConfig::PSF_Ms;
    config.lensing.gal_smooth = LensingConfig::gal_smooth;
    config.lensing.star_smooth = LensingConfig::star_smooth;
    config.lensing.pixel_size = LensingConfig::pixel_size;
    config.lensing.nmax_chip = LensingConfig::DEFAULT_CHIP_COUNT;
    config.lensing.chipnx = LensingConfig::chipnx;
    config.lensing.chipny = LensingConfig::chipny;
    return config;
}

bool findRuntimeConfigPath(int argc, char* argv[], std::string& config_path,
                           bool& help_requested, std::string& error) {
    config_path.clear();
    help_requested = false;
    error.clear();
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            help_requested = true;
        } else if (argument == "--config") {
            if (index + 1 >= argc || std::string(argv[index + 1]).rfind("--", 0) == 0) {
                error = "missing value after --config";
                return false;
            }
            config_path = argv[++index];
        } else if (argument.rfind("--config=", 0) == 0) {
            config_path = argument.substr(std::string("--config=").size());
            if (config_path.empty()) {
                error = "--config path must not be empty";
                return false;
            }
        }
    }
    return true;
}

bool applyRuntimeConfigText(const std::string& text,
                            const std::string& source_name,
                            RuntimeConfig& config,
                            std::string& error) {
    RuntimeConfig candidate = config;
    std::istringstream input(text);
    std::string section;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string cleaned = trim(stripInlineComment(line));
        if (cleaned.empty()) {
            continue;
        }
        if (cleaned.front() == '[') {
            if (cleaned.size() < 3 || cleaned.back() != ']') {
                error = formatConfigError(source_name, section, "<section>", cleaned,
                                          "malformed section header at line " + std::to_string(line_number));
                return false;
            }
            section = lowercase(trim(cleaned.substr(1, cleaned.size() - 2)));
            if (section != "process" && section != "extcat"
                && section != "init" && section != "lensing") {
                error = formatConfigError(source_name, section, "<section>", cleaned,
                                          "unknown section at line " + std::to_string(line_number));
                return false;
            }
            continue;
        }
        const std::size_t equals = cleaned.find('=');
        if (equals == std::string::npos) {
            error = formatConfigError(source_name, section, "<none>", cleaned,
                                      "expected key=value at line " + std::to_string(line_number));
            return false;
        }
        if (section.empty()) {
            error = formatConfigError(source_name, section,
                                      trim(cleaned.substr(0, equals)),
                                      trim(cleaned.substr(equals + 1)),
                                      "key appears before any section at line " + std::to_string(line_number));
            return false;
        }
        const std::string key = lowercase(trim(cleaned.substr(0, equals)));
        const std::string raw_value = cleaned.substr(equals + 1);
        std::string reason;
        if (key.empty() || !applyIniValue(section, key, raw_value, candidate, reason)) {
            if (reason.empty()) reason = "empty key";
            error = formatConfigError(source_name, section, key,
                                      trim(raw_value), reason + " at line " + std::to_string(line_number));
            return false;
        }
    }
    config = std::move(candidate);
    error.clear();
    return true;
}

bool parseRuntimeCommandLine(int argc, char* argv[], RuntimeConfig& config,
                             std::string& error) {
    std::string legacy_exposure_list;
    bool legacy_exposure_list_supplied = false;
    CommandLineState state;
    error.clear();
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            config.help_requested = true;
            continue;
        }
        if (argument == "--config") {
            if (index + 1 >= argc) {
                error = "missing value after --config";
                return false;
            }
            ++index;
            continue;
        }
        if (argument.rfind("--config=", 0) == 0) {
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
        if (!applyNamedOption(name, value, config, state, error)) {
            return false;
        }
    }
    if (!config.external_expo_list_supplied && legacy_exposure_list_supplied) {
        config.process.expo_list = legacy_exposure_list;
        config.external_expo_list_supplied = true;
    }
    return true;
}

bool validateRuntimeConfig(const RuntimeConfig& config, std::string& error) {
    error.clear();
    const ProcessRuntimeConfig& process = config.process;
    const ExtCatRuntimeConfig& extcat = config.extcat;
    const InitRuntimeConfig& init = config.init;
    const LensingRuntimeConfig& lensing = config.lensing;

    if ((process.run_process_extcat
                || (process.run_process_main && lensing.ext_cat == 1))
               && extcat.output_directory.empty()) {
        error = "external source-catalog output directory must not be empty";
    } else if (process.run_process_extcat && extcat.input_directory.empty()) {
        error = "external source-catalog input directory must not be empty";
    } else if (process.run_process_fd && lensing.ext_cat != 1) {
        error = "process_fd currently requires ext_cat=1";
    } else if ((process.run_process_init || process.run_process_main
                || process.run_process_rearr || process.run_process_fd)
               && init.datasets.empty()) {
        error = "at least one dataset must be configured";
    } else if (!valueInSet(lensing.astrometry_trivial, {0, 1})) {
        error = "lensing.astrometry_trivial must be 0 or 1";
    } else if (!valueInSet(lensing.include_flat, {0, 1})) {
        error = "lensing.include_flat must be 0 or 1";
    } else if (!valueInSet(lensing.include_mask, {0, 1, 2, 3})) {
        error = "lensing.include_mask must be 0, 1, 2, or 3";
    } else if (!valueInSet(lensing.ext_cat, {0, 1})) {
        error = "lensing.ext_cat must be 0 or 1";
    } else if (!valueInSet(lensing.ext_psf, {0, 1})) {
        error = "lensing.ext_psf must be 0 or 1";
    } else if (!valueInSet(lensing.ccd_split, {1, 2})) {
        error = "lensing.ccd_split must be 1 or 2";
    } else if (!valueInSet(lensing.deblending, {0, 1})) {
        error = "lensing.deblending must be 0 or 1";
    } else if (!valueInSet(lensing.psf_type, {1, 2})) {
        error = "lensing.psf_type must be 1 or 2";
    } else if (!valueInSet(lensing.psf_ms, {0, 1})) {
        error = "lensing.psf_ms must be 0 or 1";
    } else if (lensing.process_stage <= 0) {
        error = "lensing.process_stage must be positive";
    } else if (lensing.process_stage % 23 == 0
               && lensing.process_stage % 19 != 0) {
        error = "Stage 9 requires Stage 8: process_stage has factor 23 without 19";
    } else if (lensing.gal_smooth < 0 || lensing.star_smooth < 0) {
        error = "lensing smoothing parameters must be non-negative";
    } else if (!(lensing.pixel_size > 0.0) || !std::isfinite(lensing.pixel_size)) {
        error = "lensing.pixel_size must be finite and positive";
    } else if (lensing.nmax_chip <= 0 || lensing.chipnx <= 0
               || lensing.chipny <= 0) {
        error = "lensing.nmax_chip, chipnx, and chipny must be positive";
    } else if (extcat.chunk_mib == 0 || extcat.total_columns == 0) {
        error = "extcat.chunk_mib and total_columns must be positive";
    } else if (extcat.ra_column_one_based == 0
               || extcat.dec_column_one_based == 0
               || extcat.zp_column_one_based == 0) {
        error = "extcat RA, Dec, and ZP columns must be positive";
    } else if (extcat.delimiter != "auto" && extcat.delimiter != "whitespace"
               && extcat.delimiter != "comma" && extcat.delimiter != "tab") {
        error = "extcat.delimiter must be auto, whitespace, comma, or tab";
    } else if (extcat.header_mode != "auto" && extcat.header_mode != "present"
               && extcat.header_mode != "absent") {
        error = "extcat.header_mode must be auto, present, or absent";
    } else if (extcat.malformed_policy != "fail"
               && extcat.malformed_policy != "skip") {
        error = "extcat.malformed_policy must be fail or skip";
    } else if (extcat.existing_policy != "fail"
               && extcat.existing_policy != "overwrite") {
        error = "extcat.existing_policy must be fail or overwrite";
    } else if (extcat.use_explicit_columns
               && extcat.input_columns_one_based.empty()) {
        error = "extcat explicit input column list must not be empty";
    } else if (init.existing != "fail" && init.existing != "resume"
               && init.existing != "overwrite") {
        error = "init.existing must be fail, resume, or overwrite";
    } else if (init.f77_max_path < 0) {
        error = "init.f77_max_path must be non-negative";
    }
    if (!error.empty()) {
        return false;
    }

    std::set<std::string> targets;
    for (const InitConfig::DatasetSpec& dataset : init.datasets) {
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
    for (const std::string& token : init.contains) {
        if (token.empty()) {
            error = "init contains tokens must be non-empty";
            return false;
        }
    }
    for (const std::string& token : extcat.filename_tokens) {
        if (token.empty()) {
            error = "extcat filename tokens must be non-empty";
            return false;
        }
    }
    if ((process.run_process_main || process.run_process_rearr)
        && !process.run_process_init && init.datasets.size() > 1
        && !process.expo_list.empty()) {
        error = "one expo_list cannot serve multiple downstream-only datasets";
        return false;
    }
    return true;
}

namespace RuntimeConfigStore {

bool initialize(RuntimeConfig config, std::string& error) {
    if (stored_config) {
        error = "RuntimeConfigStore is already initialized";
        return false;
    }
    stored_config = std::make_unique<const RuntimeConfig>(std::move(config));
    error.clear();
    return true;
}

bool isInitialized() {
    return static_cast<bool>(stored_config);
}

const RuntimeConfig& get() {
    if (!stored_config) {
        throw std::logic_error("RuntimeConfigStore is not initialized");
    }
    return *stored_config;
}

}  // namespace RuntimeConfigStore
