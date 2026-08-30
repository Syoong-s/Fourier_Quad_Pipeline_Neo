#include "general/ExposureList.hpp"

#include <fstream>
#include <utility>

namespace {

// ==========================================
// Function: Remove the pipeline's matching double quotes from one path token
// Method: Preserve the historical token parser by stripping only a complete
//         pair of double quotes after formatted extraction.
// ==========================================
std::string stripDoubleQuotes(std::string value) {
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

// ==========================================
// Function: Validate the shared exposure-list cardinality contract
// Method: Reject empty input and enforce a nonzero caller limit while leaving
//         protocol-specific error prefixes to the caller.
// ==========================================
bool validateCount(std::size_t count,
                   std::size_t max_entries,
                   const std::string& filename,
                   std::string& error) {
    if (count == 0) {
        error = "exposure list contains no entries: " + filename;
        return false;
    }
    if (max_entries > 0 && count > max_entries) {
        error = "exposure list exceeds the configured entry limit: " + filename;
        return false;
    }
    error.clear();
    return true;
}

}  // namespace

namespace ExposureList {

bool loadPipelineList(const std::string& filename,
                      std::vector<Entry>& entries,
                      std::size_t max_entries,
                      std::string& error) {
    entries.clear();
    std::ifstream input(filename);
    if (!input.is_open()) {
        error = "cannot open exposure list: " + filename;
        return false;
    }

    std::string path;
    int chip_count = 0;
    while (input >> path >> chip_count) {
        entries.push_back({stripDoubleQuotes(std::move(path)), chip_count});
    }
    if (!input.eof()) {
        entries.clear();
        error = "exposure list contains an invalid record: " + filename;
        return false;
    }
    if (!validateCount(entries.size(), max_entries, filename, error)) {
        entries.clear();
        return false;
    }
    return true;
}

bool loadPathList(const std::string& filename,
                  std::vector<std::string>& paths,
                  std::size_t max_entries,
                  std::string& error) {
    paths.clear();
    std::ifstream input(filename);
    if (!input.is_open()) {
        error = "cannot open exposure list: " + filename;
        return false;
    }

    std::string path;
    while (input >> path) {
        paths.push_back(stripDoubleQuotes(std::move(path)));
    }
    if (!input.eof()) {
        paths.clear();
        error = "exposure list contains an invalid record: " + filename;
        return false;
    }
    if (!validateCount(paths.size(), max_entries, filename, error)) {
        paths.clear();
        return false;
    }
    return true;
}

}  // namespace ExposureList
