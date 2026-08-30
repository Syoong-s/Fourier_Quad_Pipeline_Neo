#ifndef GENERAL_EXPOSURE_LIST_HPP
#define GENERAL_EXPOSURE_LIST_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace ExposureList {

struct Entry {
    std::string path;
    int chip_count = 0;
};

// ==========================================
// Function: Load a top-level pipeline exposure list
// Method: Parse path-and-chip-count records in file order and enforce a
//         caller-supplied entry limit without changing chip-count policy.
// ==========================================
bool loadPipelineList(const std::string& filename,
                      std::vector<Entry>& entries,
                      std::size_t max_entries,
                      std::string& error);

// ==========================================
// Function: Load a path-only exposure list
// Method: Parse whitespace-delimited paths in file order and enforce a
//         caller-supplied entry limit.
// ==========================================
bool loadPathList(const std::string& filename,
                  std::vector<std::string>& paths,
                  std::size_t max_entries,
                  std::string& error);

}  // namespace ExposureList

#endif  // GENERAL_EXPOSURE_LIST_HPP
