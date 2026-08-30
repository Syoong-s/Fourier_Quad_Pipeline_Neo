#include "general/PathUtils.hpp"

namespace PathUtils {

std::filesystem::path normalizedAbsolute(const std::filesystem::path& path) {
    return std::filesystem::weakly_canonical(std::filesystem::absolute(path));
}

bool isPathWithin(const std::filesystem::path& candidate,
                  const std::filesystem::path& parent) {
    auto candidate_iterator = candidate.begin();
    auto parent_iterator = parent.begin();
    for (; parent_iterator != parent.end();
         ++parent_iterator, ++candidate_iterator) {
        if (candidate_iterator == candidate.end()
            || *candidate_iterator != *parent_iterator) {
            return false;
        }
    }
    return true;
}

bool parentAtLevel(const std::filesystem::path& path,
                   int levels,
                   std::filesystem::path& parent,
                   std::string& error) {
    if (levels < 0) {
        error = "parent level must be non-negative";
        parent.clear();
        return false;
    }
    parent = path;
    for (int level = 0; level < levels; ++level) {
        const std::filesystem::path next = parent.parent_path();
        if (next.empty() || next == parent) {
            error = "path has fewer parent levels than requested: " + path.string();
            parent.clear();
            return false;
        }
        parent = next;
    }
    error.clear();
    return true;
}

}  // namespace PathUtils
