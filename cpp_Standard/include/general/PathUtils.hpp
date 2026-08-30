#ifndef GENERAL_PATH_UTILS_HPP
#define GENERAL_PATH_UTILS_HPP

#include <filesystem>
#include <string>

namespace PathUtils {

// ==========================================
// Function: Normalize one filesystem path
// Method: Resolve existing components while allowing a not-yet-created tail.
// ==========================================
std::filesystem::path normalizedAbsolute(const std::filesystem::path& path);

// ==========================================
// Function: Test whether one path is equal to or below another
// Method: Compare complete normalized path components instead of string
//         prefixes.
// ==========================================
bool isPathWithin(const std::filesystem::path& candidate,
                  const std::filesystem::path& parent);

// ==========================================
// Function: Resolve an ancestor at a fixed level
// Method: Walk parent_path exactly levels times and report insufficient depth.
// ==========================================
bool parentAtLevel(const std::filesystem::path& path,
                   int levels,
                   std::filesystem::path& parent,
                   std::string& error);

}  // namespace PathUtils

#endif  // GENERAL_PATH_UTILS_HPP
