#ifndef CATALOGTILENAMING_HPP
#define CATALOGTILENAMING_HPP

#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

namespace CatalogTileNaming {

// ==========================================
// Function: Format one signed declination tile boundary
// Method: Prefix the absolute integer value with p/m and at least two digits.
// ==========================================
inline std::string formatDeclination(int value) {
    std::ostringstream output;
    output << (value >= 0 ? 'p' : 'm')
           << std::setw(2) << std::setfill('0') << std::abs(value);
    return output.str();
}

// ==========================================
// Function: Build one canonical one-degree catalog tile basename
// Method: Append the fixed RA/Dec boundary grammar to the caller-selected prefix.
// ==========================================
inline std::string tileFilename(std::string_view prefix,
                                int ra_lower,
                                int dec_lower) {
    std::ostringstream output;
    output << prefix << "RA_"
           << std::setw(3) << std::setfill('0') << ra_lower << '_'
           << std::setw(3) << std::setfill('0') << ra_lower + 1
           << "_Dec_" << formatDeclination(dec_lower) << '_'
           << formatDeclination(dec_lower + 1) << ".dat";
    return output.str();
}

// ==========================================
// Function: Recognize one canonical one-degree catalog tile basename
// Method: Match the configured prefix literally, then validate the fixed suffix by position.
// ==========================================
inline bool isTileFilename(std::string_view basename,
                           std::string_view prefix) {
    constexpr std::size_t suffix_size = 26;
    if (basename.size() != prefix.size() + suffix_size
        || basename.substr(0, prefix.size()) != prefix) {
        return false;
    }

    const std::string_view suffix = basename.substr(prefix.size());
    if (suffix.substr(0, 3) != "RA_"
        || suffix[6] != '_'
        || suffix.substr(10, 5) != "_Dec_"
        || suffix[18] != '_'
        || suffix.substr(22, 4) != ".dat"
        || (suffix[15] != 'p' && suffix[15] != 'm')
        || (suffix[19] != 'p' && suffix[19] != 'm')) {
        return false;
    }

    constexpr std::size_t digit_positions[] = {
        3, 4, 5, 7, 8, 9, 16, 17, 20, 21};
    for (const std::size_t position : digit_positions) {
        if (suffix[position] < '0' || suffix[position] > '9') {
            return false;
        }
    }
    return true;
}

}  // namespace CatalogTileNaming

#endif  // CATALOGTILENAMING_HPP
