#ifndef GENERAL_OUTPUT_LAYOUT_HPP
#define GENERAL_OUTPUT_LAYOUT_HPP

#include "pathconfig.hpp"

#include <stdexcept>
#include <string>

namespace OutputLayout {

// Extract the exposure component from the established <exposure>_<ccdnum>
// prefix contract used by process_init and UniversalUtils::getPrefix.
inline std::string exposureFromChipPrefix(const std::string& chip_prefix) {
    const std::size_t delimiter = chip_prefix.find_last_of('_');
    if (delimiter == std::string::npos || delimiter == 0
        || delimiter + 1 >= chip_prefix.size()) {
        throw std::invalid_argument("invalid chip prefix: " + chip_prefix);
    }
    return chip_prefix.substr(0, delimiter);
}

// Build <output-root>/<product-directory>/<exposure>/<ccdnum-prefix><suffix>.
inline std::string chipPath(const std::string& output_root,
                            const std::string& product_directory,
                            const std::string& chip_prefix,
                            const std::string& suffix) {
    return output_root + "/" + product_directory + "/"
           + exposureFromChipPrefix(chip_prefix) + "/"
           + chip_prefix + suffix;
}

}  // namespace OutputLayout

#endif  // GENERAL_OUTPUT_LAYOUT_HPP
