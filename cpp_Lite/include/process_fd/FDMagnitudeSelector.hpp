#ifndef FD_MAGNITUDE_SELECTOR_HPP
#define FD_MAGNITUDE_SELECTOR_HPP

#include "CatalogLayout.hpp"

#include <cstddef>
#include <string>

namespace ProcessFD {

// ==========================================
// Structure: Identify the external magnitude used by FD star selection
// Method: Carry one band label and its effective zero-based catalog position.
// ==========================================
struct MagnitudeSelection {
    std::string band;
    std::size_t column = 0;
};

// ==========================================
// Function: Select the magnitude band used by FD
// Method: Choose the first available external magnitude in i-z-r-g-y priority
//         order and fail explicitly when the layout contains no magnitude.
// ==========================================
bool selectMagnitude(const PipelineCatalog::CatalogLayout& layout,
                     MagnitudeSelection& selection,
                     std::string& error);

}  // namespace ProcessFD

#endif  // FD_MAGNITUDE_SELECTOR_HPP
