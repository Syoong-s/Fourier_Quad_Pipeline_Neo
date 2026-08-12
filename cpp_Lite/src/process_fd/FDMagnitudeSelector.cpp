#include "process_fd/FDMagnitudeSelector.hpp"

#include <array>
#include <optional>
#include <utility>

namespace ProcessFD {

// ==========================================
// Function: Select the magnitude band used by FD
// Method: Choose the first available external magnitude in i-z-r-g-y priority
//         order and fail explicitly when the layout contains no magnitude.
// ==========================================
bool selectMagnitude(const PipelineCatalog::CatalogLayout& layout,
                     MagnitudeSelection& selection,
                     std::string& error) {
    selection = MagnitudeSelection{};
    const std::array<
        std::pair<const char*, const std::optional<std::size_t>*>, 5>
        candidates = {{
            {"i", &layout.external.mag_i},
            {"z", &layout.external.mag_z},
            {"r", &layout.external.mag_r},
            {"g", &layout.external.mag_g},
            {"y", &layout.external.mag_y},
        }};

    for (const auto& candidate : candidates) {
        if (candidate.second->has_value()) {
            selection.band = candidate.first;
            selection.column = candidate.second->value();
            error.clear();
            return true;
        }
    }

    error = "FD requires at least one external magnitude column "
            "(selection priority: i -> z -> r -> g -> y)";
    return false;
}

}  // namespace ProcessFD
