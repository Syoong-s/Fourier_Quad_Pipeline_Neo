#ifndef F77_NOISE_SELECTION_HPP
#define F77_NOISE_SELECTION_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace SourceExtractor {
namespace Internal {

// ==========================================
// Structure: Describe the deterministic F77 blank-noise candidate
// Method: Preserve the zero-based origin and raw maximum used by the caller.
// ==========================================
struct F77NoiseCandidate {
    bool found = false;
    int x0 = 0;
    int y0 = 0;
    float peak = 100000.0f;
};

// ==========================================
// Function: Select the historical F77 blank-noise region
// Method: Scan the 5x5 outer ring in Fortran order and retain the first
//         candidate with nd<=nl and a strictly smaller raw maximum.
// ==========================================
inline F77NoiseCandidate selectF77NoiseCandidate(
    int nx,
    int ny,
    const std::vector<float>& image,
    const std::vector<int>& weight,
    double xp,
    double yp,
    int nl,
    int nl_half,
    int chip_edge_margin) {
    F77NoiseCandidate selected;
    const std::size_t expected = static_cast<std::size_t>(nx)
                               * static_cast<std::size_t>(ny);
    if (nx <= 0 || ny <= 0 || nl <= 0 || nl_half <= 0
        || chip_edge_margin <= 0 || image.size() != expected
        || weight.size() != expected || !std::isfinite(xp)
        || !std::isfinite(yp)) {
        return selected;
    }

    const int center_x0 = static_cast<int>(xp + 0.5) - 1;
    const int center_y0 = static_cast<int>(yp + 0.5) - 1;
    const int base_x0 = center_x0 - nl_half - 1;
    const int base_y0 = center_y0 - nl_half - 1;
    constexpr int step_count = 5;
    constexpr float sentinel = 100000.0f;

    for (int i = 0; i < step_count; ++i) {
        for (int j = 0; j < step_count; ++j) {
            if (i != 0 && i != step_count - 1
                && j != 0 && j != step_count - 1) {
                continue;
            }
            const int x0 = base_x0 + nl * (i - 2);
            const int y0 = base_y0 + nl * (j - 2);
            const int x1 = x0 + nl - 1;
            const int y1 = y0 + nl - 1;
            if (x0 < chip_edge_margin - 1
                || x1 > nx - chip_edge_margin - 1
                || y0 < chip_edge_margin - 1
                || y1 > ny - chip_edge_margin - 1) {
                continue;
            }

            float peak = -sentinel;
            int defect_count = 0;
            for (int x = x0; x <= x1; ++x) {
                for (int y = y0; y <= y1; ++y) {
                    const std::size_t index = static_cast<std::size_t>(y) * nx + x;
                    peak = std::max(peak, image[index]);
                    if (weight[index] == 0) ++defect_count;
                }
            }
            if (defect_count <= nl && peak < selected.peak) {
                selected.found = true;
                selected.x0 = x0;
                selected.y0 = y0;
                selected.peak = peak;
            }
        }
    }
    return selected;
}

}  // namespace Internal
}  // namespace SourceExtractor

#endif  // F77_NOISE_SELECTION_HPP
