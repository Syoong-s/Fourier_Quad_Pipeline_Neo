#ifndef PSF_CANDIDATE_QUALITY_HPP
#define PSF_CANDIDATE_QUALITY_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace PSFModel {
namespace Internal {

enum class CandidatePowerStatus {
    Accepted,
    InvalidShape,
    NonFinitePower,
    NegativeCoreMedian,
    NonPositiveSum
};

// ==========================================
// Function: Assess one corrected star-candidate power spectrum
// Method: Require a complete finite image, a non-negative median around DC,
//         and a finite positive signed sum for pairwise normalization.
// ==========================================
inline CandidatePowerStatus assessCandidatePower(
    int nx, int ny, const std::vector<float>& power, double& sum_power) {
    sum_power = 0.0;
    if (nx < 3 || ny < 3
        || power.size() != static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny)) {
        return CandidatePowerStatus::InvalidShape;
    }

    for (float value : power) {
        if (!std::isfinite(value)) {
            return CandidatePowerStatus::NonFinitePower;
        }
    }

    const int cx = nx / 2;
    const int cy = ny / 2;
    std::array<float, 8> center_neighbors = {
        power[(cy - 1) * nx + (cx - 1)],
        power[(cy - 1) * nx + cx],
        power[(cy - 1) * nx + (cx + 1)],
        power[cy * nx + (cx - 1)],
        power[cy * nx + (cx + 1)],
        power[(cy + 1) * nx + (cx - 1)],
        power[(cy + 1) * nx + cx],
        power[(cy + 1) * nx + (cx + 1)]
    };
    std::sort(center_neighbors.begin(), center_neighbors.end());
    const double median8 = 0.5
        * (static_cast<double>(center_neighbors[3])
           + static_cast<double>(center_neighbors[4]));
    if (median8 < 0.0) {
        return CandidatePowerStatus::NegativeCoreMedian;
    }

    for (float value : power) {
        sum_power += static_cast<double>(value);
    }
    if (!std::isfinite(sum_power) || sum_power <= 0.0) {
        return CandidatePowerStatus::NonPositiveSum;
    }
    return CandidatePowerStatus::Accepted;
}

// ==========================================
// Function: Validate diagnostics derived from an accepted candidate spectrum
// Method: Reject only non-finite numerical results so NaN cannot reach star
//         selection or Stage-8 aggregation after recovery is removed.
// ==========================================
inline bool candidateDiagnosticsAreFinite(
    double size, double e1, double e2, double fwhm) {
    return std::isfinite(size) && std::isfinite(e1)
        && std::isfinite(e2) && std::isfinite(fwhm);
}

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_CANDIDATE_QUALITY_HPP
