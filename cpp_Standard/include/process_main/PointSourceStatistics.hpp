#ifndef POINT_SOURCE_STATISTICS_HPP
#define POINT_SOURCE_STATISTICS_HPP

#include <vector>

// ==========================================
// PointSourceStatisticsResult - Fourier-power morphology measurements.
// Method: Publish the two catalog statistics together with an in-memory
//         validity flag; invalid measurements retain the -999 sentinel.
// ==========================================
struct PointSourceStatisticsResult {
    float delta_chi2 = -999.0f;
    float orth_ext = -999.0f;
    bool valid = false;
};

namespace PointSourceStatistics {

// ==========================================
// Function: Measure whether a source is PSF-like or spatially extended.
// Method: Compare the source power with the local PSF and one fixed-beta
//         extended template using one joint inner-product traversal.
// ==========================================
[[nodiscard]] PointSourceStatisticsResult measure(
    int ns,
    const std::vector<float>& source_power,
    const std::vector<float>& psf_power);

}  // namespace PointSourceStatistics

#endif  // POINT_SOURCE_STATISTICS_HPP
