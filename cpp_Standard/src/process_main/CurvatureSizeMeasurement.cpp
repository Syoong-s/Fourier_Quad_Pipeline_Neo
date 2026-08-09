#include "ShearMeasurement.hpp"
#include "LensingConfig.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace ShearMeasurement {
namespace {

constexpr float invalid_size = -999.0f;

}  // namespace

// ==========================================
// Function: Measure the low-frequency curvature size of a Fourier power spectrum.
// Method: Fit B0 + B1*x + B2*x^2 with equal pixel weights inside a centered
//         Fourier disk, then evaluate T = -2*B1/(B0*qmax^2).
// ==========================================
float measurePowerCurvatureSize(int ns, const std::vector<float>& power, int rmax) {
    if (ns <= 0 || rmax <= 0 ||
        power.size() != static_cast<std::size_t>(ns) * static_cast<std::size_t>(ns)) {
        return invalid_size;
    }

    const int center = ns / 2;
    const int complete_disk_radius = std::min(center, ns - 1 - center);
    if (rmax > complete_disk_radius) {
        return invalid_size;
    }

    const int rmax2 = rmax * rmax;
    const double dq = 2.0 * LensingConfig::pi / static_cast<double>(ns);
    const double qmax2 = dq * dq * static_cast<double>(rmax2);

    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    int pixel_count = 0;

    for (int y = 0; y < ns; ++y) {
        const int v = y - center;
        for (int xpix = 0; xpix < ns; ++xpix) {
            const int u = xpix - center;
            const int r2 = u * u + v * v;
            if (r2 > rmax2) {
                continue;
            }

            const double observed = static_cast<double>(power[y * ns + xpix]);
            if (!std::isfinite(observed)) {
                continue;
            }

            const double x = static_cast<double>(r2) / static_cast<double>(rmax2);
            const Eigen::Vector3d basis(1.0, x, x * x);
            normal.noalias() += basis * basis.transpose();
            rhs.noalias() += basis * observed;
            ++pixel_count;
        }
    }

    if (pixel_count < 6) {
        return invalid_size;
    }

    const Eigen::FullPivLU<Eigen::Matrix3d> solver(normal);
    if (!solver.isInvertible()) {
        return invalid_size;
    }

    const Eigen::Vector3d coefficients = solver.solve(rhs);
    if (!coefficients.allFinite() || coefficients[0] <= 0.0) {
        return invalid_size;
    }

    const double size = -2.0 * coefficients[1] / (coefficients[0] * qmax2);
    if (!std::isfinite(size) || size <= 0.0 ||
        size > static_cast<double>(std::numeric_limits<float>::max())) {
        return invalid_size;
    }

    return static_cast<float>(size);
}

}  // namespace ShearMeasurement
