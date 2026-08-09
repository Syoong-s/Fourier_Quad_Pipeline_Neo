#include "PointSourceStatistics.hpp"
#include "LensingConfig.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace PointSourceStatistics {
namespace {

constexpr int kMinValidPixels = 8;
constexpr double kTaperStartFraction = 0.8;

// ==========================================
// Function: Taper Fourier modes near the reliable-window boundary.
// Method: Use unit weight through 0.8*kmax and a raised-cosine taper to zero.
// ==========================================
double radialWeight(double k, double kmax) {
    if (!(k >= 0.0) || !(kmax > 0.0) || k >= kmax) {
        return 0.0;
    }

    const double ratio = k / kmax;
    if (ratio <= kTaperStartFraction) {
        return 1.0;
    }

    const double phase = LensingConfig::pi
        * (ratio - kTaperStartFraction)
        / (1.0 - kTaperStartFraction);
    return 0.5 * (1.0 + std::cos(phase));
}

// ==========================================
// Function: Find the PSF-supported Fourier radius.
// Method: Locate the nearest mode whose finite PSF power is no greater than
//         1e-4 of the global positive PSF peak, matching Stage 7's cutoff.
// ==========================================
double findReliableKWindow(int ns, const std::vector<float>& psf_power) {
    double peak = 0.0;
    for (const float value : psf_power) {
        if (std::isfinite(value) && value > peak) {
            peak = static_cast<double>(value);
        }
    }
    if (!(peak > 0.0) || !std::isfinite(peak)) {
        return 0.0;
    }

    const double threshold = 1.0e-4 * peak;
    const int center = ns / 2;
    double min_k2 = std::numeric_limits<double>::infinity();

    for (int y = 0; y < ns; ++y) {
        const double ky = static_cast<double>(y - center);
        for (int x = 0; x < ns; ++x) {
            const std::size_t idx = static_cast<std::size_t>(y)
                * static_cast<std::size_t>(ns) + static_cast<std::size_t>(x);
            const double value = static_cast<double>(psf_power[idx]);
            if (std::isfinite(value) && value > threshold) {
                continue;
            }

            const double kx = static_cast<double>(x - center);
            const double k2 = kx * kx + ky * ky;
            min_k2 = std::min(min_k2, k2);
        }
    }

    return std::isfinite(min_k2) && min_k2 > 0.0
        ? std::sqrt(min_k2)
        : 0.0;
}

// ==========================================
// Function: Check whether a statistic can be stored in the float catalog.
// Method: Reject NaN, infinity, and finite doubles outside float range.
// ==========================================
bool isRepresentableFloat(double value) {
    return std::isfinite(value)
        && std::fabs(value) <= static_cast<double>(std::numeric_limits<float>::max());
}

}  // namespace

// ==========================================
// Function: Measure whether a source is PSF-like or spatially extended.
// Method: Accumulate six tapered inner products for the PSF and the single
//         fixed-beta extended template, then evaluate both statistics
//         analytically without allocating template images.
// ==========================================
PointSourceStatisticsResult measure(
    int ns,
    const std::vector<float>& source_power,
    const std::vector<float>& psf_power) {
    PointSourceStatisticsResult result;

    if (ns <= 1) {
        return result;
    }
    const std::size_t side = static_cast<std::size_t>(ns);
    if (side > std::numeric_limits<std::size_t>::max() / side) {
        return result;
    }
    const std::size_t expected_size = side * side;
    if (source_power.size() != expected_size || psf_power.size() != expected_size) {
        return result;
    }

    const double r_win = findReliableKWindow(ns, psf_power);
    const double kmax = LensingConfig::point_stat_k_frac * r_win;
    if (!(kmax > 0.0) || !std::isfinite(kmax)) {
        return result;
    }

    double mm = 0.0;
    double mp = 0.0;
    double pp = 0.0;
    double mt = 0.0;
    double tt = 0.0;
    double pt = 0.0;
    int valid_pixels = 0;

    const int center = ns / 2;
    const double kmax2 = kmax * kmax;
    for (int y = 0; y < ns; ++y) {
        const double ky = static_cast<double>(y - center);
        for (int x = 0; x < ns; ++x) {
            const double kx = static_cast<double>(x - center);
            const double k2 = kx * kx + ky * ky;
            if (k2 <= 0.0 || k2 >= kmax2) {
                continue;
            }

            const std::size_t idx = static_cast<std::size_t>(y) * side
                + static_cast<std::size_t>(x);
            const double source = static_cast<double>(source_power[idx]);
            const double psf = static_cast<double>(psf_power[idx]);
            if (!std::isfinite(source) || !std::isfinite(psf) || psf <= 0.0) {
                continue;
            }

            const double weight = radialWeight(std::sqrt(k2), kmax);
            if (!(weight > 0.0)) {
                continue;
            }

            const double extension_factor = std::exp(
                -LensingConfig::point_stat_beta * k2 / kmax2);
            const double extended = psf * extension_factor;

            mm += weight * source * source;
            mp += weight * source * psf;
            pp += weight * psf * psf;
            mt += weight * source * extended;
            tt += weight * extended * extended;
            pt += weight * psf * extended;
            ++valid_pixels;
        }
    }

    const double eps = LensingConfig::point_stat_eps;
    if (valid_pixels < kMinValidPixels
        || !std::isfinite(mm) || !std::isfinite(mp) || !std::isfinite(pp)
        || !std::isfinite(mt) || !std::isfinite(tt) || !std::isfinite(pt)
        || mm <= eps || pp <= eps || tt <= eps) {
        return result;
    }

    const double projection_norm = std::sqrt(mm) * std::sqrt(pp);
    if (!(projection_norm > eps) || !std::isfinite(projection_norm)) {
        return result;
    }
    const double correlation = std::fabs(mp) / projection_norm;
    if (!std::isfinite(correlation)
        || correlation < LensingConfig::point_stat_min_corr) {
        return result;
    }

    const double psf_residual = std::max(0.0, mm - mp * mp / pp);
    const double extended_residual = std::max(0.0, mm - mt * mt / tt);
    const double delta_chi2 = (psf_residual - extended_residual) / (mm + eps);

    const double extension_projection = mt - pt * mp / pp;
    const double extension_norm2 = std::max(0.0, tt - pt * pt / pp);
    if (!(extension_norm2 > eps) || !std::isfinite(extension_norm2)) {
        return result;
    }
    const double orth_ext = (extension_projection / std::sqrt(extension_norm2))
        * (std::sqrt(pp) / (std::fabs(mp) + eps));

    if (!isRepresentableFloat(delta_chi2) || !isRepresentableFloat(orth_ext)) {
        return result;
    }

    result.delta_chi2 = static_cast<float>(delta_chi2);
    result.orth_ext = static_cast<float>(orth_ext);
    result.valid = true;
    return result;
}

}  // namespace PointSourceStatistics
