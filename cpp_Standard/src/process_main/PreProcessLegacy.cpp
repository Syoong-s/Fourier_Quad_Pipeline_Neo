#include "process_main/PreProcessLegacy.hpp"

#include "general/NumericalRecipes.hpp"
#include "process_main/UniversalUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace PreProcessLegacy {
namespace {

struct LegacySample {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// ==========================================
// Function: Validate one polynomial coefficient vector
// Method: Require the requested shape and finite values while allowing a
//         scientifically valid all-zero background.
// ==========================================
bool coefficientsAreFinite(
    const std::vector<double>& coefficients,
    int expected_size) {
    if (expected_size <= 0
        || coefficients.size() != static_cast<std::size_t>(expected_size)) {
        return false;
    }
    return std::all_of(
        coefficients.begin(), coefficients.end(),
        [](double value) { return std::isfinite(value); });
}

// ==========================================
// Function: Evaluate one tensor-product polynomial on a legacy coordinate
// Method: Reuse the current checked basis definition and coefficient ordering.
// ==========================================
double evaluateLegacy(
    double x,
    double y,
    int nct,
    int ncx,
    const std::vector<double>& coefficients) {
    return UniversalUtils::funcVal(x, y, nct, ncx, coefficients);
}

// ==========================================
// Function: Compute a small exact binomial coefficient
// Method: Multiply the shorter factorial ratio in double precision.
// ==========================================
double choose(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    k = std::min(k, n - k);
    double value = 1.0;
    for (int i = 1; i <= k; ++i) {
        value *= static_cast<double>(n - k + i)
               / static_cast<double>(i);
    }
    return value;
}

// ==========================================
// Function: Convert a legacy absolute-coordinate background to Stage-3 basis
// Method: Expand u=Ax+Bx*X and v=Ay+By*Y binomially in the current
//         amplifier-local normalized coordinates.
// ==========================================
bool convertLegacyBackground(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    double ratio,
    int nct,
    int ncx,
    const std::vector<double>& legacy_coefficients,
    std::vector<double>& normalized_coefficients) {
    if (!coefficientsAreFinite(legacy_coefficients, nct)
        || !std::isfinite(ratio) || ratio <= 0.0 || ncx <= 0) {
        return false;
    }
    const double x_mid = 0.5 * (
        static_cast<double>(x_start + 1) + static_cast<double>(x_end));
    const double y_mid = 0.5 * (
        static_cast<double>(y_start + 1) + static_cast<double>(y_end));
    const double x_half_inv = 2.0
        / static_cast<double>(std::max(x_end - x_start - 1, 1));
    const double y_half_inv = 2.0
        / static_cast<double>(std::max(y_end - y_start - 1, 1));
    const double ax = ratio * x_mid;
    const double bx = ratio / x_half_inv;
    const double ay = ratio * y_mid;
    const double by = ratio / y_half_inv;

    normalized_coefficients.assign(static_cast<std::size_t>(nct), 0.0);
    for (int term = 0; term < nct; ++term) {
        const int px = term % ncx;
        const int py = term / ncx;
        for (int rx = 0; rx <= px; ++rx) {
            for (int ry = 0; ry <= py; ++ry) {
                const int target = ry * ncx + rx;
                if (target >= nct) return false;
                normalized_coefficients[target] += legacy_coefficients[term]
                    * choose(px, rx) * std::pow(ax, px - rx) * std::pow(bx, rx)
                    * choose(py, ry) * std::pow(ay, py - ry) * std::pow(by, ry);
            }
        }
    }
    return coefficientsAreFinite(normalized_coefficients, nct);
}

// ==========================================
// Function: Validate the legacy sigma plane without modern quality cuts
// Method: Require only finite positive values at the four rectangle corners.
// ==========================================
bool legacySigmaPlaneIsPositive(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    double aa,
    double bb,
    double cc) {
    const double x1 = static_cast<double>(x_start + 1);
    const double x2 = static_cast<double>(x_end);
    const double y1 = static_cast<double>(y_start + 1);
    const double y2 = static_cast<double>(y_end);
    const double corners[4] = {
        aa + bb * x1 + cc * y1,
        aa + bb * x2 + cc * y1,
        aa + bb * x1 + cc * y2,
        aa + bb * x2 + cc * y2};
    for (double value : corners) {
        if (!std::isfinite(value) || value <= 0.0) return false;
    }
    return true;
}

}  // namespace

bool setBackground(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    int nx,
    int ny,
    std::vector<float>& image,
    int blocksize,
    int nct,
    int ncx,
    std::vector<double>& normalized_coefficients) {
    normalized_coefficients.clear();
    const std::size_t image_size = static_cast<std::size_t>(nx)
                                 * static_cast<std::size_t>(ny);
    if (nx <= 0 || ny <= 0 || x_start < 0 || x_end > nx
        || y_start < 0 || y_end > ny || x_end <= x_start + 1
        || y_end <= y_start + 1 || image.size() < image_size
        || blocksize <= 0 || nct <= 0 || ncx <= 0) {
        std::cerr << "Error / F77 set_background invalid geometry or configuration\n";
        return false;
    }

    constexpr int sample_count = 1000;
    const int nx1 = x_start + 1;
    const int nx2 = x_end;
    const int ny1 = y_start + 1;
    const int ny2 = y_end;
    const int span = std::max(nx2 - nx1, ny2 - ny1);
    if (span <= 0) return false;
    const double ratio = 1.0 / static_cast<double>(span);

    std::vector<LegacySample> rough_samples;
    rough_samples.reserve(sample_count);
    std::vector<double> rough_flux;
    rough_flux.reserve(sample_count);
    for (int sample = 0; sample < sample_count; ++sample) {
        const int ix = static_cast<int>(
            NumericalRecipes::ran1() * (nx2 - nx1) + nx1);
        const int iy = static_cast<int>(
            NumericalRecipes::ran1() * (ny2 - ny1) + ny1);
        const double flux = image[
            static_cast<std::size_t>(iy - 1) * nx + ix - 1];
        if (!std::isfinite(flux)) return false;
        rough_samples.push_back({ix * ratio, iy * ratio, flux});
        rough_flux.push_back(flux);
    }
    std::sort(rough_flux.begin(), rough_flux.end());
    if (rough_flux.front() == rough_flux.back()) return false;
    const double lower = rough_flux[sample_count / 3 - 1];
    const double upper = rough_flux[2 * sample_count / 3 - 1];
    std::vector<Point3D> rough_points;
    rough_points.reserve(sample_count);
    for (const LegacySample& sample : rough_samples) {
        if (sample.z >= lower && sample.z <= upper) {
            rough_points.push_back({sample.x, sample.y, sample.z});
        }
    }
    if (rough_points.size() < static_cast<std::size_t>(sample_count / 10)) {
        return false;
    }
    std::vector<double> rough_coefficients;
    UniversalUtils::fit2D(rough_points, 4, 2, rough_coefficients);
    if (!coefficientsAreFinite(rough_coefficients, 4)) return false;

    for (int x = x_start; x < x_end; ++x) {
        for (int y = y_start; y < y_end; ++y) {
            const double background = evaluateLegacy(
                (x + 1) * ratio, (y + 1) * ratio,
                4, 2, rough_coefficients);
            image[static_cast<std::size_t>(y) * nx + x] -=
                static_cast<float>(background);
        }
    }

    const int nbx = std::max((nx2 - nx1) / blocksize, 1);
    const int nby = std::max((ny2 - ny1) / blocksize, 1);
    std::vector<Point3D> block_points;
    block_points.reserve(static_cast<std::size_t>(nbx) * nby);
    for (int ib = 0; ib < nbx; ++ib) {
        const int xmin = ib * blocksize + nx1;
        const int xmax = std::min(xmin + blocksize, nx2);
        for (int jb = 0; jb < nby; ++jb) {
            const int ymin = jb * blocksize + ny1;
            const int ymax = std::min(ymin + blocksize, ny2);
            if (xmax <= xmin || ymax <= ymin) return false;
            std::vector<LegacySample> samples;
            samples.reserve(sample_count);
            for (int sample = 0; sample < sample_count; ++sample) {
                const int ix = static_cast<int>(
                    NumericalRecipes::ran1() * (xmax - xmin) + xmin);
                const int iy = static_cast<int>(
                    NumericalRecipes::ran1() * (ymax - ymin) + ymin);
                const double flux = image[
                    static_cast<std::size_t>(iy - 1) * nx + ix - 1];
                if (!std::isfinite(flux)) return false;
                samples.push_back({ix * ratio, iy * ratio, flux});
            }
            std::stable_sort(
                samples.begin(), samples.end(),
                [](const LegacySample& first, const LegacySample& second) {
                    return first.z < second.z;
                });
            const LegacySample& median = samples[sample_count / 2 - 1];
            block_points.push_back({median.x, median.y, median.z});
        }
    }

    const int minimum_points = nct * 3 / 2;
    bool changed = true;
    while (changed
           && block_points.size() >= static_cast<std::size_t>(minimum_points)) {
        double aa = 0.0;
        double bb = 0.0;
        double cc = 0.0;
        UniversalUtils::findSlope2D(block_points, aa, bb, cc);
        if (!std::isfinite(aa) || !std::isfinite(bb) || !std::isfinite(cc)) {
            return false;
        }
        std::vector<double> residuals;
        residuals.reserve(block_points.size());
        for (const Point3D& point : block_points) {
            residuals.push_back(point.z - aa - bb * point.x - cc * point.y);
        }
        std::vector<double> ordered = residuals;
        std::sort(ordered.begin(), ordered.end());
        const std::size_t count = ordered.size();
        const double center = ordered[count / 2 - 1];
        const double width = 0.5 * (
            ordered[count * 5 / 6 - 1] - ordered[count / 6 - 1]);
        std::vector<Point3D> retained;
        retained.reserve(block_points.size());
        for (std::size_t index = 0; index < block_points.size(); ++index) {
            if (std::abs(residuals[index] - center) < 3.0 * width) {
                retained.push_back(block_points[index]);
            }
        }
        changed = retained.size() != block_points.size();
        block_points.swap(retained);
    }
    if (block_points.size() < static_cast<std::size_t>(minimum_points)) {
        std::cerr << "Error / F77 set_background unstable block population "
                  << block_points.size() << '\n';
        return false;
    }

    std::vector<double> fine_coefficients;
    UniversalUtils::fit2D(block_points, nct, ncx, fine_coefficients);
    if (!coefficientsAreFinite(fine_coefficients, nct)) return false;
    for (int x = x_start; x < x_end; ++x) {
        for (int y = y_start; y < y_end; ++y) {
            const double background = evaluateLegacy(
                (x + 1) * ratio, (y + 1) * ratio,
                nct, ncx, fine_coefficients);
            image[static_cast<std::size_t>(y) * nx + x] -=
                static_cast<float>(background);
        }
    }

    std::vector<double> total_legacy = fine_coefficients;
    for (int term = 0; term < 4; ++term) {
        const int px = term % 2;
        const int py = term / 2;
        const int target = py * ncx + px;
        if (target < 0 || target >= nct) return false;
        total_legacy[target] += rough_coefficients[term];
    }
    return convertLegacyBackground(
        x_start, x_end, y_start, y_end, ratio,
        nct, ncx, total_legacy, normalized_coefficients);
}

bool setSig(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    int nx,
    int ny,
    std::vector<float>& image,
    double& aa,
    double& bb,
    double& cc) {
    aa = 0.0;
    bb = 0.0;
    cc = 0.0;
    const std::size_t image_size = static_cast<std::size_t>(nx)
                                 * static_cast<std::size_t>(ny);
    if (nx <= 0 || ny <= 0 || x_start < 0 || x_end > nx
        || y_start < 0 || y_end > ny || x_end <= x_start + 2
        || y_end <= y_start + 2 || image.size() < image_size) {
        std::cerr << "Error / F77 set_sig invalid amplifier geometry\n";
        return false;
    }

    constexpr int sample_count = 2000;
    const int nx1 = x_start + 1;
    const int nx2 = x_end;
    const int ny1 = y_start + 1;
    const int ny2 = y_end;
    const int x_range = nx2 - nx1 - 1;
    const int y_range = ny2 - ny1 - 1;
    if (x_range <= 0 || y_range <= 0) return false;

    std::vector<Point3D> samples;
    samples.reserve(sample_count);
    std::vector<double> values;
    values.reserve(sample_count);
    for (int sample = 0; sample < sample_count; ++sample) {
        const int ix = static_cast<int>(
            NumericalRecipes::ran1() * x_range + nx1);
        const int iy = static_cast<int>(
            NumericalRecipes::ran1() * y_range + ny1);
        const std::size_t index =
            static_cast<std::size_t>(iy - 1) * nx + ix - 1;
        const double center = image[index];
        const double dx = center - image[index + 1];
        const double dy = center - image[index + nx];
        const double value = 0.5 * (dx * dx + dy * dy);
        if (!std::isfinite(value)) return false;
        samples.push_back({static_cast<double>(ix), static_cast<double>(iy), value});
        values.push_back(value);
    }
    std::sort(values.begin(), values.end());
    const double lower = values[sample_count / 3 - 1];
    const double upper = values[2 * sample_count / 3 - 1];
    std::vector<Point3D> retained;
    retained.reserve(sample_count);
    for (const Point3D& sample : samples) {
        if (sample.z >= lower && sample.z <= upper) retained.push_back(sample);
    }
    if (retained.size() < static_cast<std::size_t>(sample_count / 10)) {
        return false;
    }
    UniversalUtils::findSlope2D(retained, aa, bb, cc);
    if (!std::isfinite(aa) || !std::isfinite(bb) || !std::isfinite(cc)) {
        return false;
    }

    aa = static_cast<double>(static_cast<float>(aa));
    bb = static_cast<double>(static_cast<float>(bb));
    cc = static_cast<double>(static_cast<float>(cc));
    if (!legacySigmaPlaneIsPositive(
            x_start, x_end, y_start, y_end, aa, bb, cc)) {
        std::cerr << "Error / F77 set_sig nonpositive noise plane\n";
        return false;
    }
    for (int x = x_start; x < x_end; ++x) {
        for (int y = y_start; y < y_end; ++y) {
            const std::size_t index = static_cast<std::size_t>(y) * nx + x;
            const double plane = aa + bb * (x + 1) + cc * (y + 1);
            image[index] = static_cast<float>(
                static_cast<double>(image[index]) / std::sqrt(0.5 * plane));
        }
    }
    return true;
}

}  // namespace PreProcessLegacy
