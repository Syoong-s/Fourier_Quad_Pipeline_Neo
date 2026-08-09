#include "PointSourceStatistics.hpp"
#include "LensingConfig.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double kSyntheticReliableRadius = 10.0;
constexpr double kSyntheticKmax =
    LensingConfig::point_stat_k_frac * kSyntheticReliableRadius;

// ==========================================
// Function: Stop the test program when a requirement is not met.
// Method: Print a focused failure message and return a non-zero process code.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PointSourceStatistics test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Compare two morphology statistics with absolute tolerance.
// Method: Require both operands to be finite before applying the tolerance.
// ==========================================
void requireNear(float actual, float expected, float tolerance,
                 const std::string& message) {
    require(std::isfinite(actual) && std::isfinite(expected)
                && std::fabs(actual - expected) <= tolerance,
            message);
}

// ==========================================
// Function: Build a PSF with an exactly known reliable Fourier radius.
// Method: Set unit power for k<10 and zero power outside that disk so the
//         Stage 7 threshold finder returns r_win=10 Fourier pixels.
// ==========================================
std::vector<float> makeSyntheticPsf(int ns) {
    std::vector<float> psf(static_cast<std::size_t>(ns * ns), 0.0f);
    const int center = ns / 2;
    for (int y = 0; y < ns; ++y) {
        const double ky = static_cast<double>(y - center);
        for (int x = 0; x < ns; ++x) {
            const double kx = static_cast<double>(x - center);
            if (std::sqrt(kx * kx + ky * ky) < kSyntheticReliableRadius) {
                psf[static_cast<std::size_t>(y * ns + x)] = 1.0f;
            }
        }
    }
    return psf;
}

// ==========================================
// Function: Build the fixed-beta extended source used by production code.
// Method: Apply the configured Gaussian extension factor at the known test
//         kmax and an arbitrary brightness scale.
// ==========================================
std::vector<float> makeExtendedSource(
    int ns, const std::vector<float>& psf, double amplitude) {
    std::vector<float> source(psf.size(), 0.0f);
    const int center = ns / 2;
    const double kmax2 = kSyntheticKmax * kSyntheticKmax;
    for (int y = 0; y < ns; ++y) {
        const double ky = static_cast<double>(y - center);
        for (int x = 0; x < ns; ++x) {
            const double kx = static_cast<double>(x - center);
            const double k2 = kx * kx + ky * ky;
            const std::size_t idx = static_cast<std::size_t>(y * ns + x);
            source[idx] = static_cast<float>(
                amplitude * static_cast<double>(psf[idx])
                * std::exp(-LensingConfig::point_stat_beta * k2 / kmax2));
        }
    }
    return source;
}

// ==========================================
// Function: Exercise PSF-like, extended, scaled, negative, and invalid inputs.
// Method: Use synthetic powers with a known Fourier window and compare the
//         direction and scale invariance of both catalog statistics.
// ==========================================
void runSyntheticCases() {
    const int ns = LensingConfig::ns;
    const std::vector<float> psf = makeSyntheticPsf(ns);

    std::vector<float> point_source = psf;
    const PointSourceStatisticsResult point =
        PointSourceStatistics::measure(ns, point_source, psf);
    require(point.valid, "ideal PSF source must be valid");
    require(point.delta_chi2 < 0.0f && std::fabs(point.delta_chi2) < 0.01f,
            "ideal PSF source must prefer the PSF hypothesis");
    requireNear(point.orth_ext, 0.0f, 1.0e-5f,
                "ideal PSF source must have zero orthogonal extension");

    const std::vector<float> extended_source = makeExtendedSource(ns, psf, 1.0);
    const PointSourceStatisticsResult extended =
        PointSourceStatistics::measure(ns, extended_source, psf);
    require(extended.valid, "fixed-beta extended source must be valid");
    require(extended.delta_chi2 > 0.0f,
            "fixed-beta extended source must prefer the extended hypothesis");
    require(extended.orth_ext > 0.0f,
            "fixed-beta extended source must have positive extension projection");

    for (const double amplitude : {0.1, 10.0, 100.0}) {
        const PointSourceStatisticsResult scaled = PointSourceStatistics::measure(
            ns, makeExtendedSource(ns, psf, amplitude), psf);
        require(scaled.valid, "brightness-scaled extended source must be valid");
        requireNear(scaled.delta_chi2, extended.delta_chi2, 2.0e-5f,
                    "delta_chi2 must be brightness invariant");
        requireNear(scaled.orth_ext, extended.orth_ext, 2.0e-5f,
                    "orth_ext must be brightness invariant");
    }

    std::vector<float> noisy_source = extended_source;
    noisy_source[static_cast<std::size_t>((ns / 2) * ns + ns / 2 + 1)] = -0.25f;
    const PointSourceStatisticsResult noisy =
        PointSourceStatistics::measure(ns, noisy_source, psf);
    require(noisy.valid && std::isfinite(noisy.delta_chi2)
                && std::isfinite(noisy.orth_ext),
            "negative noise-subtracted source power must remain usable");

    std::vector<float> wrong_size(psf.size() - 1, 1.0f);
    const PointSourceStatisticsResult bad_size =
        PointSourceStatistics::measure(ns, wrong_size, psf);
    require(!bad_size.valid && bad_size.delta_chi2 == -999.0f
                && bad_size.orth_ext == -999.0f,
            "dimension mismatch must return the invalid sentinel");

    std::vector<float> no_boundary(psf.size(), 1.0f);
    const PointSourceStatisticsResult bad_window =
        PointSourceStatistics::measure(ns, no_boundary, no_boundary);
    require(!bad_window.valid,
            "PSF without a reliable-window boundary must be invalid");
}

}  // namespace

// ==========================================
// Function: Run the focused point-source statistics regression suite.
// Method: Validate compile-time catalog widths and execute synthetic cases.
// ==========================================
int main() {
    static_assert(LensingConfig::shear_cat_ncols == 28,
                  "the Stage 7 catalog must contain both morphology statistics");
    static_assert(LensingConfig::expo_cat_ncols == 29,
                  "the exposure catalog must append chi2 after 28 shear columns");

    runSyntheticCases();
    std::cout << "PointSourceStatistics synthetic tests passed\n";
    return EXIT_SUCCESS;
}
