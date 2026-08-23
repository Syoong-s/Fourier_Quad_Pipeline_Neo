#include "ImageProcessing.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

    // ==========================================
    // Function: Compare two float power grids within scaled tolerance
    // Method: Apply a mixed absolute/relative bound and report the first differing Fourier pixel.
    // ==========================================
    bool vectorsNear(const std::vector<float>& first,
                     const std::vector<float>& second,
                     double tolerance,
                     const char* label) {
        if (first.size() != second.size()) {
            std::cerr << label << ": size mismatch\n";
            return false;
        }
        for (std::size_t i = 0; i < first.size(); ++i) {
            const double scale = std::max({
                1.0, std::abs(static_cast<double>(first[i])),
                std::abs(static_cast<double>(second[i]))});
            if (std::abs(static_cast<double>(first[i])
                         - static_cast<double>(second[i])) > tolerance * scale) {
                std::cerr << label << ": mismatch at " << i << "\n";
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Measure the production outer-edge mean
    // Method: Sum four boundary sides without corners and divide by 4*(n-2).
    // ==========================================
    double outerEdgeMean(int n, const std::vector<float>& power) {
        double mean = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            mean += power[static_cast<std::size_t>(i) * n]
                  + power[static_cast<std::size_t>(i) * n + (n - 1)]
                  + power[i]
                  + power[static_cast<std::size_t>(n - 1) * n + i];
        }
        return mean / (4.0 * static_cast<double>(n - 2));
    }

    // ==========================================
    // Function: Reproduce the former positive-only logarithmic smoother
    // Method: Use the legacy 1e-4 span offset around the unchanged hole smoother for regression.
    // ==========================================
    void legacyPositiveLogSmooth(int nx, int ny, std::vector<float>& power) {
        const auto extrema = std::minmax_element(power.begin(), power.end());
        const float offset = 1.0e-4f * (*extrema.second - *extrema.first);
        for (float& value : power) {
            value = std::log(value + offset);
        }
        ImageProcessing::smoothImage55Hole(nx, ny, power);
        for (float& value : power) {
            value = std::exp(value) - offset;
        }
    }

    // ==========================================
    // Function: Verify signed-safe and positive-only logarithmic smoothing
    // Method: Require finite output for mixed-sign/zero input and near-legacy output when the old
    //         positive-map offset was already valid.
    // ==========================================
    bool testSignedLogSmoothing() {
        constexpr int n = 8;
        std::vector<float> signedPower(n * n, 0.0f);
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                signedPower[static_cast<std::size_t>(y) * n + x]
                    = static_cast<float>(0.4 * std::sin(0.7 * x)
                                         - 0.6 * std::cos(0.4 * y));
            }
        }
        signedPower[0] = 0.0f;
        ImageProcessing::smoothImage55HoleLn(n, n, signedPower);
        if (!std::all_of(signedPower.begin(), signedPower.end(),
                         [](float value) { return std::isfinite(value); })) {
            return false;
        }

        std::vector<float> positive(n * n, 0.0f);
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                positive[static_cast<std::size_t>(y) * n + x]
                    = static_cast<float>(2.0 + 0.04 * x + 0.03 * y
                                         + 0.02 * std::sin(0.2 * x * y));
            }
        }
        std::vector<float> legacy = positive;
        legacyPositiveLogSmooth(n, n, legacy);
        ImageProcessing::smoothImage55HoleLn(n, n, positive);
        return vectorsNear(positive, legacy, 2.0e-6, "positive log regression");
    }

    // ==========================================
    // Function: Verify corrected-power construction for star and galaxy smoothing modes
    // Method: Compare the unified helper with explicit raw-FFT, subtract, smooth, edge order for
    //         modes 0/1/2 and prove smoothing-before-subtraction gives a different result.
    // ==========================================
    bool testCorrectedPowerOrder() {
        constexpr int n = 8;
        std::vector<float> source(n * n, 0.0f);
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                source[static_cast<std::size_t>(y) * n + x]
                    = static_cast<float>(0.2 + 0.03 * x - 0.01 * y
                                         + 0.05 * std::sin(0.5 * x + 0.2 * y));
            }
        }
        source[static_cast<std::size_t>(n / 2) * n + n / 2] += 5.0f;

        std::vector<float> rawPower;
        double rawPc = 0.0;
        ImageProcessing::getPower(n, n, source, rawPower, 0, rawPc);

        std::vector<float> storedNoise(rawPower.size(), 0.0f);
        for (std::size_t i = 0; i < rawPower.size(); ++i) {
            const int signedIndex = static_cast<int>(i % 9U) - 4;
            storedNoise[i] = 0.35f * rawPower[i]
                           + static_cast<float>(signedIndex) * 2.0e-4f;
        }

        std::vector<float> starPower;
        for (int smoothMode = 0; smoothMode <= 2; ++smoothMode) {
            std::vector<float> expected = rawPower;
            ImageProcessing::subtractNoisePower(n, expected, storedNoise);
            ImageProcessing::smoothPower(n, n, expected, smoothMode);
            ImageProcessing::subtractPowerEdgeMean(n, expected);

            std::vector<float> actual;
            double actualPc = 0.0;
            if (!ImageProcessing::buildCorrectedPower(
                    n, n, source, storedNoise, smoothMode, actual, actualPc)
                || std::abs(actualPc - rawPc) > 1.0e-12
                || !vectorsNear(actual, expected, 2.0e-7, "corrected power order")
                || std::abs(outerEdgeMean(n, actual)) > 2.0e-7) {
                return false;
            }

            if (smoothMode > 0) {
                std::vector<float> wrong = rawPower;
                ImageProcessing::smoothPower(n, n, wrong, smoothMode);
                ImageProcessing::subtractNoisePower(n, wrong, storedNoise);
                ImageProcessing::subtractPowerEdgeMean(n, wrong);
                double maximumDifference = 0.0;
                for (std::size_t i = 0; i < actual.size(); ++i) {
                    maximumDifference = std::max(
                        maximumDifference,
                        std::abs(static_cast<double>(actual[i])
                                 - static_cast<double>(wrong[i])));
                }
                if (maximumDifference <= 1.0e-7) return false;
            }

            if (smoothMode == 2) starPower = actual;
        }

        ImageProcessing::regularizePower(n, n, starPower, 2);
        const float center = starPower[static_cast<std::size_t>(n / 2) * n + n / 2];
        return std::all_of(starPower.begin(), starPower.end(),
                           [](float value) { return std::isfinite(value); })
            && std::abs(center - 1.0f) <= 2.0e-6f;
    }

}

// ==========================================
// Function: Run downstream corrected-power regression tests
// Method: Validate signed log safety, positive-only compatibility, raw source FFT behavior, exact
//         subtract-smooth-edge ordering, all galaxy modes, and final star regularization.
// ==========================================
int main() {
    if (!testSignedLogSmoothing()) {
        std::cerr << "signed-log smoothing tests failed\n";
        return 1;
    }
    if (!testCorrectedPowerOrder()) {
        std::cerr << "corrected-power order tests failed\n";
        return 1;
    }
    std::cout << "Power processing tests passed\n";
    return 0;
}

