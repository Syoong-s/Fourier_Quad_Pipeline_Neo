#include "general/NumericalRecipes.hpp"
#include "process_main/PreProcessLegacy.hpp"
#include "process_main/UniversalUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

// ==========================================
// Function: Stop the legacy-preprocessing test on one failed invariant
// Method: Emit a focused diagnostic and terminate with failure.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "legacy preprocessing test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Build a deterministic non-degenerate amplifier image
// Method: Combine a smooth polynomial background with bounded texture.
// ==========================================
std::vector<float> makeImage(int nx, int ny) {
    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 0.0f);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            image[static_cast<std::size_t>(y) * nx + x] = static_cast<float>(
                30.0 + 0.04 * x - 0.03 * y + 0.0007 * x * y
                + 1.2 * std::sin(0.41 * x + 0.17 * y)
                + 0.8 * std::cos(0.13 * x - 0.37 * y));
        }
    }
    return image;
}

// ==========================================
// Function: Verify F77 background reproducibility and metadata round trip
// Method: Repeat one seeded fit and rebuild every subtraction from the
//         normalized coefficients consumed by Stage 3.
// ==========================================
void testBackgroundRoundTrip() {
    constexpr int nx = 128;
    constexpr int ny = 256;
    constexpr int x_start = nx / 2;
    constexpr int x_end = nx;
    constexpr int y_start = 0;
    constexpr int y_end = ny;
    const std::vector<float> raw = makeImage(nx, ny);
    std::vector<float> first = raw;
    std::vector<float> second = raw;
    std::vector<double> first_coefficients;
    std::vector<double> second_coefficients;
    NumericalRecipes::seedRandom(24681357U);
    require(PreProcessLegacy::setBackground(
                x_start, x_end, y_start, y_end, nx, ny, first, 20, 12, 3,
                first_coefficients),
            "seeded F77 background fit must succeed");
    NumericalRecipes::seedRandom(24681357U);
    require(PreProcessLegacy::setBackground(
                x_start, x_end, y_start, y_end, nx, ny, second, 20, 12, 3,
                second_coefficients),
            "repeated seeded F77 background fit must succeed");
    require(first == second && first_coefficients == second_coefficients,
            "same seed and image must reproduce the exact background result");

    const double x_mid = 0.5 * (x_start + 1.0 + x_end);
    const double y_mid = 0.5 * (y_start + 1.0 + y_end);
    const double x_half_inv = 2.0 / (x_end - x_start - 1.0);
    const double y_half_inv = 2.0 / (y_end - y_start - 1.0);
    double maximum_error = 0.0;
    for (int y = y_start; y < y_end; ++y) {
        for (int x = x_start; x < x_end; ++x) {
            const double normalized_x = (x + 1.0 - x_mid) * x_half_inv;
            const double normalized_y = (y + 1.0 - y_mid) * y_half_inv;
            const double background = UniversalUtils::funcVal(
                normalized_x, normalized_y, 12, 3, first_coefficients);
            const double rebuilt = static_cast<double>(raw[
                static_cast<std::size_t>(y) * nx + x]) - background;
            maximum_error = std::max(
                maximum_error,
                std::abs(rebuilt - first[static_cast<std::size_t>(y) * nx + x]));
        }
    }
    require(maximum_error < 2.0e-5,
            "serialized total background must reproduce the applied subtraction");
}

// ==========================================
// Function: Verify F77 sigma reproducibility and coefficient round trip
// Method: Repeat one seeded fit and reconstruct the exact published plane.
// ==========================================
void testSigmaRoundTrip() {
    constexpr int nx = 128;
    constexpr int ny = 256;
    constexpr int x_start = nx / 2;
    constexpr int x_end = nx;
    constexpr int y_start = 0;
    constexpr int y_end = ny;
    const std::vector<float> input = makeImage(nx, ny);
    std::vector<float> first = input;
    std::vector<float> second = input;
    double aa1 = 0.0;
    double bb1 = 0.0;
    double cc1 = 0.0;
    double aa2 = 0.0;
    double bb2 = 0.0;
    double cc2 = 0.0;
    NumericalRecipes::seedRandom(97531U);
    require(PreProcessLegacy::setSig(
                x_start, x_end, y_start, y_end, nx, ny,
                first, aa1, bb1, cc1),
            "seeded F77 sigma fit must succeed");
    NumericalRecipes::seedRandom(97531U);
    require(PreProcessLegacy::setSig(
                x_start, x_end, y_start, y_end, nx, ny,
                second, aa2, bb2, cc2),
            "repeated seeded F77 sigma fit must succeed");
    require(first == second && aa1 == aa2 && bb1 == bb2 && cc1 == cc2,
            "same seed and image must reproduce the exact sigma result");

    double maximum_error = 0.0;
    for (int y = y_start; y < y_end; ++y) {
        for (int x = x_start; x < x_end; ++x) {
            const double plane = aa1 + bb1 * (x + 1.0) + cc1 * (y + 1.0);
            const double rebuilt = input[static_cast<std::size_t>(y) * nx + x]
                / std::sqrt(0.5 * plane);
            maximum_error = std::max(
                maximum_error,
                std::abs(rebuilt - first[static_cast<std::size_t>(y) * nx + x]));
        }
    }
    require(maximum_error < 1.0e-5,
            "published sigma coefficients must reproduce normalization; max error="
                + std::to_string(maximum_error));
}

}  // namespace

// ==========================================
// Function: Run focused historical preprocessing regressions
// Method: Verify deterministic sampling and both downstream coefficient contracts.
// ==========================================
int main() {
    testBackgroundRoundTrip();
    testSigmaRoundTrip();
    std::cout << "legacy preprocessing tests passed\n";
    return EXIT_SUCCESS;
}
