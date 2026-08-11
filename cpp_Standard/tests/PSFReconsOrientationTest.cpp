#include "PCAImageLayout.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

// ==========================================
// Function: Stop the test program when a requirement is not met
// Method: Print a focused failure message and return a non-zero process code.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFRecons orientation test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Compare two reconstructed pixels with floating-point tolerance
// Method: Require finite operands and an absolute error below tolerance.
// ==========================================
void requireNear(float actual, float expected, float tolerance,
                 const std::string& message) {
    require(std::isfinite(actual) && std::isfinite(expected)
                && std::fabs(actual - expected) <= tolerance,
            message);
}

// ==========================================
// Function: Verify production PCA helpers preserve asymmetric row-major pixels
// Method: Copy and center a 5x3 residual whose row and column contributions
//         differ by two orders of magnitude.
// ==========================================
void testPcaInputFeatureOrder() {
    constexpr int nx = 5;
    constexpr int ny = 3;
    constexpr int pixelCount = nx * ny;
    std::vector<float> residual(pixelCount);
    std::vector<double> copied(pixelCount, 0.0);
    std::vector<double> mean(pixelCount, 0.5);
    std::vector<double> centered(pixelCount, 0.0);

    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx; ++col) {
            residual[row * nx + col] = static_cast<float>(100 * row + col);
        }
    }

    PCAImageLayout::copySampleRowMajor(
        residual.data(), pixelCount, copied.data());
    PCAImageLayout::centerSampleRowMajor(
        residual.data(), mean.data(), pixelCount, centered.data());
    for (int pixel = 0; pixel < pixelCount; ++pixel) {
        requireNear(static_cast<float>(copied[pixel]), residual[pixel], 0.0f,
                    "PCA sample copy transposed a pixel");
        requireNear(static_cast<float>(centered[pixel]), residual[pixel] - 0.5f,
                    0.0f, "PCA centering transposed a pixel");
    }
}

// ==========================================
// Function: Compare row-major PCA reconstruction with the legacy permutation
// Method: Build equivalent asymmetric mean/components in both feature orders,
//         reconstruct spatial images, and compare only the final image pixels.
// ==========================================
void testNumericalEquivalenceToLegacyPermutation() {
    constexpr int nx = 5;
    constexpr int ny = 3;
    constexpr int pixelCount = nx * ny;
    constexpr int modeCount = 2;

    std::vector<double> mean(pixelCount);
    std::vector<double> components(pixelCount * modeCount);
    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx; ++col) {
            const int pixel = row * nx + col;
            mean[pixel] = 100.0 * row + col;
            components[pixel * modeCount] = 0.25 * (row + 1) + col;
            components[pixel * modeCount + 1] = row - 0.5 * col;
        }
    }
    const std::vector<float> coefficients = {1.25f, -0.75f};

    std::vector<float> rowMajor(pixelCount, 0.0f);
    PCAImageLayout::reconstructRowMajor(
        mean.data(), components.data(), pixelCount, modeCount,
        coefficients.data(), rowMajor.data());

    std::vector<double> legacyMean(pixelCount, 0.0);
    std::vector<double> legacyComponents(pixelCount * modeCount, 0.0);
    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx; ++col) {
            const int rowMajorPixel = row * nx + col;
            const int legacyFeature = col * ny + row;
            legacyMean[legacyFeature] = mean[rowMajorPixel];
            for (int mode = 0; mode < modeCount; ++mode) {
                legacyComponents[legacyFeature * modeCount + mode] =
                    components[rowMajorPixel * modeCount + mode];
            }
        }
    }

    std::vector<float> legacyImage(pixelCount, 0.0f);
    for (int row = 0; row < ny; ++row) {
        for (int col = 0; col < nx; ++col) {
            const int legacyFeature = col * ny + row;
            double value = legacyMean[legacyFeature];
            for (int mode = 0; mode < modeCount; ++mode) {
                value += legacyComponents[legacyFeature * modeCount + mode]
                    * coefficients[mode];
            }
            legacyImage[row * nx + col] = static_cast<float>(value);
        }
    }

    for (int pixel = 0; pixel < pixelCount; ++pixel) {
        requireNear(rowMajor[pixel], legacyImage[pixel], 1.0e-6f,
                    "row-major PCA changed the reconstructed image");
    }
    require(rowMajor[1 * nx + 4] != rowMajor[2 * nx + 1],
            "synthetic reconstruction must remain orientation-sensitive");
}

}  // namespace

// ==========================================
// Function: Run focused PSF PCA orientation and equivalence regressions
// Method: Exercise the same row-major helpers used by production PSFRecons.
// ==========================================
int main() {
    testPcaInputFeatureOrder();
    testNumericalEquivalenceToLegacyPermutation();
    std::cout << "PSFRecons orientation tests passed\n";
    return EXIT_SUCCESS;
}
