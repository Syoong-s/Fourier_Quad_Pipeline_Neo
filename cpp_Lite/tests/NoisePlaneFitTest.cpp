#include "LensingConfig.hpp"
#include "NoisePlaneFit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

namespace {

    constexpr double coefficientTolerance = 2.0e-6;
    constexpr int ccdSplit = 2;

    // ==========================================
    // Function: Compare recovered and injected plane coefficients
    // Method: Apply one absolute tolerance suitable for float-valued synthetic images.
    // ==========================================
    bool coefficientsMatch(
        double aa, double bb, double cc,
        double expectedA, double expectedB, double expectedC) {
        return std::abs(aa - expectedA) <= coefficientTolerance
            && std::abs(bb - expectedB) <= coefficientTolerance
            && std::abs(cc - expectedC) <= coefficientTolerance;
    }

    // ==========================================
    // Function: Populate a local image with a plane in the production coordinate convention
    // Method: Evaluate A+B*x+C*y using local coordinates shifted by sourceOffset plus one.
    // ==========================================
    void fillPlane(
        std::vector<float>& image,
        int regionSize,
        int sourceOffset,
        double aa,
        double bb,
        double cc) {
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                const double planeX = static_cast<double>(x - sourceOffset + 1);
                const double planeY = static_cast<double>(y - sourceOffset + 1);
                image[static_cast<std::size_t>(y) * static_cast<std::size_t>(regionSize)
                      + static_cast<std::size_t>(x)] =
                    static_cast<float>(aa + bb * planeX + cc * planeY);
            }
        }
    }

    // ==========================================
    // Function: Verify configurable square-shell geometry and central-source exclusion
    // Method: Recover a plane for both production and alternate sizes while replacing every
    //         central-exclusion pixel with a large non-planar source signal.
    // ==========================================
    bool testSquareShellGeometry() {
        const int regionSizes[2] = {LensingConfig::noise_region_size, 48};
        const int innerSizes[2] = {LensingConfig::noise_inner_size, 20};
        for (int testIndex = 0; testIndex < 2; ++testIndex) {
            const int regionSize = regionSizes[testIndex];
            const int innerSize = innerSizes[testIndex];
            const int sourceOffset = regionSize / 2 - 7;
            const int localStartX = 8;
            const int localStartY = 6;
            const int chipWidth = 4 * regionSize;
            const int chipHeight = 2 * regionSize;
            const int sourceChipX = localStartX + regionSize / 2;
            const std::size_t elementCount = static_cast<std::size_t>(regionSize)
                                           * static_cast<std::size_t>(regionSize);
            std::vector<float> image(elementCount, 0.0f);
            std::vector<int> weight(elementCount, 1);
            constexpr double expectedA = 7.5;
            constexpr double expectedB = 0.03125;
            constexpr double expectedC = -0.046875;
            fillPlane(image, regionSize, sourceOffset, expectedA, expectedB, expectedC);

            const int innerStart = (regionSize - innerSize) / 2;
            const int innerEnd = innerStart + innerSize;
            for (int y = innerStart; y < innerEnd; ++y) {
                for (int x = innerStart; x < innerEnd; ++x) {
                    image[static_cast<std::size_t>(y) * static_cast<std::size_t>(regionSize)
                          + static_cast<std::size_t>(x)] =
                        static_cast<float>(1.0e5 + 17 * x - 23 * y + (x * y) % 31);
                }
            }

            double aa = 0.0;
            double bb = 0.0;
            double cc = 0.0;
            if (!NoisePlaneFit::fitNoiseRegionPlane(
                    image, weight, regionSize, innerSize, sourceOffset,
                    localStartX, localStartY, chipWidth, chipHeight, ccdSplit, sourceChipX,
                    aa, bb, cc)
                || !coefficientsMatch(
                    aa, bb, cc, expectedA, expectedB, expectedC)) {
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Verify masked and non-finite shell samples cannot bias the fit
    // Method: Corrupt deterministic shell subsets with weight 0/2 or NaN and recover the plane
    //         from the remaining finite weight-one samples.
    // ==========================================
    bool testMaskedAndNonFinitePixels() {
        const int regionSize = LensingConfig::noise_region_size;
        const int innerSize = LensingConfig::noise_inner_size;
        const int sourceOffset = regionSize / 2 - LensingConfig::nl_2;
        const int localStartX = 12;
        const int localStartY = 10;
        const int chipWidth = 4 * regionSize;
        const int chipHeight = 2 * regionSize;
        const int sourceChipX = localStartX + regionSize / 2;
        const std::size_t elementCount = static_cast<std::size_t>(regionSize)
                                       * static_cast<std::size_t>(regionSize);
        std::vector<float> image(elementCount, 0.0f);
        std::vector<int> weight(elementCount, 1);
        constexpr double expectedA = -2.25;
        constexpr double expectedB = 0.0125;
        constexpr double expectedC = 0.01875;
        fillPlane(image, regionSize, sourceOffset, expectedA, expectedB, expectedC);

        const int innerStart = (regionSize - innerSize) / 2;
        const int innerEnd = innerStart + innerSize;
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                if (x >= innerStart && x < innerEnd && y >= innerStart && y < innerEnd) {
                    continue;
                }
                const std::size_t index = static_cast<std::size_t>(y)
                                        * static_cast<std::size_t>(regionSize)
                                        + static_cast<std::size_t>(x);
                if ((3 * x + 5 * y) % 11 == 0) {
                    weight[index] = 2;
                    image[index] = 1.0e6f;
                } else if ((7 * x + 2 * y) % 17 == 0) {
                    image[index] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        }

        double aa = 0.0;
        double bb = 0.0;
        double cc = 0.0;
        return NoisePlaneFit::fitNoiseRegionPlane(
                   image, weight, regionSize, innerSize, sourceOffset,
                   localStartX, localStartY, chipWidth, chipHeight, ccdSplit, sourceChipX,
                   aa, bb, cc)
            && coefficientsMatch(aa, bb, cc, expectedA, expectedB, expectedC);
    }

    // ==========================================
    // Function: Verify another amplifier cannot contaminate the source plane
    // Method: Place a source-local region across the boundary, replace all opposite-amplifier
    //         pixels by two different planes in turn, and require invariant source coefficients.
    // ==========================================
    bool testAmplifierIsolation() {
        const int regionSize = LensingConfig::noise_region_size;
        const int innerSize = LensingConfig::noise_inner_size;
        const int sourceOffset = regionSize / 2 - LensingConfig::nl_2;
        const int chipWidth = 2 * regionSize;
        const int chipHeight = regionSize + 32;
        const int sourceChipX = chipWidth / 2 - 16;
        const int localStartX = sourceChipX - regionSize / 2;
        const int localStartY = 16;
        const std::size_t elementCount = static_cast<std::size_t>(regionSize)
                                       * static_cast<std::size_t>(regionSize);
        std::vector<float> image(elementCount, 0.0f);
        std::vector<int> weight(elementCount, 1);
        constexpr double expectedA = 3.0;
        constexpr double expectedB = -0.021;
        constexpr double expectedC = 0.017;

        for (int contaminationCase = 0; contaminationCase < 2; ++contaminationCase) {
            fillPlane(image, regionSize, sourceOffset, expectedA, expectedB, expectedC);
            for (int y = 0; y < regionSize; ++y) {
                for (int x = 0; x < regionSize; ++x) {
                    const int chipX = localStartX + x;
                    if (chipX < chipWidth / 2) {
                        continue;
                    }
                    const double scale = contaminationCase == 0 ? 1.0 : -3.0;
                    image[static_cast<std::size_t>(y) * static_cast<std::size_t>(regionSize)
                          + static_cast<std::size_t>(x)] =
                        static_cast<float>(scale * (1.0e4 + 9.0 * x - 13.0 * y));
                }
            }

            double aa = 0.0;
            double bb = 0.0;
            double cc = 0.0;
            if (!NoisePlaneFit::fitNoiseRegionPlane(
                    image, weight, regionSize, innerSize, sourceOffset,
                    localStartX, localStartY, chipWidth, chipHeight, ccdSplit, sourceChipX,
                    aa, bb, cc)
                || !coefficientsMatch(
                    aa, bb, cc, expectedA, expectedB, expectedC)) {
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Verify the runtime single-amplifier mode and split validation
    // Method: Recover one full-chip shell plane with ccdSplit=1 and reject an unsupported value.
    // ==========================================
    bool testRuntimeSplitModes() {
        const int regionSize = LensingConfig::noise_region_size;
        const int innerSize = LensingConfig::noise_inner_size;
        const int sourceOffset = regionSize / 2 - LensingConfig::nl_2;
        const int chipWidth = regionSize;
        const int chipHeight = regionSize;
        const int sourceChipX = regionSize / 2;
        const std::size_t elementCount = static_cast<std::size_t>(regionSize)
                                       * static_cast<std::size_t>(regionSize);
        std::vector<float> image(elementCount, 0.0f);
        std::vector<int> weight(elementCount, 1);
        constexpr double expectedA = 2.5;
        constexpr double expectedB = 0.015;
        constexpr double expectedC = -0.025;
        fillPlane(image, regionSize, sourceOffset, expectedA, expectedB, expectedC);

        double aa = 0.0;
        double bb = 0.0;
        double cc = 0.0;
        if (!NoisePlaneFit::fitNoiseRegionPlane(
                image, weight, regionSize, innerSize, sourceOffset,
                0, 0, chipWidth, chipHeight, 1, sourceChipX, aa, bb, cc)
            || !coefficientsMatch(aa, bb, cc, expectedA, expectedB, expectedC)) {
            return false;
        }
        return !NoisePlaneFit::fitNoiseRegionPlane(
            image, weight, regionSize, innerSize, sourceOffset,
            0, 0, chipWidth, chipHeight, 3, sourceChipX, aa, bb, cc);
    }

    // ==========================================
    // Function: Verify partial chip regions and minimum-valid-fraction rejection
    // Method: Recover from an edge-truncated shell while ignoring padded samples, then mask all
    //         but four geometric candidates and require a stable false return.
    // ==========================================
    bool testChipEdgeAndValidFraction() {
        const int regionSize = LensingConfig::noise_region_size;
        const int innerSize = LensingConfig::noise_inner_size;
        const int sourceOffset = regionSize / 2 - LensingConfig::nl_2;
        const int chipWidth = 2 * regionSize;
        const int chipHeight = regionSize + 24;
        const int sourceChipX = 50;
        const int sourceChipY = 50;
        const int localStartX = sourceChipX - regionSize / 2;
        const int localStartY = sourceChipY - regionSize / 2;
        const std::size_t elementCount = static_cast<std::size_t>(regionSize)
                                       * static_cast<std::size_t>(regionSize);
        std::vector<float> image(elementCount, 0.0f);
        std::vector<int> weight(elementCount, 1);
        constexpr double expectedA = 1.75;
        constexpr double expectedB = 0.009;
        constexpr double expectedC = -0.014;
        fillPlane(image, regionSize, sourceOffset, expectedA, expectedB, expectedC);

        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                const int chipX = localStartX + x;
                const int chipY = localStartY + y;
                if (chipX >= 0 && chipX < chipWidth && chipY >= 0 && chipY < chipHeight) {
                    continue;
                }
                image[static_cast<std::size_t>(y) * static_cast<std::size_t>(regionSize)
                      + static_cast<std::size_t>(x)] = 1.0e7f;
            }
        }

        double aa = 0.0;
        double bb = 0.0;
        double cc = 0.0;
        if (!NoisePlaneFit::fitNoiseRegionPlane(
                image, weight, regionSize, innerSize, sourceOffset,
                localStartX, localStartY, chipWidth, chipHeight, ccdSplit, sourceChipX,
                aa, bb, cc)
            || !coefficientsMatch(aa, bb, cc, expectedA, expectedB, expectedC)) {
            return false;
        }

        std::fill(weight.begin(), weight.end(), 0);
        const int shellLow = (regionSize - innerSize) / 2 + 2;
        const int shellHigh = (regionSize + innerSize) / 2 + 6;
        const int validCoordinates[4][2] = {{shellLow, shellHigh},
                                            {regionSize / 2, shellHigh},
                                            {shellHigh, shellLow},
                                            {shellHigh, regionSize / 2}};
        for (const auto& coordinate : validCoordinates) {
            const int x = coordinate[0];
            const int y = coordinate[1];
            weight[static_cast<std::size_t>(y) * static_cast<std::size_t>(regionSize)
                   + static_cast<std::size_t>(x)] = 1;
        }
        return !NoisePlaneFit::fitNoiseRegionPlane(
            image, weight, regionSize, innerSize, sourceOffset,
            localStartX, localStartY, chipWidth, chipHeight, ccdSplit, sourceChipX,
            aa, bb, cc);
    }

}

// ==========================================
// Function: Run Stage-3 outer-noise plane-fit regression tests
// Method: Exercise configurable shell geometry, masking, amplifier isolation, and edge handling.
// ==========================================
int main() {
    if (!testSquareShellGeometry()) {
        std::cerr << "square-shell geometry test failed\n";
        return 1;
    }
    if (!testMaskedAndNonFinitePixels()) {
        std::cerr << "masked/non-finite pixel test failed\n";
        return 1;
    }
    if (!testAmplifierIsolation()) {
        std::cerr << "amplifier-isolation test failed\n";
        return 1;
    }
    if (!testRuntimeSplitModes()) {
        std::cerr << "runtime split-mode test failed\n";
        return 1;
    }
    if (!testChipEdgeAndValidFraction()) {
        std::cerr << "chip-edge/valid-fraction test failed\n";
        return 1;
    }
    std::cout << "Noise plane-fit tests passed\n";
    return 0;
}
