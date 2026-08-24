#include "ImageProcessing.hpp"
#include "NoiseCovariance.hpp"
#include "NumericalRecipes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

    constexpr double pi = 3.14159265358979323846;

    // ==========================================
    // Function: Map a signed lag into the compact covariance-test grid
    // Method: Mirror the production (2*L+1)-side row-major layout.
    // ==========================================
    std::size_t lagIndex(int dx, int dy, int maxLag) {
        const int side = 2 * maxLag + 1;
        return static_cast<std::size_t>(dy + maxLag) * side + dx + maxLag;
    }

    // ==========================================
    // Function: Build a shifted signed target power grid for decoration tests
    // Method: Create one Hermitian anisotropic spectrum, optionally insert a symmetric negative
    //         pair, and place DC at the production center coordinate.
    // ==========================================
    std::vector<float> makeShiftedPower(int n, bool includeNegativePair) {
        std::vector<float> shifted(static_cast<std::size_t>(n) * n, 0.0f);
        const int half = n / 2;
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                double value = 1.0
                             + 0.55 * std::cos(2.0 * pi * x / n)
                             + 0.10 * std::cos(2.0 * pi * y / n);
                if (includeNegativePair
                    && ((x == 1 && y == 2) || (x == n - 1 && y == n - 2))) {
                    value = -0.4;
                }
                const int shiftedX = (x + half) % n;
                const int shiftedY = (y + half) % n;
                shifted[static_cast<std::size_t>(shiftedY) * n + shiftedX]
                    = static_cast<float>(value);
            }
        }
        return shifted;
    }

    // ==========================================
    // Function: Reconstruct the synthesis-only target PSD used by production
    // Method: Clip signed modes and renormalize the shifted grid so its sum equals C(0,0).
    // ==========================================
    std::vector<double> expectedFillPower(
        const std::vector<float>& storedPower, double zeroLagCovariance) {
        double positiveSum = 0.0;
        for (float value : storedPower) {
            positiveSum += std::max(0.0, static_cast<double>(value));
        }
        std::vector<double> expected(storedPower.size(), 0.0);
        for (std::size_t i = 0; i < storedPower.size(); ++i) {
            expected[i] = std::max(0.0, static_cast<double>(storedPower[i]))
                        * zeroLagCovariance / positiveSum;
        }
        return expected;
    }

    // ==========================================
    // Function: Verify the clean-stamp fast path
    // Method: Require bitwise stamp preservation and compare the next Gaussian draw across a
    //         reseed to prove that no-mask decoration consumes no random numbers.
    // ==========================================
    bool testNoMaskBypass() {
        constexpr int n = 8;
        std::vector<int> weights(n * n, 1);
        std::vector<float> stamp(n * n, 0.0f);
        for (std::size_t i = 0; i < stamp.size(); ++i) {
            stamp[i] = static_cast<float>(0.25 + 0.01 * i);
        }
        const std::vector<float> original = stamp;

        NumericalRecipes::seedRandom(73129U);
        const double expectedNextGaussian = NumericalRecipes::gasdev();
        NumericalRecipes::seedRandom(73129U);
        if (!ImageProcessing::decorateStampCorrelated(
                n, 0, {}, 1.7, weights, stamp)) {
            return false;
        }
        const double actualNextGaussian = NumericalRecipes::gasdev();
        return stamp == original && actualNextGaussian == expectedNextGaussian;
    }

    // ==========================================
    // Function: Verify same-coordinate masked-only replacement and failure policy
    // Method: Preserve signed power and every unmasked value exactly, require finite replacements,
    //         and reject a non-positive synthesis PSD without white-noise fallback.
    // ==========================================
    bool testMaskedOnlyReplacement() {
        constexpr int n = 8;
        constexpr int synthesisSize = 12;
        std::vector<float> synthesisPower = makeShiftedPower(synthesisSize, true);
        const std::vector<float> originalPower = synthesisPower;
        std::vector<int> weights(n * n, 1);
        weights[3] = 0;
        weights[37] = 0;
        std::vector<float> stamp(n * n, 0.0f);
        for (std::size_t i = 0; i < stamp.size(); ++i) {
            stamp[i] = static_cast<float>(1000.0 + i);
        }
        const std::vector<float> originalStamp = stamp;

        NumericalRecipes::seedRandom(991U);
        if (!ImageProcessing::decorateStampCorrelated(
                n, synthesisSize, synthesisPower, 2.3, weights, stamp)
            || synthesisPower != originalPower) {
            return false;
        }
        for (std::size_t i = 0; i < stamp.size(); ++i) {
            if (weights[i] != 0 && stamp[i] != originalStamp[i]) return false;
            if (weights[i] == 0
                && (!std::isfinite(stamp[i]) || stamp[i] == originalStamp[i])) {
                return false;
            }
        }

        std::vector<float> invalidPower(synthesisSize * synthesisSize, -1.0f);
        std::vector<float> rejectedStamp = originalStamp;
        return !ImageProcessing::decorateStampCorrelated(
                   n, synthesisSize, invalidPower, 2.3, weights, rejectedStamp)
            && rejectedStamp == originalStamp;
    }

    // ==========================================
    // Function: Exercise the existing primary-source mask protection contract
    // Method: Place a bad pixel at the core, next to the connected source, or only in the
    //         periphery and require respectively negative, negative, and non-negative outcomes.
    // ==========================================
    bool testPrimarySourceMaskProtection() {
        constexpr int n = 32;
        constexpr int center = n / 2;
        auto runCase = [](int badX, int badY, int& outputFlag,
                          std::vector<int>& outputWeight) {
            std::vector<float> stamp(n * n, 0.0f);
            outputWeight.assign(n * n, 1);
            stamp[static_cast<std::size_t>(center) * n + center] = 10.0f;
            outputWeight[static_cast<std::size_t>(badY) * n + badX] = 0;
            int boundx[2] = {0, 0};
            int boundy[2] = {0, 0};
            double totalFlux = 0.0;
            int totalArea = 0;
            double peak = 0.0;
            double halfLightFlux = 0.0;
            int halfLightArea = 0;
            double radius = 0.0;
            int xp = 0;
            int yp = 0;
            outputFlag = 0;
            ImageProcessing::markSource(
                n, stamp, outputWeight, 1.0, 2.0, 4.0,
                boundx, boundy, totalFlux, totalArea, peak,
                halfLightFlux, halfLightArea, outputFlag, radius, xp, yp);
        };

        int coreFlag = 0;
        std::vector<int> coreWeight;
        runCase(center, center, coreFlag, coreWeight);
        int adjacentFlag = 0;
        std::vector<int> adjacentWeight;
        runCase(center + 1, center, adjacentFlag, adjacentWeight);
        int peripheralFlag = 0;
        std::vector<int> peripheralWeight;
        runCase(center + 8, center, peripheralFlag, peripheralWeight);

        return coreFlag < 0 && adjacentFlag < 0 && peripheralFlag >= 0
            && std::any_of(peripheralWeight.begin(), peripheralWeight.end(),
                           [](int value) { return value == 0; });
    }

    // ==========================================
    // Function: Verify Fourier-power, variance, and anisotropy closure
    // Method: Generate a deterministic ensemble from a signed anisotropic stored spectrum and
    //         compare mean shifted power and real-space second moments with the clipped target.
    // ==========================================
    bool testEnsembleClosure() {
        constexpr int n = 8;
        constexpr int samples = 3000;
        constexpr double zeroLagCovariance = 1.7;
        const std::vector<float> storedPower = makeShiftedPower(n, true);
        const std::vector<float> originalPower = storedPower;
        const std::vector<double> expected = expectedFillPower(
            storedPower, zeroLagCovariance);
        std::vector<double> meanPower(n * n, 0.0);
        std::vector<int> weights(n * n, 0);
        double meanSecondMoment = 0.0;
        double meanCovarianceX = 0.0;
        double meanCovarianceY = 0.0;

        NumericalRecipes::seedRandom(481516U);
        for (int sample = 0; sample < samples; ++sample) {
            std::vector<float> realization(n * n, 0.0f);
            if (!ImageProcessing::decorateStampCorrelated(
                    n, n, storedPower, zeroLagCovariance, weights, realization)) {
                return false;
            }
            std::vector<float> power;
            double pc = 0.0;
            ImageProcessing::getPower(n, n, realization, power, 0, pc);
            for (std::size_t i = 0; i < power.size(); ++i) meanPower[i] += power[i];

            double secondMoment = 0.0;
            double covarianceX = 0.0;
            double covarianceY = 0.0;
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const double value = realization[static_cast<std::size_t>(y) * n + x];
                    secondMoment += value * value;
                    covarianceX += value * realization[
                        static_cast<std::size_t>(y) * n + (x + 1) % n];
                    covarianceY += value * realization[
                        static_cast<std::size_t>((y + 1) % n) * n + x];
                }
            }
            const double inverseElements = 1.0 / static_cast<double>(n * n);
            meanSecondMoment += secondMoment * inverseElements;
            meanCovarianceX += covarianceX * inverseElements;
            meanCovarianceY += covarianceY * inverseElements;
        }

        for (double& value : meanPower) value /= samples;
        meanSecondMoment /= samples;
        meanCovarianceX /= samples;
        meanCovarianceY /= samples;

        double expectedCovarianceX = 0.0;
        double expectedCovarianceY = 0.0;
        double meanPowerSum = 0.0;
        double maximumRelativeError = 0.0;
        const int half = n / 2;
        for (int shiftedY = 0; shiftedY < n; ++shiftedY) {
            const int y = (shiftedY + half) % n;
            for (int shiftedX = 0; shiftedX < n; ++shiftedX) {
                const int x = (shiftedX + half) % n;
                const std::size_t index = static_cast<std::size_t>(shiftedY) * n + shiftedX;
                meanPowerSum += meanPower[index];
                expectedCovarianceX += expected[index] * std::cos(2.0 * pi * x / n);
                expectedCovarianceY += expected[index] * std::cos(2.0 * pi * y / n);
                if (expected[index] > 1.0e-12) {
                    maximumRelativeError = std::max(
                        maximumRelativeError,
                        std::abs(meanPower[index] - expected[index]) / expected[index]);
                } else if (std::abs(meanPower[index]) > 1.0e-10) {
                    return false;
                }
            }
        }

        return storedPower == originalPower
            && maximumRelativeError < 0.16
            && std::abs(meanPowerSum - zeroLagCovariance) < 0.04
            && std::abs(meanSecondMoment - zeroLagCovariance) < 0.04
            && std::abs(meanCovarianceX - expectedCovarianceX) < 0.04
            && std::abs(meanCovarianceY - expectedCovarianceY) < 0.04
            && meanCovarianceX > meanCovarianceY + 0.10;
    }

    // ==========================================
    // Function: Verify dynamic synthesis, central-crop covariance, and finite-periodogram closure
    // Method: Transform one positive anisotropic compact covariance, generate cached 20-side fields,
    //         crop 16-side stamps, and compare ensemble covariance/power with production expectations,
    //         including a far-edge pair that would expose periodic wrap-around.
    // ==========================================
    bool testDynamicCropClosure() {
        constexpr int stampSize = 16;
        constexpr int maxLag = 3;
        constexpr int samples = 2500;
        constexpr double zeroLagCovariance = 1.7;
        constexpr double covarianceX = 0.18;
        constexpr double covarianceY = 0.07;
        const int lagSide = 2 * maxLag + 1;
        std::vector<double> covariance(
            static_cast<std::size_t>(lagSide) * lagSide, 0.0);
        covariance[lagIndex(0, 0, maxLag)] = zeroLagCovariance;
        covariance[lagIndex(1, 0, maxLag)] = covarianceX;
        covariance[lagIndex(-1, 0, maxLag)] = covarianceX;
        covariance[lagIndex(0, 1, maxLag)] = covarianceY;
        covariance[lagIndex(0, -1, maxLag)] = covarianceY;

        std::vector<float> finitePower;
        double finiteMaxImaginary = 0.0;
        double negativeFraction = 0.0;
        if (!NoiseCovariance::covarianceToFiniteStampNoisePower(
                stampSize, maxLag, covariance, finitePower,
                finiteMaxImaginary, negativeFraction)) {
            return false;
        }
        int synthesisSize = 0;
        std::vector<float> synthesisPower;
        double synthesisMaxImaginary = 0.0;
        if (!NoiseCovariance::covarianceToSynthesisPower(
                stampSize, maxLag, covariance, synthesisSize,
                synthesisPower, synthesisMaxImaginary)
            || synthesisSize != 20 || finiteMaxImaginary > 1.0e-12
            || synthesisMaxImaginary > 1.0e-12 || negativeFraction != 0.0) {
            return false;
        }
        const std::vector<float> originalFinitePower = finitePower;
        const std::vector<float> originalSynthesisPower = synthesisPower;

        std::vector<double> meanPower(finitePower.size(), 0.0);
        std::vector<int> weights(stampSize * stampSize, 0);
        double meanVariance = 0.0;
        double meanCovarianceX = 0.0;
        double meanCovarianceY = 0.0;
        double meanCovarianceTwoX = 0.0;
        double meanFarEdgeCovariance = 0.0;

        NumericalRecipes::seedRandom(8675309U);
        for (int sample = 0; sample < samples; ++sample) {
            std::vector<float> realization(stampSize * stampSize, 0.0f);
            if (!ImageProcessing::decorateStampCorrelated(
                    stampSize, synthesisSize, synthesisPower,
                    zeroLagCovariance, weights, realization)) {
                return false;
            }

            std::vector<float> measuredPower;
            double centralPower = 0.0;
            ImageProcessing::getPower(
                stampSize, stampSize, realization, measuredPower, 0, centralPower);
            for (std::size_t i = 0; i < meanPower.size(); ++i) {
                meanPower[i] += measuredPower[i];
            }

            double variance = 0.0;
            double covarianceOneX = 0.0;
            double covarianceOneY = 0.0;
            double covarianceTwoX = 0.0;
            double farEdge = 0.0;
            for (int y = 0; y < stampSize; ++y) {
                for (int x = 0; x < stampSize; ++x) {
                    const double value = realization[
                        static_cast<std::size_t>(y) * stampSize + x];
                    variance += value * value;
                    if (x + 1 < stampSize) {
                        covarianceOneX += value * realization[
                            static_cast<std::size_t>(y) * stampSize + x + 1];
                    }
                    if (y + 1 < stampSize) {
                        covarianceOneY += value * realization[
                            static_cast<std::size_t>(y + 1) * stampSize + x];
                    }
                    if (x + 2 < stampSize) {
                        covarianceTwoX += value * realization[
                            static_cast<std::size_t>(y) * stampSize + x + 2];
                    }
                }
                farEdge += realization[static_cast<std::size_t>(y) * stampSize]
                         * realization[static_cast<std::size_t>(y) * stampSize
                                       + stampSize - 1];
            }
            meanVariance += variance / static_cast<double>(stampSize * stampSize);
            meanCovarianceX += covarianceOneX
                / static_cast<double>(stampSize * (stampSize - 1));
            meanCovarianceY += covarianceOneY
                / static_cast<double>(stampSize * (stampSize - 1));
            meanCovarianceTwoX += covarianceTwoX
                / static_cast<double>(stampSize * (stampSize - 2));
            meanFarEdgeCovariance += farEdge / stampSize;
        }

        for (double& value : meanPower) value /= samples;
        meanVariance /= samples;
        meanCovarianceX /= samples;
        meanCovarianceY /= samples;
        meanCovarianceTwoX /= samples;
        meanFarEdgeCovariance /= samples;

        double maximumRelativePowerError = 0.0;
        for (std::size_t i = 0; i < meanPower.size(); ++i) {
            if (!(finitePower[i] > 0.0f)) return false;
            maximumRelativePowerError = std::max(
                maximumRelativePowerError,
                std::abs(meanPower[i] - finitePower[i]) / finitePower[i]);
        }

        return finitePower == originalFinitePower
            && synthesisPower == originalSynthesisPower
            && maximumRelativePowerError < 0.12
            && std::abs(meanVariance - zeroLagCovariance) < 0.035
            && std::abs(meanCovarianceX - covarianceX) < 0.035
            && std::abs(meanCovarianceY - covarianceY) < 0.035
            && std::abs(meanCovarianceTwoX) < 0.035
            && std::abs(meanFarEdgeCovariance) < 0.035
            && meanCovarianceX > meanCovarianceY + 0.07;
    }

}

// ==========================================
// Function: Run correlated mask-decoration regression tests
// Method: Validate fast-path RNG preservation, masked-only writes, hard failure, source protection,
//         signed-power immutability, FFT normalization, variance closure, and anisotropy recovery.
// ==========================================
int main() {
    if (!testNoMaskBypass()) {
        std::cerr << "no-mask bypass test failed\n";
        return 1;
    }
    if (!testMaskedOnlyReplacement()) {
        std::cerr << "masked-only replacement test failed\n";
        return 1;
    }
    if (!testPrimarySourceMaskProtection()) {
        std::cerr << "primary-source mask protection test failed\n";
        return 1;
    }
    if (!testEnsembleClosure()) {
        std::cerr << "correlated ensemble closure test failed\n";
        return 1;
    }
    if (!testDynamicCropClosure()) {
        std::cerr << "dynamic crop closure test failed\n";
        return 1;
    }
    std::cout << "Correlated decoration tests passed\n";
    return 0;
}
