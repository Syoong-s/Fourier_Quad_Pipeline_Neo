#ifndef NOISE_COVARIANCE_HPP
#define NOISE_COVARIANCE_HPP

#include <cstddef>
#include <vector>

namespace NoiseCovariance {

    bool computeMaskedCovarianceFFT(
        int regionSize, int paddedSize, int maxLag, double minPairFraction,
        const std::vector<double>& residual,
        const std::vector<unsigned char>& mask,
        std::vector<double>& covariance,
        std::vector<double>* pairCounts = nullptr);

    // ==========================================
    // Function: Detect a source stamp that straddles a two-amplifier boundary
    // Method: Compare the half-open stamp x interval with the chip midpoint only when split mode 2
    //         is active; single-amplifier mode never rejects on this geometry.
    // ==========================================
    bool sourceStampCrossesAmplifier(
        int chipWidth, int ccdSplit, int stampStartX, int stampSize);

    // ==========================================
    // Function: Identify an even FFT side composed only of factors 2, 3, and 5
    // Method: Expose the synthesis-grid performance predicate for deterministic tests.
    // ==========================================
    bool isFastEvenFFTSize(int size);

    // ==========================================
    // Function: Select the smallest FFT-friendly even side at least as large as requested
    // Method: Scan upward with overflow protection for dynamic correlated-fill grids.
    // ==========================================
    int nextFastEvenFFTSize(int requiredSize);

    // ==========================================
    // Function: Convert covariance into finite-stamp signed subtraction power
    // Method: Apply pair-window weighting and modulo accumulation before the normalized FFT.
    // ==========================================
    bool covarianceToFiniteStampNoisePower(
        int outputSize, int maxLag,
        const std::vector<double>& covariance,
        std::vector<float>& noisePower,
        double& maxImaginary,
        double& negativeFraction);

    // ==========================================
    // Function: Convert unwindowed covariance into dynamic-grid signed synthesis power
    // Method: Derive a no-wrap FFT-friendly grid and use a synthesis-only cached transform.
    // ==========================================
    bool covarianceToSynthesisPower(
        int stampSize, int maxLag,
        const std::vector<double>& covariance,
        int& synthesisSize,
        std::vector<float>& synthesisPower,
        double& maxImaginary);

    bool copyStoredNoisePower(
        const std::vector<float>& storedNoisePower,
        std::size_t offset,
        int stampSize,
        std::vector<float>& noisePower);

}

#endif // NOISE_COVARIANCE_HPP
