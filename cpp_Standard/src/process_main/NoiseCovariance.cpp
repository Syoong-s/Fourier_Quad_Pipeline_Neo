#include "process_main/NoiseCovariance.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include <fftw3.h>

namespace NoiseCovariance {

    namespace {

        // ==========================================
        // Class: Reusable zero-padded autocorrelation workspace
        // Method: Bind fixed-size FFTW plans to two complex buffers and reuse them for every
        //         residual/mask pair evaluated by one execution thread.
        // ==========================================
        class AutocorrelationWorkspace {
        public:
            // ==========================================
            // Function: Construct fixed-size forward/backward autocorrelation plans
            // Method: Allocate caller-validated linear-correlation padding and bind in-place FFTW
            //         transforms for the residual and binary mask buffers.
            // ==========================================
            AutocorrelationWorkspace(int regionSize, int paddedSize)
                : regionSize_(regionSize),
                  paddedSize_(paddedSize),
                  residualBuffer_(elementCount()),
                  maskBuffer_(elementCount()) {
                residualForward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, residualData(), residualData(),
                    FFTW_FORWARD, FFTW_ESTIMATE);
                residualBackward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, residualData(), residualData(),
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                maskForward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, maskData(), maskData(),
                    FFTW_FORWARD, FFTW_ESTIMATE);
                maskBackward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, maskData(), maskData(),
                    FFTW_BACKWARD, FFTW_ESTIMATE);
            }

            // ==========================================
            // Function: Destroy reusable autocorrelation plans
            // Method: Release every successfully constructed FFTW plan while vector storage is
            //         still alive.
            // ==========================================
            ~AutocorrelationWorkspace() {
                if (residualForward_ != nullptr) fftw_destroy_plan(residualForward_);
                if (residualBackward_ != nullptr) fftw_destroy_plan(residualBackward_);
                if (maskForward_ != nullptr) fftw_destroy_plan(maskForward_);
                if (maskBackward_ != nullptr) fftw_destroy_plan(maskBackward_);
            }

            AutocorrelationWorkspace(const AutocorrelationWorkspace&) = delete;
            AutocorrelationWorkspace& operator=(const AutocorrelationWorkspace&) = delete;

            // ==========================================
            // Function: Report whether all autocorrelation plans were created
            // Method: Require the four bound FFTW plan handles used by execute().
            // ==========================================
            bool valid() const {
                return residualForward_ != nullptr && residualBackward_ != nullptr
                    && maskForward_ != nullptr && maskBackward_ != nullptr;
            }

            // ==========================================
            // Function: Compute linear residual and mask autocorrelations
            // Method: Zero-pad at the origin, transform both fields, replace each spectrum by
            //         its squared magnitude, and inverse transform in place.
            // ==========================================
            void execute(const std::vector<double>& residual,
                         const std::vector<unsigned char>& mask) {
                std::fill(residualBuffer_.begin(), residualBuffer_.end(),
                          std::complex<double>(0.0, 0.0));
                std::fill(maskBuffer_.begin(), maskBuffer_.end(),
                          std::complex<double>(0.0, 0.0));

                for (int y = 0; y < regionSize_; ++y) {
                    for (int x = 0; x < regionSize_; ++x) {
                        const std::size_t inputIndex = static_cast<std::size_t>(y)
                                                     * static_cast<std::size_t>(regionSize_)
                                                     + static_cast<std::size_t>(x);
                        const std::size_t paddedIndex = static_cast<std::size_t>(y)
                                                      * static_cast<std::size_t>(paddedSize_)
                                                      + static_cast<std::size_t>(x);
                        if (mask[inputIndex] != 0U) {
                            residualBuffer_[paddedIndex] = residual[inputIndex];
                            maskBuffer_[paddedIndex] = 1.0;
                        }
                    }
                }

                fftw_execute(residualForward_);
                fftw_execute(maskForward_);
                for (std::size_t i = 0; i < elementCount(); ++i) {
                    residualBuffer_[i] = std::complex<double>(
                        std::norm(residualBuffer_[i]), 0.0);
                    maskBuffer_[i] = std::complex<double>(
                        std::norm(maskBuffer_[i]), 0.0);
                }
                fftw_execute(residualBackward_);
                fftw_execute(maskBackward_);
            }

            // ==========================================
            // Function: Read one normalized autocorrelation numerator
            // Method: Map a signed linear lag to the periodic zero-padded FFT index and remove
            //         FFTW's unnormalized inverse-transform scale.
            // ==========================================
            double residualCorrelation(int dx, int dy) const {
                return residualBuffer_[lagIndex(dx, dy)].real()
                     / static_cast<double>(elementCount());
            }

            // ==========================================
            // Function: Read one normalized valid-pair count
            // Method: Apply the same lag mapping and inverse scale as the residual correlation.
            // ==========================================
            double maskCorrelation(int dx, int dy) const {
                return maskBuffer_[lagIndex(dx, dy)].real()
                     / static_cast<double>(elementCount());
            }

            // ==========================================
            // Function: Return the configured source-region side length
            // Method: Expose the immutable plan key for thread-local cache reuse.
            // ==========================================
            int regionSize() const {
                return regionSize_;
            }

            // ==========================================
            // Function: Return the configured autocorrelation FFT side length
            // Method: Expose the second immutable plan key for thread-local cache reuse.
            // ==========================================
            int paddedSize() const {
                return paddedSize_;
            }

        private:
            // ==========================================
            // Function: Count complex elements in one padded FFT buffer
            // Method: Square the fixed padded side in size_t arithmetic.
            // ==========================================
            std::size_t elementCount() const {
                return static_cast<std::size_t>(paddedSize_)
                     * static_cast<std::size_t>(paddedSize_);
            }

            // ==========================================
            // Function: Convert a signed lag to a padded row-major index
            // Method: Wrap negative x/y lags to the upper periodic indices of the full
            //         2*N-1 linear-correlation grid.
            // ==========================================
            std::size_t lagIndex(int dx, int dy) const {
                const int x = dx >= 0 ? dx : paddedSize_ + dx;
                const int y = dy >= 0 ? dy : paddedSize_ + dy;
                return static_cast<std::size_t>(y) * static_cast<std::size_t>(paddedSize_)
                     + static_cast<std::size_t>(x);
            }

            // ==========================================
            // Function: Access residual storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* residualData() {
                return reinterpret_cast<fftw_complex*>(residualBuffer_.data());
            }

            // ==========================================
            // Function: Access mask storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* maskData() {
                return reinterpret_cast<fftw_complex*>(maskBuffer_.data());
            }

            int regionSize_ = 0;
            int paddedSize_ = 0;
            std::vector<std::complex<double>> residualBuffer_;
            std::vector<std::complex<double>> maskBuffer_;
            fftw_plan residualForward_ = nullptr;
            fftw_plan residualBackward_ = nullptr;
            fftw_plan maskForward_ = nullptr;
            fftw_plan maskBackward_ = nullptr;
        };

        // ==========================================
        // Class: Reusable covariance-to-power workspace
        // Method: Bind one fixed-size forward FFTW plan to the output covariance grid.
        // ==========================================
        class PowerWorkspace {
        public:
            // ==========================================
            // Function: Construct a fixed-size covariance transform
            // Method: Allocate one square complex grid and bind an in-place forward FFTW plan.
            // ==========================================
            explicit PowerWorkspace(int outputSize)
                : outputSize_(outputSize),
                  buffer_(static_cast<std::size_t>(outputSize)
                          * static_cast<std::size_t>(outputSize)) {
                forward_ = fftw_plan_dft_2d(
                    outputSize_, outputSize_, data(), data(), FFTW_FORWARD, FFTW_ESTIMATE);
            }

            // ==========================================
            // Function: Destroy the reusable covariance transform
            // Method: Release the FFTW plan before its bound vector storage is destroyed.
            // ==========================================
            ~PowerWorkspace() {
                if (forward_ != nullptr) fftw_destroy_plan(forward_);
            }

            PowerWorkspace(const PowerWorkspace&) = delete;
            PowerWorkspace& operator=(const PowerWorkspace&) = delete;

            // ==========================================
            // Function: Report whether the covariance transform plan exists
            // Method: Validate the one forward FFTW plan required by execute().
            // ==========================================
            bool valid() const {
                return forward_ != nullptr;
            }

            // ==========================================
            // Function: Transform an unshifted covariance grid
            // Method: Copy the caller grid into bound storage and execute the reusable forward
            //         transform without altering signed Fourier modes.
            // ==========================================
            void execute(const std::vector<std::complex<double>>& covarianceGrid) {
                std::copy(covarianceGrid.begin(), covarianceGrid.end(), buffer_.begin());
                fftw_execute(forward_);
            }

            // ==========================================
            // Function: Access one transformed Fourier coefficient
            // Method: Read the immutable result buffer by row-major index.
            // ==========================================
            const std::complex<double>& value(std::size_t index) const {
                return buffer_[index];
            }

            // ==========================================
            // Function: Return the configured Fourier-grid side length
            // Method: Expose the immutable plan key for thread-local cache reuse.
            // ==========================================
            int outputSize() const {
                return outputSize_;
            }

        private:
            // ==========================================
            // Function: Access power-grid storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* data() {
                return reinterpret_cast<fftw_complex*>(buffer_.data());
            }

            int outputSize_ = 0;
            std::vector<std::complex<double>> buffer_;
            fftw_plan forward_ = nullptr;
        };

        // ==========================================
        // Function: Reuse the current thread's autocorrelation workspace
        // Method: Rebuild the cached workspace when either the region or padded FFT side changes.
        // ==========================================
        AutocorrelationWorkspace& autocorrelationWorkspace(
            int regionSize, int paddedSize) {
            thread_local std::unique_ptr<AutocorrelationWorkspace> workspace;
            if (!workspace || workspace->regionSize() != regionSize
                || workspace->paddedSize() != paddedSize) {
                workspace = std::make_unique<AutocorrelationWorkspace>(
                    regionSize, paddedSize);
            }
            return *workspace;
        }

        // ==========================================
        // Function: Reuse the current thread's finite-stamp covariance-to-power workspace
        // Method: Keep subtraction-plan ownership separate from dynamic synthesis transforms.
        // ==========================================
        PowerWorkspace& finitePowerWorkspace(int outputSize) {
            thread_local std::unique_ptr<PowerWorkspace> workspace;
            if (!workspace || workspace->outputSize() != outputSize) {
                workspace = std::make_unique<PowerWorkspace>(outputSize);
            }
            return *workspace;
        }

        // ==========================================
        // Function: Reuse the current thread's synthesis covariance-to-power workspace
        // Method: Cache the dynamic fill-grid plan independently so alternating ns/G transforms
        //         never destroy and rebuild one another.
        // ==========================================
        PowerWorkspace& synthesisPowerWorkspace(int synthesisSize) {
            thread_local std::unique_ptr<PowerWorkspace> workspace;
            if (!workspace || workspace->outputSize() != synthesisSize) {
                workspace = std::make_unique<PowerWorkspace>(synthesisSize);
            }
            return *workspace;
        }

        // ==========================================
        // Function: Map a signed lag to the compact covariance-vector index
        // Method: Offset both axes by maxLag in a (2*maxLag+1)^2 row-major grid.
        // ==========================================
        std::size_t compactLagIndex(int dx, int dy, int maxLag) {
            const int side = 2 * maxLag + 1;
            return static_cast<std::size_t>(dy + maxLag) * static_cast<std::size_t>(side)
                 + static_cast<std::size_t>(dx + maxLag);
        }

        // ==========================================
        // Function: Map a signed lag onto one periodic Fourier-grid coordinate
        // Method: Normalize the C++ remainder into [0, period) for both positive and negative lags.
        // ==========================================
        int positiveModulo(int value, int period) {
            const int remainder = value % period;
            return remainder >= 0 ? remainder : remainder + period;
        }

    }

    // ==========================================
    // Function: Identify an even FFT side composed only of factors 2, 3, and 5
    // Method: Reject non-positive/odd sides, divide out FFTW-friendly factors, and require unity.
    // ==========================================
    bool isFastEvenFFTSize(int size) {
        if (size <= 0 || size % 2 != 0) {
            return false;
        }
        int remainder = size;
        for (int factor : {2, 3, 5}) {
            while (remainder % factor == 0) {
                remainder /= factor;
            }
        }
        return remainder == 1;
    }

    // ==========================================
    // Function: Select the smallest even 2/3/5-smooth FFT side at least as large as requested
    // Method: Round upward to even and scan by two with explicit signed-integer overflow guards.
    // ==========================================
    int nextFastEvenFFTSize(int requiredSize) {
        if (requiredSize <= 0) {
            return 0;
        }
        if (requiredSize == std::numeric_limits<int>::max()) {
            return 0;
        }
        int candidate = requiredSize + requiredSize % 2;
        while (!isFastEvenFFTSize(candidate)) {
            if (candidate > std::numeric_limits<int>::max() - 2) {
                return 0;
            }
            candidate += 2;
        }
        return candidate;
    }

    // ==========================================
    // Function: Estimate a masked two-dimensional covariance with FFT acceleration
    // Method: Compute full zero-padded linear autocorrelations for M*R and M, divide by valid
    //         pair counts, retain qualified short lags, and enforce C(d)=C(-d) symmetry.
    // ==========================================
    bool computeMaskedCovarianceFFT(
        int regionSize, int paddedSize, int maxLag, double minPairFraction,
        const std::vector<double>& residual,
        const std::vector<unsigned char>& mask,
        std::vector<double>& covariance,
        std::vector<double>* pairCounts) {
        if (regionSize <= 0 || paddedSize < 2 * regionSize - 1
            || maxLag < 0 || maxLag >= regionSize
            || !std::isfinite(minPairFraction)
            || minPairFraction <= 0.0 || minPairFraction > 1.0) {
            return false;
        }

        const std::size_t regionElements = static_cast<std::size_t>(regionSize)
                                         * static_cast<std::size_t>(regionSize);
        if (residual.size() != regionElements || mask.size() != regionElements) {
            return false;
        }
        for (std::size_t i = 0; i < regionElements; ++i) {
            if (mask[i] != 0U && !std::isfinite(residual[i])) {
                return false;
            }
        }

        AutocorrelationWorkspace& workspace = autocorrelationWorkspace(
            regionSize, paddedSize);
        if (!workspace.valid()) {
            return false;
        }
        workspace.execute(residual, mask);

        const int lagSide = 2 * maxLag + 1;
        const std::size_t lagElements = static_cast<std::size_t>(lagSide)
                                      * static_cast<std::size_t>(lagSide);
        covariance.assign(lagElements, 0.0);
        std::vector<double> localPairCounts(lagElements, 0.0);

        const double zeroLagPairs = workspace.maskCorrelation(0, 0);
        if (!std::isfinite(zeroLagPairs) || zeroLagPairs < 1.0) {
            return false;
        }
        const double minimumPairs = minPairFraction * zeroLagPairs;

        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                const std::size_t index = compactLagIndex(dx, dy, maxLag);
                const double pairs = std::max(0.0, workspace.maskCorrelation(dx, dy));
                localPairCounts[index] = pairs;
                if (pairs + 1.0e-8 < minimumPairs) {
                    continue;
                }
                const double numerator = workspace.residualCorrelation(dx, dy);
                const double value = numerator / pairs;
                if (!std::isfinite(value)) {
                    return false;
                }
                covariance[index] = value;
            }
        }

        for (int dy = 0; dy <= maxLag; ++dy) {
            const int dxStart = dy == 0 ? 0 : -maxLag;
            for (int dx = dxStart; dx <= maxLag; ++dx) {
                const std::size_t forwardIndex = compactLagIndex(dx, dy, maxLag);
                const std::size_t reverseIndex = compactLagIndex(-dx, -dy, maxLag);
                const bool forwardValid = localPairCounts[forwardIndex] + 1.0e-8 >= minimumPairs;
                const bool reverseValid = localPairCounts[reverseIndex] + 1.0e-8 >= minimumPairs;
                if (forwardValid && reverseValid) {
                    const double symmetricValue = 0.5 * (
                        covariance[forwardIndex] + covariance[reverseIndex]);
                    covariance[forwardIndex] = symmetricValue;
                    covariance[reverseIndex] = symmetricValue;
                } else {
                    covariance[forwardIndex] = 0.0;
                    covariance[reverseIndex] = 0.0;
                }
            }
        }

        if (pairCounts != nullptr) {
            *pairCounts = std::move(localPairCounts);
        }
        return true;
    }

    // ==========================================
    // Function: Detect a source stamp that straddles a two-amplifier boundary
    // Method: Treat the stamp x range as half-open and compare it with the chip midpoint only in
    //         split mode 2; invalid geometry is conservatively rejected.
    // ==========================================
    bool sourceStampCrossesAmplifier(
        int chipWidth, int ccdSplit, int stampStartX, int stampSize) {
        if (ccdSplit != 2) {
            return false;
        }
        if (chipWidth <= 0 || stampSize <= 0 || stampStartX < 0
            || stampStartX > chipWidth - stampSize) {
            return true;
        }
        const int amplifierBoundary = chipWidth / 2;
        const int stampEndX = stampStartX + stampSize;
        return stampStartX < amplifierBoundary && stampEndX > amplifierBoundary;
    }

    // ==========================================
    // Function: Convert covariance into the expected power of one finite source stamp
    // Method: Apply the finite pair window, accumulate signed lags modulo the stamp side, run the
    //         normalized forward transform, fftshift, and retain signed modes for subtraction.
    // ==========================================
    bool covarianceToFiniteStampNoisePower(
        int outputSize, int maxLag,
        const std::vector<double>& covariance,
        std::vector<float>& noisePower,
        double& maxImaginary,
        double& negativeFraction) {
        if (outputSize <= 0 || maxLag < 0
            || maxLag > (std::numeric_limits<int>::max() - 1) / 2) {
            return false;
        }
        const int lagSide = 2 * maxLag + 1;
        const std::size_t expectedCovariance = static_cast<std::size_t>(lagSide)
                                               * static_cast<std::size_t>(lagSide);
        if (covariance.size() != expectedCovariance
            || !std::all_of(covariance.begin(), covariance.end(),
                            [](double value) { return std::isfinite(value); })) {
            return false;
        }

        const std::size_t outputElements = static_cast<std::size_t>(outputSize)
                                         * static_cast<std::size_t>(outputSize);
        std::vector<std::complex<double>> covarianceGrid(
            outputElements, std::complex<double>(0.0, 0.0));
        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            const int absoluteY = std::abs(dy);
            if (absoluteY >= outputSize) continue;
            const double windowY = 1.0 - static_cast<double>(absoluteY) / outputSize;
            const int y = positiveModulo(dy, outputSize);
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                const int absoluteX = std::abs(dx);
                if (absoluteX >= outputSize) continue;
                const double windowX = 1.0 - static_cast<double>(absoluteX) / outputSize;
                const int x = positiveModulo(dx, outputSize);
                covarianceGrid[static_cast<std::size_t>(y)
                               * static_cast<std::size_t>(outputSize)
                               + static_cast<std::size_t>(x)]
                    += covariance[compactLagIndex(dx, dy, maxLag)] * windowX * windowY;
            }
        }

        PowerWorkspace& workspace = finitePowerWorkspace(outputSize);
        if (!workspace.valid()) {
            return false;
        }
        workspace.execute(covarianceGrid);

        noisePower.assign(outputElements, 0.0f);
        maxImaginary = 0.0;
        negativeFraction = 0.0;
        double negativeMagnitude = 0.0;
        double totalMagnitude = 0.0;
        const double normalization = 1.0 / static_cast<double>(outputElements);
        const int half = outputSize / 2;

        for (int y = 0; y < outputSize; ++y) {
            const int shiftedY = (y + half) % outputSize;
            for (int x = 0; x < outputSize; ++x) {
                const int shiftedX = (x + half) % outputSize;
                const std::size_t inputIndex = static_cast<std::size_t>(y)
                                             * static_cast<std::size_t>(outputSize)
                                             + static_cast<std::size_t>(x);
                const std::size_t outputIndex = static_cast<std::size_t>(shiftedY)
                                              * static_cast<std::size_t>(outputSize)
                                              + static_cast<std::size_t>(shiftedX);
                const std::complex<double> value = workspace.value(inputIndex) * normalization;
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                    return false;
                }
                maxImaginary = std::max(maxImaginary, std::abs(value.imag()));
                totalMagnitude += std::abs(value.real());
                if (value.real() < 0.0) {
                    negativeMagnitude += -value.real();
                }
                noisePower[outputIndex] = static_cast<float>(value.real());
            }
        }

        if (totalMagnitude > 0.0) {
            negativeFraction = negativeMagnitude / totalMagnitude;
        }
        return std::isfinite(negativeFraction);
    }

    // ==========================================
    // Function: Convert unwindowed covariance into signed dynamic-grid synthesis power
    // Method: Select a no-wrap even 2/3/5-smooth grid, embed each retained lag uniquely, transform
    //         with a synthesis-only cached plan, normalize by G squared, and fftshift the result.
    // ==========================================
    bool covarianceToSynthesisPower(
        int stampSize, int maxLag,
        const std::vector<double>& covariance,
        int& synthesisSize,
        std::vector<float>& synthesisPower,
        double& maxImaginary) {
        synthesisSize = 0;
        synthesisPower.clear();
        maxImaginary = 0.0;
        if (stampSize <= 0 || maxLag < 0
            || maxLag > (std::numeric_limits<int>::max() - 1) / 2) {
            return false;
        }
        const int lagSide = 2 * maxLag + 1;
        const std::size_t expectedCovariance = static_cast<std::size_t>(lagSide)
                                               * static_cast<std::size_t>(lagSide);
        if (covariance.size() != expectedCovariance
            || !std::all_of(covariance.begin(), covariance.end(),
                            [](double value) { return std::isfinite(value); })) {
            return false;
        }

        const long long required = std::max(
            2LL * maxLag + 1LL,
            static_cast<long long>(stampSize) + maxLag);
        if (required <= 0 || required > std::numeric_limits<int>::max()) {
            return false;
        }
        synthesisSize = nextFastEvenFFTSize(static_cast<int>(required));
        if (synthesisSize <= 0) {
            return false;
        }

        const std::size_t synthesisElements = static_cast<std::size_t>(synthesisSize)
                                            * static_cast<std::size_t>(synthesisSize);
        std::vector<std::complex<double>> covarianceGrid(
            synthesisElements, std::complex<double>(0.0, 0.0));
        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            const int y = positiveModulo(dy, synthesisSize);
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                const int x = positiveModulo(dx, synthesisSize);
                covarianceGrid[static_cast<std::size_t>(y)
                               * static_cast<std::size_t>(synthesisSize)
                               + static_cast<std::size_t>(x)]
                    = covariance[compactLagIndex(dx, dy, maxLag)];
            }
        }

        PowerWorkspace& workspace = synthesisPowerWorkspace(synthesisSize);
        if (!workspace.valid()) {
            synthesisSize = 0;
            return false;
        }
        workspace.execute(covarianceGrid);

        synthesisPower.assign(synthesisElements, 0.0f);
        const double normalization = 1.0 / static_cast<double>(synthesisElements);
        const int half = synthesisSize / 2;
        for (int y = 0; y < synthesisSize; ++y) {
            const int shiftedY = (y + half) % synthesisSize;
            for (int x = 0; x < synthesisSize; ++x) {
                const int shiftedX = (x + half) % synthesisSize;
                const std::size_t inputIndex = static_cast<std::size_t>(y)
                                             * static_cast<std::size_t>(synthesisSize)
                                             + static_cast<std::size_t>(x);
                const std::size_t outputIndex = static_cast<std::size_t>(shiftedY)
                                              * static_cast<std::size_t>(synthesisSize)
                                              + static_cast<std::size_t>(shiftedX);
                const std::complex<double> value = workspace.value(inputIndex) * normalization;
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                    synthesisSize = 0;
                    synthesisPower.clear();
                    return false;
                }
                maxImaginary = std::max(maxImaginary, std::abs(value.imag()));
                synthesisPower[outputIndex] = static_cast<float>(value.real());
            }
        }
        return true;
    }

    // ==========================================
    // Function: Load one stored signed noise-power stamp without re-transforming it
    // Method: Validate the requested contiguous ns-by-ns slice and copy its values verbatim,
    //         preserving negative Fourier modes for linear subtraction.
    // ==========================================
    bool copyStoredNoisePower(
        const std::vector<float>& storedNoisePower,
        std::size_t offset,
        int stampSize,
        std::vector<float>& noisePower) {
        if (stampSize <= 0) {
            return false;
        }
        const std::size_t stampElements = static_cast<std::size_t>(stampSize)
                                        * static_cast<std::size_t>(stampSize);
        if (offset > storedNoisePower.size()
            || stampElements > storedNoisePower.size() - offset) {
            return false;
        }
        noisePower.assign(storedNoisePower.begin() + static_cast<std::ptrdiff_t>(offset),
                          storedNoisePower.begin()
                              + static_cast<std::ptrdiff_t>(offset + stampElements));
        return true;
    }

}
