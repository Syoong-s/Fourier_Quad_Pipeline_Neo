#ifndef PCA_IMAGE_LAYOUT_HPP
#define PCA_IMAGE_LAYOUT_HPP

#include <algorithm>
#include <cstddef>

namespace PCAImageLayout {

// ==========================================
// Function: Copy one residual stamp into a PCA sample row
// Method: Preserve the image's contiguous row-major pixel order exactly.
// ==========================================
inline void copySampleRowMajor(const float* source, int pixelCount,
                               double* destination) {
    std::transform(source, source + pixelCount, destination,
                   [](float value) { return static_cast<double>(value); });
}

// ==========================================
// Function: Center one residual stamp for PCA projection
// Method: Subtract the matching row-major PCA mean from each image pixel.
// ==========================================
inline void centerSampleRowMajor(const float* source, const double* mean,
                                 int pixelCount, double* centered) {
    for (int pixel = 0; pixel < pixelCount; ++pixel) {
        centered[pixel] = static_cast<double>(source[pixel]) - mean[pixel];
    }
}

// ==========================================
// Function: Reconstruct one row-major residual image from PCA coefficients
// Method: Store components as [pixel][mode] and apply every retained mode
//         to the mean at the same image pixel.
// ==========================================
inline void reconstructRowMajor(const double* mean, const double* components,
                                int pixelCount, int modeCount,
                                const float* coefficients, float* output) {
    for (int pixel = 0; pixel < pixelCount; ++pixel) {
        double value = mean[pixel];
        for (int mode = 0; mode < modeCount; ++mode) {
            const std::size_t componentIndex =
                static_cast<std::size_t>(pixel) * modeCount + mode;
            value += components[componentIndex]
                * coefficients[mode];
        }
        output[pixel] = static_cast<float>(value);
    }
}

}  // namespace PCAImageLayout

#endif  // PCA_IMAGE_LAYOUT_HPP
