#include "Universalblock.hpp"

#include "FitsIO.hpp"
#include "OutputLayout.hpp"
#include "UniversalUtils.hpp"

#include <cmath>
#include <iostream>

namespace Universalblock {

// ==========================================
// Function: Build the normalized-image path for one chip
// Method: Apply the shared chip prefix and sharded process-main output layout.
// ==========================================
std::string normFilename(const std::string& imageFile,
                         const std::string& dirOutput) {
    const std::string prefix = UniversalUtils::getPrefix(imageFile);
    return OutputLayout::chipPath(
        dirOutput, "stamps/Norm", prefix, "_norm.fits");
}

// ==========================================
// Function: Classify one chip from its normalized-image sentinel
// Method: Read only the first FITS pixel, preserving missing/read-error states while treating
//         the established non-finite, non-negative, and failure sentinels as invalid data.
// ==========================================
NormStatus checkNorm(const std::string& imageFile,
                     const std::string& dirOutput,
                     float* sentinel) {
    float norm0 = 0.0f;
    const FitsIO::PixelReadStatus readStatus = FitsIO::readFirstPixel(
        normFilename(imageFile, dirOutput), norm0);

    if (readStatus == FitsIO::PixelReadStatus::Missing) {
        return NormStatus::Missing;
    }
    if (readStatus != FitsIO::PixelReadStatus::Ok) {
        return NormStatus::ReadError;
    }

    if (sentinel != nullptr) {
        *sentinel = norm0;
    }

    if (!std::isfinite(norm0) || norm0 >= 0.0f || norm0 < -99990.0f) {
        return NormStatus::Invalid;
    }
    return NormStatus::Valid;
}

// ==========================================
// Function: Report a real normalized-image input failure
// Method: Emit one shared diagnostic for missing/read-error states and remain silent otherwise.
// ==========================================
void reportNormError(NormStatus status,
                     const std::string& imageFile,
                     const std::string& dirOutput) {
    if (status == NormStatus::Missing) {
        std::cerr << "Error / norm FITS file missing: "
                  << normFilename(imageFile, dirOutput) << std::endl;
    } else if (status == NormStatus::ReadError) {
        std::cerr << "Error / norm FITS read failure: "
                  << normFilename(imageFile, dirOutput) << std::endl;
    }
}

}  // namespace Universalblock
