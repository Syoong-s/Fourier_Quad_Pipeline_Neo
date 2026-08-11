#ifndef UNIVERSAL_BLOCK_HPP
#define UNIVERSAL_BLOCK_HPP

#include <string>

namespace Universalblock {

// ==========================================
// Enum: Classify the authoritative Stage-1 norm sentinel
// Method: Keep ordinary invalid-chip filtering distinct from missing or unreadable FITS input.
// ==========================================
enum class NormStatus {
    Valid,
    Invalid,
    Missing,
    ReadError
};

// ==========================================
// Function: Build the normalized-image path for one chip
// Method: Apply the shared chip prefix and sharded process-main output layout.
// ==========================================
std::string normFilename(const std::string& imageFile,
                         const std::string& dirOutput);

// ==========================================
// Function: Classify one chip from its normalized-image sentinel
// Method: Read only the first FITS pixel and distinguish invalid data from real I/O failures.
// ==========================================
NormStatus checkNorm(const std::string& imageFile,
                     const std::string& dirOutput,
                     float* sentinel = nullptr);

// ==========================================
// Function: Report a real normalized-image input failure
// Method: Emit one shared diagnostic for missing/read-error states and remain silent otherwise.
// ==========================================
void reportNormError(NormStatus status,
                     const std::string& imageFile,
                     const std::string& dirOutput);

}  // namespace Universalblock

#endif  // UNIVERSAL_BLOCK_HPP
