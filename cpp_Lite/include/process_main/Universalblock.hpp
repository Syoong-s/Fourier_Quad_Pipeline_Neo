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
// Function: Describe a normalized-image input failure
// Method: Return a stable missing/read-error diagnostic for MPI-wide fatal reporting.
// ==========================================
std::string normErrorDetail(NormStatus status,
                            const std::string& imageFile,
                            const std::string& dirOutput);

}  // namespace Universalblock

#endif  // UNIVERSAL_BLOCK_HPP
