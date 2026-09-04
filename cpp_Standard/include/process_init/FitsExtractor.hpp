#ifndef FQ_INIT_FITS_EXTRACTOR_HPP
#define FQ_INIT_FITS_EXTRACTOR_HPP

#include <filesystem>
#include <string>
#include <vector>

namespace fqinit {

enum class ProductKind {
    Science,
    DqMask,
};

enum class ExistingPolicy {
    Fail,
    Resume,
    Overwrite,
};

struct ExtractionResult {
    bool success = false;
    bool resumed = false;
    std::vector<std::filesystem::path> output_paths;
    std::string error;
    int skipped_hdus = 0;
};

// ==========================================
// Function: Extract one multi-HDU FITS/FZ archive into pipeline chip images
// Method: Read the source archive in place with CFITSIO, stage uncompressed
//         two-dimensional HDUs, then commit the completed output set.
// ==========================================
ExtractionResult extractArchive(const std::filesystem::path& source,
                                ProductKind kind,
                                const std::filesystem::path& final_directory,
                                const std::filesystem::path& staging_directory,
                                ExistingPolicy policy);

// ==========================================
// Function: Return the exposure stem used by science chip output names
// Method: Remove the exact configured archive suffix from the source basename.
// ==========================================
std::string archiveStem(const std::filesystem::path& source);

// ==========================================
// Function: Return the DQ exposure stem expected by the pipeline
// Method: Remove the archive suffix and apply the configured DQ stem replacement.
// ==========================================
std::string dqOutputStem(const std::filesystem::path& source);

}  // namespace fqinit

#endif
