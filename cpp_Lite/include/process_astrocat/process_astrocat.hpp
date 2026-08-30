#ifndef PROCESS_ASTROCAT_PROCESS_ASTROCAT_HPP
#define PROCESS_ASTROCAT_PROCESS_ASTROCAT_HPP

#include "RuntimeConfig.hpp"

#include <filesystem>

namespace ProcessAstrocat {

enum class ExistingPolicy {
    Fail,
    Overwrite,
};

// ==========================================
// Configuration: Runtime contract for process_astrocat
// Method: Keep the raw input, generated-tile output, header behavior, and
//         replacement policy independent from process_main catalog paths.
// ==========================================
struct Config {
    std::filesystem::path input_directory;
    std::filesystem::path output_directory;
    bool add_header = true;
    ExistingPolicy existing_policy = ExistingPolicy::Fail;
};

// ==========================================
// Function: Normalize and validate process_astrocat paths
// Method: Resolve absolute paths and reject empty, identical, or nested input
//         and output roots before directory discovery begins.
// ==========================================
void normalizeAndValidateConfig(Config& config);

}  // namespace ProcessAstrocat

// ==========================================
// Function: Repartition raw two-column Gaia catalogs
// Method: Read whole files through dynamic MPI scheduling, redistribute rows
//         to deterministic tile owners, de-duplicate, and publish final tiles.
// ==========================================
int process_astrocat(ProcessAstrocat::Config config);

// ==========================================
// Function: Run process_astrocat from immutable runtime configuration
// Method: Translate the dedicated astrocat domain without changing the
//         process_main lensing astrometry catalog path.
// ==========================================
int process_astrocat(const RuntimeConfig& runtime_config);

#endif  // PROCESS_ASTROCAT_PROCESS_ASTROCAT_HPP
