#ifndef FQ_INIT_INITIALIZER_HPP
#define FQ_INIT_INITIALIZER_HPP

#include "InitConfig.hpp"
#include "process_init/FitsExtractor.hpp"

#include <mpi.h>

#include <filesystem>
#include <string>
#include <vector>

namespace fqinit {

struct Config {
    std::filesystem::path science_root;
    std::filesystem::path dq_root;
    std::filesystem::path output_root;
    std::string target;
    std::string filename_prefix;
    std::vector<std::string> filename_tokens = {"v1"};
    ExistingPolicy existing_policy = ExistingPolicy::Fail;
    int f77_max_path = InitConfig::F77_MAX_PATH;
    int max_chip = 62;
};

// ==========================================
// Function: Normalize and validate one initializer configuration
// Method: Resolve all paths and enforce the same required-field and target-name
//         contract for integrated workflow execution.
// ==========================================
void normalizeAndValidateConfig(Config& config);

// ==========================================
// Function: Build the pipeline input layout from the original FITS/FZ archives
// Method: Discover on rank zero, broadcast in-memory paths, extract archives
//         directly in parallel, and publish corrected deterministic lists.
// ==========================================
int runInitializer(const Config& config, MPI_Comm communicator);

}  // namespace fqinit

#endif
