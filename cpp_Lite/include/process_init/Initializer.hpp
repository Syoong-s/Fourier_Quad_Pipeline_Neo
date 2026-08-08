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
// Function: Parse the standalone initializer command line
// Method: Require explicit archive/output roots and matching keys while
//         retaining safe defaults for existing outputs and F77 path limits.
// ==========================================
Config parseArguments(int argc, char** argv);

// ==========================================
// Function: Normalize and validate one initializer configuration
// Method: Resolve all paths and enforce the same required-field and target-name
//         contract for standalone parsing and integrated workflow execution.
// ==========================================
void normalizeAndValidateConfig(Config& config);

// ==========================================
// Function: Print the portable initializer command-line contract
// Method: Describe required paths, filters, existing-output policy, and limits.
// ==========================================
void printUsage(const char* program_name);

// ==========================================
// Function: Build the pipeline input layout from the original FITS/FZ archives
// Method: Discover on rank zero, broadcast in-memory paths, extract archives
//         directly in parallel, and publish corrected deterministic lists.
// ==========================================
int runInitializer(const Config& config, MPI_Comm communicator);

}  // namespace fqinit

#endif
