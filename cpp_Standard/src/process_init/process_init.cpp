#include "process_init/process_init.hpp"

#include "process_init/Initializer.hpp"

#include <mpi.h>

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

// ==========================================
// Function: Run the integrated MPI pipeline initializer
// Method: Convert unified runtime options into the preserved initializer modules
//         and return the generated absolute exposure-list path on success.
// ==========================================
int process_init(const ProcessConfig::RuntimeOptions& options,
                 const InitConfig::DatasetSpec& dataset,
                 std::string& generated_expo_list) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    generated_expo_list.clear();

    fqinit::Config config;
    config.science_root = options.science_root;
    config.dq_root = options.dq_root;
    config.output_root = options.output_root;
    config.target = dataset.target;
    config.filename_prefix = dataset.prefix;
    config.filename_tokens = options.contains;
    config.f77_max_path = options.f77_max_path;
    config.max_chip = LensingConfig::NMAX_CHIP;

    if (options.existing == "fail") {
        config.existing_policy = fqinit::ExistingPolicy::Fail;
    } else if (options.existing == "resume") {
        config.existing_policy = fqinit::ExistingPolicy::Resume;
    } else if (options.existing == "overwrite") {
        config.existing_policy = fqinit::ExistingPolicy::Overwrite;
    } else {
        if (rank == 0) {
            std::cerr << "Initializer argument error: --existing must be fail, resume, or overwrite"
                      << std::endl;
        }
        return 2;
    }

    int local_config_ok = 1;
    std::string config_error;
    try {
        fqinit::normalizeAndValidateConfig(config);
    } catch (const std::exception& exception) {
        local_config_ok = 0;
        config_error = exception.what();
    }

    int global_config_ok = 0;
    MPI_Allreduce(&local_config_ok, &global_config_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (global_config_ok == 0) {
        if (rank == 0) {
            std::cerr << "Initializer argument error: "
                      << (config_error.empty() ? "configuration validation failed on another rank"
                                               : config_error)
                      << std::endl;
        }
        return 2;
    }

    const int return_code = fqinit::runInitializer(config, MPI_COMM_WORLD);
    if (return_code != 0) {
        return return_code;
    }

    const std::filesystem::path list_path = config.output_root
                                            / ("expo_" + config.target + ".list");
    int local_list_ok = std::filesystem::is_regular_file(list_path) ? 1 : 0;
    int global_list_ok = 0;
    MPI_Allreduce(&local_list_ok, &global_list_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (global_list_ok == 0) {
        if (rank == 0) {
            std::cerr << "Initializer completed without a readable exposure list: "
                      << list_path << std::endl;
        }
        return 1;
    }

    generated_expo_list = std::filesystem::weakly_canonical(list_path).string();
    return 0;
}
