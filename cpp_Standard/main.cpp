#include "CatalogLayout.hpp"
#include "ExtCatConfig.hpp"
#include "InitConfig.hpp"
#include "ProcessConfig.hpp"
#include "RuntimeConfig.hpp"
#include "process_extcat/process_extcat.hpp"
#include "process_fd/process_fd.hpp"
#include "process_init/process_init.hpp"
#include "process_main/MPIScheduler.hpp"
#include "process_main/NumericalRecipes.hpp"
#include "process_main/process_main.hpp"
#include "process_rearr/process_rearr.hpp"

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

// ==========================================
// Function: broadcastString
// Method: Broadcast one length-prefixed startup string from rank zero while
//         rejecting text that cannot fit an MPI int count.
// ==========================================
bool broadcastString(std::string& value, int rank, std::string& error) {
    int length = 0;
    if (rank == 0) {
        if (value.size() > static_cast<std::size_t>(
                               std::numeric_limits<int>::max())) {
            error = "runtime config text exceeds the MPI broadcast count range";
            length = -1;
        } else {
            length = static_cast<int>(value.size());
        }
    }
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (length < 0) {
        return false;
    }
    if (rank != 0) {
        value.resize(static_cast<std::size_t>(length));
    }
    if (length > 0) {
        MPI_Bcast(value.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    return true;
}

// ==========================================
// Function: loadAndBroadcastConfigText
// Method: Read the selected INI only on rank zero, broadcast one success flag,
//         then distribute the exact text to every rank for identical parsing.
// ==========================================
bool loadAndBroadcastConfigText(const std::string& path, int rank,
                                std::string& text, std::string& error) {
    int read_ok = 1;
    if (rank == 0) {
        std::ifstream input(path, std::ios::binary);
        if (!input.is_open()) {
            error = "cannot open runtime config file: " + path;
            read_ok = 0;
        } else {
            std::ostringstream buffer;
            buffer << input.rdbuf();
            if (!input.good() && !input.eof()) {
                error = "I/O failure while reading runtime config file: " + path;
                read_ok = 0;
            } else {
                text = buffer.str();
            }
        }
    }
    MPI_Bcast(&read_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (read_ok == 0) {
        return false;
    }
    return broadcastString(text, rank, error);
}

// ==========================================
// Function: configuredDatasetsText
// Method: Join compiled fallback dataset pairs for help output.
// ==========================================
std::string configuredDatasetsText() {
    if (InitConfig::DATASETS.empty()) return "none";
    std::string text;
    for (const InitConfig::DatasetSpec& dataset : InitConfig::DATASETS) {
        if (!text.empty()) text += ", ";
        text += dataset.target + ":" + dataset.prefix;
    }
    return text;
}

// ==========================================
// Function: configuredContainsText
// Method: Join compiled initializer tokens or describe an empty filter.
// ==========================================
std::string configuredContainsText() {
    if (InitConfig::CONTAINS.empty()) return "none (no token filter)";
    std::string text;
    for (const std::string& token : InitConfig::CONTAINS) {
        if (!text.empty()) text += ", ";
        text += token;
    }
    return text;
}

// ==========================================
// Function: configuredExtcatContainsText
// Method: Join compiled raw-catalog tokens or describe an all-file match.
// ==========================================
std::string configuredExtcatContainsText() {
    if (ExtCatConfig::EXTCAT_FILENAME_TOKENS.empty()) return "none (all files)";
    std::string text;
    for (const std::string& token : ExtCatConfig::EXTCAT_FILENAME_TOKENS) {
        if (!text.empty()) text += ", ";
        text += token;
    }
    return text;
}

// ==========================================
// Function: printUsage
// Method: Document config-file precedence, retained CLI overrides, and the
//         config-only Standard lensing controls without duplicating the INI example.
// ==========================================
void printUsage(const char* program_name) {
    std::cout
        << "Usage: " << program_name << " [options] [LEGACY_EXPO_LIST]\n"
        << "  --config PATH         Load an INI runtime config before CLI overrides\n"
        << "  --run-extcat BOOL     Repartition raw external catalogs first (default: "
        << (ProcessConfig::RUN_PROCESS_EXTCAT ? "true" : "false") << ")\n"
        << "  --run-init BOOL       Run initializer (default: "
        << (ProcessConfig::RUN_PROCESS_INIT ? "true" : "false") << ")\n"
        << "  --run-main BOOL       Run numerical pipeline (default: "
        << (ProcessConfig::RUN_PROCESS_MAIN ? "true" : "false") << ")\n"
        << "  --run-rearr BOOL      Rearrange _all.cat (default: "
        << (ProcessConfig::RUN_PROCESS_REARR ? "true" : "false") << ")\n"
        << "  --run-fd BOOL         Run field-distortion test (default: "
        << (ProcessConfig::RUN_PROCESS_FD ? "true" : "false") << ")\n"
        << "  --extcat-input PATH   Directory containing raw external catalogs\n"
        << "  --extcat-output PATH  External tile directory and effective SOURCE_CAT\n"
        << "  --extcat-contains T   Repeatable raw basename token (default: "
        << configuredExtcatContainsText() << ")\n"
        << "  --extcat-recursive B  Recurse below extcat input\n"
        << "  --extcat-delimiter M  auto, whitespace, comma, or tab\n"
        << "  --extcat-header M     auto, present, or absent\n"
        << "  --extcat-columns LIST Ordered one-based projection indices\n"
        << "  --extcat-ra-column N  Raw one-based RA column\n"
        << "  --extcat-dec-column N Raw one-based Dec column\n"
        << "  --extcat-zp-column N  Raw one-based ZP column\n"
        << "  --extcat-chunk-mib N  MPI task size in MiB (default: "
        << ExtCatConfig::EXTCAT_CHUNK_MIB << ")\n"
        << "  --extcat-malformed P  fail or skip malformed rows\n"
        << "  --extcat-existing P   fail or overwrite generated tiles\n"
        << "  --science-root PATH   Original Science FITS/FZ repository\n"
        << "  --dq-root PATH        Original DQ FITS/FZ repository\n"
        << "  --output-root PATH    Parent of targets and generated lists\n"
        << "  --dataset T:P         Repeatable dataset pair (default: "
        << configuredDatasetsText() << ")\n"
        << "  --target NAME         Legacy single-dataset target\n"
        << "  --prefix TEXT         Legacy single-dataset prefix\n"
        << "  --contains TEXT       Repeatable archive token (default: "
        << configuredContainsText() << ")\n"
        << "  --existing MODE       fail, resume, or overwrite\n"
        << "  --f77-max-path N      Generated path limit; zero disables (default: "
        << InitConfig::F77_MAX_PATH << ")\n"
        << "  --expo-list PATH      Exposure list for downstream-only mode\n"
        << "  --rearr-output-dir D  Rearrangement output directory\n"
        << "  --rearr-output-base P Rearrangement output base\n"
        << "  --rearr-list-name F   Rearranged list filename (default: "
        << ProcessConfig::REARRANGED_EXPO_LIST_FILENAME << ")\n"
        << "  --rearr-list-dir P    Rearranged list directory\n"
        << "  --fd-expo-list PATH   FD exposure-list override\n"
        << "  --fd-output-dir D     FD output directory\n"
        << "  --fd-output-base P    FD output base\n"
        << "  --help                Show this help\n"
        << "Precedence: compiled defaults < --config INI < CLI. Options accept both "
           "--name value and --name=value; duplicate scalars use the last value.\n"
        << "The first CLI --dataset, --contains, or --extcat-contains replaces its "
           "configured list; repeats append. Standard lensing controls, geometry, "
           "and scientific paths are set in the INI [lensing] section.\n";
}

// ==========================================
// Function: resolveExposureList
// Method: Prefer the configured list or derive output_root/expo_TARGET.list for
//         the current downstream-only dataset.
// ==========================================
std::string resolveExposureList(const RuntimeConfig& config,
                                const InitConfig::DatasetSpec& dataset) {
    std::filesystem::path path;
    if (!config.process.expo_list.empty()) {
        path = config.process.expo_list;
    } else {
        if (config.init.output_root.empty()) {
            throw std::runtime_error(
                "expo_list is absent and output_root cannot derive its default");
        }
        path = std::filesystem::path(config.init.output_root)
               / ("expo_" + dataset.target + ".list");
    }
    return std::filesystem::weakly_canonical(
               std::filesystem::absolute(path)).string();
}

// ==========================================
// Function: deriveDatasetRootFromExpoList
// Method: Read the first usable image path and return its great-grandparent,
//         matching process_main's dataset-root convention.
// ==========================================
std::string deriveDatasetRootFromExpoList(const std::string& exposure_list) {
    std::ifstream expo_input(exposure_list);
    if (!expo_input.is_open()) return "";
    std::string path;
    int chip_count = 0;
    while (expo_input >> path >> chip_count) {
        if (path.size() >= 2 && path.front() == '"' && path.back() == '"') {
            path = path.substr(1, path.size() - 2);
        }
        std::ifstream list_input(path);
        if (!list_input.is_open()) continue;
        std::string line;
        while (std::getline(list_input, line)) {
            const std::size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) continue;
            const std::size_t last = line.find_last_not_of(" \t\r\n");
            line = line.substr(first, last - first + 1);
            if (line.size() >= 2 && line.front() == '"' && line.back() == '"') {
                line = line.substr(1, line.size() - 2);
            }
            const std::filesystem::path image_path(line);
            const std::filesystem::path root = image_path.parent_path()
                                                   .parent_path()
                                                   .parent_path();
            if (root.empty()) return "";
            return std::filesystem::absolute(root).lexically_normal().string();
        }
    }
    return "";
}

}  // namespace

// ==========================================
// Function: main
// Method: Read one rank-zero INI, apply identical text and CLI on all ranks,
//         freeze RuntimeConfig, resolve only required schemas, then dispatch the
//         five pipeline phases in their fixed order.
// ==========================================
int main(int argc, char* argv[]) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::my_id;
    if (rank == 0) std::cout << "MPI Init Done..." << std::endl;

    int return_code = 0;
    std::string config_path;
    bool phase_a_help = false;
    std::string startup_error;
    const int local_scan_ok = findRuntimeConfigPath(
        argc, argv, config_path, phase_a_help, startup_error) ? 1 : 0;
    int global_scan_ok = 0;
    MPI_Allreduce(&local_scan_ok, &global_scan_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    RuntimeConfig config = makeDefaultRuntimeConfig();
    std::string config_text;
    int global_config_ok = global_scan_ok;
    if (global_scan_ok != 0 && !config_path.empty() && !phase_a_help) {
        const int local_load_ok = loadAndBroadcastConfigText(
            config_path, rank, config_text, startup_error) ? 1 : 0;
        MPI_Allreduce(&local_load_ok, &global_config_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if (global_config_ok != 0) {
            const int local_apply_ok = applyRuntimeConfigText(
                config_text, config_path, config, startup_error) ? 1 : 0;
            MPI_Allreduce(&local_apply_ok, &global_config_ok, 1, MPI_INT,
                          MPI_MIN, MPI_COMM_WORLD);
        }
    }

    int global_parse_ok = global_config_ok;
    if (global_config_ok != 0) {
        const int local_parse_ok = parseRuntimeCommandLine(
            argc, argv, config, startup_error) ? 1 : 0;
        MPI_Allreduce(&local_parse_ok, &global_parse_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
    }

    int global_validation_ok = global_parse_ok;
    if (global_parse_ok != 0 && !config.help_requested) {
        const int local_validation_ok = validateRuntimeConfig(
            config, startup_error) ? 1 : 0;
        MPI_Allreduce(&local_validation_ok, &global_validation_ok, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);
    }

    if (global_scan_ok == 0 || global_config_ok == 0
        || global_parse_ok == 0 || global_validation_ok == 0) {
        if (rank == 0) {
            std::cerr << "Config error:\n  "
                      << (startup_error.empty()
                              ? "startup failed on another MPI rank"
                              : startup_error)
                      << std::endl;
            printUsage(argv[0]);
        }
        return_code = 2;
    } else if (config.help_requested) {
        if (rank == 0) printUsage(argv[0]);
    } else {
        std::string store_error;
        const int local_store_ok = RuntimeConfigStore::initialize(
            config, store_error) ? 1 : 0;
        int global_store_ok = 0;
        MPI_Allreduce(&local_store_ok, &global_store_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if (global_store_ok == 0) {
            if (rank == 0) {
                std::cerr << "Runtime config error: "
                          << (store_error.empty()
                                  ? "store initialization failed on another rank"
                                  : store_error)
                          << std::endl;
            }
            return_code = 2;
        } else {
            const RuntimeConfig& runtime_config = RuntimeConfigStore::get();
            std::optional<PipelineCatalog::CatalogLayout> external_layout;
            std::string schema_error;
            int global_schema_ok = 1;
            if (runtime_config.lensing.ext_cat == 1) {
                PipelineCatalog::CatalogLayout layout;
                const int local_schema_ok = PipelineCatalog::resolveCatalogLayout(
                    runtime_config, layout, schema_error) ? 1 : 0;
                MPI_Allreduce(&local_schema_ok, &global_schema_ok, 1, MPI_INT,
                              MPI_MIN, MPI_COMM_WORLD);
                if (global_schema_ok != 0) external_layout = std::move(layout);
            }

            PipelineCatalog::RearrCatalogSchema rearr_schema;
            if (global_schema_ok != 0
                && runtime_config.process.run_process_rearr) {
                const int local_rearr_schema_ok =
                    PipelineCatalog::resolveRearrCatalogSchema(
                        runtime_config,
                        external_layout ? &*external_layout : nullptr,
                        rearr_schema, schema_error) ? 1 : 0;
                MPI_Allreduce(&local_rearr_schema_ok, &global_schema_ok, 1,
                              MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            }

            if (global_schema_ok == 0) {
                if (rank == 0) {
                    std::cerr << "Catalog schema error: "
                              << (schema_error.empty()
                                      ? "resolution failed on another rank"
                                      : schema_error)
                              << std::endl;
                }
                return_code = 2;
            } else {
                if (rank == 0 && external_layout) {
                    std::cout << PipelineCatalog::describeCatalogLayout(
                                     *external_layout)
                              << std::endl;
                }
                if (runtime_config.process.run_process_extcat) {
                    if (rank == 0) {
                        std::cout << "Running process_extcat before all dataset phases"
                                  << std::endl;
                    }
                    return_code = process_extcat(runtime_config,
                                                 MPI_COMM_WORLD);
                    if (return_code == 0) MPIScheduler::barrier();
                }

                bool rng_initialized = false;
                for (std::size_t index = 0;
                     index < runtime_config.init.datasets.size()
                         && return_code == 0
                         && (runtime_config.process.run_process_init
                             || runtime_config.process.run_process_main
                             || runtime_config.process.run_process_rearr
                             || runtime_config.process.run_process_fd);
                     ++index) {
                    const InitConfig::DatasetSpec& dataset =
                        runtime_config.init.datasets[index];
                    if (rank == 0) {
                        std::cout << "Dataset " << (index + 1) << "/"
                                  << runtime_config.init.datasets.size()
                                  << ": target=" << dataset.target
                                  << " prefix=" << dataset.prefix << std::endl;
                    }

                    std::string generated_exposure_list;
                    if (runtime_config.process.run_process_init) {
                        return_code = process_init(runtime_config, dataset,
                                                   generated_exposure_list);
                    }

                    std::string selected_exposure_list;
                    if (return_code == 0
                        && (runtime_config.process.run_process_main
                            || runtime_config.process.run_process_rearr
                            || runtime_config.process.run_process_fd)) {
                        if (runtime_config.process.run_process_init) {
                            selected_exposure_list = generated_exposure_list;
                            if (rank == 0) {
                                std::cout << "Downstream phases will use process_init output: "
                                          << selected_exposure_list << std::endl;
                            }
                            MPIScheduler::barrier();
                        } else {
                            int local_path_ok = 1;
                            std::string path_error;
                            try {
                                selected_exposure_list = resolveExposureList(
                                    runtime_config, dataset);
                            } catch (const std::exception& exception) {
                                local_path_ok = 0;
                                path_error = exception.what();
                            }
                            int global_path_ok = 0;
                            MPI_Allreduce(&local_path_ok, &global_path_ok, 1,
                                          MPI_INT, MPI_MIN, MPI_COMM_WORLD);
                            if (global_path_ok == 0) {
                                if (rank == 0) {
                                    std::cerr << "Exposure-list error: "
                                              << (path_error.empty()
                                                      ? "resolution failed on another rank"
                                                      : path_error)
                                              << std::endl;
                                }
                                return_code = 2;
                            }
                        }
                    }

                    if (return_code == 0
                        && runtime_config.process.run_process_main) {
                        if (!rng_initialized) {
                            const unsigned int seed =
                                NumericalRecipes::initializeRan1Seed(
                                    rank, MPIScheduler::num_procs);
                            std::cout << "RNG_SEED rank seed: " << rank << " "
                                      << seed << std::endl;
                            MPIScheduler::barrier();
                            rng_initialized = true;
                        }
                        return_code = process_main(
                            selected_exposure_list, runtime_config,
                            external_layout ? &*external_layout : nullptr);
                    }

                    if (return_code == 0
                        && runtime_config.process.run_process_rearr) {
                        MPIScheduler::barrier();
                        if (rank == 0) std::cout << "Running process_rearr" << std::endl;
                        return_code = process_rearr(
                            selected_exposure_list, runtime_config,
                            rearr_schema, MPI_COMM_WORLD);
                    }

                    if (return_code == 0
                        && runtime_config.process.run_process_fd) {
                        MPIScheduler::barrier();
                        if (!rng_initialized) {
                            NumericalRecipes::initializeRan1Seed(
                                rank, MPIScheduler::num_procs);
                            MPIScheduler::barrier();
                            rng_initialized = true;
                        }
                        std::string fd_expo_list;
                        if (runtime_config.process.fd_expo_list.empty()) {
                            const std::filesystem::path list_dir =
                                runtime_config.process
                                        .rearranged_expo_list_directory.empty()
                                    ? std::filesystem::path(selected_exposure_list)
                                          .parent_path()
                                    : std::filesystem::path(
                                          runtime_config.process
                                              .rearranged_expo_list_directory);
                            fd_expo_list = std::filesystem::absolute(
                                list_dir / runtime_config.process
                                               .rearranged_expo_list_filename)
                                               .lexically_normal().string();
                        } else {
                            fd_expo_list = runtime_config.process.fd_expo_list;
                        }
                        const std::string dataset_root =
                            deriveDatasetRootFromExpoList(selected_exposure_list);
                        if (rank == 0) {
                            std::cout << "Running process_fd: expo_list="
                                      << fd_expo_list
                                      << " dataset_root=" << dataset_root
                                      << std::endl;
                        }
                        return_code = process_fd(fd_expo_list, runtime_config,
                                                 dataset_root,
                                                 *external_layout);
                    }
                    if (return_code == 0
                        && index + 1 < runtime_config.init.datasets.size()) {
                        MPIScheduler::barrier();
                    }
                }
            }
        }
    }

    MPIScheduler::finalize();
    return return_code;
}
