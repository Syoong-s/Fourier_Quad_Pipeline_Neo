#include "process_main/process_main.hpp"

#include "Astrometry.hpp"
#include "CatalogCombiner.hpp"
#include "ExposureInfo.hpp"
#include "ExternalCatalogReader.hpp"
#include "FourierTransformSt1.hpp"
#include "FourierTransformSt2.hpp"
#include "LensingConfig.hpp"
#include "MPIScheduler.hpp"
#include "OutputFile.hpp"
#include "PSFModel.hpp"
#include "PreProcess.hpp"
#include "ShearMeasurement.hpp"
#include "SourceExtractor.hpp"
#include "UniversalUtils.hpp"

#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

std::vector<std::string> EXPO_FILE;
int N_EXPO = 0;

namespace {

// ==========================================
// Function: Load and validate the top-level exposure list on rank zero
// Method: Parse quoted per-exposure list paths and reject empty or malformed
//         input before any rank enters the numerical stage scheduler.
// ==========================================
bool loadExposureList(const std::string& exposure_list, std::string& error) {
    EXPO_FILE.clear();
    N_EXPO = 0;

    std::ifstream input(exposure_list);
    if (!input.is_open()) {
        error = "EXPO_LIST reading error: " + exposure_list;
        return false;
    }

    std::string exposure_name;
    int chip_count = 0;
    while (input >> exposure_name >> chip_count) {
        if (exposure_name.size() >= 2 && exposure_name.front() == '"'
            && exposure_name.back() == '"') {
            exposure_name = exposure_name.substr(1, exposure_name.size() - 2);
        }
        EXPO_FILE.push_back(exposure_name);
    }
    if (!input.eof()) {
        error = "EXPO_LIST contains an invalid record: " + exposure_list;
        EXPO_FILE.clear();
        return false;
    }
    if (EXPO_FILE.empty()) {
        error = "EXPO_LIST contains no exposures: " + exposure_list;
        return false;
    }
    N_EXPO = static_cast<int>(EXPO_FILE.size());
    std::cout << "Total number of EXPOSURE: " << N_EXPO << std::endl;
    return true;
}

// ==========================================
// Function: Broadcast the validated exposure list to every MPI rank
// Method: Send the count, then length-prefix each mutable C++ string.
// ==========================================
void broadcastExposureList(int rank) {
    MPI_Bcast(&N_EXPO, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        EXPO_FILE.resize(static_cast<std::size_t>(N_EXPO));
    }

    for (int index = 0; index < N_EXPO; ++index) {
        int length = rank == 0 ? static_cast<int>(EXPO_FILE[index].size()) : 0;
        MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            EXPO_FILE[index].resize(static_cast<std::size_t>(length));
        }
        if (length > 0) {
            MPI_Bcast(EXPO_FILE[index].data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
    }
}

}  // namespace

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline with compiled defaults
// Method: Preserve the historical one-argument API by installing and forwarding
//         a default RuntimeConfig with its default external layout.
// ==========================================
int process_main(const std::string& exposure_list) {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    std::string store_error;
    if (!RuntimeConfigStore::isInitialized()
        && !RuntimeConfigStore::initialize(config, store_error)) {
        if (MPIScheduler::my_id == 0) {
            std::cerr << "Runtime config error: " << store_error << std::endl;
        }
        return 1;
    }
    return process_main(exposure_list, RuntimeConfigStore::get());
}

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline with unified runtime options
// Method: Resolve one compatibility layout, then forward to the startup-layout implementation.
// ==========================================
int process_main(const std::string& exposure_list,
                 const RuntimeConfig& runtime_config) {
    PipelineCatalog::CatalogLayout layout;
    std::string layout_error;
    const int local_layout_ok =
        PipelineCatalog::resolveCatalogLayout(runtime_config, layout, layout_error)
            ? 1
            : 0;
    int global_layout_ok = 0;
    MPI_Allreduce(&local_layout_ok, &global_layout_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (global_layout_ok == 0) {
        if (MPIScheduler::my_id == 0) {
            std::cerr << "Catalog layout error: "
                      << (layout_error.empty()
                              ? "resolution failed on another MPI rank"
                              : layout_error)
                      << std::endl;
        }
        return 1;
    }
    return process_main(exposure_list, runtime_config, layout);
}

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline with one shared catalog layout
// Method: Configure readers from the startup schema, then execute MPI stages without re-resolution.
// ==========================================
int process_main(const std::string& exposure_list,
                 const RuntimeConfig& runtime_config,
                 const PipelineCatalog::CatalogLayout& layout) {
    const int rank = MPIScheduler::my_id;
    if (!RuntimeConfigStore::isInitialized()) {
        std::string store_error;
        if (!RuntimeConfigStore::initialize(runtime_config, store_error)) {
            if (rank == 0) {
                std::cerr << "Runtime config error: " << store_error
                          << std::endl;
            }
            return 1;
        }
    }

    std::string column_error;
    const int local_columns_ok =
        ExternalCatalogReader::configure(layout, column_error) ? 1 : 0;
    int global_columns_ok = 0;
    MPI_Allreduce(&local_columns_ok, &global_columns_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (global_columns_ok == 0) {
        if (rank == 0) {
            std::cerr << "External-catalog column error: "
                      << (column_error.empty()
                              ? "validation failed on another MPI rank"
                              : column_error)
                      << std::endl;
        }
        return 1;
    }

    int load_ok = 1;
    std::string load_error;
    if (rank == 0 && !loadExposureList(exposure_list, load_error)) {
        load_ok = 0;
    }
    MPI_Bcast(&load_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (load_ok == 0) {
        if (rank == 0) {
            std::cerr << load_error << std::endl;
        }
        return 1;
    }

    broadcastExposureList(rank);
    MPIScheduler::barrier();

    // ==========================================
    // Function: Validate stage dependency before execution
    // Method: Stage 9 consumes Stage 8 exposure chi2, so reject PROCESS_stage with 23 but without 19.
    // ==========================================
    const LensingRuntimeConfig& lensing = runtime_config.lensing;
    if (lensing.process_stage % 23 == 0 && lensing.process_stage % 19 != 0) {
        if (rank == 0) {
            std::cerr << "Error: Stage 9 requires Stage 8. PROCESS_stage enables "
                         "CatalogCombiner without ExposureInfo."
                      << std::endl;
        }
        return 1;
    }

    if (lensing.process_stage % 2 == 0) {
        MPIScheduler::distribute(N_EXPO, PreProcess::preProcess, "Pre-process...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 3 == 0) {
        MPIScheduler::distribute(N_EXPO, Astrometry::procAstrometry, "Astrometry...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 5 == 0) {
        MPIScheduler::distribute(N_EXPO, SourceExtractor::procSource, "Sources ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 7 == 0) {
        MPIScheduler::distribute(N_EXPO, FourierTransformSt1::procFourierTSt1, "FFT st1...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 11 == 0) {
        MPIScheduler::distribute(N_EXPO, PSFModel::procPSF, "PSF ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 13 == 0) {
        MPIScheduler::distribute(N_EXPO, FourierTransformSt2::procFourierTSt2, "FFT st2...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 17 == 0) {
        MPIScheduler::distribute(N_EXPO, ShearMeasurement::procShear, "Shear ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 19 == 0) {
        const int mpi_parameter_count = N_EXPO * 6;
        const std::size_t parameter_count =
            static_cast<std::size_t>(mpi_parameter_count);

        ExposureInfo::expo_para.assign(parameter_count, 0.0f);
        MPIScheduler::distribute(N_EXPO, ExposureInfo::procInfo, "Info ...");
        MPIScheduler::barrier();

        std::vector<float> reduced_exposure_parameters(parameter_count, 0.0f);
        MPI_Allreduce(ExposureInfo::expo_para.data(),
                      reduced_exposure_parameters.data(),
                      mpi_parameter_count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        ExposureInfo::expo_para = std::move(reduced_exposure_parameters);

        if (rank == 0) {
            const std::string root_directory = UniversalUtils::getDir(exposure_list, 1);
            const std::string filename = root_directory + "/expo_info.dat";
            MainIO::OutputFile output(filename);
            output << std::setprecision(10);
            output << "N-valid-chip PSF-FWHM(arcsec) chi_d-stars nstar-per-chip "
                      "cRVAL1 cRVAL2 expo_name\n";
            for (int exposure = 0; exposure < N_EXPO; ++exposure) {
                for (int parameter = 0; parameter < 6; ++parameter) {
                    output << ExposureInfo::expo_para[exposure * 6 + parameter] << " ";
                }
                output << EXPO_FILE[exposure] << "\n";
            }
            output.close();
        }
        MPIScheduler::barrier();
    }

    if (lensing.process_stage % 23 == 0) {
        MPIScheduler::distribute(N_EXPO, CatalogCombiner::procComb, "combine ...");
    }
    MPIScheduler::barrier();
    return 0;
}
