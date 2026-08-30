#include "process_main/process_main.hpp"

#include "process_main/Astrometry.hpp"
#include "process_main/CatalogCombiner.hpp"
#include "process_main/ExposureInfo.hpp"
#include "process_main/ExternalCatalogReader.hpp"
#include "process_main/FourierTransformSt1.hpp"
#include "process_main/FourierTransformSt2.hpp"
#include "LensingConfig.hpp"
#include "general/MPIScheduler.hpp"
#include "process_main/OutputFile.hpp"
#include "process_main/PSFModel.hpp"
#include "process_main/PreProcess.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/ShearMeasurement.hpp"
#include "process_main/SourceExtractor.hpp"
#include "process_main/UniversalUtils.hpp"
#include "general/ExposureList.hpp"
#include "general/MPIUtils.hpp"

#include <mpi.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

ProcessMain::State ProcessMain::state;

namespace {

// ==========================================
// Function: Load and validate the top-level exposure list on rank zero
// Method: Parse quoted per-exposure list paths and reject empty or malformed
//         input before any rank enters the numerical stage scheduler.
// ==========================================
bool loadExposureList(const std::string& exposure_list, std::string& error) {
    ProcessMain::state.clear();
    std::vector<ExposureList::Entry> entries;
    if (!ExposureList::loadPipelineList(
            exposure_list, entries, 0, error)) {
        return false;
    }

    ProcessMain::state.exposure_files.reserve(entries.size());
    for (const ExposureList::Entry& entry : entries) {
        ProcessMain::state.exposure_files.push_back(entry.path);
    }

    std::cout << "Total number of EXPOSURE: "
              << ProcessMain::state.exposureCount() << std::endl;
    return true;
}

// ==========================================
// Function: Broadcast the validated exposure list to every MPI rank
// Method: Send the count, then length-prefix each mutable C++ string.
// ==========================================
bool broadcastExposureList(std::string& error) {
    return MPIUtils::broadcastStrings(
        ProcessMain::state.exposure_files, 0, error);
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
        if (MPIScheduler::state.rank == 0) {
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
                  MPIScheduler::state.communicator);
    if (global_layout_ok == 0) {
            if (MPIScheduler::state.rank == 0) {
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
    const int rank = MPIScheduler::state.rank;
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
                  MPIScheduler::state.communicator);
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
    MPI_Bcast(&load_ok, 1, MPI_INT, 0, MPIScheduler::state.communicator);
    if (load_ok == 0) {
        if (rank == 0) {
            std::cerr << load_error << std::endl;
        }
        return 1;
    }

    std::string broadcast_error;
    if (!broadcastExposureList(broadcast_error)) {
        std::cerr << "Exposure-list broadcast error on rank " << rank << ": "
                  << broadcast_error << std::endl;
        return 1;
    }
    const int exposure_count = ProcessMain::state.exposureCount();
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
        MPIScheduler::distribute(exposure_count, PreProcess::preProcess, "Pre-process...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 3 == 0) {
        MPIScheduler::distribute(exposure_count, Astrometry::procAstrometry, "Astrometry...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 5 == 0) {
        MPIScheduler::distribute(exposure_count, SourceExtractor::procSource, "Sources ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 7 == 0) {
        MPIScheduler::distribute(exposure_count, FourierTransformSt1::procFourierTSt1, "FFT st1...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 11 == 0) {
        MPIScheduler::distribute(exposure_count, PSFModel::procPSF, "PSF ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 13 == 0) {
        MPIScheduler::distribute(exposure_count, FourierTransformSt2::procFourierTSt2, "FFT st2...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 17 == 0) {
        MPIScheduler::distribute(exposure_count, ShearMeasurement::procShear, "Shear ...");
    }
    MPIScheduler::barrier();

    if (lensing.process_stage % 19 == 0) {
        const int mpi_parameter_count = exposure_count * 6;
        const std::size_t parameter_count =
            static_cast<std::size_t>(mpi_parameter_count);

        ExposureInfo::state.parameters.assign(parameter_count, 0.0f);
        MPIScheduler::distribute(exposure_count, ExposureInfo::procInfo, "Info ...");
        MPIScheduler::barrier();

        std::vector<float> reduced_exposure_parameters(parameter_count, 0.0f);
        MPI_Allreduce(ExposureInfo::state.parameters.data(),
                      reduced_exposure_parameters.data(),
                      mpi_parameter_count, MPI_FLOAT, MPI_SUM,
                      MPIScheduler::state.communicator);
        ExposureInfo::state.parameters = std::move(reduced_exposure_parameters);

        if (rank == 0) {
            const std::string root_directory = UniversalUtils::getDir(exposure_list, 1);
            const std::string filename = root_directory + "/expo_info.dat";
            MainIO::OutputFile output(filename);
            output << std::setprecision(10);
            output << "N-valid-chip PSF-FWHM(arcsec) chi_d-stars nstar-per-chip "
                      "cRVAL1 cRVAL2 expo_name\n";
            for (int exposure = 0; exposure < exposure_count; ++exposure) {
                for (int parameter = 0; parameter < 6; ++parameter) {
                    output << ExposureInfo::state.parameters[exposure * 6 + parameter] << " ";
                }
                output << ProcessMain::state.exposure_files[exposure] << "\n";
            }
            output.close();
        }
        MPIScheduler::barrier();
    }

    if (lensing.process_stage % 23 == 0) {
        MPIScheduler::distribute(exposure_count, CatalogCombiner::procComb, "combine ...");
    }
    MPIScheduler::barrier();
    return 0;
}
