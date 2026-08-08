#include "process_fd/process_fd.hpp"
#include "FDConfig.hpp"
#include "process_fd/FDData.hpp"
#include "process_fd/ShearCatalogReader.hpp"
#include "process_fd/StarCutCalculator.hpp"
#include "process_fd/KMeansClusterer.hpp"
#include "process_fd/FDMeasurement.hpp"
#include "process_main/MPIScheduler.hpp"
#include "process_main/NumericalRecipes.hpp"

#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fc = FDConfig;
namespace lc = LensingConfig;

// ==========================================
// Function: Load and broadcast the exposure list on rank zero
// ==========================================
namespace {

bool loadExposureList(const std::string& path, std::vector<std::string>& files,
                      int rank) {
    files.clear();
    if (rank == 0) {
        std::ifstream input(path);
        if (!input.is_open()) {
            std::cerr << "EXPO_LIST reading error: " << path << std::endl;
            return false;
        }
        std::string name;
        while (input >> name) {
            if (name.size() >= 2 && name.front() == '"'
                && name.back() == '"') {
                name = name.substr(1, name.size() - 2);
            }
            files.push_back(name);
        }
        if (files.empty()) {
            std::cerr << "EXPO_LIST contains no exposures: " << path << std::endl;
            return false;
        }
        if (static_cast<int>(files.size()) > LensingConfig::NMAX_EXPO) {
            std::cerr << "EXPO_LIST exceeds NMAX_EXPO: " << path << std::endl;
            return false;
        }
        std::cout << "Total number of EXPOSURE: " << files.size() << std::endl;
    }
    return true;
}

void broadcastExposureList(std::vector<std::string>& files, int rank) {
    int n_expo = static_cast<int>(files.size());
    MPI_Bcast(&n_expo, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) files.resize(n_expo);
    for (int i = 0; i < n_expo; ++i) {
        int len = (rank == 0) ? static_cast<int>(files[i].size()) : 0;
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0) files[i].resize(len);
        if (len > 0)
            MPI_Bcast(files[i].data(), len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
}

}  // namespace

// ==========================================
// Function: process_fd
// Method: Run the FD (field-distortion) shear test. Read per-exposure shear
//         catalogs, remove point sources via star-bar fitting (single or
//         per-exposure), recover mean shear per spatial bin (PDF or
//         jackknife), and write FD_test_comb.dat.
// ==========================================
int process_fd(const std::string& exposure_list,
               const ProcessConfig::RuntimeOptions& options,
               const std::string& dataset_root,
               MPI_Comm communicator) {
    const int rank = MPIScheduler::my_id;
    const int num_procs = MPIScheduler::num_procs;

    // 1. Load and broadcast exposure list
    std::vector<std::string> expo_files;
    bool ok = loadExposureList(exposure_list, expo_files, rank);
    int global_ok = 0;
    int local_ok = ok ? 1 : 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (global_ok == 0) return 1;
    broadcastExposureList(expo_files, rank);
    int n_expo = static_cast<int>(expo_files.size());

    // Initialize RNG
    NumericalRecipes::initializeRan1Seed(rank, num_procs);
    MPIScheduler::barrier();

    // 2. Allocate data arrays and read catalogs
    FDData data;
    data.reserve(fc::nmax_per_core);
    data.ng = 0;

    if (rank == 0)
        std::cout << "Reading " << n_expo << " shear catalogs..." << std::endl;

    MPIScheduler::distribute(n_expo,
        [&](int iexpo) {
            ShearCatalogReader::readExposure(iexpo, data, expo_files, rank);
        },
        "reading FD cat");

    MPIScheduler::barrier();
    if (rank == 0)
        std::cout << "Catalog reading complete. Local ng=" << data.ng << std::endl;

    // 3. Calculate star cut
    if (fc::FD_PER_EXPOSURE_STAR_BAR) {
        std::vector<float> S_mean_arr, S_std_arr, S_cut_arr;
        StarCutCalculator::calculateGlobalStarCutAuto(data, rank, num_procs,
                                                      S_mean_arr, S_std_arr,
                                                      S_cut_arr);
        // Apply advanced cuts (per-exposure star cut + SNR cuts)
        StarCutCalculator::applyAdvancedCuts(data, S_cut_arr);
    } else {
        float S_mean = 0.0, S_std = 0.0, S_cut = 0.0;
        StarCutCalculator::calculateGlobalStarCut(data, rank, num_procs,
                                                   S_mean, S_std, S_cut);
        StarCutCalculator::applySingleStarCut(data, S_cut);
    }

    MPIScheduler::barrier();
    if (rank == 0)
        std::cout << "After star cut: local ng=" << data.ng << std::endl;

    // 4-5. K-means jackknife region clustering + label assignment
    //      (skipped in PDF_SIGMA mode where jackknife is not used)
    if constexpr (fc::FD_USE_JACKKNIFE) {
        std::vector<float> centers;
        KMeansClusterer::runMPI(data.ng, data.rra, data.ddec, rank, num_procs,
                                centers);
        MPIScheduler::barrier();

        for (int idx = 0; idx < data.ng; ++idx) {
            float ra_rad = data.rra[idx] * LensingConfig::arc_convert;
            float dec_rad = data.ddec[idx] * LensingConfig::arc_convert;
            float pos[3] = {
                std::cos(dec_rad) * std::cos(ra_rad),
                std::cos(dec_rad) * std::sin(ra_rad),
                std::sin(dec_rad)
            };
            float max_dot = -2.0;
            int best_k = 0;
            for (int k = 0; k < fc::N_jack; ++k) {
                float dot = pos[0] * centers[0 * fc::N_jack + k]
                          + pos[1] * centers[1 * fc::N_jack + k]
                          + pos[2] * centers[2 * fc::N_jack + k];
                if (dot > max_dot) { max_dot = dot; best_k = k; }
            }
            data.labels[idx] = best_k;
        }
    }

    // 6. Run shear measurement for g1 and g2
    FDMeasurement measurer(rank, num_procs);
    int nbin = fc::fd_num;

    // g1 component
    std::vector<float> s_arr(nbin * 4, 0.0);
    measurer.plotComparison(data.ng, fc::PDF_BINS, data.x1, data.y1, data.de1,
                            data.labels, nbin, fc::N_jack, s_arr);

    std::vector<float> comb_data(nbin * 8, 0.0);
    if (rank == 0)
        for (int i = 0; i < nbin; ++i) {
            comb_data[i * 8 + 0] = s_arr[i * 4 + 3];  // N1
            comb_data[i * 8 + 1] = s_arr[i * 4 + 0];  // gf1 center
            comb_data[i * 8 + 2] = s_arr[i * 4 + 1];  // g1(GAL)
            comb_data[i * 8 + 3] = s_arr[i * 4 + 2];  // sigma1
        }
    MPIScheduler::barrier();

    // g2 component
    std::fill(s_arr.begin(), s_arr.end(), 0.0);
    measurer.plotComparison(data.ng, fc::PDF_BINS, data.x2, data.y2, data.de2,
                            data.labels, nbin, fc::N_jack, s_arr);
    if (rank == 0)
        for (int i = 0; i < nbin; ++i) {
            comb_data[i * 8 + 4] = s_arr[i * 4 + 3];  // N2
            comb_data[i * 8 + 5] = s_arr[i * 4 + 0];  // gf2 center
            comb_data[i * 8 + 6] = s_arr[i * 4 + 1];  // g2(GAL)
            comb_data[i * 8 + 7] = s_arr[i * 4 + 2];  // sigma2
        }
    MPIScheduler::barrier();

    // 7. Write output
    if (rank == 0) {
        const std::string base_dir_str(options.fd_output_base_directory);
        const std::filesystem::path base_dir =
            base_dir_str.empty() ? std::filesystem::path(dataset_root)
                                 : std::filesystem::path(base_dir_str);
        std::filesystem::path configured_output(options.fd_output_directory);
        if (configured_output.empty()) {
            configured_output = base_dir;
        } else if (configured_output.is_relative()) {
            configured_output = base_dir / configured_output;
        }
        const std::string output_directory =
            std::filesystem::absolute(configured_output)
                .lexically_normal().string();
        std::filesystem::create_directories(output_directory);
        std::string filename =
            (std::filesystem::path(output_directory) / "FD_test_comb.dat").string();
        std::ofstream out(filename, std::ios::trunc);
        if (!out.is_open()) {
            std::cerr << "Error opening output: " << filename << std::endl;
            return 1;
        }
        out << "Selected_NUM1  g1(FD)  g1(GAL)  sigma1"
            << "  Selected_NUM2  g2(FD)  g2(GAL)  sigma2\n";
        for (int i = 0; i < nbin; ++i) {
            for (int j = 0; j < 8; ++j)
                out << comb_data[i * 8 + j] << "  ";
            out << "\n";
        }
        out.close();
        std::cout << "FD test results written to: " << filename << std::endl;
    }

    MPIScheduler::barrier();
    return 0;
}
