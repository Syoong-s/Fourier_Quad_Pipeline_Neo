#include "PSFModel.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "FitsIO.hpp"
#include "Astrometry.hpp"
#include "NumericalRecipes.hpp"
#include "UniversalUtils.hpp"
#include "ImageProcessing.hpp"
#include "ExStar.hpp"
#include "LinearSolve.hpp"
#include <mpi.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <complex>
#include <memory>

// Extern variables defined elsewhere (e.g. main.cpp)
extern std::vector<std::string> EXPO_FILE;

namespace PSFModel {

    // Global PSF PCA components and coefficients loaded from disk (replacing psf_storage_mod)
    std::vector<double> global_components;
    std::vector<double> global_mean_psf;
    std::vector<float> global_poly_coefs;
    bool is_data_loaded = false;

    // Heap-allocated structure representing the exposure-wide star catalog & chi^2 grid to prevent stack overflow
    struct ExposurePSFState {
        std::vector<int> nstar;
        std::vector<double> star_para; // Size: NMAX_CHIP * nstar_max * npara
        std::vector<float> chi_d;      // Size: NMAX_CHIP * nstar_max * nstar_max

        ExposurePSFState() {
            nstar.assign(LensingConfig::NMAX_CHIP, 0);
            star_para.assign(static_cast<size_t>(LensingConfig::NMAX_CHIP) * LensingConfig::nstar_max * LensingConfig::npara, 0.0);
            chi_d.assign(static_cast<size_t>(LensingConfig::NMAX_CHIP) * LensingConfig::nstar_max * LensingConfig::nstar_max, 0.0f);
        }

        double& getStarPara(int chip, int star, int para) {
            return star_para[((static_cast<size_t>(chip) * LensingConfig::nstar_max) + star) * LensingConfig::npara + para];
        }

        const double& getStarPara(int chip, int star, int para) const {
            return star_para[((static_cast<size_t>(chip) * LensingConfig::nstar_max) + star) * LensingConfig::npara + para];
        }

        float& getChiD(int chip, int star1, int star2) {
            return chi_d[((static_cast<size_t>(chip) * LensingConfig::nstar_max) + star1) * LensingConfig::nstar_max + star2];
        }
    };

    // Forward declarations of local helper functions
    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void starSelection(int nchip, ExposurePSFState& state);
    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void saveRescaleFactor(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void makePSFHybrid(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);

    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics);
    LinearSolve::SolveStatus interpolatePSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);
   void getPSFModelVeryLocal(const std::vector<float>& psfmap, double x, double y, std::vector<float>& model, double& dmax, int stride);
    void genPSFFits(std::vector<float>& psfmap, int nums, int nx, int ny, const std::vector<float>& star, const std::vector<std::array<double, 2>>& posi);

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio);
    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio);
    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio);
    void getPSFFWHM(const std::vector<float>& power, double& FWHM);

    // ==========================================
    // Function: Validate one PSF fitting sample
    // Method: Require finite position, shape diagnostics, and every power-spectrum pixel before retaining the star.
    // ==========================================
    static bool isFinitePSFStar(const std::vector<float>& star, std::size_t offset, int pixel_count,
                                double x, double y, double size, double e1, double e2) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(size) ||
            !std::isfinite(e1) || !std::isfinite(e2) ||
            offset + static_cast<std::size_t>(pixel_count) > star.size()) {
            return false;
        }
        for (int idx = 0; idx < pixel_count; ++idx) {
            if (!std::isfinite(star[offset + idx])) {
                return false;
            }
        }
        return true;
    }

    // Load and broadcast PSF PCA components
    void initAndLoadAllPSF(const std::string& dirOutput, int myRank) {
        if (is_data_loaded) return;

        if (myRank == 0) {
            std::cout << "Allocating memory on all ranks..." << std::endl;
        }

        global_components.assign(static_cast<size_t>(LensingConfig::NMAX_CHIP) * LensingConfig::nsns * LensingConfig::n_pcs, 0.0);
        global_mean_psf.assign(static_cast<size_t>(LensingConfig::NMAX_CHIP) * LensingConfig::nsns, 0.0);
        global_poly_coefs.assign(static_cast<size_t>(LensingConfig::NMAX_CHIP) * 2 * 2 * LensingConfig::n_pcs * LensingConfig::npp6th, 0.0f);

        if (myRank == 0) {
            std::cout << "Rank 0 is reading files from disk..." << std::endl;

            for (int i_ccd = 1; i_ccd <= LensingConfig::NMAX_CHIP; ++i_ccd) {
                if (i_ccd == 2 || i_ccd == 61) continue;

                std::ostringstream ss_ccd;
                ss_ccd << std::setfill('0') << std::setw(2) << i_ccd;
                std::string filename = dirOutput + "/stamps/dat_Pcs/pcs_ccd" + ss_ccd.str() + ".dat";

                std::ifstream infile(filename);
                if (infile.is_open()) {
                    for (int k = 0; k < LensingConfig::nsns; ++k) {
                        for (int j = 0; j < LensingConfig::n_pcs; ++j) {
                            infile >> global_components[getCompIndex(i_ccd - 1, k, j)];
                        }
                        infile >> global_mean_psf[getMeanIndex(i_ccd - 1, k)];
                    }
                    infile.close();

                    if (global_components[getCompIndex(i_ccd - 1, 0, 0)] < -1.0e20) {
                        std::cout << "CCD " << i_ccd << " has bad PCS data." << std::endl;
                    }
                } else {
                    global_components[getCompIndex(i_ccd - 1, 0, 0)] = -1.0e30;
                    std::cout << "CCD " << i_ccd << " PCS data read error: " << filename << std::endl;
                    std::exit(1);
                }

                for (int bx = 1; bx <= LensingConfig::nblocks; ++bx) {
                    for (int by = 1; by <= LensingConfig::nblocks; ++by) {
                        std::string filename_coeff = dirOutput + "/stamps/dat_Pcs/coeff_ccd" + ss_ccd.str() + "_" + std::to_string(bx) + std::to_string(by) + ".dat";
                        std::ifstream coeff_file(filename_coeff);
                        if (coeff_file.is_open()) {
                            for (int u = 0; u < LensingConfig::n_pcs; ++u) {
                                for (int j = 0; j < LensingConfig::npp6th; ++j) {
                                    coeff_file >> global_poly_coefs[getPolyIndex(i_ccd - 1, bx - 1, by - 1, u, j)];
                                }
                                if (std::isnan(global_poly_coefs[getPolyIndex(i_ccd - 1, bx - 1, by - 1, u, 0)])) {
                                    global_poly_coefs[getPolyIndex(i_ccd - 1, bx - 1, by - 1, 0, 0)] = -1.0e30f;
                                }
                            }
                            coeff_file.close();

                            if (global_poly_coefs[getPolyIndex(i_ccd - 1, bx - 1, by - 1, 0, 0)] < -1.0e20f) {
                                std::cout << "CCD " << i_ccd << " field " << bx << " " << by << " has bad polynomial data." << std::endl;
                            }
                        } else {
                            global_poly_coefs[getPolyIndex(i_ccd - 1, bx - 1, by - 1, 0, 0)] = -1.0e30f;
                            std::cout << "CCD " << i_ccd << " field " << bx << " " << by << " polynomial data read error: " << filename_coeff << std::endl;
                            std::exit(1);
                        }
                    }
                }
            }
            std::cout << "Rank 0 finished reading. Starting Broadcast..." << std::endl;
        }

        int total_count = LensingConfig::NMAX_CHIP * LensingConfig::nsns * LensingConfig::n_pcs;
        MPI_Bcast(global_components.data(), total_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        total_count = LensingConfig::NMAX_CHIP * LensingConfig::nsns;
        MPI_Bcast(global_mean_psf.data(), total_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        total_count = LensingConfig::NMAX_CHIP * 2 * 2 * LensingConfig::n_pcs * LensingConfig::npp6th;
        MPI_Bcast(global_poly_coefs.data(), total_count, MPI_FLOAT, 0, MPI_COMM_WORLD);

        is_data_loaded = true;
        MPI_Barrier(MPI_COMM_WORLD);

        if (myRank == 0) {
            std::cout << "Broadcast finished. All ready." << std::endl;
        }
    }

    void freePSFMemory() {
        global_components.clear();
        global_mean_psf.clear();
        global_poly_coefs.clear();
        is_data_loaded = false;
    }

    // Stage 5 main entry
    void procPSF(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = EXPO_FILE[iexpo - 1];
        std::vector<std::string> imageFiles;
        std::string dirOutput;
        UniversalUtils::getImageList(expo_file_path, imageFiles, dirOutput);

        int nchip = static_cast<int>(imageFiles.size());

        auto state_ptr = std::make_unique<ExposurePSFState>();
        ExposurePSFState& state = *state_ptr;

        int nc = 0;
        std::vector<std::array<double, 4>> p_chip(LensingConfig::NMAX_CHIP, {0.0, 0.0, 0.0, 0.0});

        readInCandidates(nchip, imageFiles, dirOutput, nc, p_chip, state);

        starSelection(nchip, state);

        plotStarExpo(nchip, imageFiles, dirOutput, state);

        plotStars(nchip, imageFiles, dirOutput, nc, p_chip, state);

        if (LensingConfig::PSF_Ms == 1) {
            saveRescaleFactor(nchip, imageFiles, dirOutput, state);
        }

        if (LensingConfig::PSF_type == 1) {
            makePSFLocalFit(nchip, imageFiles, dirOutput, state);
        } else if (LensingConfig::PSF_type == 2) {
            makePSFHybrid(nchip, imageFiles, dirOutput, state);
        } else {
            std::cerr << "Error / PSF Invalid PSF fitting method!" << std::endl;
            std::exit(1);
        }
    }

    // Coordinate scaling helpers (used by both PSFModel and PSFRecons)
    void PSF_rescale(std::vector<float>& resi, float refactor) {
        int ns = LensingConfig::ns;
        double xc = ns / 2.0;
        double yc = ns / 2.0;

        std::vector<float> temp(ns * ns, 0.0f);

        for (int row = 0; row < ns; ++row) {
            for (int col = 0; col < ns; ++col) {
                double rx = (col - xc) * refactor + xc;
                double ry = (row - yc) * refactor + yc;
                int lx = static_cast<int>(std::floor(rx));
                int ly = static_cast<int>(std::floor(ry));

                float val1 = 0.0f;
                float val2 = 0.0f;
                float val3 = 0.0f;
                float val4 = 0.0f;

                if (lx >= 0 && lx < ns && ly >= 0 && ly < ns) {
                    val1 = resi[ly * ns + lx];
                }
                if (lx + 1 >= 0 && lx + 1 < ns && ly >= 0 && ly < ns) {
                    val2 = resi[ly * ns + (lx + 1)];
                }
                if (lx >= 0 && lx < ns && ly + 1 >= 0 && ly + 1 < ns) {
                    val3 = resi[(ly + 1) * ns + lx];
                }
                if (lx + 1 >= 0 && lx + 1 < ns && ly + 1 >= 0 && ly + 1 < ns) {
                    val4 = resi[(ly + 1) * ns + (lx + 1)];
                }

                double d1 = (lx + 1) - rx;
                double d2 = rx - lx;
                double d3 = (ly + 1) - ry;
                double d4 = ry - ly;

                float w1 = static_cast<float>(val1 * d1 * d3);
                float w2 = static_cast<float>(val2 * d2 * d3);
                float w3 = static_cast<float>(val3 * d1 * d4);
                float w4 = static_cast<float>(val4 * d2 * d4);

                temp[row * ns + col] = w1 + w2 + w3 + w4;
            }
        }
        resi = temp;
    }

    void PSF_unscale(std::vector<float>& model_in, float scale_factor) {
        int ns = LensingConfig::ns;
        double xc = ns / 2.0;
        double yc = ns / 2.0;

        float inv_scale = 1.0f / scale_factor;
        std::vector<float> psf_out(ns * ns, 0.0f);

        for (int row = 0; row < ns; ++row) {
            for (int col = 0; col < ns; ++col) {
                double rx = (col - xc) * inv_scale + xc;
                double ry = (row - yc) * inv_scale + yc;
                int lx = static_cast<int>(std::floor(rx));
                int ly = static_cast<int>(std::floor(ry));

                float val1 = 0.0f;
                float val2 = 0.0f;
                float val3 = 0.0f;
                float val4 = 0.0f;

                if (lx >= 0 && lx < ns && ly >= 0 && ly < ns) {
                    val1 = model_in[ly * ns + lx];
                }
                if (lx + 1 >= 0 && lx + 1 < ns && ly >= 0 && ly < ns) {
                    val2 = model_in[ly * ns + (lx + 1)];
                }
                if (lx >= 0 && lx < ns && ly + 1 >= 0 && ly + 1 < ns) {
                    val3 = model_in[(ly + 1) * ns + lx];
                }
                if (lx + 1 >= 0 && lx + 1 < ns && ly + 1 >= 0 && ly + 1 < ns) {
                    val4 = model_in[(ly + 1) * ns + (lx + 1)];
                }

                double d1 = (lx + 1) - rx;
                double d2 = rx - lx;
                double d3 = (ly + 1) - ry;
                double d4 = ry - ly;

                float w1 = static_cast<float>(val1 * d1 * d3);
                float w2 = static_cast<float>(val2 * d2 * d3);
                float w3 = static_cast<float>(val3 * d1 * d4);
                float w4 = static_cast<float>(val4 * d2 * d4);

                psf_out[row * ns + col] = w1 + w2 + w3 + w4;
            }
        }
        model_in = psf_out;
    }

    // Local Helper Routines

    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int len_s = LensingConfig::len_s;
        int nstar_max = LensingConfig::nstar_max;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string headname = dirOutput + "/astrometry/Head/" + prefix_e + ".head";

        nc = 0;

        for (int k = 0; k < nchip; ++k) {
            state.nstar[k] = 0;

            double cRPIX[2] = {0.0, 0.0};
            double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
            double cRVAL[2] = {0.0, 0.0};
            double PU[2][LensingConfig::npd];
            int ierror = 0;

            Astrometry::readAstrometryPara(headname, k + 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd, ierror);

            if (ierror == 1) continue;

            nc++;
            double x = 1.0;
            double y = 1.0;
            double xx = 0.0, yy = 0.0;
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][0] = xx;
            p_chip[nc - 1][1] = yy;

            x = 2046.0;
            y = 4094.0;
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][2] = xx;
            p_chip[nc - 1][3] = yy;

            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/dat_StarCanInfo", prefix, "_star_can_info.dat");

            std::ifstream infile(filepath);
            if (!infile.is_open()) {
                std::cerr << filepath << "\n";
                std::cerr << "Error / PSF star_can_info catalog file error!!" << std::endl;
                std::exit(1);
            }

            // Skip header line
            std::string header;
            std::getline(infile, header);

            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                float aa[4] = {0.0f, 0.0f, 0.0f, 0.0f};
                if (iss >> aa[0] >> aa[1] >> aa[2] >> aa[3]) {
                    int star_idx = state.nstar[k];
                    if (star_idx < nstar_max) {
                        state.getStarPara(k, star_idx, 0) = aa[0];
                        state.getStarPara(k, star_idx, 1) = aa[1];
                        state.getStarPara(k, star_idx, 2) = aa[2];
                        state.getStarPara(k, star_idx, 3) = aa[3];
                        state.nstar[k]++;
                    }
                }
            }
            infile.close();

            if (state.nstar[k] > 0) {
                int nn1 = ns * len_s;
                int nn2 = ns * ((state.nstar[k] / len_s) + 1);
                std::string stampPath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                std::vector<float> star(static_cast<size_t>(nstar_max) * ns * ns, 0.0f);
                if (!FitsIO::readStamps(nstar_max, 1, state.nstar[k], ns, ns, star, nn1, nn2, stampPath)) {
                    std::cerr << "readInCandidates: readStamps failed for " << stampPath << std::endl;
                    std::exit(1);
                }

                for (int i = 0; i < state.nstar[k]; ++i) {
                    state.getStarPara(k, i, 4) = 1.0;

                    std::vector<float> source_p(ns * ns);
                    for (int v = 0; v < ns; ++v) {
                        for (int u = 0; u < ns; ++u) {
                            source_p[v * ns + u] = star[static_cast<size_t>(i) * ns * ns + v * ns + u];
                        }
                    }

                    x = state.getStarPara(k, i, 1);
                    y = state.getStarPara(k, i, 2);
                    Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
                    state.getStarPara(k, i, 5) = xx;
                    state.getStarPara(k, i, 6) = yy;

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, source_p, ee, size, 0.02f);
                    state.getStarPara(k, i, 7) = size;
                    state.getStarPara(k, i, 8) = ee[0];
                    state.getStarPara(k, i, 9) = ee[1];

                    double FWHM = 0.0;
                    getPSFFWHM(source_p, FWHM);
                    state.getStarPara(k, i, 10) = FWHM;

                    double temp = 0.0;
                    for (int idx = 0; idx < ns * ns; ++idx) {
                        temp += source_p[idx];
                    }
                    state.getStarPara(k, i, 11) = 1.0 / temp;
                }

                for (int i = 0; i < nstar_max; ++i) {
                    for (int j = 0; j < nstar_max; ++j) {
                        state.getChiD(k, i, j) = 0.0f;
                    }
                }

                for (int i = 0; i < state.nstar[k] - 1; ++i) {
                    for (int j = i + 1; j < state.nstar[k]; ++j) {
                        std::vector<float> map1(ns * ns);
                        std::vector<float> map2(ns * ns);
                        double sp_i_12 = state.getStarPara(k, i, 11);
                        double sp_j_12 = state.getStarPara(k, j, 11);
                        for (int idx = 0; idx < ns * ns; ++idx) {
                            map1[idx] = static_cast<float>(star[static_cast<size_t>(i) * ns * ns + idx] * sp_i_12);
                            map2[idx] = static_cast<float>(star[static_cast<size_t>(j) * ns * ns + idx] * sp_j_12);
                        }
                        double temp_chi = 0.0;
                        UniversalUtils::anaChi2(ns, map1, map2, temp_chi);
                        float temp_val = static_cast<float>(std::sqrt(temp_chi));
                        state.getChiD(k, i, j) = temp_val;
                        state.getChiD(k, j, i) = temp_val;
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Select PSF stars from candidates
    // Method: Follow F77 star_selection thresholding and connected chi-group selection.
    // ==========================================
    void starSelection(int nchip, ExposurePSFState& state) {
        int nstar_min = LensingConfig::nstar_min;
        int nstar_max = LensingConfig::nstar_max;
        int nstar_min_local = LensingConfig::nstar_min_local;

        int ntot = 0;
        for (int k = 0; k < nchip; ++k) {
            ntot += state.nstar[k];
        }

        if (ntot < nstar_min * 2) {
            for (int k = 0; k < nchip; ++k) {
                for (int i = 0; i < state.nstar[k]; ++i) {
                    state.getStarPara(k, i, 4) = -1.0;
                }
            }
            return;
        }

        // Determine size threshold
        std::vector<float> tmp_size;
        tmp_size.reserve(ntot);
        for (int k = 0; k < nchip; ++k) {
            for (int i = 0; i < state.nstar[k]; ++i) {
                tmp_size.push_back(static_cast<float>(state.getStarPara(k, i, 7))); // F77 index 8 size
            }
        }
        NumericalRecipes::sort(tmp_size);
        int thresh_size_idx = static_cast<int>((tmp_size.size() * 2) / 3) - 1;
        thresh_size_idx = std::max(0, std::min(thresh_size_idx, static_cast<int>(tmp_size.size()) - 1));
        float thresh_size = tmp_size[thresh_size_idx];

        // Determine chi^2 threshold
        std::vector<float> chimin(static_cast<size_t>(LensingConfig::NMAX_CHIP) * nstar_max, 1000.0f);
        std::vector<float> tmp_chi;

        for (int k = 0; k < nchip; ++k) {
            for (int i = 0; i < state.nstar[k] - 1; ++i) {
                for (int j = i + 1; j < state.nstar[k]; ++j) {
                    float temp = state.getChiD(k, i, j);
                    chimin[k * nstar_max + i] = std::min(temp, chimin[k * nstar_max + i]);
                    chimin[k * nstar_max + j] = std::min(temp, chimin[k * nstar_max + j]);

                    if (state.getStarPara(k, i, 7) < thresh_size) continue;
                    if (state.getStarPara(k, j, 7) < thresh_size) continue;

                    tmp_chi.push_back(temp);
                }
            }
        }

        float peak = 0.0f, sig = 0.0f;
        NumericalRecipes::getPeakWidthLowSide(tmp_chi, peak, sig);
        float thresh_chi = peak + 4.0f * sig;

        for (int k = 0; k < nchip; ++k) {
            int n_valid_local = 0;
            for (int i = 0; i < state.nstar[k]; ++i) {
                if (chimin[k * nstar_max + i] > thresh_chi) {
                    state.getStarPara(k, i, 4) = -1.0;
                    continue;
                }
                n_valid_local++;
            }

            if (n_valid_local < nstar_min_local) {
                for (int i = 0; i < state.nstar[k]; ++i) {
                    state.getStarPara(k, i, 4) = -1.0;
                }
                continue;
            }

            std::vector<int> id;
            for (int i = 0; i < state.nstar[k]; ++i) {
                if (state.getStarPara(k, i, 4) < 0.0) continue;
                id.push_back(i);
            }

            int ntot_local = static_cast<int>(id.size());
            std::vector<int> group_id(ntot_local, 0);
            int max_group_id = 0;

            for (int i = 0; i < ntot_local; ++i) {
                if (group_id[i] == 0) {
                    max_group_id++;
                    group_id[i] = max_group_id;
                }
                for (int j = i + 1; j < ntot_local; ++j) {
                    if (group_id[j] == group_id[i]) continue;
                    if (state.getChiD(k, id[i], id[j]) > thresh_chi) continue;
                    if (group_id[j] == 0) {
                        group_id[j] = group_id[i];
                    } else {
                        int u = group_id[j];
                        int v = group_id[i];
                        for (int w = 0; w < ntot_local; ++w) {
                            if (group_id[w] == u) {
                                group_id[w] = v;
                            }
                        }
                    }
                }
            }

            std::vector<int> group_size(max_group_id + 1, 0);
            for (int i = 0; i < ntot_local; ++i) {
                int u = group_id[i];
                group_size[u]++;
            }

            int max_gsize = 0;
            int max_gid = 0;
            int max2_gsize = 0;
            int max2_gid = 0;

            if (state.nstar[k] > 0 && max_group_id > 0) {
                max_gsize = group_size[1];
                max_gid = 1;
                for (int i = 2; i <= max_group_id; ++i) {
                    if (group_size[i] > max_gsize) {
                        max2_gsize = max_gsize;
                        max2_gid = max_gid;
                        max_gsize = group_size[i];
                        max_gid = i;
                    } else if (group_size[i] > max2_gsize) {
                        max2_gsize = group_size[i];
                        max2_gid = i;
                    }
                }
            }

            for (int i = 0; i < ntot_local; ++i) {
                int j = id[i];
                if (group_id[i] != max_gid) {
                    state.getStarPara(k, j, 4) = -1.0;
                }
            }
        }

        // Final local star counts pass
        for (int k = 0; k < nchip; ++k) {
            int n = 0;
            for (int i = 0; i < state.nstar[k]; ++i) {
                if (state.getStarPara(k, i, 4) < 0.0) continue;
                n++;
            }
            if (n < nstar_min_local) {
                for (int i = 0; i < state.nstar[k]; ++i) {
                    state.getStarPara(k, i, 4) = -1.0;
                }
            }
        }
    }

    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int len_s = LensingConfig::len_s;
        int nstar_max = LensingConfig::nstar_max;
        int nmax_stamp = 5000;

        size_t total_stars_limit = static_cast<size_t>(LensingConfig::NMAX_CHIP) * nstar_max;
        std::vector<int> opt(total_stars_limit, 0);

        // star_test size: 62 * 2000 * 64 * 64 floats (Heap allocated)
        std::vector<float> star_test(total_stars_limit * ns * ns, 0.0f);
        std::vector<float> star(static_cast<size_t>(nstar_max) * ns * ns, 0.0f);

        int ntot = 0;
        int w = 0;
        int start = 0;

        for (int ichip = 0; ichip < nchip; ++ichip) {
            if (state.nstar[ichip] == 0) continue;
            int nn1 = ns * len_s;
            int nn2 = ns * ((state.nstar[ichip] / len_s) + 1);

            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");

            if (!FitsIO::readStamps(nstar_max, 1, state.nstar[ichip], ns, ns, star, nn1, nn2, filepath)) {
                std::cerr << "plotStarExpo: readStamps failed for " << filepath << std::endl;
                std::exit(1);
            }

            for (int i = 0; i < state.nstar[ichip]; ++i) {
                for (int u = 0; u < ns; ++u) {
                    for (int v = 0; v < ns; ++v) {
                        size_t srcIdx = static_cast<size_t>(i) * ns * ns + v * ns + u;
                        size_t destIdx = static_cast<size_t>(start + i) * ns * ns + v * ns + u;
                        star_test[destIdx] = star[srcIdx];
                    }
                }
                w = start + i;
                if (state.getStarPara(ichip, i, 4) <= 0.0) continue;
                ntot++;
                if (ntot < nmax_stamp) {
                    opt[w] = 1;
                }
            }
            start += state.nstar[ichip];
        }

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string out_filename = dirOutput + "/stamps/fits_StarP/" + prefix_e + "_star_power_expo.fits";

        if (ntot > 0) {
            int len_sam = LensingConfig::len_sam;
            int nn1 = ns * len_sam;
            int nn2 = ns * ((std::min(ntot, nmax_stamp) / len_sam) + 1);
            FitsIO::writeStamps2(total_stars_limit, w + 1, ns, ns, star_test, opt, 1, nn1, nn2, out_filename);
        }
    }

    // ==========================================
    // Function: Publish exposure-level star diagnostics and PSF visualization
    // Method: Aggregate selected stars and route text/FITS products through
    //         checked main-process writers.
    // ==========================================
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        int nm = 1000;
        int nstar_max = LensingConfig::nstar_max;
        int nstar_min_local = LensingConfig::nstar_min_local;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        // prefix_dir inlined: per-type stamps/ subdirs (reorganized layout)

        std::string info_filename = dirOutput + "/stamps/dat_StarInfo/" + prefix_e + "_star_info_expo.dat";
        MainIO::OutputFile outfile(info_filename);
        if (!outfile.is_open()) {
            std::cerr << "plotStars: Error opening " << info_filename << std::endl;
            return;
        }
        outfile << std::setprecision(10);

        outfile << "# ichip nstar FWHM e1 e2 chi_d\n";

        size_t total_stars_limit = static_cast<size_t>(LensingConfig::NMAX_CHIP) * nstar_max;
        std::vector<std::array<double, 5>> sk;
        sk.reserve(total_stars_limit);

        for (int k = 0; k < nchip; ++k) {
            double FWHM_ave = 0.0;
            double e1_ave = 0.0;
            double e2_ave = 0.0;
            double chi_d_ave = 0.0;
            int nums = 0;
            int prev_idx = -1;

            for (int i = 0; i < state.nstar[k]; ++i) {
                if (state.getStarPara(k, i, 4) <= 0.0) continue;
                nums++;
                sk.push_back({
                    state.getStarPara(k, i, 5),
                    state.getStarPara(k, i, 6),
                    state.getStarPara(k, i, 7),
                    state.getStarPara(k, i, 8),
                    state.getStarPara(k, i, 9)
                });

                FWHM_ave += state.getStarPara(k, i, 10);
                e1_ave += state.getStarPara(k, i, 8);
                e2_ave += state.getStarPara(k, i, 9);
                if (nums >= 2 && prev_idx != -1) {
                    chi_d_ave += state.getChiD(k, i, prev_idx);
                }
                prev_idx = i;
            }

            if (nums >= nstar_min_local) {
                FWHM_ave /= nums;
                e1_ave /= nums;
                e2_ave /= nums;
                chi_d_ave /= (nums - 1.0);
                outfile << (k + 1) << " " << nums << " "
                        << std::scientific << std::setprecision(10)
                        << FWHM_ave << " " << e1_ave << " " << e2_ave << " " << chi_d_ave << "\n";
            } else {
                outfile << (k + 1) << " 0 -99.0 -99.0 -99.0 -99.0\n";
            }
        }

        // std::cout << prefix_dir << " total no. of stars: " << sk.size() << std::endl;
        outfile.close();

        std::vector<float> PSFmap;
        ImageProcessing::drawShearExpo(nm, PSFmap, p_chip, sk, 200.0, 1.0);

        std::string fits_filename = dirOutput + "/stamps/fits_PsfSrc/" + prefix_e + "_PSF_source.fits";
        FitsIO::writeImage(fits_filename, nm, nm, PSFmap);
    }

    // ==========================================
    // Function: Publish one exposure's PSF residual rescale factor
    // Method: Average valid stellar sizes and serialize the derived factor with
    //         the checked main-process text writer.
    // ==========================================
    void saveRescaleFactor(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        double summ = 0.0;
        double ntot = 0.0;
        for (int ichip = 0; ichip < nchip; ++ichip) {
            if (state.nstar[ichip] == 0) continue;
            for (int i = 0; i < state.nstar[ichip]; ++i) {
                if (state.getStarPara(ichip, i, 4) <= 0.0) continue;
                summ += state.getStarPara(ichip, i, 10);
                ntot += 1.0;
            }
        }
        if (ntot > 0.0) {
            summ = summ / ntot;
            summ = LensingConfig::rescale_size / summ;
        } else {
            summ = 0.0;
        }

        std::string prefix = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string filename = dirOutput + "/stamps/dat_Rescale/" + prefix + "_factor.dat";
        MainIO::OutputFile outfile(filename);
        if (outfile.is_open()) {
            outfile << std::scientific << std::setprecision(10) << summ << std::endl;
            outfile.close();
        } else {
            std::cerr << "saveRescaleFactor: Error opening file " << filename << std::endl;
        }
    }

    // ==========================================
    // Function: Fit and serialize local PSF models.
    // Method: Preserve F77 model layout with 17-digit double serialization.
    // ==========================================
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int len_s = LensingConfig::len_s;
        int nstar_max = LensingConfig::nstar_max;
        int npl = LensingConfig::npl;
        int nplx = LensingConfig::nplx;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string comp_filename = dirOutput + "/stamps/dat_StarComp/" + prefix_e + "_star_comp_expo.dat";
        MainIO::OutputFile file90(comp_filename);
        if (!file90.is_open()) {
            std::cerr << "makePSFLocalFit: Error opening " << comp_filename << std::endl;
            return;
        }
        file90 << std::setprecision(17);

        float res_factor = 1.0f;
        if (LensingConfig::PSF_Ms == 1) {
            std::string rescale_filename = dirOutput + "/stamps/dat_Rescale/" + prefix_e + "_factor.dat";
            std::ifstream rescale_file(rescale_filename);
            if (rescale_file.is_open()) {
                rescale_file >> res_factor;
                rescale_file.close();
            } else {
                std::cerr << "cannot find rescale factor file" << std::endl;
                std::exit(1);
            }
        }

        for (int k = 0; k < nchip; ++k) {
            int nums = 0;
            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);

            std::vector<float> star(static_cast<size_t>(nstar_max) * ns * ns, 0.0f);
            if (state.nstar[k] > 0) {
                int nn1 = ns * len_s;
                int nn2 = ns * ((state.nstar[k] / len_s) + 1);
                std::string filepath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                if (!FitsIO::readStamps(nstar_max, 1, state.nstar[k], ns, ns, star, nn1, nn2, filepath)) {
                    std::cerr << "makePSFLocalFit: readStamps failed for " << filepath << std::endl;
                    std::exit(1);
                }
            }

            std::string coe_filename = OutputLayout::chipPath(
                dirOutput, "stamps/dat_PsfFit", prefix, "_PSF_coe_local.dat");
            MainIO::OutputFile file10(coe_filename);
            if (!file10.is_open()) {
                std::cerr << "makePSFLocalFit: Error opening " << coe_filename << std::endl;
                continue;
            }
            file10 << std::setprecision(17);

            std::vector<std::array<double, 2>> posi;
            std::vector<std::array<double, 3>> sshape;
            std::vector<float> star_local;
            int removed_non_finite = 0;

            for (int i = 0; i < state.nstar[k]; ++i) {
                if (state.getStarPara(k, i, 4) < 0.0) continue;
                double px = state.getStarPara(k, i, 1);
                double py = state.getStarPara(k, i, 2);
                double shape_size = state.getStarPara(k, i, 7);
                double shape_e1 = state.getStarPara(k, i, 8);
                double shape_e2 = state.getStarPara(k, i, 9);
                std::size_t star_offset = static_cast<std::size_t>(i) * ns * ns;
                if (!isFinitePSFStar(star, star_offset, ns * ns, px, py,
                                     shape_size, shape_e1, shape_e2)) {
                    removed_non_finite++;
                    continue;
                }
                nums++;
                posi.push_back({px, py});
                sshape.push_back({shape_size, shape_e1, shape_e2});
                for (int idx = 0; idx < ns * ns; ++idx) {
                    star_local.push_back(star[star_offset + idx]);
                }
            }
            MainIO::OutputFile file20;
            std::string ccd_id;
            if (LensingConfig::PSF_Ms == 1) {
                int chip_index = UniversalUtils::getChipId(imageFiles[k]);
                ccd_id = std::to_string(chip_index);
                std::string star_xy_filename = OutputLayout::chipPath(
                    dirOutput, "stamps/dat_StarXY", prefix_e + "_" + ccd_id,
                    "_star_xy.dat");
                file20.open(star_xy_filename);
                file20 << std::setprecision(10);
            }

            std::vector<double> PSF_coe_l;
            LinearSolve::SolveDiagnostics fit_diagnostics;
            LinearSolve::SolveStatus fit_status = LinearSolve::SolveStatus::FailedRankDeficient;
            if (nums >= LensingConfig::nstar_min_local) {
                fit_status = itpNormPSF(
                    nums, star_local, posi, ns, npl,
                    LensingConfig::chipnx, LensingConfig::chipny,
                    PSF_coe_l, &fit_diagnostics);
            }

            if (nums >= LensingConfig::nstar_min_local &&
                fit_status == LinearSolve::SolveStatus::Normal) {

                file90 << (k + 1) << " " << nums << " 1\n";

                std::vector<float> poly_cochi2(nums);
                float poly_ave = 0.0f, poly_std = 0.0f;

                if (LensingConfig::PSF_Ms == 1) {
                    file20 << nums << " 1\n";
                }

                std::vector<float> psf_residual(static_cast<size_t>(nstar_max) * ns * ns, 0.0f);
                for (int i = 0; i < nums; ++i) {
                    double xx = 2.0 * (posi[i][0] / static_cast<double>(LensingConfig::chipnx)) - 1.0;
                    double yy = 2.0 * (posi[i][1] / static_cast<double>(LensingConfig::chipny)) - 1.0;
                    std::vector<float> model, model0;
                    getPSFModel(ns, npl, PSF_coe_l, xx, yy, model, model0);
                    ExStar::anaChi2Simple(ns, model.data(), model0.data(), poly_cochi2[i]);

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, model, ee, size, 0.02f);

                    double msshape_size = size;
                    double msshape_e1 = ee[0];
                    double msshape_e2 = ee[1];

                    float px = static_cast<float>(posi[i][0]);
                    float py = static_cast<float>(posi[i][1]);

                    file90 << px << " " << py << " "
                           << sshape[i][0] << " " << sshape[i][1] << " " << sshape[i][2] << " "
                           << msshape_size << " " << msshape_e1 << " " << msshape_e2 << "\n";

                    if (LensingConfig::PSF_Ms == 1) {
                        if (size < 0.1 || std::isnan(model[0])) {
                            file20 << "-1 -1\n";
                        } else {
                            file20 << px << " " << py << "\n";
                        }

                        std::vector<float> temp_res(ns * ns);
                        for (int idx = 0; idx < ns * ns; ++idx) {
                            temp_res[idx] = star_local[static_cast<size_t>(i) * ns * ns + idx] - model[idx];
                        }
                        PSF_rescale(temp_res, res_factor);
                        for (int idx = 0; idx < ns * ns; ++idx) {
                            psf_residual[static_cast<size_t>(i) * ns * ns + idx] = temp_res[idx];
                        }
                    }
                }

                ExStar::getArrayAveStd(poly_cochi2, poly_ave, poly_std);
                file10 << nums << " 1 " << poly_ave << " " << poly_std << "\n";
                for (int i = 0; i < ns; ++i) {
                    for (int j = 0; j < ns; ++j) {
                        for (int u = 0; u < npl + 1; ++u) {
                            file10 << PSF_coe_l[(j * ns + i) * (npl + 1) + u] << (u == npl ? "" : " ");
                        }
                        file10 << "\n";
                    }
                }

                if (LensingConfig::PSF_Ms == 1) {
                    int nn1 = ns * len_s;
                    int nn2 = ns * ((nums / len_s) + 1);
                    std::string resi_fits = OutputLayout::chipPath(
                        dirOutput, "stamps/fits_PsfResi", prefix_e + "_" + ccd_id,
                        "_psf_p_resi.fits");
                    FitsIO::writeStamps(nstar_max, 1, nums, ns, ns, psf_residual, nn1, nn2, resi_fits);
                }
            } else {
                if (nums < LensingConfig::nstar_min_local) {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", LinearSolve::SolveStatus::FailedRankDeficient,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " valid_samples=" + std::to_string(nums) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " required=" + std::to_string(LensingConfig::nstar_min_local) +
                            " action=MARK_CHIP_INVALID");
                } else {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", fit_status,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " " + LinearSolve::diagnosticsContext(fit_diagnostics) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " action=MARK_CHIP_INVALID");
                }
                file10 << nums << " -1 -1 -1\n";
                file90 << (k + 1) << " " << nums << " -1\n";
                if (LensingConfig::PSF_Ms == 1) {
                    file20 << nums << " -1\n";
                }
            }
            file10.close();
            if (LensingConfig::PSF_Ms == 1) {
                file20.close();
            }
        }
        file90.close();
    }

    // ==========================================
    // Function: Fit and serialize hybrid PSF models.
    // Method: Preserve F77 model layout with 17-digit double serialization.
    // ==========================================
    void makePSFHybrid(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int len_s = LensingConfig::len_s;
        int nstar_max = LensingConfig::nstar_max;
        int npl = LensingConfig::npl;
        int nplx = LensingConfig::nplx;
        int npx = LensingConfig::npx;
        int npy = LensingConfig::npy;
        int step_psf = LensingConfig::step_psf;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string comp_filename = dirOutput + "/stamps/dat_StarComp/" + prefix_e + "_star_comp_expo.dat";
        MainIO::OutputFile file90(comp_filename);
        if (!file90.is_open()) {
            std::cerr << "makePSFHybrid: Error opening " << comp_filename << std::endl;
            return;
        }
        file90 << std::setprecision(17);

        for (int k = 0; k < nchip; ++k) {
            int nums = 0;
            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);

            std::vector<float> star(static_cast<size_t>(nstar_max) * ns * ns, 0.0f);
            int nn1 = ns * len_s;
            int nn2 = ns * ((state.nstar[k] / len_s) + 1);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
            if (!FitsIO::readStamps(nstar_max, 1, state.nstar[k], ns, ns, star, nn1, nn2, filepath)) {
                std::cerr << "makePSFHybrid: readStamps failed for " << filepath << std::endl;
                std::exit(1);
            }

            std::string coe_filename = OutputLayout::chipPath(
                dirOutput, "stamps/dat_PsfFit", prefix, "_PSF_coe_local.dat");
            MainIO::OutputFile file10(coe_filename);
            if (!file10.is_open()) {
                std::cerr << "makePSFHybrid: Error opening " << coe_filename << std::endl;
                continue;
            }
            file10 << std::setprecision(17);

            std::vector<std::array<double, 2>> posi;
            std::vector<std::array<double, 3>> sshape;
            std::vector<float> star_local;
            int removed_non_finite = 0;

            for (int i = 0; i < state.nstar[k]; ++i) {
                if (state.getStarPara(k, i, 4) < 0.0) continue;
                double px = state.getStarPara(k, i, 2 - 1);
                double py = state.getStarPara(k, i, 3 - 1);
                double shape_size = state.getStarPara(k, i, 8 - 1);
                double shape_e1 = state.getStarPara(k, i, 9 - 1);
                double shape_e2 = state.getStarPara(k, i, 10 - 1);
                std::size_t star_offset = static_cast<std::size_t>(i) * ns * ns;
                if (!isFinitePSFStar(star, star_offset, ns * ns, px, py,
                                     shape_size, shape_e1, shape_e2)) {
                    removed_non_finite++;
                    continue;
                }
                nums++;
                posi.push_back({px, py}); // F77 index 2 & 3
                sshape.push_back({shape_size, shape_e1, shape_e2}); // F77 index 8, 9, 10
                for (int idx = 0; idx < ns * ns; ++idx) {
                    star_local.push_back(star[star_offset + idx]);
                }
            }
            std::vector<float> psfmap(static_cast<size_t>(npx) * npy, 0.0f);

            std::vector<double> PSF_coe_l;
            LinearSolve::SolveDiagnostics fit_diagnostics;
            LinearSolve::SolveStatus fit_status = LinearSolve::SolveStatus::FailedRankDeficient;
            if (nums >= LensingConfig::nstar_min_local) {
                fit_status = interpolatePSF(
                    nums, star_local, posi, ns, npl, PSF_coe_l, &fit_diagnostics);
            }

            if (nums >= LensingConfig::nstar_min_local &&
                fit_status == LinearSolve::SolveStatus::Normal) {

                file10 << nums << " 1\n";
                for (int i = 0; i < ns; ++i) {
                    for (int j = 0; j < ns; ++j) {
                        for (int u = 0; u < npl + 1; ++u) {
                            file10 << PSF_coe_l[(j * ns + i) * (npl + 1) + u] << (u == npl ? "" : " ");
                        }
                        file10 << "\n";
                    }
                }

                // Residual computation
                std::vector<float> star_residual = star_local;
                for (int i = 0; i < nums; ++i) {
                    double xx = posi[i][0];
                    double yy = posi[i][1];
                    std::vector<float> model, model0;
                    getPSFModel(ns, npl, PSF_coe_l, xx, yy, model, model0);
                    for (int idx = 0; idx < ns * ns; ++idx) {
                        star_residual[static_cast<size_t>(i) * ns * ns + idx] -= model[idx];
                    }
                }

                int nx = 0, ny = 0;
                if (!FitsIO::readPara(imageFiles[k], nx, ny)) {
                    std::cerr << "makePSFHybrid: readPara failed for " << imageFiles[k] << std::endl;
                    std::exit(1);
                }
                genPSFFits(psfmap, nums, nx, ny, star_residual, posi);

                int ixt = nx / step_psf + 1;
                int iyt = ny / step_psf + 1;

                for (int ix = 1; ix <= ixt; ++ix) {
                    int xs = (ix - 1) * step_psf;
                    double xx = (ix - 0.5) * step_psf;
                    for (int iy = 1; iy <= iyt; ++iy) {
                        int ys = (iy - 1) * step_psf;
                        double yy = (iy - 0.5) * step_psf;
                        std::vector<float> model, model0;
                        getPSFModel(ns, npl, PSF_coe_l, xx, yy, model, model0);
                        for (int v = 0; v < ns; ++v) {
                            for (int u = 0; u < ns; ++u) {
                                psfmap[(ys + v) * npx + (xs + u)] += model[v * ns + u];
                            }
                        }
                    }
                }

                nx = (nx / step_psf + 1) * step_psf;
                ny = (ny / step_psf + 1) * step_psf;

                std::string out_fits = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_PsfLocal", prefix, "_PSF_local.fits");
                std::vector<float> sub_psfmap(static_cast<size_t>(nx) * ny);
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        sub_psfmap[y * nx + x] = psfmap[y * npx + x];
                    }
                }
                FitsIO::writeImage(out_fits, nx, ny, sub_psfmap);

                file90 << (k + 1) << " " << nums << " 1\n";
                for (int i = 0; i < nums; ++i) {
                    float px = static_cast<float>(posi[i][0]);
                    float py = static_cast<float>(posi[i][1]);
                    std::vector<float> model;
                    double dmax = 0.0;
                    getPSFModelVeryLocal(psfmap, px, py, model, dmax);

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, model, ee, size, 0.02f);

                    double mshape_size = size;
                    double mshape_e1 = ee[0];
                    double mshape_e2 = ee[1];

                    file90 << px << " " << py << " "
                           << sshape[i][0] << " " << sshape[i][1] << " " << sshape[i][2] << " "
                           << mshape_size << " " << mshape_e1 << " " << mshape_e2 << " " << dmax << "\n";
                }
            } else {
                if (nums < LensingConfig::nstar_min_local) {
                    LinearSolve::reportFailure(
                        "PSFModel::interpolatePSF", LinearSolve::SolveStatus::FailedRankDeficient,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " valid_samples=" + std::to_string(nums) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " required=" + std::to_string(LensingConfig::nstar_min_local) +
                            " action=MARK_CHIP_INVALID");
                } else {
                    LinearSolve::reportFailure(
                        "PSFModel::interpolatePSF", fit_status,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " " + LinearSolve::diagnosticsContext(fit_diagnostics) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " action=MARK_CHIP_INVALID");
                }
                for (int i = 0; i < npx; ++i) {
                    for (int j = 0; j < npy; ++j) {
                        psfmap[j * npx + i] = 0.0f;
                    }
                }
                psfmap[(step_psf - 1) * npx + (step_psf - 1)] = -100.0f;
                int nx = 3 * step_psf;
                int ny = 3 * step_psf;

                std::string out_fits = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_PsfLocal", prefix, "_PSF_local.fits");
                std::vector<float> sub_psfmap(static_cast<size_t>(nx) * ny);
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        sub_psfmap[y * nx + x] = psfmap[y * npx + x];
                    }
                }
                FitsIO::writeImage(out_fits, nx, ny, sub_psfmap);

                file10 << nums << " -1\n";
                file90 << (k + 1) << " " << nums << " -1\n";
            }
            file10.close();
        }
        file90.close();
    }

    // Mathematical Interpolation Helpers

    // ==========================================
    // Function: Fit all PSF-frequency pixels with a shared spatial design
    // Method: Factor the constant and polynomial designs once, then reuse both QR objects for every RHS pixel.
    // ==========================================
    static LinearSolve::SolveStatus fitPSFCoefficients(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, bool normalize_positions,
        std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics) {
        const int pixel_count = ns * ns;
        PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);

        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = nsam;
        diag.cols = npp;
        diag.required_rank = npp;

        if (nsam < npp || ns <= 0 || npp <= 0 ||
            static_cast<int>(posi.size()) < nsam ||
            image.size() < static_cast<size_t>(nsam) * pixel_count ||
            (normalize_positions && (nx <= 0 || ny <= 0))) {
            diag.rank = std::min(nsam, npp);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd polynomial_design(nsam, npp);
        Eigen::MatrixXd constant_design = Eigen::MatrixXd::Ones(nsam, 1);
        for (int sample = 0; sample < nsam; ++sample) {
            double xx = posi[sample][0];
            double yy = posi[sample][1];
            if (normalize_positions) {
                xx = 2.0 * (xx / static_cast<double>(nx)) - 1.0;
                yy = 2.0 * (yy / static_cast<double>(ny)) - 1.0;
            }
            if (!std::isfinite(xx) || !std::isfinite(yy)) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
            for (int term = 0; term < npp; ++term) {
                polynomial_design(sample, term) = UniversalUtils::fitFunc2(xx, yy, term);
            }
        }
        if (!polynomial_design.allFinite()) {
            return LinearSolve::SolveStatus::FailedSolver;
        }
        for (std::size_t idx = 0; idx < static_cast<std::size_t>(nsam) * pixel_count; ++idx) {
            if (!std::isfinite(image[idx])) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
        }

        LinearSolve::LeastSquaresQR constant_solver;
        LinearSolve::SolveDiagnostics constant_diagnostics;
        LinearSolve::SolveStatus status = constant_solver.factorize(constant_design, constant_diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            diag = constant_diagnostics;
            return status;
        }

        LinearSolve::LeastSquaresQR polynomial_solver;
        status = polynomial_solver.factorize(polynomial_design, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        Eigen::VectorXd rhs(nsam);
        Eigen::VectorXd constant_solution;
        Eigen::VectorXd polynomial_solution;
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                const int pixel = j * ns + i;
                for (int sample = 0; sample < nsam; ++sample) {
                    rhs(sample) = image[static_cast<size_t>(sample) * pixel_count + pixel];
                }

                status = constant_solver.solve(rhs, constant_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }
                status = polynomial_solver.solve(rhs, polynomial_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }

                for (int term = 0; term < npp; ++term) {
                    PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + term] = polynomial_solution(term);
                }
                PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + npp] = constant_solution(0);
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // ==========================================
    // Function: Fit normalized local PSF coefficients
    // Method: Normalize chip coordinates to [-1,1] and delegate to the shared-design batch fitter.
    // ==========================================
    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics) {
        return fitPSFCoefficients(
            nsam, image, posi, ns, npp, nx, ny, true, PSF_coe, diagnostics);
    }

    // ==========================================
    // Function: Fit hybrid PSF coefficients in native coordinates
    // Method: Preserve the native-coordinate polynomial model and reuse one design factorization for all pixels.
    // ==========================================
    LinearSolve::SolveStatus interpolatePSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics) {
        return fitPSFCoefficients(
            nsam, image, posi, ns, npp, 1, 1, false, PSF_coe, diagnostics);
    }

    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0) {
        modelp.assign(static_cast<size_t>(ns) * ns, 0.0f);
        model0.assign(static_cast<size_t>(ns) * ns, 0.0f);
        std::vector<double> coep(npp);
        std::vector<double> coe0(1);
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                for (int k = 0; k < npp; ++k) {
                    coep[k] = PSF_coe[(j * ns + i) * (npp + 1) + k];
                }
                coe0[0] = PSF_coe[(j * ns + i) * (npp + 1) + npp];
                model0[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, 1, coe0));
                modelp[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, npp, coep));
                if (std::isnan(modelp[j * ns + i])) {
                    modelp[0] = modelp[j * ns + i];
                    return;
                }
            }
        }
    }

    void getPSFModelVeryLocal(const std::vector<float>& psfmap, double x, double y, std::vector<float>& model, double& dmax, int stride) {
        int ns = LensingConfig::ns;
        int step_psf = LensingConfig::step_psf;

        int ixs = static_cast<int>(x / step_psf) * step_psf;
        int iys = static_cast<int>(y / step_psf) * step_psf;

        model.assign(static_cast<size_t>(ns) * ns, 0.0f);
        for (int v = 0; v < ns; ++v) {
            for (int u = 0; u < ns; ++u) {
                model[v * ns + u] = psfmap[(iys + v) * stride + (ixs + u)];
            }
        }
        dmax = psfmap[iys * stride + (ixs + ns)];
    }

    void genPSFFits(std::vector<float>& psfmap, int nums, int nx, int ny, const std::vector<float>& star, const std::vector<std::array<double, 2>>& posi) {
        int npx = LensingConfig::npx;
        int npy = LensingConfig::npy;
        int ns = LensingConfig::ns;
        int step_psf = LensingConfig::step_psf;

        psfmap.assign(static_cast<size_t>(npx) * npy, 0.0f);

        int ixt = nx / step_psf + 1;
        int iyt = ny / step_psf + 1;

        if (nums < LensingConfig::nstar_min_local) {
            psfmap[(step_psf - 1) * npx + (step_psf - 1)] = -100.0f;
            psfmap[(step_psf - 2) * npx + (step_psf - 2)] = static_cast<float>(nums);
            return;
        }

        psfmap[(step_psf - 1) * npx + (step_psf - 1)] = 1.0f;
        psfmap[(step_psf - 2) * npx + (step_psf - 2)] = static_cast<float>(nums);

        double d2min = std::pow(step_psf / 4.0, 2); // 625.0

        for (int ix = 1; ix <= ixt; ++ix) {
            int xs = (ix - 1) * step_psf;
            double x = (ix - 0.5) * step_psf;
            for (int iy = 1; iy <= iyt; ++iy) {
                int ys = (iy - 1) * step_psf;
                double y = (iy - 0.5) * step_psf;

                std::vector<float> d2(nums);
                for (int i = 0; i < nums; ++i) {
                    d2[i] = static_cast<float>(std::pow(x - posi[i][0], 2) + std::pow(y - posi[i][1], 2));
                }

                std::vector<int> indx;
                NumericalRecipes::indexx(d2, indx);

                double dmax = std::sqrt(d2[indx[LensingConfig::nstar_min_local - 1]]);
                double norm = 0.0;
                std::vector<float> model(ns * ns, 0.0f);
                std::vector<float> weight(LensingConfig::n_neighbor);

                for (int i = 0; i < LensingConfig::n_neighbor; ++i) {
                    int idx = indx[i];
                    weight[i] = 1.0f / std::max(static_cast<float>(d2min), d2[idx]);
                    norm += weight[i];
                    for (int v = 0; v < ns; ++v) {
                        for (int u = 0; u < ns; ++u) {
                            model[v * ns + u] += weight[i] * star[static_cast<size_t>(idx) * ns * ns + v * ns + u];
                        }
                    }
                }

                for (int v = 0; v < ns; ++v) {
                    for (int u = 0; u < ns; ++u) {
                        psfmap[(ys + v) * npx + (xs + u)] = model[v * ns + u] / static_cast<float>(norm);
                    }
                }
                psfmap[ys * npx + (xs + ns)] = static_cast<float>(dmax);
            }
        }
    }

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;
                        }
                    }
                }
            }
            area0 = tempi;
        }

        area = (area_cnt - 1) / 2;
    }

    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        e[0] = 0.0;
        e[1] = 0.0;
        double norm = 0.0;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;

                            double kx = u - cx;
                            double ky = v - cy;
                            double p_val = power[idx];
                            e[0] += p_val * (kx * kx - ky * ky);
                            e[1] += p_val * 2.0 * kx * ky;
                            norm += p_val * (kx * kx + ky * ky);
                        }
                    }
                }
            }
            area0 = tempi;
        }

        e[0] /= norm;
        e[1] /= norm;
    }

    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio) {
        int area = 0;
        getPowerArea(nx, ny, power, area, thresh_ratio);
        size = area;
        getPowerE(nx, ny, power, e, thresh_ratio);
    }

    void getPSFFWHM(const std::vector<float>& power, double& FWHM) {
        int ns = LensingConfig::ns;
        float thresh = power[(ns / 2) * ns + (ns / 2)] * std::exp(-1.0f);
        double area = -1e-5;
        for (int idx = 0; idx < ns * ns; ++idx) {
            if (power[idx] >= thresh) {
                area += 1.0;
            }
        }
        if (area <= 0.0) {
            FWHM = 0.0;
            return;
        }
        double beta = ns / (2.0 * LensingConfig::pi) / std::sqrt(area / LensingConfig::pi);
        FWHM = beta * 2.0 * std::sqrt(2.0 * std::log(2.0)) * LensingConfig::pixel_size;
    }
}
