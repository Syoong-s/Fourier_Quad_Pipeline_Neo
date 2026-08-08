#include "PSFRecons.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "PSFModel.hpp"
#include "FitsIO.hpp"
#include "UniversalUtils.hpp"
#include "MPIScheduler.hpp"
#include "NumericalRecipes.hpp"
#include "LinearSolve.hpp"
#include <Eigen/Dense>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>


// Extern exposures defined in main
extern std::vector<std::string> EXPO_FILE;

namespace PSFRecons {

    // ==========================================
    // Function: Validate one residual stamp
    // Method: Require a complete finite ns*ns sample before adding it to PCA statistics or projections.
    // ==========================================
    static bool isFiniteResidualSample(const std::vector<float>& residual,
                                       std::size_t offset, int pixel_count) {
        if (offset + static_cast<std::size_t>(pixel_count) > residual.size()) {
            return false;
        }
        for (int idx = 0; idx < pixel_count; ++idx) {
            if (!std::isfinite(residual[offset + idx])) {
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Fit all PCA-mode coefficient surfaces for one CCD block
    // Method: Check the shared design rank once, invert the F77 single-precision normal matrix once, and reuse it for every active mode.
    // ==========================================
    static LinearSolve::SolveStatus fitPcaBlockCoefficients(
        const std::vector<float>& x, const std::vector<float>& y,
        const std::vector<int>& sample_indices, const std::vector<double>& projected_coefficients,
        int effective_pcs, int nc, std::vector<float>& fitted_coefficients,
        LinearSolve::SolveDiagnostics& diagnostics) {
        const int nsam = static_cast<int>(x.size());
        fitted_coefficients.assign(
            static_cast<size_t>(LensingConfig::n_pcs) * nc, 0.0f);
        diagnostics = {};
        diagnostics.rows = nsam;
        diagnostics.cols = nc;
        diagnostics.required_rank = nc;

        if (nsam < nc || nc <= 0 || y.size() != x.size() ||
            sample_indices.size() != x.size() || effective_pcs < 0 ||
            effective_pcs > LensingConfig::n_pcs) {
            diagnostics.rank = std::min(nsam, nc);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd design(nsam, nc);
        std::vector<float> basis(static_cast<size_t>(nsam) * nc, 0.0f);
        for (int sample = 0; sample < nsam; ++sample) {
            if (!std::isfinite(x[sample]) || !std::isfinite(y[sample])) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
            for (int term = 0; term < nc; ++term) {
                double value = UniversalUtils::fitFunc2(x[sample], y[sample], term);
                if (!std::isfinite(value)) {
                    return LinearSolve::SolveStatus::FailedSolver;
                }
                design(sample, term) = value;
                basis[static_cast<size_t>(sample) * nc + term] = static_cast<float>(value);
            }
        }

        LinearSolve::LeastSquaresQR rank_solver;
        LinearSolve::SolveStatus status = rank_solver.factorize(design, diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        std::vector<float> normal(static_cast<size_t>(nc) * nc, 0.0f);
        for (int sample = 0; sample < nsam; ++sample) {
            const float* sample_basis = basis.data() + static_cast<size_t>(sample) * nc;
            for (int u = 0; u < nc; ++u) {
                for (int v = 0; v < nc; ++v) {
                    normal[static_cast<size_t>(u) * nc + v] += sample_basis[u] * sample_basis[v];
                }
            }
        }
        for (float value : normal) {
            if (!std::isfinite(value)) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
        }

        std::vector<float> inverse;
        if (!UniversalUtils::invertMatrixF77(normal, nc, inverse, false)) {
            return LinearSolve::SolveStatus::FailedSolver;
        }
        for (float value : inverse) {
            if (!std::isfinite(value)) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
        }

        for (int mode = 0; mode < effective_pcs; ++mode) {
            std::vector<float> rhs(nc, 0.0f);
            for (int sample = 0; sample < nsam; ++sample) {
                int source_index = sample_indices[sample];
                if (source_index < 0 ||
                    static_cast<size_t>(source_index) * LensingConfig::n_pcs + mode >=
                        projected_coefficients.size()) {
                    return LinearSolve::SolveStatus::FailedSolver;
                }
                double coefficient = projected_coefficients[
                    static_cast<size_t>(source_index) * LensingConfig::n_pcs + mode];
                if (!std::isfinite(coefficient)) {
                    return LinearSolve::SolveStatus::FailedSolver;
                }
                float coefficient_f = static_cast<float>(coefficient);
                if (!std::isfinite(coefficient_f)) {
                    return LinearSolve::SolveStatus::FailedSolver;
                }
                const float* sample_basis = basis.data() + static_cast<size_t>(sample) * nc;
                for (int term = 0; term < nc; ++term) {
                    rhs[term] += coefficient_f * sample_basis[term];
                }
            }

            for (int output_term = 0; output_term < nc; ++output_term) {
                float value = 0.0f;
                for (int input_term = 0; input_term < nc; ++input_term) {
                    value += inverse[static_cast<size_t>(output_term) * nc + input_term] *
                             rhs[input_term];
                }
                if (!std::isfinite(value)) {
                    return LinearSolve::SolveStatus::FailedSolver;
                }
                fitted_coefficients[static_cast<size_t>(mode) * nc + output_term] = value;
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // Static helper function to accumulate covariance and mean block by block
    static void accumulateBlock(int n, int m, const std::vector<double>& blk, std::vector<double>& sum_x, std::vector<double>& sum_xxt) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                sum_x[j] += blk[static_cast<size_t>(i) * m + j];
            }
        }
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < m; ++k) {
                double value = 0.0;
                for (int i = 0; i < n; ++i) {
                    value += blk[static_cast<size_t>(i) * m + j] * blk[static_cast<size_t>(i) * m + k];
                }
                sum_xxt[static_cast<size_t>(j) * m + k] += value;
            }
        }
    }

    // Stage 6 main entry: coordinates PSF fitting and reconstruction across chips and exposures
    void chipPSFRecons(int nexpo) {
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(EXPO_FILE[0], image_files, dir_output);

        // Call forcecov to run PCA fitting on CCDs in parallel
        MPIScheduler::forcecov(
            LensingConfig::procs_pn,
            LensingConfig::work_pn,
            LensingConfig::NMAX_CHIP,
            [](int ichip, int nexpo_inner) {
                chipResPCAFit(ichip, nexpo_inner);
            },
            "fitting residual...",
            nexpo
        );

        if (MPIScheduler::my_id == 0) {
            std::cout << "PSF PCA fitting completed for all chips." << std::endl;
        }

        MPIScheduler::barrier();
        PSFModel::initAndLoadAllPSF(dir_output, MPIScheduler::my_id);
        MPIScheduler::barrier();

        // Map modified residuals distributed across exposures
        MPIScheduler::distribute(
            nexpo,
            [](int iexpo) {
                plotResidualsV2(iexpo);
            },
            "Mapping Modified Residuals..."
        );
        MPIScheduler::barrier();
    }

    // ==========================================
    // Function: Fit residual PCA and its spatial coefficient surfaces for one CCD
    // Method: Filter non-finite stars, retain only positive PCA modes, and route numerical failures to polynomial-only output.
    // ==========================================
    void chipResPCAFit(int ichip, int nexpo) {
        if (ichip == 2 || ichip == 61) {
            return;
        }

        int ns = LensingConfig::ns;
        int nsns = LensingConfig::nsns;
        int block_size = 2000;

        std::vector<double> block_dble(static_cast<size_t>(block_size) * nsns, 0.0);
        std::vector<double> cov_arr(static_cast<size_t>(nsns) * nsns, 0.0);
        std::vector<double> mean_arr(nsns, 0.0);
        int buf_cnt = 0;
        int ntot = 0;
        int removed_non_finite = 0;

        // 1. Read star residuals block by block and compute mean and covariance
        for (int i = 1; i <= nexpo; ++i) {
            std::vector<std::string> image_files;
            std::string dir_out;
            UniversalUtils::getImageList(EXPO_FILE[i - 1], image_files, dir_out);
            std::string prefix_e = UniversalUtils::getPrefixExpo(image_files[0]);
            
            std::string filename_xy = OutputLayout::chipPath(
                dir_out, "stamps/dat_StarXY",
                prefix_e + "_" + std::to_string(ichip), "_star_xy.dat");
            std::ifstream xy_file(filename_xy);
            if (!xy_file.is_open()) {
                continue;
            }

            int nstar_file = 0, valid_file = 0;
            if (xy_file >> nstar_file >> valid_file) {
                if (nstar_file > 0 && valid_file >= 0) {
                    std::vector<float> psf_residual;
                    int nn1 = ns * LensingConfig::len_s;
                    int nn2 = ns * ((nstar_file / LensingConfig::len_s) + 1);
                    std::string fits_filename = OutputLayout::chipPath(
                        dir_out, "stamps/fits_PsfResi",
                        prefix_e + "_" + std::to_string(ichip), "_psf_p_resi.fits");
                    
                    if (FitsIO::readStamps(nstar_file, 1, nstar_file, ns, ns, psf_residual, nn1, nn2, fits_filename)) {
                        double px = 0.0, py = 0.0;
                        int valid_num = 0;
                        while (xy_file >> px >> py) {
                            valid_num++;
                            if (valid_num > nstar_file) break;
                            size_t src_offset = static_cast<size_t>(valid_num - 1) * nsns;
                            if (!std::isfinite(px) || !std::isfinite(py) ||
                                !isFiniteResidualSample(psf_residual, src_offset, nsns)) {
                                removed_non_finite++;
                                continue;
                            }
                            if (px < 0.0 || py < 0.0) continue;

                            ntot++;
                            buf_cnt++;

                            size_t dest_offset = static_cast<size_t>(buf_cnt - 1) * nsns;
                            for (int row = 0; row < ns; ++row) {
                                for (int col = 0; col < ns; ++col) {
                                    block_dble[dest_offset + col * ns + row] = psf_residual[src_offset + row * ns + col];
                                }
                            }

                            if (buf_cnt == block_size) {
                                accumulateBlock(block_size, nsns, block_dble, mean_arr, cov_arr);
                                buf_cnt = 0;
                            }
                        }
                    } else {
                        std::exit(1);
                    }
                }
            }
            xy_file.close();
        }

        if (buf_cnt > 0) {
            accumulateBlock(buf_cnt, nsns, block_dble, mean_arr, cov_arr);
        }

        std::cout << "chip " << ichip << " total stars: " << ntot << std::endl;

        std::ostringstream ss_ccd;
        ss_ccd << std::setfill('0') << std::setw(2) << ichip;
        std::string c_chip_2digit = ss_ccd.str();

        std::vector<std::string> dummy_image_files;
        std::string dirOutput;
        UniversalUtils::getImageList(EXPO_FILE[0], dummy_image_files, dirOutput);

        std::vector<double> components(static_cast<size_t>(nsns) * LensingConfig::n_pcs, 0.0);
        bool pca_failed = false;
        int effective_pcs = 0;

        if (ntot < 2) {
            pca_failed = true;
            components[0] = -1.0e30;
            LinearSolve::reportFailure(
                "PSFRecons::chipResPCAFit", LinearSolve::SolveStatus::FailedRankDeficient,
                "ccd=" + std::to_string(ichip) + " valid_samples=" + std::to_string(ntot) +
                    " required=2 removed_samples=" + std::to_string(removed_non_finite) +
                    " action=POLYNOMIAL_ONLY");
        } else {
            for (int j = 0; j < nsns; ++j) {
                mean_arr[j] /= static_cast<double>(ntot);
            }
            for (int j = 0; j < nsns; ++j) {
                for (int k = 0; k < nsns; ++k) {
                    cov_arr[static_cast<size_t>(j) * nsns + k] =
                        (cov_arr[static_cast<size_t>(j) * nsns + k]
                         - static_cast<double>(ntot) * mean_arr[j] * mean_arr[k])
                        / static_cast<double>(ntot - 1);
                }
            }

            for (int j = 0; j < nsns; ++j) {
                for (int k = j + 1; k < nsns; ++k) {
                    double symmetric_value = 0.5 * (
                        cov_arr[static_cast<size_t>(j) * nsns + k] +
                        cov_arr[static_cast<size_t>(k) * nsns + j]);
                    cov_arr[static_cast<size_t>(j) * nsns + k] = symmetric_value;
                    cov_arr[static_cast<size_t>(k) * nsns + j] = symmetric_value;
                }
            }

            using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic,
                                                   Eigen::Dynamic, Eigen::RowMajor>;
            const Eigen::Map<const RowMajorMatrixXd> cov_mat(
                cov_arr.data(), nsns, nsns);

            if (!cov_mat.allFinite()) {
                pca_failed = true;
                components[0] = -1.0e30;
                LinearSolve::reportFailure(
                    "PSFRecons::chipResPCAFit", LinearSolve::SolveStatus::FailedSolver,
                    "ccd=" + std::to_string(ichip) + " action=POLYNOMIAL_ONLY");
            } else {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;
                eig_solver.compute(cov_mat, Eigen::ComputeEigenvectors);

                if (eig_solver.info() != Eigen::Success) {
                    pca_failed = true;
                    components[0] = -1.0e30;
                    LinearSolve::reportFailure(
                        "PSFRecons::chipResPCAFit", LinearSolve::SolveStatus::FailedSolver,
                        "ccd=" + std::to_string(ichip) + " action=POLYNOMIAL_ONLY");
                } else if (!eig_solver.eigenvalues().allFinite() ||
                           !eig_solver.eigenvectors().allFinite()) {
                    pca_failed = true;
                    components[0] = -1.0e30;
                    LinearSolve::reportFailure(
                        "PSFRecons::chipResPCAFit", LinearSolve::SolveStatus::FailedSolver,
                        "ccd=" + std::to_string(ichip) + " action=POLYNOMIAL_ONLY");
                } else {
                    const Eigen::VectorXd& eigenvalues = eig_solver.eigenvalues();
                    LinearSolve::EigenSpectrumDiagnostics spectrum_diagnostics;
                    LinearSolve::SolveStatus spectrum_status =
                        LinearSolve::analyzeCovarianceSpectrum(
                            eigenvalues, ntot, LensingConfig::n_pcs,
                            LensingConfig::pca_negative_eigenvalue_threshold,
                            spectrum_diagnostics);
                    if (spectrum_status != LinearSolve::SolveStatus::Normal) {
                        pca_failed = true;
                        components[0] = -1.0e30;
                        std::ostringstream context;
                        context << std::scientific << std::setprecision(17)
                                << "ccd=" << ichip
                                << " lambda_min=" << spectrum_diagnostics.lambda_min
                                << " threshold=" << LensingConfig::pca_negative_eigenvalue_threshold
                                << " action=POLYNOMIAL_ONLY";
                        LinearSolve::reportFailure(
                            "PSFRecons::chipResPCAFit",
                            spectrum_status,
                            context.str());
                    } else {
                        effective_pcs = spectrum_diagnostics.effective_modes;

                        // Eigen returns ascending eigenpairs. Retain only the numerically positive top subspace.
                        const Eigen::MatrixXd& eigvec = eig_solver.eigenvectors();
                        for (int j = 0; j < effective_pcs; ++j) {
                            const int eig_col = nsns - 1 - j;
                            for (int k = 0; k < nsns; ++k) {
                                components[static_cast<size_t>(k) * LensingConfig::n_pcs + j] =
                                    eigvec(k, eig_col);
                            }
                        }
                    }
                }
            }
        }

        // Write PCA results to disk
        std::string filename_pcs = dirOutput + "/stamps/dat_Pcs/pcs_ccd" + c_chip_2digit + ".dat";
        MainIO::OutputFile pcs_file(filename_pcs);
        if (pcs_file.is_open()) {
            for (int k = 0; k < nsns; ++k) {
                for (int j = 0; j < LensingConfig::n_pcs; ++j) {
                    pcs_file << std::scientific << std::setprecision(17)
                             << components[static_cast<size_t>(k) * LensingConfig::n_pcs + j] << " ";
                }
                pcs_file << std::scientific << std::setprecision(17) << mean_arr[k] << "\n";
            }
            pcs_file.close();
            std::cout << "PCA finished chip ..." << ichip << std::endl;
        } else {
            std::cerr << "Error writing PCA file: " << filename_pcs << std::endl;
        }

        // 2. Project PCA coefficients
        std::vector<double> coeff(static_cast<size_t>(ntot) * LensingConfig::n_pcs, 0.0);
        std::vector<double> xc(ntot, 0.0);
        std::vector<double> yc(ntot, 0.0);

        if (!pca_failed) {
            int curr_tot = 0;
            for (int i = 1; i <= nexpo; ++i) {
                std::vector<std::string> image_files;
                std::string dir_out;
                UniversalUtils::getImageList(EXPO_FILE[i - 1], image_files, dir_out);
                std::string prefix_e = UniversalUtils::getPrefixExpo(image_files[0]);
                std::string filename_xy = OutputLayout::chipPath(
                    dir_out, "stamps/dat_StarXY",
                    prefix_e + "_" + std::to_string(ichip), "_star_xy.dat");

                std::ifstream xy_file(filename_xy);
                if (xy_file.is_open()) {
                    int nstar_file = 0, valid_file = 0;
                    if (xy_file >> nstar_file >> valid_file) {
                        if (nstar_file > 0 && valid_file >= 0) {
                            std::vector<float> psf_residual;
                            int nn1 = ns * LensingConfig::len_s;
                            int nn2 = ns * ((nstar_file / LensingConfig::len_s) + 1);
                            std::string fits_filename = OutputLayout::chipPath(
                                dir_out, "stamps/fits_PsfResi",
                                prefix_e + "_" + std::to_string(ichip),
                                "_psf_p_resi.fits");
                            
                            if (FitsIO::readStamps(nstar_file, 1, nstar_file, ns, ns, psf_residual, nn1, nn2, fits_filename)) {
                                double px = 0.0, py = 0.0;
                                int valid_num = 0;
                                while (xy_file >> px >> py) {
                                    valid_num++;
                                    if (valid_num > nstar_file) break;
                                    int src_offset = (valid_num - 1) * nsns;
                                    if (!std::isfinite(px) || !std::isfinite(py) ||
                                        !isFiniteResidualSample(
                                            psf_residual, static_cast<std::size_t>(src_offset), nsns)) {
                                        continue;
                                    }
                                    if (px < 0.0 || py < 0.0) continue;
                                    if (curr_tot >= ntot) break;

                                    xc[curr_tot] = px;
                                    yc[curr_tot] = py;

                                    std::vector<double> res_slice(nsns, 0.0);
                                    for (int row = 0; row < ns; ++row) {
                                        for (int col = 0; col < ns; ++col) {
                                            int idx = col * ns + row;
                                            res_slice[idx] = psf_residual[src_offset + row * ns + col] - mean_arr[idx];
                                        }
                                    }
                                    for (int j = 0; j < effective_pcs; ++j) {
                                        double value = 0.0;
                                        for (int idx = 0; idx < nsns; ++idx) {
                                            value += res_slice[idx] * components[static_cast<size_t>(idx) * LensingConfig::n_pcs + j];
                                        }
                                        coeff[static_cast<size_t>(curr_tot) * LensingConfig::n_pcs + j] = value;
                                    }
                                    curr_tot++;
                                }
                }
                    }
                    xy_file.close();
                }
            }
        }
        }

        // 3. Fit PCA coefficients as functions of position
        double nxc = static_cast<double>(LensingConfig::chipnx) / LensingConfig::nblocks;
        double nyc = static_cast<double>(LensingConfig::chipny) / LensingConfig::nblocks;

        for (int j = 1; j <= LensingConfig::nblocks; ++j) {
            double lxd = (j - 1.0) * nxc;
            double uxd = j * nxc;
            for (int k = 1; k <= LensingConfig::nblocks; ++k) {
                double lyd = (k - 1.0) * nyc;
                double uyd = k * nyc;

                int fit_num = 0;
                std::vector<int> idx_list;
                std::vector<float> xfit, yfit;

                if (!pca_failed) {
                    for (int ibstar = 0; ibstar < ntot; ++ibstar) {
                        if (xc[ibstar] < lxd || xc[ibstar] > uxd) continue;
                        if (yc[ibstar] < lyd || yc[ibstar] > uyd) continue;
                        fit_num++;
                        idx_list.push_back(ibstar);
                        xfit.push_back(static_cast<float>(2.0 * (xc[ibstar] - lxd) / (uxd - lxd) - 1.0));
                        yfit.push_back(static_cast<float>(2.0 * (yc[ibstar] - lyd) / (uyd - lyd) - 1.0));
                    }
                }

                std::string filename_coeff = dirOutput + "/stamps/dat_Pcs/coeff_ccd" + c_chip_2digit + "_" + std::to_string(j) + std::to_string(k) + ".dat";
                MainIO::OutputFile coeff_file(filename_coeff);
                if (!coeff_file.is_open()) {
                    std::cerr << "Error writing coeff file: " << filename_coeff << std::endl;
                    continue;
                }

                if (pca_failed || fit_num <= (LensingConfig::npp6th + 10)) {
                    if (!pca_failed) {
                        LinearSolve::reportFailure(
                            "PSFRecons::interpolate_6th",
                            LinearSolve::SolveStatus::FailedRankDeficient,
                            "ccd=" + std::to_string(ichip) +
                                " block=" + std::to_string(j) + std::to_string(k) +
                                " valid_samples=" + std::to_string(fit_num) +
                                " required=" + std::to_string(LensingConfig::npp6th + 11) +
                                " action=POLYNOMIAL_ONLY");
                    }
                    for (int i = 0; i < LensingConfig::npp6th; ++i) {
                        coeff_file << std::scientific << std::setprecision(10)
                                   << -1.0e30f << (i == LensingConfig::npp6th - 1 ? "" : " ");
                    }
                    coeff_file << "\n";
                    coeff_file.close();
                    continue;
                }

                std::vector<float> block_coefficients;
                LinearSolve::SolveDiagnostics fit_diagnostics;
                LinearSolve::SolveStatus fit_status = fitPcaBlockCoefficients(
                    xfit, yfit, idx_list, coeff, effective_pcs,
                    LensingConfig::npp6th, block_coefficients, fit_diagnostics);
                if (fit_status != LinearSolve::SolveStatus::Normal) {
                    LinearSolve::reportFailure(
                        "PSFRecons::interpolate_6th", fit_status,
                        "ccd=" + std::to_string(ichip) +
                            " block=" + std::to_string(j) + std::to_string(k) +
                            " " + LinearSolve::diagnosticsContext(fit_diagnostics) +
                            " action=POLYNOMIAL_ONLY");
                    for (int i = 0; i < LensingConfig::npp6th; ++i) {
                        coeff_file << std::scientific << std::setprecision(10)
                                   << -1.0e30f
                                   << (i == LensingConfig::npp6th - 1 ? "" : " ");
                    }
                    coeff_file << "\n";
                    coeff_file.close();
                    continue;
                }

                for (int u = 0; u < LensingConfig::n_pcs; ++u) {
                    for (int i = 0; i < LensingConfig::npp6th; ++i) {
                        coeff_file << std::scientific << std::setprecision(10)
                                   << block_coefficients[
                                          static_cast<size_t>(u) * LensingConfig::npp6th + i]
                                   << (i == LensingConfig::npp6th - 1 ? "" : " ");
                    }
                    coeff_file << "\n";
                }
                coeff_file.close();
            }
        }
    }

    // ==========================================
    // Function: Fit one sixth-order PCA coefficient surface
    // Method: Reuse the checked block fitter while preserving the F77 single-precision normal-equation solution.
    // ==========================================
    LinearSolve::SolveStatus interpolate_6th(
        int nsam, const std::vector<float>& x, const std::vector<float>& y,
        const std::vector<float>& z, int nc, std::vector<float>& coef,
        LinearSolve::SolveDiagnostics* diagnostics) {
        coef.assign(nc, 0.0f);
        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        if (nsam < 0 || static_cast<int>(x.size()) < nsam ||
            static_cast<int>(y.size()) < nsam || static_cast<int>(z.size()) < nsam) {
            diag = {};
            diag.rows = nsam;
            diag.cols = nc;
            diag.rank = std::min(nsam, nc);
            diag.required_rank = nc;
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        std::vector<float> x_used(x.begin(), x.begin() + nsam);
        std::vector<float> y_used(y.begin(), y.begin() + nsam);
        std::vector<int> sample_indices(nsam, 0);
        std::vector<double> projected_coefficients(
            static_cast<size_t>(nsam) * LensingConfig::n_pcs, 0.0);
        for (int sample = 0; sample < nsam; ++sample) {
            sample_indices[sample] = sample;
            projected_coefficients[
                static_cast<size_t>(sample) * LensingConfig::n_pcs] = z[sample];
        }

        std::vector<float> all_coefficients;
        LinearSolve::SolveStatus status = fitPcaBlockCoefficients(
            x_used, y_used, sample_indices, projected_coefficients,
            1, nc, all_coefficients, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }
        std::copy_n(all_coefficients.begin(), nc, coef.begin());
        return LinearSolve::SolveStatus::Normal;
    }

    // Plot residuals and map modified residuals for a specific exposure
    // ==========================================
    // Function: Reconstruct and write exposure-wide PSF residual diagnostics.
    // Method: Match F77 invalid-chip fallback by writing -999 rows and continuing.
    // ==========================================
    void plotResidualsV2(int iexpo) {
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(EXPO_FILE[iexpo - 1], image_files, dir_output);
        std::string prefix_e = UniversalUtils::getPrefixExpo(image_files[0]);

        float res_factor = 1.0f;
        std::string rescale_filename = dir_output + "/stamps/dat_Rescale/" + prefix_e + "_factor.dat";
        std::ifstream rescale_file(rescale_filename);
        if (rescale_file.is_open() && (rescale_file >> res_factor)) {
            rescale_file.close();
        } else {
            std::cerr << "cannot find rescale factor file" << std::endl;
            std::exit(1);
        }

        std::string out_filename = dir_output + "/stamps/dat_StarCompV2/" + prefix_e + "_star_comp_expo_v2.dat";
        MainIO::OutputFile file11(out_filename);
        if (!file11.is_open()) {
            std::cerr << "plotResidualsV2: Error opening " << out_filename << std::endl;
            return;
        }

        std::string in_filename = dir_output + "/stamps/dat_StarComp/" + prefix_e + "_star_comp_expo.dat";
        std::ifstream file10(in_filename);
        if (!file10.is_open()) {
            std::cerr << in_filename << " does not exist!!" << std::endl;
            return;
        }

        auto writeInvalidChip = [&](int ichip, int nstar, int valid) {
            file11 << ichip << " " << nstar << " " << valid << "\n";
            for (int chip_circle = 0; chip_circle < nstar; ++chip_circle) {
                double dummy = 0.0;
                for (int u = 0; u < 8; ++u) {
                    file10 >> dummy;
                }
                for (int u = 0; u < 8; ++u) {
                    file11 << -999.0 << (u == 7 ? "" : " ");
                }
                file11 << "\n";
            }
        };

        int ichip = 0, nstar = 0, valid = 0;
        while (file10 >> ichip >> nstar >> valid) {
            int proc_error = 0;
            
            if (valid < 0) {
                writeInvalidChip(ichip, nstar, valid);
                continue;
            }

            std::string prefix_c = UniversalUtils::getPrefix(image_files[ichip - 1]);
            std::string coe_filename = OutputLayout::chipPath(
                dir_output, "stamps/dat_PsfFit", prefix_c, "_PSF_coe_local.dat");
            std::ifstream file13(coe_filename);
            int nstar_coe = 0, status = 0;
            int ns = LensingConfig::ns;
            int npl = LensingConfig::npl;
            std::vector<double> local_coe(static_cast<size_t>(ns) * ns * (npl + 1), 0.0);

            if (file13.is_open()) {
                if (file13 >> nstar_coe >> status) {
                    if (status == -1) {
                        proc_error = 1;
                    } else {
                        for (int i = 0; i < ns; ++i) {
                            for (int j = 0; j < ns; ++j) {
                                for (int k = 0; k < npl + 1; ++k) {
                                    if (!(file13 >> local_coe[(j * ns + i) * (npl + 1) + k])) {
                                        proc_error = 1;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    proc_error = 1;
                }
                file13.close();
            } else {
                proc_error = 1;
            }

            if (proc_error == 1) {
                valid = -1;
                writeInvalidChip(ichip, nstar, valid);
                continue;
            }

            int chip_index = UniversalUtils::getChipId(image_files[ichip - 1]);
            if (chip_index < 0) {
                std::exit(1);
            }
            file11 << ichip << " " << nstar << " " << valid << "\n";

            for (int chip_circle = 0; chip_circle < nstar; ++chip_circle) {
                std::vector<double> psf_para(8);
                for (int u = 0; u < 8; ++u) {
                    if (!(file10 >> psf_para[u])) {
                        std::cerr << "Error in reading file: " << in_filename << std::endl;
                        return;
                    }
                }

                bool undefined_flag = false;
                for (int u = 0; u < 8; ++u) {
                    if (std::isnan(psf_para[u]) || std::isinf(psf_para[u])) {
                        undefined_flag = true;
                        break;
                    }
                }

                if (undefined_flag) {
                    for (int u = 0; u < 8; ++u) {
                        file11 << std::scientific << std::setprecision(17) << psf_para[u] << (u == 7 ? "" : " ");
                    }
                    file11 << "\n";
                } else {
                    double px = psf_para[0];
                    double py = psf_para[1];
                    std::vector<float> psf_model(ns * ns, 0.0f);
                    getPSFModelHierarchical(chip_index, px, py, res_factor, local_coe, psf_model);

                    if (std::isnan(psf_model[0])) {
                        file11 << -999.0 << " " << -999.0;
                        for (int u = 0; u < 6; ++u) {
                            file11 << " " << std::numeric_limits<float>::quiet_NaN();
                        }
                        file11 << "\n";
                    } else {
                        std::array<double, 2> ee = {0.0, 0.0};
                        double size = 0.0;
                        PSFModel::getPowerAll(ns, ns, psf_model, ee, size, 0.02f);
                        for (int u = 0; u < 5; ++u) {
                            file11 << std::scientific << std::setprecision(17) << psf_para[u] << " ";
                        }
                        file11 << std::scientific << std::setprecision(17) << size << " " << ee[0] << " " << ee[1] << "\n";
                    }
                }
            }
        }
        file10.close();
        file11.close();
    }

    // Hierarchical PSF model reconstruction at a specific position (x, y) on a CCD chip
    void getPSFModelHierarchical(int i_ccd, double x, double y, float refactor, 
                                 const std::vector<double>& local_coe, std::vector<float>& psf_model) {
        int ns = LensingConfig::ns;
        int npl = LensingConfig::npl;
        
        psf_model.assign(ns * ns, 0.0f);

        if (!PSFModel::is_data_loaded) {
            std::cerr << "Error: PSF data not loaded! Call init first." << std::endl;
            std::exit(1);
        }

        double xx_norm = 2.0 * (x / static_cast<double>(LensingConfig::chipnx)) - 1.0;
        double yy_norm = 2.0 * (y / static_cast<double>(LensingConfig::chipny)) - 1.0;

        std::vector<float> psf_layer1(ns * ns, 0.0f);
        std::vector<float> psf_layer2(ns * ns, 0.0f);
        PSFModel::getPSFModel(ns, npl, local_coe, xx_norm, yy_norm, psf_layer1, psf_layer2);

        if (std::isnan(psf_layer1[0])) {
            std::cerr << "Error in getting PSF model layer1 for chip " << i_ccd << std::endl;
            psf_model = psf_layer1;
            return;
        }

        int bx = 0;
        int by = 0;
        double nxc = static_cast<double>(LensingConfig::chipnx) / LensingConfig::nblocks;
        double nyc = static_cast<double>(LensingConfig::chipny) / LensingConfig::nblocks;
        double x_norm = 0.0;
        double y_norm = 0.0;

        for (int j = 1; j <= LensingConfig::nblocks; ++j) {
            double lxd = (j - 1.0) * nxc;
            double uxd = j * nxc;
            for (int k = 1; k <= LensingConfig::nblocks; ++k) {
                double lyd = (k - 1.0) * nyc;
                double uyd = k * nyc;
                if (x < lxd || x > uxd) continue;
                if (y < lyd || y > uyd) continue;
                bx = j;
                by = k;
                x_norm = 2.0 * (x - lxd) / (uxd - lxd) - 1.0;
                y_norm = 2.0 * (y - lyd) / (uyd - lyd) - 1.0;
                break;
            }
            if (bx != 0) break;
        }

        if (bx == 0 || by == 0) {
            std::cerr << "Error: cannot find the block for source " << i_ccd << std::endl;
            psf_model = psf_layer1;
            return;
        }

        if (PSFModel::global_components[PSFModel::getCompIndex(i_ccd - 1, 0, 0)] < -1.0e20) {
            psf_model = psf_layer1;
            return;
        }

        std::vector<double> vec_b(LensingConfig::npp6th);
        for (int j = 0; j < LensingConfig::npp6th; ++j) {
            vec_b[j] = UniversalUtils::fitFunc2(x_norm, y_norm, j);
        }

        if (PSFModel::global_poly_coefs[PSFModel::getPolyIndex(i_ccd - 1, bx - 1, by - 1, 0, 0)] < -1.0e20f) {
            psf_model = psf_layer1;
            return;
        }

        std::vector<float> coeff_val(LensingConfig::n_pcs, 0.0f);
        for (int u = 0; u < LensingConfig::n_pcs; ++u) {
            double val = 0.0;
            for (int j = 0; j < LensingConfig::npp6th; ++j) {
                val += PSFModel::global_poly_coefs[PSFModel::getPolyIndex(i_ccd - 1, bx - 1, by - 1, u, j)] * vec_b[j];
            }
            coeff_val[u] = static_cast<float>(val);
        }
        for (int row = 0; row < ns; ++row) {
            for (int col = 0; col < ns; ++col) {
                int k_idx = col * ns + row;
                double mean_val = PSFModel::global_mean_psf[PSFModel::getMeanIndex(i_ccd - 1, k_idx)];
                double val = mean_val;
                for (int j = 0; j < LensingConfig::n_pcs; ++j) {
                    val += PSFModel::global_components[PSFModel::getCompIndex(i_ccd - 1, k_idx, j)] * coeff_val[j];
                }
                psf_layer2[row * ns + col] = static_cast<float>(val);
            }
        }

        PSFModel::PSF_unscale(psf_layer2, refactor);

        for (int idx = 0; idx < ns * ns; ++idx) {
            psf_model[idx] = psf_layer1[idx] + psf_layer2[idx];
        }
    }

    // ==========================================
    // Function: Fit normalized PSF coefficients and covariance with failure status
    // Method: Propagate the first non-finite, rank, or solver failure without emitting partial covariance output.
    // ==========================================
    LinearSolve::SolveStatus itpNormPSFCov(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        std::vector<float>& sigmarr, std::vector<float>& comat,
        LinearSolve::SolveDiagnostics* diagnostics) {
        PSF_coe.assign(static_cast<size_t>(ns) * ns * (npp + 1), 0.0);
        sigmarr.assign(static_cast<size_t>(ns) * ns, 0.0f);
        comat.assign(static_cast<size_t>(npp) * npp, 0.0f);

        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        if (nsam <= npp || ns <= 0 || npp <= 0 || nx <= 0 || ny <= 0 ||
            static_cast<int>(posi.size()) < nsam ||
            image.size() < static_cast<size_t>(nsam) * ns * ns) {
            diag.rows = nsam;
            diag.cols = npp;
            diag.rank = std::min(nsam, npp);
            diag.required_rank = npp;
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        std::vector<Point3D> points(nsam);

        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                for (int k = 0; k < nsam; ++k) {
                    points[k].x = 2.0 * (posi[k][0] / static_cast<double>(nx)) - 1.0;
                    points[k].y = 2.0 * (posi[k][1] / static_cast<double>(ny)) - 1.0;
                    points[k].z = image[static_cast<size_t>(k) * ns * ns + j * ns + i];
                }

                std::vector<double> coe0;
                LinearSolve::SolveStatus status =
                    UniversalUtils::fit2D2(points, 1, coe0, &diag);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(ns) * ns * (npp + 1), 0.0);
                    sigmarr.assign(static_cast<size_t>(ns) * ns, 0.0f);
                    comat.assign(static_cast<size_t>(npp) * npp, 0.0f);
                    return status;
                }

                std::vector<double> coep;
                std::vector<double> cov_mat;
                double sigma = 0.0;
                status = UniversalUtils::fit2D2Cov(
                    points, npp, coep, sigma, cov_mat, &diag);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(ns) * ns * (npp + 1), 0.0);
                    sigmarr.assign(static_cast<size_t>(ns) * ns, 0.0f);
                    comat.assign(static_cast<size_t>(npp) * npp, 0.0f);
                    return status;
                }

                for (int k = 0; k < npp; ++k) {
                    PSF_coe[(j * ns + i) * (npp + 1) + k] = coep[k];
                }
                PSF_coe[(j * ns + i) * (npp + 1) + npp] = coe0[0];
                sigmarr[j * ns + i] = static_cast<float>(sigma);
                for (int u = 0; u < npp; ++u) {
                    for (int v = 0; v < npp; ++v) {
                        comat[u * npp + v] = static_cast<float>(cov_mat[u * npp + v]);
                    }
                }
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }
}
