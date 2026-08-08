#include "PreProcess.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "Astrometry.hpp"
#include "ImageProcessing.hpp"
#include "NumericalRecipes.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

// Global/extern variables representing exposure filenames list (defined in main.cpp)
extern std::vector<std::string> EXPO_FILE;
extern int N_EXPO;

namespace PreProcess {

    namespace {
        // ==========================================
        // Function: Normalize a sub-region onto [-1,1]x[-1,1] for polynomial fitting
        // Method: Build the affine map (p - mid) / half from the F77 one-based inclusive
        //         bounds [x_start+1, x_end] x [y_start+1, y_end], so the monomial basis stays
        //         well conditioned on every amplifier instead of only on the one nearest the origin.
        // Note:   numerical_fix F1/F2. The legacy code multiplied the ABSOLUTE pixel index by a
        //         single ratio 1/max(dx,dy); on the second amplifier that placed x in [0.275,0.470]
        //         where {1,x,x^2} are nearly collinear (cond(A^T A)=2.0e9, eps32*cond=235).
        //         The tensor-product monomial span is invariant under an affine change of x and of
        //         y, so the fitted surface is unchanged in exact arithmetic; only the conditioning
        //         improves (cond(A^T A) drops to 5.3e2).
        // ==========================================
        struct SubRegionFrame {
            double x_mid = 0.0;
            double y_mid = 0.0;
            double x_half_inv = 0.0;
            double y_half_inv = 0.0;

            SubRegionFrame(int x_start, int x_end, int y_start, int y_end) {
                x_mid = 0.5 * (static_cast<double>(x_start + 1) + static_cast<double>(x_end));
                y_mid = 0.5 * (static_cast<double>(y_start + 1) + static_cast<double>(y_end));
                x_half_inv = 2.0 / static_cast<double>(std::max(x_end - x_start - 1, 1));
                y_half_inv = 2.0 / static_cast<double>(std::max(y_end - y_start - 1, 1));
            }

            // Map a zero-based pixel index to the normalized frame (the +1 restores F77 indexing).
            double normX(int i) const { return (static_cast<double>(i + 1) - x_mid) * x_half_inv; }
            double normY(int j) const { return (static_cast<double>(j + 1) - y_mid) * y_half_inv; }
        };
        struct SigPlane {
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
        };

        // ==========================================
        // Function: Compute a deterministic median of a finite sample
        // Method: Sort in ascending order and average the two central values for an even count.
        // ==========================================
        double medianValue(std::vector<double>& values) {
            std::sort(values.begin(), values.end());
            const size_t n = values.size();
            if ((n & 1U) == 0U) return 0.5 * (values[n / 2 - 1] + values[n / 2]);
            return values[n / 2];
        }

        // ==========================================
        // Function: Estimate stripe median and width from valid blocks
        // Method: Preserve the F77 1000-sample robust estimator while drawing only from
        //         blocks that contain at least one usable pixel.
        // ==========================================
        void getStripeSigMed(const std::vector<float>& map,
                             const std::vector<unsigned char>& valid,
                             float& sig, float& med) {
            std::vector<float> values;
            values.reserve(map.size());
            for (size_t i = 0; i < map.size(); ++i) {
                if (valid[i] != 0U) values.push_back(map[i]);
            }

            if (values.empty()) {
                sig = -1.0f;
                med = 0.0f;
                return;
            }

            constexpr int sample_count = 1000;
            std::vector<float> samples(sample_count);
            for (int i = 0; i < sample_count; ++i) {
                size_t index = static_cast<size_t>(
                    NumericalRecipes::ran1() * static_cast<double>(values.size()));
                if (index >= values.size()) index = values.size() - 1;
                samples[static_cast<size_t>(i)] = values[index];
            }

            std::sort(samples.begin(), samples.end());
            sig = 0.5f * (samples[5 * sample_count / 6 - 1]
                          - samples[sample_count / 6 - 1]);
            med = samples[sample_count / 2 - 1];
        }

        // ==========================================
        // Function: Validate a linear F6 plane on an amplifier
        // Method: Check the four one-based rectangle corners for finite positive extrema and the
        //         configured maximum-to-minimum ratio.
        // ==========================================
        bool validateSigPlane(int x_start, int x_end, int y_start, int y_end,
                              const SigPlane& plane, double& plane_min, double& plane_max) {
            const double x1 = static_cast<double>(x_start + 1);
            const double x2 = static_cast<double>(x_end);
            const double y1 = static_cast<double>(y_start + 1);
            const double y2 = static_cast<double>(y_end);
            const double corners[4] = {
                plane.a + plane.b * x1 + plane.c * y1,
                plane.a + plane.b * x2 + plane.c * y1,
                plane.a + plane.b * x1 + plane.c * y2,
                plane.a + plane.b * x2 + plane.c * y2
            };

            plane_min = corners[0];
            plane_max = corners[0];
            for (double value : corners) {
                if (!std::isfinite(value)) {
                    std::cerr << "Error / set_sig nonfinite plane corner" << std::endl;
                    return false;
                }
                plane_min = std::min(plane_min, value);
                plane_max = std::max(plane_max, value);
            }
            if (plane_min <= LensingConfig::sig_plane_min) {
                std::cerr << "Error / set_sig nonpositive noise plane " << plane_min << " "
                          << plane_max << std::endl;
                return false;
            }
            if (plane_max / plane_min > LensingConfig::sig_max_plane_ratio) {
                std::cerr << "Error / set_sig excessive plane variation " << plane_min << " "
                          << plane_max << std::endl;
                return false;
            }
            return true;
        }

        // ==========================================
        // Function: Build the image-only private F6 brightness mask
        // Method: Symmetrically clip around the amplifier mode, optionally scale the width by a
        //         provisional plane, and dilate without modifying the caller's weight map.
        // ==========================================
        bool buildSigPrivateMask(int x_start, int x_end, int y_start, int y_end, int nx,
                                 const std::vector<float>& image, const std::vector<int>& weight,
                                 double center, double sigma0, bool use_local,
                                 const SigPlane& plane, std::vector<unsigned char>& src,
                                 long& nmask) {
            std::fill(src.begin(), src.end(), static_cast<unsigned char>(0));
            double plane_center = 1.0;
            if (use_local) {
                const double xmid = 0.5 * static_cast<double>(x_start + 1 + x_end);
                const double ymid = 0.5 * static_cast<double>(y_start + 1 + y_end);
                plane_center = plane.a + plane.b * xmid + plane.c * ymid;
                if (!std::isfinite(plane_center) || plane_center <= LensingConfig::sig_plane_min) {
                    std::cerr << "Error / set_sig invalid local mask plane" << std::endl;
                    return false;
                }
            }

            long nbase = 0;
            for (int i = x_start; i < x_end; ++i) {
                for (int j = y_start; j < y_end; ++j) {
                    const size_t idx = static_cast<size_t>(j) * nx + i;
                    if (weight[idx] <= 0) continue;
                    ++nbase;
                    const double value = image[idx];
                    if (!std::isfinite(value)) {
                        std::cerr << "Error / set_sig nonfinite mask pixel" << std::endl;
                        return false;
                    }
                    double sigma_local = sigma0;
                    if (use_local) {
                        const double local_plane = plane.a + plane.b * static_cast<double>(i + 1)
                                                 + plane.c * static_cast<double>(j + 1);
                        if (!std::isfinite(local_plane)
                            || local_plane <= LensingConfig::sig_plane_min) {
                            std::cerr << "Error / set_sig invalid local mask value" << std::endl;
                            return false;
                        }
                        sigma_local *= std::sqrt(local_plane / plane_center);
                    }
                    if (std::abs(value - center) <= LensingConfig::sig_clip_k * sigma_local) continue;

                    const int i1 = std::max(x_start, i - LensingConfig::sig_rdil);
                    const int i2 = std::min(x_end - 1, i + LensingConfig::sig_rdil);
                    const int j1 = std::max(y_start, j - LensingConfig::sig_rdil);
                    const int j2 = std::min(y_end - 1, j + LensingConfig::sig_rdil);
                    for (int ii = i1; ii <= i2; ++ii) {
                        for (int jj = j1; jj <= j2; ++jj) {
                            src[static_cast<size_t>(jj) * nx + ii] = 1;
                        }
                    }
                }
            }
            if (nbase < LensingConfig::sig_min_fit_triples) {
                std::cerr << "Error / set_sig too few base pixels " << nbase << std::endl;
                return false;
            }

            nmask = 0;
            for (int i = x_start; i < x_end; ++i) {
                for (int j = y_start; j < y_end; ++j) {
                    const size_t idx = static_cast<size_t>(j) * nx + i;
                    if (weight[idx] > 0 && src[idx] != 0) ++nmask;
                }
            }
            return true;
        }

        // ==========================================
        // Function: Solve the scaled three-by-three F6 normal equations
        // Method: Use checked Cholesky factorization and reject singular or nonfinite pivots.
        // ==========================================
        bool solveSigPlane3(const Eigen::Matrix3d& matrix, const Eigen::Vector3d& rhs,
                            Eigen::Vector3d& solution) {
            Eigen::Matrix3d lower = Eigen::Matrix3d::Zero();
            Eigen::Vector3d work = Eigen::Vector3d::Zero();
            solution.setZero();

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double value = matrix(i, j);
                    for (int k = 0; k < j; ++k) value -= lower(i, k) * lower(j, k);
                    if (!std::isfinite(value)) {
                        std::cerr << "Error / set_sig nonfinite Cholesky term" << std::endl;
                        return false;
                    }
                    if (i == j) {
                        if (value <= LensingConfig::sig_pivot_min) {
                            std::cerr << "Error / set_sig singular normal matrix" << std::endl;
                            return false;
                        }
                        lower(i, j) = std::sqrt(value);
                    } else {
                        if (lower(j, j) <= LensingConfig::sig_pivot_min) {
                            std::cerr << "Error / set_sig invalid Cholesky pivot" << std::endl;
                            return false;
                        }
                        lower(i, j) = value / lower(j, j);
                    }
                }
            }

            for (int i = 0; i < 3; ++i) {
                double value = rhs(i);
                for (int k = 0; k < i; ++k) value -= lower(i, k) * work(k);
                work(i) = value / lower(i, i);
                if (!std::isfinite(work(i))) {
                    std::cerr << "Error / set_sig nonfinite forward solve" << std::endl;
                    return false;
                }
            }
            for (int i = 2; i >= 0; --i) {
                double value = work(i);
                for (int k = i + 1; k < 3; ++k) value -= lower(k, i) * solution(k);
                solution(i) = value / lower(i, i);
                if (!std::isfinite(solution(i))) {
                    std::cerr << "Error / set_sig nonfinite backward solve" << std::endl;
                    return false;
                }
            }
            return true;
        }

        // ==========================================
        // Function: Fit the raw F6 2*sigma^2 plane from private-mask survivors
        // Method: Accumulate every eligible triple in double precision, Jacobi-scale the normal
        //         equations, solve by checked Cholesky, and validate the raw plane.
        // ==========================================
        bool fitSigMaskedPlane(int x_start, int x_end, int y_start, int y_end, int nx,
                               const std::vector<float>& image, const std::vector<int>& weight,
                               const std::vector<unsigned char>& src, SigPlane& plane,
                               long& npossible, long& nuse) {
            const SubRegionFrame frame(x_start, x_end, y_start, y_end);
            Eigen::Matrix3d matrix = Eigen::Matrix3d::Zero();
            Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
            npossible = 0;
            nuse = 0;

            for (int i = x_start; i < x_end - 1; ++i) {
                for (int j = y_start; j < y_end - 1; ++j) {
                    const size_t idx = static_cast<size_t>(j) * nx + i;
                    const size_t idx_x = idx + 1;
                    const size_t idx_y = idx + nx;
                    if (weight[idx] <= 0 || weight[idx_x] <= 0 || weight[idx_y] <= 0) continue;
                    ++npossible;
                    if (src[idx] != 0 || src[idx_x] != 0 || src[idx_y] != 0) continue;

                    const double value = image[idx];
                    const double dx = value - image[idx_x];
                    const double dy = value - image[idx_y];
                    const double pix = 0.5 * (dx * dx + dy * dy);
                    if (!std::isfinite(pix)) {
                        std::cerr << "Error / set_sig nonfinite final pix" << std::endl;
                        return false;
                    }
                    const Eigen::Vector3d basis(1.0, frame.normX(i), frame.normY(j));
                    matrix.noalias() += basis * basis.transpose();
                    rhs.noalias() += pix * basis;
                    ++nuse;
                }
            }

            if (npossible < LensingConfig::sig_min_fit_triples
                || nuse < LensingConfig::sig_min_fit_triples) {
                std::cerr << "Error / set_sig too few final triples " << nuse << std::endl;
                return false;
            }
            if (static_cast<double>(nuse) / static_cast<double>(npossible)
                < LensingConfig::sig_min_fit_frac) {
                std::cerr << "Error / set_sig private mask too large " << nuse << std::endl;
                return false;
            }

            Eigen::Vector3d scale;
            for (int i = 0; i < 3; ++i) {
                if (!std::isfinite(matrix(i, i)) || matrix(i, i) <= 0.0) {
                    std::cerr << "Error / set_sig invalid normal matrix" << std::endl;
                    return false;
                }
                scale(i) = 1.0 / std::sqrt(matrix(i, i));
            }
            for (int i = 0; i < 3; ++i) {
                rhs(i) *= scale(i);
                for (int j = 0; j < 3; ++j) {
                    matrix(i, j) *= scale(i) * scale(j);
                    if (!std::isfinite(matrix(i, j))) {
                        std::cerr << "Error / set_sig nonfinite normal matrix" << std::endl;
                        return false;
                    }
                }
            }

            Eigen::Vector3d solution;
            if (!solveSigPlane3(matrix, rhs, solution)) return false;
            plane.b = solution(1) * scale(1) * frame.x_half_inv;
            plane.c = solution(2) * scale(2) * frame.y_half_inv;
            plane.a = solution(0) * scale(0) - plane.b * frame.x_mid - plane.c * frame.y_mid;
            if (!std::isfinite(plane.a) || !std::isfinite(plane.b) || !std::isfinite(plane.c)) {
                std::cerr << "Error / set_sig nonfinite raw coefficients" << std::endl;
                return false;
            }
            double plane_min = 0.0;
            double plane_max = 0.0;
            return validateSigPlane(x_start, x_end, y_start, y_end, plane,
                                    plane_min, plane_max);
        }

        // ==========================================
        // Function: Apply a validated F6 noise plane to one amplifier
        // Method: Divide the image immediately by sqrt(plane/2) in double precision.
        // ==========================================
        void applySigPlane(int x_start, int x_end, int y_start, int y_end, int nx,
                           std::vector<float>& image, const SigPlane& plane) {
            for (int i = x_start; i < x_end; ++i) {
                for (int j = y_start; j < y_end; ++j) {
                    const size_t idx = static_cast<size_t>(j) * nx + i;
                    const double value = plane.a + plane.b * static_cast<double>(i + 1)
                                               + plane.c * static_cast<double>(j + 1);
                    image[idx] = static_cast<float>(static_cast<double>(image[idx])
                                                    / std::sqrt(0.5 * value));
                }
            }
        }
    }

    // Stage 1 driver
    void preProcess(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = EXPO_FILE[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);
        
        for (const auto& image_file : image_files) {
            int cid = UniversalUtils::getChipId(image_file);

            chipPreProcess(image_file, dir_output, cid);
        }
    }

    // ==========================================
    // Function: Individual chip preprocessing
    // Method: Match the Fortran Stage 1 flow while keeping diagnostics side-effect free.
    // ==========================================
    void chipPreProcess(const std::string& imageFile, const std::string& dirOutput, int cid) {
        int proc_error = 0;
        int nx = 0, ny = 0;
        std::vector<float> array;
        WCSParams wcs;
        std::string prefix = UniversalUtils::getPrefix(imageFile);

        bool image_read_ok = FitsIO::readImagePara(imageFile, nx, ny, array, wcs);
        if (!image_read_ok || (!array.empty() && array[0] < -99990.0f)) {
            std::cerr << "Error reading image parameters of: " << imageFile << std::endl;
            proc_error = 1;
            if (nx <= 0 || ny <= 0 || array.empty()) {
                nx = LensingConfig::npx;
                ny = LensingConfig::npy;
                array.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), -99999.0f);
            } else if (array.size() < static_cast<size_t>(nx) * static_cast<size_t>(ny)) {
                array.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), -99999.0f);
            }
            array[0] = -99999.0f;
        }

        std::vector<int> weight(nx * ny, 1);
        std::vector<float> normap(nx * ny);

        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int idx = y * nx + x;
                if (array[idx] > LensingConfig::saturation_thresh) {
                    weight[idx] = 0;
                }
                normap[idx] = array[idx];
            }
        }

        double sigabc[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        int nxc = nx / 2;
        double aa = 0.0, bb = 0.0, cc = 0.0;

        if (LensingConfig::CCD_split == 2) {
            setBackground(0, nxc, 0, ny, nx, ny, normap, LensingConfig::blocksize, LensingConfig::nct, LensingConfig::ncx, proc_error);
            setBackground(nxc, nx, 0, ny, nx, ny, normap, LensingConfig::blocksize, LensingConfig::nct, LensingConfig::ncx, proc_error);
            
            setSig(0, nxc, 0, ny, nx, ny, normap, weight, aa, bb, cc, proc_error,
                   LensingConfig::sig_scale);
            if (proc_error == 0) {
                sigabc[0][0] = aa; sigabc[0][1] = bb; sigabc[0][2] = cc;
            }
            if (proc_error == 0) {
                setSig(nxc, nx, 0, ny, nx, ny, normap, weight, aa, bb, cc, proc_error,
                       LensingConfig::sig_scale);
            }
            if (proc_error == 0) {
                sigabc[1][0] = aa; sigabc[1][1] = bb; sigabc[1][2] = cc;
            }
        } else {
            setBackground(0, nx, 0, ny, nx, ny, normap, LensingConfig::blocksize, LensingConfig::nct, LensingConfig::ncx, proc_error);
            setSig(0, nx, 0, ny, nx, ny, normap, weight, aa, bb, cc, proc_error,
                   LensingConfig::sig_scale);
            if (proc_error == 0) {
                sigabc[0][0] = aa; sigabc[0][1] = bb; sigabc[0][2] = cc;
            }
        }

        // ==========================================
        // Function: Apply the DQ mask before astrometry and defect merging
        // Method: Keep setSig independent of DQ, then reject every nonzero DQ pixel from all
        //         later preprocessing stages.
        // ==========================================
        if (proc_error == 0) {
            std::string prefix_e = UniversalUtils::getPrefixExpo(imageFile);
            std::string local_mask_file = dirOutput + "/dqmask/" + prefix_e + "/"
                                          + prefix_e + "_" + std::to_string(cid) + ".fits";
            int nxx = 0, nyy = 0;
            std::vector<float> flat_weight;
            if (!FitsIO::readImage(local_mask_file, nxx, nyy, flat_weight)) {
                std::cerr << "Error / cant find mask file: " << local_mask_file << std::endl;
                proc_error = 1;
            } else if (nxx != nx || nyy != ny) {
                std::cerr << "Error / wrong size of DQ file!" << std::endl;
                proc_error = 1;
            } else {
                for (int i = 0; i < nx * ny; ++i) {
                    if (std::abs(flat_weight[static_cast<size_t>(i)]) > 1e-7f) {
                        weight[static_cast<size_t>(i)] = 0;
                        normap[static_cast<size_t>(i)] = -1000.0f;
                    }
                }
            }
        }

        std::string astroFilename = OutputLayout::chipPath(
            dirOutput, "astrometry/dat_Astro", prefix, "_astro.dat");

        std::string catfile = UniversalUtils::generateGaiaFileName(LensingConfig::ASTROMETRY_CAT, wcs.crval, proc_error);
        Astrometry::genAstrometryData(catfile, nx, ny, normap, weight, wcs, astroFilename, proc_error);

        locateDefects(nx, ny, array, normap, weight, LensingConfig::area_max, LensingConfig::area_thresh, proc_error);
        mergeDefects(nx, ny, weight, normap, LensingConfig::area_max, LensingConfig::source_thresh, LensingConfig::area_thresh, proc_error);

        if (proc_error == 0) {
            // Serialize the final combined saturation, DQ, and detected-defect mask.
            for (int i = 0; i < nx * ny; ++i) {
                if (weight[static_cast<size_t>(i)] == 0) {
                    normap[static_cast<size_t>(i)] = -1000.0f;
                }
            }
        } else {
            std::fill(normap.begin(), normap.end(), -1000.0f);
        }
        if (proc_error == 0) {
            normap[0] = -1.0f;
        } else {
            normap[0] = 1.0f;
        }

        for (int i = 1; i <= LensingConfig::CCD_split; ++i) {
            for (int j = 1; j <= 3; ++j) {
                normap[(j - 1) * nx + i] = static_cast<float>(sigabc[i - 1][j - 1]);
            }
        }

        std::string normFilename = OutputLayout::chipPath(
            dirOutput, "stamps/Norm", prefix, "_norm.fits");
        if (!FitsIO::writeImageCopyHDU(imageFile, normFilename, nx, ny, normap)) {
            std::cerr << "Error writing normalized image: " << normFilename << std::endl;
            proc_error = 1;
        }

        // if (proc_error == 0) {
        //     std::cout << "Status of processing " << imageFile << ": OK." << std::endl;
        // } else {
        //     std::cout << "Status of processing " << imageFile << ": ERROR!" << std::endl;
        // }
    }

    // Helper functions for preprocessing
    void setBackground(int x_start, int x_end, int y_start, int y_end, int nx, int ny, std::vector<float>& image,
                       int blocksize, int nct, int ncx, int& ierror) {
        if (ierror == 1) return;

        // A rough flattening of the field first
        flattenChip(x_start, x_end, y_start, y_end, nx, ny, image, 4, 2, ierror);

        if (ierror == 1) return;

        // numerical_fix F1: replaces the legacy absolute-coordinate scaling by 1/max(dx,dy).
        const SubRegionFrame frame(x_start, x_end, y_start, y_end);
        int nbx = std::max((x_end - x_start) / blocksize, 1);
        int nby = std::max((y_end - y_start) / blocksize, 1);

        constexpr int npp = 1000;
        int nfit = nbx * nby;
        std::vector<Point3D> mean_sam;
        mean_sam.reserve(nfit);

        for (int i = 0; i < nbx; ++i) {
            int xmin = i * blocksize + x_start;
            int xmax = std::min(xmin + blocksize, x_end);
            for (int j = 0; j < nby; ++j) {
                int ymin = j * blocksize + y_start;
                int ymax = std::min(ymin + blocksize, y_end);

                std::vector<Point3D> arr(npp);
                std::vector<float> pix(npp);
                for (int k = 0; k < npp; ++k) {
                    int ix = xmin + static_cast<int>(NumericalRecipes::ran1() * (xmax - xmin));
                    int iy = ymin + static_cast<int>(NumericalRecipes::ran1() * (ymax - ymin));
                    arr[k].x = frame.normX(ix);
                    arr[k].y = frame.normY(iy);
                    arr[k].z = image[iy * nx + ix];
                    pix[k] = static_cast<float>(arr[k].z);
                }

                std::vector<int> indx;
                NumericalRecipes::indexx(pix, indx);

                int median_idx = indx[npp / 2 - 1];
                mean_sam.push_back(arr[median_idx]);
            }
        }

        int nct_min = nct * 3 / 2;
        int nsam1 = mean_sam.size();
        bool changed = true;

        double aa = 0.0, bb = 0.0, cc = 0.0;
        while (changed && nsam1 >= nct_min) {
            std::vector<Point3D> current_points(nsam1);
            for (int i = 0; i < nsam1; ++i) {
                current_points[i] = mean_sam[i];
            }

            UniversalUtils::findSlope2D(current_points, aa, bb, cc);

            std::vector<double> sam(nsam1);
            std::vector<double> tempo(nsam1);
            for (int i = 0; i < nsam1; ++i) {
                sam[i] = mean_sam[i].z - aa - bb * mean_sam[i].x - cc * mean_sam[i].y;
                tempo[i] = sam[i];
            }

            std::sort(tempo.begin(), tempo.end());
            double mean = tempo[nsam1 / 2 - 1];
            double sig = 0.5 * (tempo[(nsam1 * 5) / 6 - 1] - tempo[nsam1 / 6 - 1]);

            int nsam = 0;
            changed = false;
            std::vector<Point3D> filtered_mean_sam;
            filtered_mean_sam.reserve(nsam1);

            for (int i = 0; i < nsam1; ++i) {
                if (std::abs(sam[i] - mean) < 3.0 * sig) {
                    nsam++;
                    filtered_mean_sam.push_back(mean_sam[i]);
                } else {
                    changed = true;
                }
            }
            mean_sam = filtered_mean_sam;
            nsam1 = nsam;
        }

        if (nsam1 < nct_min) {
            std::cerr << "Background not stable enough! nsam1 = " << nsam1 << std::endl;
            ierror = 1;
            return;
        }

        std::vector<double> c_coeffs;
        UniversalUtils::fit2D(mean_sam, nct, ncx, c_coeffs);

        // numerical_fix F1: evaluate on the same normalized frame used to fit.
        for (int i = x_start; i < x_end; ++i) {
            for (int j = y_start; j < y_end; ++j) {
                double x = frame.normX(i);
                double y = frame.normY(j);
                image[j * nx + i] -= UniversalUtils::funcVal(x, y, nct, ncx, c_coeffs);
            }
        }
    }

    void flattenChip(int x_start, int x_end, int y_start, int y_end, int nx, int ny, std::vector<float>& array,
                     int nct, int ncx, int& ierror) {
        if (ierror == 1) return;

        // numerical_fix F1: same defect and same cure as in setBackground.
        const SubRegionFrame frame(x_start, x_end, y_start, y_end);
        constexpr int npp = 1000;

        std::vector<float> pix(npp);
        std::vector<Point3D> arr2(npp);

      for (int i = 0; i < npp; ++i) {
          int ix = x_start + static_cast<int>(NumericalRecipes::ran1() * (x_end - x_start - 1));
          int iy = y_start + static_cast<int>(NumericalRecipes::ran1() * (y_end - y_start - 1));
          pix[i] = array[iy * nx + ix];
           arr2[i].x = frame.normX(ix);
           arr2[i].y = frame.normY(iy);
           arr2[i].z = pix[i];
       }

        std::sort(pix.begin(), pix.end());

        if (pix[0] == pix[npp - 1]) {
            ierror = 1;
            return;
        }

        double arr_min = pix[npp / 3 - 1];
        double arr_max = pix[(2 * npp) / 3 - 1];

        std::vector<Point3D> arr;
        arr.reserve(npp);
        for (int i = 0; i < npp; ++i) {
            if (arr2[i].z >= arr_min && arr2[i].z <= arr_max) {
                arr.push_back(arr2[i]);
            }
        }

        int nr = arr.size();
        if (nr < npp / 10) {
            ierror = 1;
            return;
        }

        std::vector<double> c_coeffs;
        UniversalUtils::fit2D(arr, nct, ncx, c_coeffs);

        // numerical_fix F1: evaluate on the same normalized frame used to fit.
        for (int i = x_start; i < x_end; ++i) {
            for (int j = y_start; j < y_end; ++j) {
                double x = frame.normX(i);
                double y = frame.normY(j);
                array[j * nx + i] -= UniversalUtils::funcVal(x, y, nct, ncx, c_coeffs);
            }
        }
    }

    // ==========================================
    // Function: Estimate, validate, and apply one amplifier's noise-sigma plane
    // Method: Match the F77 mode-bar estimator, use a private symmetric clip mask, fit every
    //         surviving base-valid triple, and mutate the amplifier only after final validation.
    // ==========================================
    void setSig(int x_start, int x_end, int y_start, int y_end, int nx, int ny,
                std::vector<float>& image, const std::vector<int>& weight,
                double& aa, double& bb, double& cc, int& ierror, double sig_scale) {
        aa = 0.0;
        bb = 0.0;
        cc = 0.0;
        if (ierror != 0) return;

        const size_t image_size = static_cast<size_t>(nx) * static_cast<size_t>(ny);
        if (x_start < 0 || x_end > nx || y_start < 0 || y_end > ny
            || x_end <= x_start || y_end <= y_start
            || image.size() < image_size || weight.size() < image_size) {
            std::cerr << "Error / set_sig invalid amplifier geometry" << std::endl;
            ierror = 1;
            return;
        }
        if (!std::isfinite(sig_scale) || sig_scale <= 0.0) {
            std::cerr << "Error / set_sig invalid sig_scale" << std::endl;
            ierror = 1;
            return;
        }
        if (LensingConfig::sig_hist_nbin < 4
            || (LensingConfig::sig_hist_nbin % 2) != 0
            || LensingConfig::sig_clip_niter < 1) {
            std::cerr << "Error / set_sig invalid configured iteration" << std::endl;
            ierror = 1;
            return;
        }

        const int width = x_end - x_start;
        const int height = y_end - y_start;
        const int nbx = (width + LensingConfig::sig_blocksize - 1)
                      / LensingConfig::sig_blocksize;
        const int nby = (height + LensingConfig::sig_blocksize - 1)
                      / LensingConfig::sig_blocksize;
        if (nbx * nby > LensingConfig::sig_max_blocks) {
            std::cerr << "Error / set_sig block table overflow " << nbx * nby << std::endl;
            ierror = 1;
            return;
        }

        std::vector<double> block_medians;
        std::vector<double> block_sigmas;
        block_medians.reserve(static_cast<size_t>(nbx * nby));
        block_sigmas.reserve(static_cast<size_t>(nbx * nby));
        std::vector<double> block_values;
        block_values.reserve(LensingConfig::sig_block_max);

        for (int ib = 0; ib < nbx; ++ib) {
            const int xmin = x_start + ib * width / nbx;
            const int xmax = x_start + (ib + 1) * width / nbx;
            for (int jb = 0; jb < nby; ++jb) {
                const int ymin = y_start + jb * height / nby;
                const int ymax = y_start + (jb + 1) * height / nby;

                block_values.clear();
                for (int i = xmin; i < xmax; ++i) {
                    for (int j = ymin; j < ymax; ++j) {
                        const size_t idx = static_cast<size_t>(j) * nx + i;
                        if (weight[idx] <= 0) continue;
                        const double value = image[idx];
                        if (!std::isfinite(value)) {
                            std::cerr << "Error / set_sig nonfinite image pixel" << std::endl;
                            ierror = 1;
                            return;
                        }
                        block_values.push_back(value);
                    }
                }
                if (block_values.size() < static_cast<size_t>(LensingConfig::sig_min_block_pixels)) {
                    std::cerr << "Error / set_sig too few block pixels " << block_values.size()
                              << std::endl;
                    ierror = 1;
                    return;
                }
                if (block_values.size() > static_cast<size_t>(LensingConfig::sig_block_max)) {
                    std::cerr << "Error / set_sig block buffer overflow" << std::endl;
                    ierror = 1;
                    return;
                }
                const double block_median = medianValue(block_values);

                block_values.clear();
                for (int i = xmin; i < xmax - 1; ++i) {
                    for (int j = ymin; j < ymax - 1; ++j) {
                        const size_t idx = static_cast<size_t>(j) * nx + i;
                        const size_t idx_x = idx + 1;
                        const size_t idx_y = idx + nx;
                        if (weight[idx] <= 0 || weight[idx_x] <= 0 || weight[idx_y] <= 0) continue;
                        const double value = image[idx];
                        const double dx = value - image[idx_x];
                        const double dy = value - image[idx_y];
                        const double pix = 0.5 * (dx * dx + dy * dy);
                        if (!std::isfinite(pix)) {
                            std::cerr << "Error / set_sig nonfinite block pix" << std::endl;
                            ierror = 1;
                            return;
                        }
                        block_values.push_back(pix);
                    }
                }
                if (block_values.size() < static_cast<size_t>(LensingConfig::sig_min_block_triples)) {
                    std::cerr << "Error / set_sig too few block triples " << block_values.size()
                              << std::endl;
                    ierror = 1;
                    return;
                }
                if (block_values.size() > static_cast<size_t>(LensingConfig::sig_block_max)) {
                    std::cerr << "Error / set_sig triple buffer overflow" << std::endl;
                    ierror = 1;
                    return;
                }
                const double pix_median = medianValue(block_values);
                if (!std::isfinite(pix_median) || pix_median <= LensingConfig::sig_plane_min) {
                    std::cerr << "Error / set_sig invalid block median" << std::endl;
                    ierror = 1;
                    return;
                }
                block_medians.push_back(block_median);
                block_sigmas.push_back(std::sqrt(pix_median / LensingConfig::sig_median_ratio));
            }
        }

        if (block_medians.size() < static_cast<size_t>(LensingConfig::sig_min_blocks)) {
            std::cerr << "Error / set_sig too few valid blocks " << block_medians.size() << std::endl;
            ierror = 1;
            return;
        }
        const double center_seed = medianValue(block_medians);
        const double sigma_seed = medianValue(block_sigmas);
        if (!std::isfinite(center_seed) || !std::isfinite(sigma_seed) || sigma_seed <= 0.0) {
            std::cerr << "Error / set_sig invalid mode seeds" << std::endl;
            ierror = 1;
            return;
        }

        std::vector<long> histogram(static_cast<size_t>(LensingConfig::sig_hist_nbin), 0);
        double hist_lo = center_seed - LensingConfig::sig_hist_range * sigma_seed;
        double hist_hi = center_seed + LensingConfig::sig_hist_range * sigma_seed;
        double hist_step = (hist_hi - hist_lo) / LensingConfig::sig_hist_nbin;
        if (!std::isfinite(hist_step) || hist_step <= 0.0) {
            std::cerr << "Error / set_sig invalid histogram range" << std::endl;
            ierror = 1;
            return;
        }

        long nhist = 0;
        for (int i = x_start; i < x_end; ++i) {
            for (int j = y_start; j < y_end; ++j) {
                const size_t idx = static_cast<size_t>(j) * nx + i;
                if (weight[idx] <= 0) continue;
                const double value = image[idx];
                if (!std::isfinite(value)) {
                    std::cerr << "Error / set_sig nonfinite mode pixel" << std::endl;
                    ierror = 1;
                    return;
                }
                if (value < hist_lo || value >= hist_hi) continue;
                int bin = static_cast<int>((value - hist_lo) / hist_step);
                bin = std::max(0, std::min(LensingConfig::sig_hist_nbin - 1, bin));
                ++histogram[static_cast<size_t>(bin)];
                ++nhist;
            }
        }

        int mode_bin = 2;
        long score_max = -1;
        for (int bin = 2; bin < LensingConfig::sig_hist_nbin - 2; ++bin) {
            const long score = histogram[bin - 2] + 2 * histogram[bin - 1]
                             + 3 * histogram[bin] + 2 * histogram[bin + 1]
                             + histogram[bin + 2];
            if (score > score_max) {
                score_max = score;
                mode_bin = bin;
            }
        }
        if (nhist < LensingConfig::sig_min_lower_count
            || histogram[mode_bin] < LensingConfig::sig_min_mode_count
            || mode_bin <= 2 || mode_bin >= LensingConfig::sig_hist_nbin - 3) {
            std::cerr << "Error / set_sig unresolved brightness mode" << std::endl;
            ierror = 1;
            return;
        }

        const long score_minus = histogram[mode_bin - 3] + 2 * histogram[mode_bin - 2]
                               + 3 * histogram[mode_bin - 1] + 2 * histogram[mode_bin]
                               + histogram[mode_bin + 1];
        const long score_center = histogram[mode_bin - 2] + 2 * histogram[mode_bin - 1]
                                + 3 * histogram[mode_bin] + 2 * histogram[mode_bin + 1]
                                + histogram[mode_bin + 2];
        const long score_plus = histogram[mode_bin - 1] + 2 * histogram[mode_bin]
                              + 3 * histogram[mode_bin + 1] + 2 * histogram[mode_bin + 2]
                              + histogram[mode_bin + 3];
        const double denominator = static_cast<double>(score_minus - 2 * score_center + score_plus);
        double delta = 0.0;
        if (denominator < 0.0) {
            delta = 0.5 * static_cast<double>(score_minus - score_plus) / denominator;
            delta = std::max(-0.5, std::min(0.5, delta));
        }
        const double center = hist_lo + (static_cast<double>(mode_bin) + 0.5 + delta) * hist_step;

        std::fill(histogram.begin(), histogram.end(), 0);
        hist_lo = center - LensingConfig::sig_hist_range * sigma_seed;
        hist_hi = center + LensingConfig::sig_hist_range * sigma_seed;
        hist_step = (hist_hi - hist_lo) / LensingConfig::sig_hist_nbin;
        for (int i = x_start; i < x_end; ++i) {
            for (int j = y_start; j < y_end; ++j) {
                const size_t idx = static_cast<size_t>(j) * nx + i;
                if (weight[idx] <= 0) continue;
                const double value = image[idx];
                if (value < hist_lo || value >= hist_hi) continue;
                int bin = static_cast<int>((value - hist_lo) / hist_step);
                bin = std::max(0, std::min(LensingConfig::sig_hist_nbin - 1, bin));
                ++histogram[static_cast<size_t>(bin)];
            }
        }

        long nbelow = 0;
        for (int bin = 0; bin < LensingConfig::sig_hist_nbin / 2; ++bin) {
            nbelow += histogram[static_cast<size_t>(bin)];
        }
        if (nbelow < LensingConfig::sig_min_lower_count) {
            std::cerr << "Error / set_sig too few lower-half pixels " << nbelow << std::endl;
            ierror = 1;
            return;
        }

        const double target = LensingConfig::sig_lower_quantile * static_cast<double>(nbelow);
        double cumulative = 0.0;
        int quantile_bin = -1;
        double quantile_value = 0.0;
        for (int bin = 0; bin < LensingConfig::sig_hist_nbin / 2; ++bin) {
            const long count = histogram[static_cast<size_t>(bin)];
            if (quantile_bin < 0 && count > 0 && cumulative + static_cast<double>(count) >= target) {
                double fraction = (target - cumulative) / static_cast<double>(count);
                fraction = std::max(0.0, std::min(1.0, fraction));
                quantile_value = hist_lo + (static_cast<double>(bin) + fraction) * hist_step;
                quantile_bin = bin;
            }
            cumulative += static_cast<double>(count);
        }
        const double sigma0 = center - quantile_value;
        if (quantile_bin < 0 || !std::isfinite(center) || !std::isfinite(sigma0)
            || sigma0 <= std::sqrt(LensingConfig::sig_plane_min)) {
            std::cerr << "Error / set_sig invalid mode width" << std::endl;
            ierror = 1;
            return;
        }

        SigPlane raw_plane;
        std::vector<unsigned char> src(image_size, 0);
        for (int iter = 0; iter < LensingConfig::sig_clip_niter; ++iter) {
            long nmask = 0;
            if (!buildSigPrivateMask(x_start, x_end, y_start, y_end, nx, image, weight,
                                    center, sigma0, iter > 0, raw_plane, src, nmask)) {
                std::cerr << "F6_DIAG mask " << x_start + 1 << " " << x_end << " " << iter + 1
                          << " " << center << " " << sigma0 << " " << nmask << std::endl;
                ierror = 1;
                return;
            }
            long npossible = 0;
            long nuse = 0;
            if (!fitSigMaskedPlane(x_start, x_end, y_start, y_end, nx, image, weight, src,
                                   raw_plane, npossible, nuse)) {
                std::cerr << "F6_DIAG fit " << x_start + 1 << " " << x_end << " " << iter + 1
                          << " " << center << " " << sigma0 << " " << nmask << " " << nuse
                          << " " << npossible << std::endl;
                ierror = 1;
                return;
            }
        }

        SigPlane scaled_plane{raw_plane.a * sig_scale,
                              raw_plane.b * sig_scale,
                              raw_plane.c * sig_scale};
        if (!std::isfinite(scaled_plane.a) || !std::isfinite(scaled_plane.b)
            || !std::isfinite(scaled_plane.c)) {
            std::cerr << "Error / set_sig nonfinite coefficients" << std::endl;
            ierror = 1;
            return;
        }
        double plane_min = 0.0;
        double plane_max = 0.0;
        if (!validateSigPlane(x_start, x_end, y_start, y_end, scaled_plane,
                              plane_min, plane_max)) {
            ierror = 1;
            return;
        }

        // Match the F77 real-coefficient inter-stage contract before mutating the image: validate
        // once in double precision, convert to the stored precision, then validate and apply that
        // exact plane so normap and the published sigabc coefficients cannot diverge.
        const SigPlane stored_plane{static_cast<double>(static_cast<float>(scaled_plane.a)),
                                    static_cast<double>(static_cast<float>(scaled_plane.b)),
                                    static_cast<double>(static_cast<float>(scaled_plane.c))};
        if (!std::isfinite(stored_plane.a) || !std::isfinite(stored_plane.b)
            || !std::isfinite(stored_plane.c)
            || !validateSigPlane(x_start, x_end, y_start, y_end, stored_plane,
                                 plane_min, plane_max)) {
            ierror = 1;
            return;
        }
        aa = stored_plane.a;
        bb = stored_plane.b;
        cc = stored_plane.c;
        applySigPlane(x_start, x_end, y_start, y_end, nx, image, stored_plane);
    }

    void locateDefects(int nx, int ny, const std::vector<float>& array, std::vector<float>& normap,
                       std::vector<int>& weight, int area_max, int area_thresh, int& ierror) {
        constexpr int margin = 10;
        constexpr double defect_halo_thresh = 1.0;
        constexpr int y_smooth = 200;
        constexpr int x_smooth = 100;

        for (int y = 0; y < ny; ++y) {
            int mid_start = nx / 2 - margin - 1;
            int mid_end = nx / 2 + margin;
            for (int x = mid_start; x < mid_end; ++x) {
                if (x >= 0 && x < nx) {
                    weight[y * nx + x] = 0;
                }
            }
            for (int x = 0; x < margin; ++x) {
                weight[y * nx + x] = 0;
            }
            for (int x = nx - margin; x < nx; ++x) {
                weight[y * nx + x] = 0;
            }
        }
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < margin; ++y) {
                weight[y * nx + x] = 0;
            }
            for (int y = ny - margin; y < ny; ++y) {
                weight[y * nx + x] = 0;
            }
        }

       if (ierror == 1) return;

       std::vector<float> map(nx * ny);

       for (int i = 0; i < nx * ny; ++i) {
           map[i] = UniversalUtils::loga(array[i], 1);
       }

        ImageProcessing::removeContinuous(nx, ny, nx, ny, map, UniversalUtils::iden, 4);

        std::vector<float> diffx(nx * ny);
        std::vector<float> diffy(nx * ny);
        for (int y = 0; y < ny; ++y) {
            int next_y = (y + 1) % ny;
            for (int x = 0; x < nx; ++x) {
                int next_x = (x + 1) % nx;
                diffx[y * nx + x] = map[y * nx + x] - map[y * nx + next_x];
                diffy[y * nx + x] = map[y * nx + x] - map[next_y * nx + x];
            }
        }

       float sigx = 0.0f, medx = 0.0f;
       ImageProcessing::getSigMed(nx, ny, diffx, sigx, medx);
       for (int i = 0; i < nx * ny; ++i) {
           if (std::abs(diffx[i]) > 8.0f * sigx) {
               weight[i] = 0;
           }
       }

       float sigy = 0.0f, medy = 0.0f;
       ImageProcessing::getSigMed(nx, ny, diffy, sigy, medy);
       for (int i = 0; i < nx * ny; ++i) {
           if (std::abs(diffy[i]) > 8.0f * sigy) {
               weight[i] = 0;
           }
       }

       maskSourceRegions(nx, ny, weight, normap, area_max, defect_halo_thresh * 2.0, area_thresh);

       detectStripes(nx, ny, normap, weight, x_smooth, y_smooth);

       detectArtificialStripes(nx, ny, weight, diffx, diffy, sigx, sigy, medx, medy);

       for (int i = 0; i < nx * ny; ++i) {
           if (weight[i] > 1) {
               weight[i] = 1;
           }
       }

       detectStellarHalo(nx, ny, normap, weight, area_max, defect_halo_thresh);

       detectDent(nx, ny, normap, weight, area_max, defect_halo_thresh);
    }

    // ==========================================
    // Function: Merge defect-connected source regions
    // Method: Traverse connected components in F77 x-major/u-major order and mask oversized regions.
    // ==========================================
   void mergeDefects(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                     int area_max, double source_thresh, int area_thresh, int& ierror) {
       if (ierror == 1) return;

       std::vector<int> mark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] >= source_thresh && weight[i] == 1) {
                mark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (mark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    mark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (mark[nidx] == 1) {
                                    component.push_back(nidx);
                                    mark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(area_max)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (mark[nidx] > 1 || (mark[nidx] == 0 && weight[nidx] == 0)) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            mark[p] = area_max;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            mark[p] = nb;
                        }
                    }
               }
           }
       }
   }

   void detectArtificialStripes(int nx, int ny, std::vector<int>& weight,
                                 const std::vector<float>& diffx, const std::vector<float>& diffy,
                                 float sigx, float sigy, float medx, float medy) {
        std::vector<float> entropy(nx * ny);
        ImageProcessing::getEntropy(nx, ny, diffx, sigx, medx, 2, entropy);
        
        float sig = 0.0f, med = 0.0f;
        ImageProcessing::getSigMed(nx, ny, entropy, sig, med);
        for (int i = 0; i < nx * ny; ++i) {
            if (weight[i] == 1 && std::abs(entropy[i] - med) > 10.0f * sig) {
                weight[i] = 0;
            }
        }

        ImageProcessing::getEntropy(nx, ny, diffy, sigy, medy, 2, entropy);
        ImageProcessing::getSigMed(nx, ny, entropy, sig, med);
        for (int i = 0; i < nx * ny; ++i) {
            if (weight[i] == 1 && std::abs(entropy[i] - med) > 10.0f * sig) {
                weight[i] = 0;
            }
        }
    }

    // ==========================================
    // Function: Mark bright source regions before stripe and halo detection
    // Method: Match F77 mask_source_regions connected-component traversal and area_max handling.
    // ==========================================
    void maskSourceRegions(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                           int area_max, double source_thresh, int area_thresh) {
        std::vector<int> mark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] >= source_thresh && weight[i] == 1) {
                mark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (mark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    mark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (mark[nidx] == 1) {
                                    component.push_back(nidx);
                                    mark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(area_max)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (mark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            mark[p] = area_max;
                            weight[p] = 2;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        if (nb >= area_thresh) {
                            for (int p : component) {
                                mark[p] = nb;
                                weight[p] = 2;
                            }
                        } else {
                            for (int p : component) {
                                mark[p] = nb;
                            }
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect stripe outliers including partial edge blocks
    // Method: Use ceiling block counts and sum/sqrt(nvalid), so full, partial, and DQ-masked
    //         blocks share one noise scale.
    // ==========================================
    void detectStripes(int nx, int ny, const std::vector<float>& normap,
                       std::vector<int>& weight, int x_smooth, int y_smooth) {
        int numy = (ny + y_smooth - 1) / y_smooth;
        std::vector<float> ymap(nx * numy, 0.0f);
        std::vector<unsigned char> yvalid(nx * numy, 0U);
        for (int x = 0; x < nx; ++x) {
            for (int y_block = 0; y_block < numy; ++y_block) {
                int starty = y_block * y_smooth;
                int endy = std::min(starty + y_smooth, ny);
                float sum = 0.0f;
                int nvalid = 0;
                for (int y = starty; y < endy; ++y) {
                    if (weight[y * nx + x] == 1) {
                        sum += normap[y * nx + x];
                        ++nvalid;
                    }
                }
                if (nvalid > 0) {
                    size_t index = static_cast<size_t>(y_block * nx + x);
                    ymap[index] = sum / std::sqrt(static_cast<float>(nvalid));
                    yvalid[index] = 1U;
                }
            }
        }

        float sig_y = 0.0f, med_y = 0.0f;
        getStripeSigMed(ymap, yvalid, sig_y, med_y);
        for (int x = 0; x < nx; ++x) {
            for (int y_block = 0; y_block < numy; ++y_block) {
                size_t index = static_cast<size_t>(y_block * nx + x);
                if (yvalid[index] != 0U && sig_y > 0.0f
                    && std::abs(ymap[index] - med_y) > sig_y * 4.0f) {
                    int starty = y_block * y_smooth;
                    int endy = std::min(starty + y_smooth, ny);
                    for (int y = starty; y < endy; ++y) {
                        weight[y * nx + x] = 0;
                    }
                }
            }
        }

        int numx = (nx + x_smooth - 1) / x_smooth;
        std::vector<float> xmap(numx * ny, 0.0f);
        std::vector<unsigned char> xvalid(numx * ny, 0U);
        for (int y = 0; y < ny; ++y) {
            for (int x_block = 0; x_block < numx; ++x_block) {
                int startx = x_block * x_smooth;
                int endx = std::min(startx + x_smooth, nx);
                float sum = 0.0f;
                int nvalid = 0;
                for (int x = startx; x < endx; ++x) {
                    if (weight[y * nx + x] == 1) {
                        sum += normap[y * nx + x];
                        ++nvalid;
                    }
                }
                if (nvalid > 0) {
                    size_t index = static_cast<size_t>(y * numx + x_block);
                    xmap[index] = sum / std::sqrt(static_cast<float>(nvalid));
                    xvalid[index] = 1U;
                }
            }
        }

        float sig_x = 0.0f, med_x = 0.0f;
        getStripeSigMed(xmap, xvalid, sig_x, med_x);
        for (int y = 0; y < ny; ++y) {
            for (int x_block = 0; x_block < numx; ++x_block) {
                size_t index = static_cast<size_t>(y * numx + x_block);
                if (xvalid[index] != 0U && sig_x > 0.0f
                    && std::abs(xmap[index] - med_x) > sig_x * 4.0f) {
                    int startx = x_block * x_smooth;
                    int endx = std::min(startx + x_smooth, nx);
                    for (int x = startx; x < endx; ++x) {
                        weight[y * nx + x] = 0;
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect stellar halo affected regions
    // Method: Match F77 detect_stellar_halo by smoothing the full normalized map before threshold-connected component masking.
    // ==========================================
    void detectStellarHalo(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                           int npmax, double defect_halo_thresh) {
        std::vector<float> smoothed = normap;
        ImageProcessing::smoothImage55(nx, ny, smoothed, 1);

        std::vector<int> dmark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (smoothed[i] >= defect_halo_thresh) {
                dmark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (dmark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    dmark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (dmark[nidx] == 1) {
                                    component.push_back(nidx);
                                    dmark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(npmax)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (dmark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            dmark[p] = npmax;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            dmark[p] = nb;
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect large negative dent regions
    // Method: Match F77 detect_dent connected-component traversal and area_max masking.
    // ==========================================
    void detectDent(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                    int npmax, double defect_halo_thresh) {
        std::vector<int> dmark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] <= -defect_halo_thresh) {
                dmark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (dmark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    dmark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (dmark[nidx] == 1) {
                                    component.push_back(nidx);
                                    dmark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(npmax)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (dmark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            dmark[p] = npmax;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            dmark[p] = nb;
                        }
                    }
                }
            }
        }
    }

}
