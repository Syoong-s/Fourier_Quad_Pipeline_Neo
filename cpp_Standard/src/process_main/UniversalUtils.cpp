#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "LinearSolve.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <Eigen/Dense>

namespace UniversalUtils {
    namespace {
        // ==========================================
        // Function: Invert a single-precision matrix
        // Method: Match F77 matrix_inverse using ludcmp/lubksb pivoting rules with optional legacy stderr output.
        // ==========================================
        bool invertMatrixF77Impl(std::vector<float> ma, int n, std::vector<float>& ma_inv,
                                 bool print_internal_error) {
            constexpr float tiny = 1.0e-20f;
            std::vector<int> indx(n, 0);
            std::vector<float> vv(n, 0.0f);
            float d = 1.0f;

            for (int i = 0; i < n; ++i) {
                float aamax = 0.0f;
                for (int j = 0; j < n; ++j) {
                    aamax = std::max(aamax, std::abs(ma[i * n + j]));
                }
                if (aamax == 0.0f) {
                    if (print_internal_error) {
                        std::cerr << "singular matrix in ludcmp" << std::endl;
                    }
                    return false;
                }
                vv[i] = 1.0f / aamax;
            }

            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < j; ++i) {
                    float sum = ma[i * n + j];
                    for (int k = 0; k < i; ++k) {
                        sum -= ma[i * n + k] * ma[k * n + j];
                    }
                    ma[i * n + j] = sum;
                }

                float aamax = 0.0f;
                int imax = j;
                for (int i = j; i < n; ++i) {
                    float sum = ma[i * n + j];
                    for (int k = 0; k < j; ++k) {
                        sum -= ma[i * n + k] * ma[k * n + j];
                    }
                    ma[i * n + j] = sum;
                    float dum = vv[i] * std::abs(sum);
                    if (dum >= aamax) {
                        imax = i;
                        aamax = dum;
                    }
                }

                if (j != imax) {
                    for (int k = 0; k < n; ++k) {
                        std::swap(ma[imax * n + k], ma[j * n + k]);
                    }
                    d = -d;
                    vv[imax] = vv[j];
                }
                indx[j] = imax;
                if (ma[j * n + j] == 0.0f) {
                    ma[j * n + j] = tiny;
                }
                if (j != n - 1) {
                    float dum = 1.0f / ma[j * n + j];
                    for (int i = j + 1; i < n; ++i) {
                        ma[i * n + j] *= dum;
                    }
                }
            }

            ma_inv.assign(static_cast<size_t>(n) * n, 0.0f);
            std::vector<float> b(n, 0.0f);
            for (int col = 0; col < n; ++col) {
                std::fill(b.begin(), b.end(), 0.0f);
                b[col] = 1.0f;

                int ii = -1;
                for (int i = 0; i < n; ++i) {
                    int ll = indx[i];
                    float sum = b[ll];
                    b[ll] = b[i];
                    if (ii >= 0) {
                        for (int j = ii; j < i; ++j) {
                            sum -= ma[i * n + j] * b[j];
                        }
                    } else if (sum != 0.0f) {
                        ii = i;
                    }
                    b[i] = sum;
                }

                for (int i = n - 1; i >= 0; --i) {
                    float sum = b[i];
                    for (int j = i + 1; j < n; ++j) {
                        sum -= ma[i * n + j] * b[j];
                    }
                    b[i] = sum / ma[i * n + i];
                }

                for (int row = 0; row < n; ++row) {
                    ma_inv[row * n + col] = b[row];
                }
            }
            (void)d;
            return true;
        }
    }

    // ==========================================
    // Function: Invert a single-precision matrix with F77 LU rules
    // Method: Preserve legacy output by default while allowing checked callers to emit only the unified Error record.
    // ==========================================
    bool invertMatrixF77(std::vector<float> ma, int n, std::vector<float>& ma_inv,
                         bool print_internal_error) {
        return invertMatrixF77Impl(
            std::move(ma), n, ma_inv, print_internal_error);
    }

    void getChi2(int ns, int n1, int n2, const std::vector<float>& image, const std::vector<float>& model, double& chi2) {
        chi2 = 0.0;
        double cc = ns / 2.0;
        for (int i = n1; i <= n2; ++i) {
            for (int j = n1; j <= n2; ++j) {
                double diff = image[i * ns + j] - model[i * ns + j];
                double dist2 = (i - cc) * (i - cc) + (j - cc) * (j - cc);
                double term = diff * dist2;
                chi2 += term * term;
            }
        }
    }

    void anaChi2(int n, const std::vector<float>& map1, const std::vector<float>& map2, double& p) {
        int n1 = n / 4 - 1;
        int n2 = (n / 4) * 3 - 1;
        if (n1 < 0) n1 = 0;
        if (n2 >= n) n2 = n - 1;
        
        double flux = 0.0;
        double p_sum = 0.0;
        for (int i = n1; i <= n2; ++i) {
            for (int j = n1; j <= n2; ++j) {
                double val1 = map1[i * n + j];
                double val2 = map2[i * n + j];
                flux += (val1 + val2) * 0.5;
                p_sum += (val1 - val2) * (val1 - val2);
            }
        }
        if (flux != 0.0) {
            p = p_sum / flux;
        } else {
            p = 0.0;
        }
    }

    double fitFunc(double x, double y, int nx, int term_idx) {
        int px = term_idx % nx;
        int py = term_idx / nx;
        return std::pow(x, px) * std::pow(y, py);
    }

    double funcVal(double x, double y, int nc, int nx, const std::vector<double>& c) {
        double val = 0.0;
        for (int i = 0; i < nc; ++i) {
            val += fitFunc(x, y, nx, i) * c[i];
        }
        return val;
    }

    // ==========================================
    // Function: Fit 2D polynomial coefficients
    // Method: Factor the raw design matrix with column-scaled rank-revealing QR and solve for the
    //         coefficients directly; zero-fill and report the diagnostics on any non-Normal status.
    // Note:   numerical_fix F3. The legacy version accumulated the normal matrix in float and then
    //         formed an EXPLICIT inverse with the single-precision NR LU before multiplying it into
    //         the right-hand side, deliberately reproducing F77 fit_2D. Squaring the condition
    //         number and inverting are two independent error amplifiers: with eps32*cond(A^T A)
    //         reaching 2.4e2 the fitted surface carried a median error of 2.1 ADU and a p99 of
    //         2.3e2 ADU on the second amplifier. Working on the design matrix keeps cond(A) at its
    //         un-squared value, and LeastSquaresQR column-scales before factoring.
    // ==========================================
    void fit2D(const std::vector<Point3D>& points, int nc, int nx, std::vector<double>& c) {
        c.assign(nc, 0.0);
        int n = static_cast<int>(points.size());
        if (n < nc) return;

        Eigen::MatrixXd design(n, nc);
        Eigen::VectorXd rhs(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < nc; ++j) {
                design(i, j) = fitFunc(points[i].x, points[i].y, nx, j);
            }
            rhs(i) = points[i].z;
        }

        if (!design.allFinite() || !rhs.allFinite()) {
            LinearSolve::reportFailure("UniversalUtils::fit2D", LinearSolve::SolveStatus::FailedSolver);
            return;
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveDiagnostics diagnostics;
        LinearSolve::SolveStatus status = solver.factorize(design, diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure("UniversalUtils::fit2D", status,
                                       LinearSolve::diagnosticsContext(diagnostics));
            return;
        }

        Eigen::VectorXd solution;
        status = solver.solve(rhs, solution);
        if (status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure("UniversalUtils::fit2D", status,
                                       LinearSolve::diagnosticsContext(diagnostics));
            return;
        }

        for (int j = 0; j < nc; ++j) {
            c[j] = solution(j);
        }
    }

    void getPolynomialPowers(int term_idx, int& px, int& py) {
        int nn = term_idx + 1;
        px = 0;
        py = -1;
        int order = -1;
        for (int k = 1; k <= nn; ++k) {
            if (py == order) {
                order++;
                px = order;
                py = 0;
            } else {
                px--;
                py++;
            }
        }
    }

    double fitFunc2(double x, double y, int term_idx) {
        int px, py;
        getPolynomialPowers(term_idx, px, py);
        return std::pow(x, px) * std::pow(y, py);
    }

    double funcVal2(double x, double y, int nc, const std::vector<double>& c) {
        double val = 0.0;
        for (int i = 0; i < nc; ++i) {
            val += fitFunc2(x, y, i) * c[i];
        }
        return val;
    }

    // ==========================================
    // Function: Fit an ordered 2D polynomial with failure status
    // Method: Build the direct design matrix and use one rank-revealing QR without a fallback solver.
    // ==========================================
    LinearSolve::SolveStatus fit2D2(const std::vector<Point3D>& points, int nc,
                                    std::vector<double>& c,
                                    LinearSolve::SolveDiagnostics* diagnostics) {
        c.assign(nc, 0.0);
        int n = static_cast<int>(points.size());
        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = n;
        diag.cols = nc;
        diag.required_rank = nc;
        if (n < nc || nc <= 0) {
            diag.rank = std::min(n, nc);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd A(n, nc);
        Eigen::VectorXd b(n);
        for (int i = 0; i < n; ++i) {
            b(i) = points[i].z;
            for (int j = 0; j < nc; ++j) {
                A(i, j) = fitFunc2(points[i].x, points[i].y, j);
            }
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveStatus status = solver.factorize(A, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        Eigen::VectorXd sol;
        status = solver.solve(b, sol);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }
        for (int j = 0; j < nc; ++j) {
            c[j] = sol(j);
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // ==========================================
    // Function: Fit an ordered 2D polynomial and coefficient covariance
    // Method: Solve with pivoted QR and derive covariance from R without explicitly inverting A^T A.
    // ==========================================
    LinearSolve::SolveStatus fit2D2Cov(const std::vector<Point3D>& points, int nc,
                                       std::vector<double>& c, double& sigma2,
                                       std::vector<double>& cov_matrix,
                                       LinearSolve::SolveDiagnostics* diagnostics) {
        c.assign(nc, 0.0);
        cov_matrix.assign(static_cast<size_t>(nc) * nc, 0.0);
        sigma2 = 0.0;
        int n = static_cast<int>(points.size());
        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = n;
        diag.cols = nc;
        diag.required_rank = nc;
        if (n <= nc || nc <= 0) {
            diag.rank = std::min(n, nc);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd A(n, nc);
        Eigen::VectorXd b(n);
        for (int i = 0; i < n; ++i) {
            b(i) = points[i].z;
            for (int j = 0; j < nc; ++j) {
                A(i, j) = fitFunc2(points[i].x, points[i].y, j);
            }
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveStatus status = solver.factorize(A, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        Eigen::VectorXd sol;
        status = solver.solve(b, sol);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }
        for (int j = 0; j < nc; ++j) {
            c[j] = sol(j);
        }

        Eigen::VectorXd residual = b - A * sol;
        double ssr = residual.dot(residual);
        sigma2 = ssr / static_cast<double>(n - nc);
        if (!residual.allFinite() || !std::isfinite(sigma2)) {
            c.assign(nc, 0.0);
            sigma2 = 0.0;
            return LinearSolve::SolveStatus::FailedSolver;
        }

        Eigen::MatrixXd covariance_base;
        status = solver.unscaledCovariance(covariance_base);
        if (status != LinearSolve::SolveStatus::Normal) {
            c.assign(nc, 0.0);
            sigma2 = 0.0;
            return status;
        }
        covariance_base *= sigma2;
        if (!covariance_base.allFinite()) {
            c.assign(nc, 0.0);
            sigma2 = 0.0;
            return LinearSolve::SolveStatus::FailedSolver;
        }

        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nc; ++j) {
                cov_matrix[i * nc + j] = covariance_base(i, j);
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // ==========================================
    // Function: Fit a 2D plane
    // Method: Factor the [1, x, y] design matrix with column-scaled rank-revealing QR and solve
    //         directly; zero the coefficients and report the diagnostics on any non-Normal status.
    // Note:   numerical_fix F3. Same cure as fit2D. The legacy float normal equations plus an
    //         explicit single-precision inverse destroyed the plane fit whenever the sample sat far
    //         from the origin -- on raw pixel coordinates cond(A^T A) reached 2.3e8, which is what
    //         made the published sigabc differ by 70% in bb and change sign in cc.
    // ==========================================
    void findSlope2D(const std::vector<Point3D>& points, double& aa, double& bb, double& cc) {
        aa = 0.0; bb = 0.0; cc = 0.0;
        int n = static_cast<int>(points.size());
        if (n < 3) {
            return;
        }

        Eigen::MatrixXd design(n, 3);
        Eigen::VectorXd rhs(n);
        for (int i = 0; i < n; ++i) {
            design(i, 0) = 1.0;
            design(i, 1) = points[i].x;
            design(i, 2) = points[i].y;
            rhs(i) = points[i].z;
        }

        if (!design.allFinite() || !rhs.allFinite()) {
            LinearSolve::reportFailure("UniversalUtils::findSlope2D",
                                       LinearSolve::SolveStatus::FailedSolver);
            return;
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveDiagnostics diagnostics;
        LinearSolve::SolveStatus status = solver.factorize(design, diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure("UniversalUtils::findSlope2D", status,
                                       LinearSolve::diagnosticsContext(diagnostics));
            return;
        }

        Eigen::VectorXd solution;
        status = solver.solve(rhs, solution);
        if (status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure("UniversalUtils::findSlope2D", status,
                                       LinearSolve::diagnosticsContext(diagnostics));
            return;
        }

        aa = solution(0);
        bb = solution(1);
        cc = solution(2);
    }

    // ==========================================
    // Function: Fit a two-coordinate affine mapping with failure status
    // Method: Factor the shared three-column design once and solve both coordinate RHS columns together.
    // ==========================================
    LinearSolve::SolveStatus fitLinear2D(const std::vector<std::array<double, 2>>& x,
                                         const std::vector<std::array<double, 2>>& xt,
                                         double coe[2][3],
                                         LinearSolve::SolveDiagnostics* diagnostics) {
        for (int row = 0; row < 2; ++row) {
            for (int col = 0; col < 3; ++col) {
                coe[row][col] = 0.0;
            }
        }
        int n = static_cast<int>(x.size());
        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = n;
        diag.cols = 3;
        diag.required_rank = 3;
        if (n < 3 || xt.size() != x.size()) {
            diag.rank = std::min(n, 3);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd A(n, 3);
        Eigen::MatrixXd rhs(n, 2);
        for (int i = 0; i < n; ++i) {
            A(i, 0) = x[i][0];
            A(i, 1) = x[i][1];
            A(i, 2) = 1.0;
            rhs(i, 0) = xt[i][0];
            rhs(i, 1) = xt[i][1];
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveStatus status = solver.factorize(A, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        Eigen::MatrixXd solution;
        status = solver.solve(rhs, solution);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        coe[0][0] = solution(0, 0); coe[0][1] = solution(1, 0); coe[0][2] = solution(2, 0);
        coe[1][0] = solution(0, 1); coe[1][1] = solution(1, 1); coe[1][2] = solution(2, 1);
        return LinearSolve::SolveStatus::Normal;
    }

    std::string generateGaiaFileName(const std::string& baseDir, const double cRVAL[2], int& proc_error) {
        if (proc_error == 1) return "";
        if (cRVAL[0] < 0) {
            proc_error = 1;
            return "";
        }
        double ra_val = cRVAL[0];
        double dec_val = cRVAL[1];
        int dec = std::min(static_cast<int>(std::abs(dec_val) / 10.0) + 1, 9);
        int ra = std::min(static_cast<int>(ra_val / 10.0), 35);

        std::ostringstream oss;
        oss << baseDir;
        if (dec_val >= 0.0) {
            oss << "/gaia_p";
        } else {
            oss << "/gaia_m";
        }
        oss << dec;
        if (dec == 9) {
            oss << ".cat";
        } else {
            oss << "_" << std::setw(2) << std::setfill('0') << ra << ".cat";
        }
        return oss.str();
    }

    std::string generateGaiaFileName(const std::string& baseDir, const double cRVAL[2]) {
        int dummy_proc_error = 0;
        return generateGaiaFileName(baseDir, cRVAL, dummy_proc_error);
    }

    std::vector<std::string> generateGalCatFileNames(const std::string& baseDir, const double cRVAL[2]) {
        double ra_val = cRVAL[0];
        double dec_val = cRVAL[1];
        double m_dec = 1.0;
        double m_ra = 1.0;
        double abs_dec = std::abs(dec_val);
        if (abs_dec < 30.0) {
            m_ra = 1.4;
        } else if (abs_dec < 40.0) {
            m_ra = 1.6;
        } else if (abs_dec < 50.0) {
            m_ra = 1.8;
        } else if (abs_dec < 60.0) {
            m_ra = 2.5;
        } else {
            m_ra = 3.0;
        }

        int dec1 = static_cast<int>(std::floor(dec_val - m_dec));
        int dec2 = static_cast<int>(std::floor(dec_val + m_dec));
        int ra1 = static_cast<int>(std::floor(ra_val - m_ra));
        int ra2 = static_cast<int>(std::floor(ra_val + m_ra));

        std::vector<std::string> filenames;
        std::string prefix = baseDir + "/des_y6_";

        for (int dec = dec1; dec <= dec2; ++dec) {
            std::string c_dec1, c_dec2;
            {
                std::ostringstream oss;
                if (dec >= 0) {
                    oss << "p" << std::setw(2) << std::setfill('0') << dec;
                } else {
                    oss << "m" << std::setw(2) << std::setfill('0') << std::abs(dec);
                }
                c_dec1 = oss.str();
            }
            {
                std::ostringstream oss;
                int dec_next = dec + 1;
                if (dec_next >= 0) {
                    oss << "p" << std::setw(2) << std::setfill('0') << dec_next;
                } else {
                    oss << "m" << std::setw(2) << std::setfill('0') << std::abs(dec_next);
                }
                c_dec2 = oss.str();
            }

            for (int ra = ra1; ra <= ra2; ++ra) {
                int mra = ra;
                if (ra < 0) {
                    mra += 360;
                } else if (ra >= 360) {
                    mra -= 360;
                }

                std::string c_ra1, c_ra2;
                {
                    std::ostringstream oss;
                    oss << std::setw(3) << std::setfill('0') << mra;
                    c_ra1 = oss.str();
                }
                {
                    std::ostringstream oss;
                    int mra_next = mra + 1;
                    oss << std::setw(3) << std::setfill('0') << mra_next;
                    c_ra2 = oss.str();
                }

                std::string fname = prefix + "RA_" + c_ra1 + "_" + c_ra2 + "_Dec_" + c_dec1 + "_" + c_dec2 + ".dat";
                filenames.push_back(fname);
            }
        }
        return filenames;
    }

    // ==========================================
    // Function: Extract chip-level file prefix
    // Method: Match F77 get_PREFIX and stop immediately on malformed paths.
    // ==========================================
    std::string getPrefix(const std::string& imagefile) {
        size_t p_slash = imagefile.find_last_of('/');
        size_t p_dot = imagefile.find_last_of('.');
        if (p_dot == std::string::npos || p_slash == std::string::npos || p_dot <= p_slash + 1) {
            std::cerr << "Image_file name is NOT normal: " << imagefile << std::endl;
            std::exit(1);
        }
        return imagefile.substr(p_slash + 1, p_dot - p_slash - 1);
    }

    // ==========================================
    // Function: Extract parent directory by slash depth
    // Method: Match F77 get_dir and stop when slash count is insufficient.
    // ==========================================
    std::string getDir(const std::string& imagefile, int level) {
        int u = 0;
        size_t i = imagefile.length();
        while (i > 0) {
            --i;
            if (imagefile[i] == '/') {
                u++;
                if (u == level) {
                    return imagefile.substr(0, i);
                }
            }
        }
        std::cerr << "Image_file name is NOT normal for level " << level << ": " << imagefile << std::endl;
        std::exit(1);
    }

    // ==========================================
    // Function: Extract exposure-level file prefix
    // Method: Match F77 get_PREFIX_expo and stop immediately on malformed paths.
    // ==========================================
    std::string getPrefixExpo(const std::string& imagefile) {
        size_t p_slash = imagefile.find_last_of('/');
        size_t p_under = imagefile.find_last_of('_');
        if (p_under == std::string::npos || p_slash == std::string::npos || p_under <= p_slash + 1) {
            std::cerr << "Image_file name is NOT normal: " << imagefile << std::endl;
            std::exit(1);
        }
        return imagefile.substr(p_slash + 1, p_under - p_slash - 1);
    }

    int getChipId(const std::string& imagefile) {
        int id = 0;
        if (!FitsIO::readCCDNUM(imagefile, id)) {
            id = -99;
            std::cout << "Error / GetChipId Failed" << std::endl;
        }
        return id;
    }

    void getMedSig(const std::vector<float>& dat, float& med, float& sig) {
        int n = static_cast<int>(dat.size());
        if (n == 0) {
            med = 0.0f; sig = 0.0f;
            return;
        }
        std::vector<float> arr = dat;
        std::sort(arr.begin(), arr.end());
        int idx_med = std::max(0, n / 2 - 1);
        int idx_high = std::max(0, n * 5 / 6 - 1);
        int idx_low = std::max(0, n / 6 - 1);
        med = arr[idx_med];
        sig = 0.5f * (arr[idx_high] - arr[idx_low]);
    }

    void getMedSig(const std::vector<double>& dat, double& med, double& sig) {
        int n = static_cast<int>(dat.size());
        if (n == 0) {
            med = 0.0; sig = 0.0;
            return;
        }
        std::vector<double> arr = dat;
        std::sort(arr.begin(), arr.end());
        int idx_med = std::max(0, n / 2 - 1);
        int idx_high = std::max(0, n * 5 / 6 - 1);
        int idx_low = std::max(0, n / 6 - 1);
        med = arr[idx_med];
        sig = 0.5 * (arr[idx_high] - arr[idx_low]);
    }

    double iden(double x, int n) {
        return x;
    }

    double loga(double x, int n) {
        if (n == 1) {
            if (x < -30.0) {
                return -std::log(-2.0 * x);
            } else {
                return std::log(x + std::sqrt(x * x + 1.0));
            }
        } else {
            return 0.5 * (std::exp(x) - std::exp(-x));
        }
    }

    double boundRA(double ra_input) {
        double res = ra_input;
        if (res >= 360.0) {
            res -= 360.0;
        } else if (res < 0.0) {
            res += 360.0;
        }
        return res;
    }

    void getImageList(const std::string& expo_file_path, std::vector<std::string>& image_files, std::string& dir_output) {
        image_files.clear();
        std::ifstream infile(expo_file_path);
        if (!infile.is_open()) {
            std::cerr << "EXPO_FILE reading error: " << expo_file_path << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(infile, line)) {
           // 1. 去除行尾的 \r, 空格, 或制表符 \t
            while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
                line.pop_back();
            }
        
            // 只有非空行才放进去
            if (!line.empty()) {
                image_files.push_back(line);
            }
            // Original Version：
            // image_files.push_back(line);
        }
        infile.close();

        if (image_files.empty()) {
            std::cerr << "EXPO_FILE contains no image entries: " << expo_file_path << std::endl;
            std::exit(1);
        }

        dir_output = getDir(image_files[0], 3);
    }
}
