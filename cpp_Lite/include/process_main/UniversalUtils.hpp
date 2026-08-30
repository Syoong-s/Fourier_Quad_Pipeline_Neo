#ifndef UNIVERSAL_UTILS_HPP
#define UNIVERSAL_UTILS_HPP

#include <string>
#include <vector>
#include <array>
#include "process_main/LinearSolve.hpp"

struct Point3D {
    double x;
    double y;
    double z;
};

namespace UniversalUtils {
    // Chi-squared utilities
    void getChi2(int ns, int n1, int n2, const std::vector<float>& image, const std::vector<float>& model, double& chi2);
    void anaChi2(int n, const std::vector<float>& map1, const std::vector<float>& map2, double& p);

    // 2D Polynomial functions (Type 1: f(x,y) = sum x^px * y^py where px=term%nx, py=term/nx)
    double fitFunc(double x, double y, int nx, int term_idx);
    double funcVal(double x, double y, int nc, int nx, const std::vector<double>& c);
    void fit2D(const std::vector<Point3D>& points, int nc, int nx, std::vector<double>& c);

    // 2D Polynomial functions (Type 2: standard ordered polynomial terms: 1, x, y, x^2, xy, y^2, ...)
    void getPolynomialPowers(int term_idx, int& px, int& py);
    double fitFunc2(double x, double y, int term_idx);
    double funcVal2(double x, double y, int nc, const std::vector<double>& c);
    LinearSolve::SolveStatus fit2D2(const std::vector<Point3D>& points, int nc,
                                    std::vector<double>& c,
                                    LinearSolve::SolveDiagnostics* diagnostics = nullptr);
    bool invertMatrixF77(std::vector<float> ma, int n, std::vector<float>& ma_inv,
                         bool print_internal_error = true);

    // 2D Polynomial fitting with covariance matrix output
    LinearSolve::SolveStatus fit2D2Cov(const std::vector<Point3D>& points, int nc,
                                       std::vector<double>& c, double& sigma2,
                                       std::vector<double>& cov_matrix,
                                       LinearSolve::SolveDiagnostics* diagnostics = nullptr);

    // 2D Plane fitting
    void findSlope2D(const std::vector<Point3D>& points, double& aa, double& bb, double& cc);

    // 2D Linear Coordinate Fitting (xt = A * x)
    LinearSolve::SolveStatus fitLinear2D(const std::vector<std::array<double, 2>>& x,
                                         const std::vector<std::array<double, 2>>& xt,
                                         double coe[2][3],
                                         LinearSolve::SolveDiagnostics* diagnostics = nullptr);

    // Path & Filename processing
    std::string generateGaiaFileName(const std::string& baseDir, const double cRVAL[2], int& proc_error);
    std::string generateGaiaFileName(const std::string& baseDir, const double cRVAL[2]);
    std::vector<std::string> generateGalCatFileNames(const std::string& baseDir, const double cRVAL[2]);
    std::string getPrefix(const std::string& imagefile);
    std::string getDir(const std::string& imagefile, int level);
    std::string getPrefixExpo(const std::string& imagefile);
    int getChipId(const std::string& imagefile);

    // Median & Sigma computation
    void getMedSig(const std::vector<float>& dat, float& med, float& sig);
    void getMedSig(const std::vector<double>& dat, double& med, double& sig);

    // Math conversions
    double iden(double x, int n);
    double loga(double x, int n);
    double boundRA(double ra_input);

    // Image list reading
    void getImageList(const std::string& expo_file_path, std::vector<std::string>& image_files, std::string& dir_output);
}

#endif // UNIVERSAL_UTILS_HPP
