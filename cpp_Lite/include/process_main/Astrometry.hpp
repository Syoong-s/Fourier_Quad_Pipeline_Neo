#ifndef ASTROMETRY_HPP
#define ASTROMETRY_HPP

#include <string>
#include <vector>
#include <array>
#include "FitsIO.hpp"
#include "LensingConfig.hpp"
#include "LinearSolve.hpp"

namespace Astrometry {
    // Stage 2 drivers
    void procAstrometry(int iexpo);
    void chipProcessAstrometry(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput);

    // Coordinate systems and distortion
    void coordinateTransferPU(double& a, double& d, double& x, double& y, int direc,
                              const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                              const double PU[2][LensingConfig::npd], int npd);

    void coordinateTransferSimple(double& a, double& d, double& x, double& y, int direc,
                                  const double cRPIX[2], const double cD[2][2], const double cRVAL[2]);

    void fieldDistortionPU(double x, double y, int npd, const double PU[2][LensingConfig::npd],
                           const double cD[2][2], const double cRPIX[2],
                           double& g1, double& g2, double& cos2, double& sin2, int& parity);

    void mappingPU(double& xx, double& yy, double& xi, double& eta, int npd,
                   const double PU[2][LensingConfig::npd], int direc);

    void raDecToXiEta(double ra, double dec, double& xi, double& eta, double cRVAL1, double cRVAL2);
    double diffra(double ra1, double ra2);
    double sumra(double dra, double ra);

    void xyToXxyy(double x, double y, double& xx, double& yy, const double cRPIX[2], const double cD[2][2]);

    void getRaDecBound(int np, int n, const std::vector<double>& a, const std::vector<double>& d,
                       double& ra, double& dra, double& dec, double& ddec);

    void getRaDecRangeFine(int nx, int ny, double& ra, double dec[2], double& dra,
                           const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                           const double PU[2][LensingConfig::npd], int npd, double astrometryShiftRatio);

    void getRaDecRange(int nx, int ny, double& ra, double dec[2], double& dra,
                       const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                       double astrometryShiftRatio);

    // Calibration & Matching
    void genAstrometryData(const std::string& catStandard, int nx, int ny,
                           const std::vector<float>& map, const std::vector<int>& weight,
                           WCSParams& wcs, const std::string& filename, int& procError);

    void getAstrometryCatalog(int nx, int ny, int npx, int npy, const std::vector<float>& image,
                              const std::vector<int>& weight, int n_user_max, int ns_max,
                              int& ns, std::vector<double>& xs, std::vector<double>& ys);

    void patternMatching(int np0, int n0, const std::vector<double>& x0, const std::vector<double>& y0,
                         int np1, int n1, const std::vector<double>& x1, const std::vector<double>& y1,
                         int shift_range, std::vector<int>& box_final);

    void checkAstrometryGlobal(int np, const std::vector<int>& n, int nc,
                               const std::vector<std::vector<double>>& ra,
                               const std::vector<std::vector<double>>& dec,
                               const std::vector<std::vector<double>>& x,
                               const std::vector<std::vector<double>>& y,
                               std::vector<std::array<double, 2>>& cRPIX,
                               std::vector<std::array<std::array<double, 2>, 2>>& cD,
                               const double cRVAL[2], double PU[2][LensingConfig::npd],
                               int npd, std::vector<int>& valid);

    LinearSolve::SolveStatus measureAstrometryGlobal(
        int np, const std::vector<int>& n, int nc,
        const std::vector<std::vector<double>>& ra,
        const std::vector<std::vector<double>>& dec,
        const std::vector<std::vector<double>>& x,
        const std::vector<std::vector<double>>& y,
        std::vector<std::array<double, 2>>& cRPIX,
        std::vector<std::array<std::array<double, 2>, 2>>& cD,
        const double cRVAL[2], double PU[2][LensingConfig::npd],
        int npd, std::vector<int>& valid,
        LinearSolve::SolveDiagnostics* diagnostics = nullptr);

    // Helpers inside proc_astrometry.f
    void getAstrometry(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput);
    void readAstrometryPara(const std::string& filename, int ichip,
                            double cRPIX[2], double cD[2][2], double cRVAL[2],
                            double PU[2][LensingConfig::npd], int npd, int& procError);
}

#endif // ASTROMETRY_HPP
