#include "Astrometry.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "LinearSolve.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <utility>

// Global/extern variables representing exposure filenames list (defined in main.cpp)
extern std::vector<std::string> EXPO_FILE;

namespace Astrometry {

    // ==========================================
    // Function: Emulate Fortran INT(value + 0.5)
    // Method: Fortran INT truncates toward zero, unlike floor for negative values.
    // ==========================================
    int fortranShiftBin(double value) {
        return static_cast<int>(value + 0.5);
    }

    double diffra(double ra1, double ra2) {
        double diff = ra1 - ra2;
        if (diff < -180.0) {
            diff += 360.0;
        } else if (diff > 180.0) {
            diff -= 360.0;
        }
        return diff;
    }

    double sumra(double dra, double ra) {
        double sum = ra + dra;
        if (sum >= 360.0) {
            sum -= 360.0;
        } else if (sum < 0.0) {
            sum += 360.0;
        }
        return sum;
    }

    void xyToXxyy(double x, double y, double& xx, double& yy, const double cRPIX[2], const double cD[2][2]) {
        xx = cD[0][0] * (x - cRPIX[0]) + cD[0][1] * (y - cRPIX[1]);
        yy = cD[1][0] * (x - cRPIX[0]) + cD[1][1] * (y - cRPIX[1]);
    }

    void raDecToXiEta(double ra, double dec, double& xi, double& eta, double cRVAL1, double cRVAL2) {
        constexpr double pi = 3.1415926;
        double const1 = pi / 180.0;
        double tandc = std::tan(cRVAL2 * const1);

        double da = diffra(ra, cRVAL1) * const1;
        double dd = (dec - cRVAL2) * const1;
        double tandd = std::tan(dd);
        double cosda = std::cos(da);
        double tanda = std::tan(da);

        double num_y = tandc * (cosda - 1.0) - (1.0 + cosda * tandc * tandc) * tandd;
        double den_y = tandd * tandc * (cosda - 1.0) - (cosda + tandc * tandc);
        double y = num_y / den_y;
        double x = tanda * (std::cos(cRVAL2 * const1) * (1.0 - y * tandc));

        xi = x / const1;
        eta = y / const1;
    }

    void mappingPU(double& xx, double& yy, double& xi, double& eta, int npd,
                   const double PU[2][LensingConfig::npd], int direc) {
        if (direc == 1) {
            xi = xx;
            eta = yy;
            for (int i = 0; i < 3; ++i) {
                double dxi = 0.0;
                double deta = 0.0;
                int px = 0;
                int py = 1;
                int order = 1;
                int n = 0;
                while (n < npd) {
                    if (py == order) {
                        order++;
                        px = order;
                        py = 0;
                    } else {
                        px--;
                        py++;
                    }
                    n++;
                    dxi += PU[0][n - 1] * std::pow(xi, px) * std::pow(eta, py);
                    deta += PU[1][n - 1] * std::pow(eta, px) * std::pow(xi, py);
                }
                xi = xx + dxi;
                eta = yy + deta;
            }
        } else {
            xx = xi;
            yy = eta;
            int px = 0;
            int py = 1;
            int order = 1;
            int n = 0;
            while (n < npd) {
                if (py == order) {
                    order++;
                    px = order;
                    py = 0;
                } else {
                    px--;
                    py++;
                }
                n++;
                xx -= PU[0][n - 1] * std::pow(xi, px) * std::pow(eta, py);
                yy -= PU[1][n - 1] * std::pow(eta, px) * std::pow(xi, py);
            }
        }
    }

    // ==========================================
    // Function: Convert between pixel and sky coordinates with PU distortion.
    // Method: Mirrors F77 coordinate_transfer_PU; x/y are in-out for direc=-1.
    // ==========================================
    void coordinateTransferPU(double& a, double& d, double& x, double& y, int direc,
                              const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                              const double PU[2][LensingConfig::npd], int npd) {
        constexpr double pi = 3.1415926;
        double const1 = pi / 180.0;
        double tandc = std::tan(cRVAL[1] * const1);

        if (direc == 1) {
            double xx = cD[0][0] * (x - cRPIX[0]) + cD[0][1] * (y - cRPIX[1]);
            double yy = cD[1][0] * (x - cRPIX[0]) + cD[1][1] * (y - cRPIX[1]);

            double xi = 0.0, eta = 0.0;
            mappingPU(xx, yy, xi, eta, npd, PU, 1);

            double xxx = xi * const1;
            double yyy = eta * const1;

            double da = xxx / (std::cos(cRVAL[1] * const1) * (1.0 - yyy * tandc));
            da = da - da * da * da * 0.3333333333333333;
            a = sumra(da / const1, cRVAL[0]);
            double cosda = 1.0 - da * da * 0.5 + da * da * da * da / 24.0;

            double dd = (yyy * (cosda + tandc * tandc) + tandc * (cosda - 1.0)) /
                        (yyy * tandc * (cosda - 1.0) + 1.0 + cosda * tandc * tandc);
            dd = dd - dd * dd * dd * 0.3333333333333333;
            d = dd / const1 + cRVAL[1];
        } else {
            double da = diffra(a, cRVAL[0]) * const1;
            double dd = (d - cRVAL[1]) * const1;
            double tandd = std::tan(dd);
            double cosda = std::cos(da);
            double tanda = std::tan(da);

            double yy = (tandc * (cosda - 1.0) - (1.0 + cosda * tandc * tandc) * tandd) /
                        (tandd * tandc * (cosda - 1.0) - (cosda + tandc * tandc));
            double xx = tanda * (std::cos(cRVAL[1] * const1) * (1.0 - yy * tandc));

            double xi = xx / const1;
            double eta = yy / const1;

            mappingPU(xx, yy, xi, eta, npd, PU, 2);

            double temp = cD[0][0] * cD[1][1] - cD[0][1] * cD[1][0];
            temp = 1.0 / temp;

            double cD_1[2][2];
            cD_1[0][0] = cD[1][1] * temp;
            cD_1[1][1] = cD[0][0] * temp;
            cD_1[0][1] = -cD[0][1] * temp;
            cD_1[1][0] = -cD[1][0] * temp;

            x = xx * cD_1[0][0] + yy * cD_1[0][1] + cRPIX[0];
            y = xx * cD_1[1][0] + yy * cD_1[1][1] + cRPIX[1];
        }
    }

    // ==========================================
    // Function: Convert between pixel and sky coordinates without PU distortion.
    // Method: Mirrors F77 coordinate_transfer_simple; x/y are in-out for direc=-1.
    // ==========================================
    void coordinateTransferSimple(double& a, double& d, double& x, double& y, int direc,
                                  const double cRPIX[2], const double cD[2][2], const double cRVAL[2]) {
        constexpr double pi = 3.1415926;
        double const1 = pi / 180.0;
        double tandc = std::tan(cRVAL[1] * const1);

        if (direc == 1) {
            double xx = cD[0][0] * (x - cRPIX[0]) + cD[0][1] * (y - cRPIX[1]);
            double yy = cD[1][0] * (x - cRPIX[0]) + cD[1][1] * (y - cRPIX[1]);

            double xxx = xx * const1;
            double yyy = yy * const1;

            double da = xxx / (std::cos(cRVAL[1] * const1) * (1.0 - yyy * tandc));
            da = da - da * da * da * 0.3333333333333333;
            a = sumra(da / const1, cRVAL[0]);
            double cosda = 1.0 - da * da * 0.5 + da * da * da * da / 24.0;

            double dd = (yyy * (cosda + tandc * tandc) + tandc * (cosda - 1.0)) /
                        (yyy * tandc * (cosda - 1.0) + 1.0 + cosda * tandc * tandc);
            dd = dd - dd * dd * dd * 0.3333333333333333;
            d = dd / const1 + cRVAL[1];
        } else {
            double da = diffra(a, cRVAL[0]) * const1;
            double dd = (d - cRVAL[1]) * const1;
            double tandd = std::tan(dd);
            double cosda = std::cos(da);
            double tanda = std::tan(da);

            double yy = (tandc * (cosda - 1.0) - (1.0 + cosda * tandc * tandc) * tandd) /
                        (tandd * tandc * (cosda - 1.0) - (cosda + tandc * tandc));
            double xx = tanda * (std::cos(cRVAL[1] * const1) * (1.0 - yy * tandc));

            xx /= const1;
            yy /= const1;

            double temp = cD[0][0] * cD[1][1] - cD[0][1] * cD[1][0];
            temp = 1.0 / temp;

            double cD_1[2][2];
            cD_1[0][0] = cD[1][1] * temp;
            cD_1[1][1] = cD[0][0] * temp;
            cD_1[0][1] = -cD[0][1] * temp;
            cD_1[1][0] = -cD[1][0] * temp;

            x = xx * cD_1[0][0] + yy * cD_1[0][1] + cRPIX[0];
            y = xx * cD_1[1][0] + yy * cD_1[1][1] + cRPIX[1];
        }
    }

    void fieldDistortionPU(double x, double y, int npd, const double PU[2][LensingConfig::npd],
                           const double cD[2][2], const double cRPIX[2],
                           double& g1, double& g2, double& cos2, double& sin2, int& parity) {
        double xx = cD[0][0] * (x - cRPIX[0]) + cD[0][1] * (y - cRPIX[1]);
        double yy = cD[1][0] * (x - cRPIX[0]) + cD[1][1] * (y - cRPIX[1]);

        double xi = 0.0, eta = 0.0;
        mappingPU(xx, yy, xi, eta, npd, PU, 1);

        double temp = 1.0 / (cD[0][0] * cD[1][1] - cD[0][1] * cD[1][0]);

        double dx_dxx = cD[1][1] * temp;
        double dx_dyy = -cD[0][1] * temp;
        double dy_dxx = -cD[1][0] * temp;
        double dy_dyy = cD[0][0] * temp;

        double dxx_dxi = 1.0;
        double dyy_deta = 1.0;
        double dxx_deta = 0.0;
        double dyy_dxi = 0.0;

        int px = 0;
        int py = 1;
        int order = 1;
        int n = 0;
        while (n < npd) {
            if (py == order) {
                order++;
                px = order;
                py = 0;
            } else {
                px--;
                py++;
            }
            n++;

            // Crucial Pitfalls Guard: guard terms when px or py are 0 to prevent std::pow error
            double term1 = (px > 0) ? PU[0][n - 1] * px * std::pow(xi, px - 1) * std::pow(eta, py) : 0.0;
            double term2 = (py > 0) ? PU[0][n - 1] * std::pow(xi, px) * py * std::pow(eta, py - 1) : 0.0;
            double term3 = (px > 0) ? PU[1][n - 1] * px * std::pow(eta, px - 1) * std::pow(xi, py) : 0.0;
            double term4 = (py > 0) ? PU[1][n - 1] * std::pow(eta, px) * py * std::pow(xi, py - 1) : 0.0;

            dxx_dxi -= term1;
            dxx_deta -= term2;
            dyy_deta -= term3;
            dyy_dxi -= term4;
        }

        double mat[2][2];
        mat[0][0] = dx_dxx * dxx_dxi + dx_dyy * dyy_dxi;
        mat[0][1] = dx_dxx * dxx_deta + dx_dyy * dyy_deta;
        mat[1][0] = dy_dxx * dxx_dxi + dy_dyy * dyy_dxi;
        mat[1][1] = dy_dxx * dxx_deta + dy_dyy * dyy_deta;

        double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        double temp2 = 1.0 / det;

        double dm[2][2];
        dm[0][0] = mat[1][1] * temp2;
        dm[0][1] = -mat[0][1] * temp2;
        dm[1][0] = -mat[1][0] * temp2;
        dm[1][1] = mat[0][0] * temp2;

        parity = 1;
        if (det < 0) {
            dm[0][0] = -dm[0][0];
            dm[0][1] = -dm[0][1];
            parity = -1;
        }

        double sqrt_det = std::sqrt(dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0]);
        dm[0][0] /= sqrt_det;
        dm[0][1] /= sqrt_det;
        dm[1][0] /= sqrt_det;
        dm[1][1] /= sqrt_det;

        double cos1 = 0.5 * (dm[0][0] + dm[1][1]);
        double sin1 = 0.5 * (dm[0][1] - dm[1][0]);

        double aa_val = -0.5 * (dm[0][1] + dm[1][0]);
        double bb_val = 0.5 * (dm[1][1] - dm[0][0]);

        g1 = aa_val * sin1 + bb_val * cos1;
        g2 = aa_val * cos1 - bb_val * sin1;

        cos2 = cos1 * cos1 - sin1 * sin1;
        sin2 = 2.0 * sin1 * cos1;

        if (parity == -1) {
            g2 = -g2;
        }
    }

    void getRaDecBound(int np, int n, const std::vector<double>& a, const std::vector<double>& d,
                       double& ra, double& dra, double& dec, double& ddec) {
        if (n <= 0) return;
        double dec1 = d[0];
        double dec2 = d[0];
        for (int i = 0; i < n; ++i) {
            dec1 = std::min(dec1, d[i]);
            dec2 = std::max(dec2, d[i]);
        }
        dec = (dec1 + dec2) * 0.5;
        ddec = dec2 - dec1;

        double ra1 = a[0];
        double ra2 = a[0];
        for (int i = 0; i < n; ++i) {
            ra1 = std::min(ra1, a[i]);
            ra2 = std::max(ra2, a[i]);
        }

        ra = 0.5 * (ra1 + ra2);
        dra = ra2 - ra1;

        if (dra > 180.0) {
            std::vector<double> tmp(n);
            for (int i = 0; i < n; ++i) {
                if (a[i] < 180.0) {
                    tmp[i] = a[i] + 360.0;
                } else {
                    tmp[i] = a[i];
                }
            }

            ra1 = tmp[0];
            ra2 = tmp[0];
            for (int i = 0; i < n; ++i) {
                ra1 = std::min(ra1, tmp[i]);
                ra2 = std::max(ra2, tmp[i]);
            }

            dra = ra2 - ra1;
            ra = 0.5 * (ra1 + ra2);
            if (ra >= 360.0) {
                ra -= 360.0;
            }
        }
    }

    void getRaDecRangeFine(int nx, int ny, double& ra, double dec[2], double& dra,
                           const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                           const double PU[2][LensingConfig::npd], int npd, double astrometryShiftRatio) {
        std::vector<double> a(4), d(4);
        double xcorn[4] = {1.0, static_cast<double>(nx), static_cast<double>(nx), 1.0};
        double ycorn[4] = {1.0, 1.0, static_cast<double>(ny), static_cast<double>(ny)};
        for (int i = 0; i < 4; ++i) {
            coordinateTransferPU(a[i], d[i], xcorn[i], ycorn[i], 1, cRPIX, cD, cRVAL, PU, npd);
        }

        double ra_c = 0.0, dec_c = 0.0, ddec = 0.0;
        getRaDecBound(4, 4, a, d, ra_c, dra, dec_c, ddec);

        dec[0] = dec_c - 0.5 * ddec;
        dec[1] = dec_c + 0.5 * ddec;

        double shift_ddec = ddec * astrometryShiftRatio;
        dec[0] -= shift_ddec;
        dec[1] += shift_ddec;

        dra = dra * (1.0 + 2.0 * astrometryShiftRatio);
        ra = ra_c;
    }

    void getRaDecRange(int nx, int ny, double& ra, double dec[2], double& dra,
                       const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                       double astrometryShiftRatio) {
        std::vector<double> a(4), d(4);
        double xcorn[4] = {1.0, static_cast<double>(nx), static_cast<double>(nx), 1.0};
        double ycorn[4] = {1.0, 1.0, static_cast<double>(ny), static_cast<double>(ny)};
        for (int i = 0; i < 4; ++i) {
            coordinateTransferSimple(a[i], d[i], xcorn[i], ycorn[i], 1, cRPIX, cD, cRVAL);
        }

        double ra_c = 0.0, dec_c = 0.0, ddec = 0.0;
        getRaDecBound(4, 4, a, d, ra_c, dra, dec_c, ddec);

        dec[0] = dec_c - 0.5 * ddec;
        dec[1] = dec_c + 0.5 * ddec;

        double shift_ddec = ddec * astrometryShiftRatio;
        dec[0] -= shift_ddec;
        dec[1] += shift_ddec;

        dra = dra * (1.0 + 2.0 * astrometryShiftRatio);
        ra = ra_c;
    }

    void getAstrometryCatalog(int nx, int ny, int npx, int npy, const std::vector<float>& image,
                              const std::vector<int>& weight, int n_user_max, int ns_max,
                              int& ns, std::vector<double>& xs, std::vector<double>& ys) {
        std::vector<int> mark(nx * ny, 0);
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int idx = y * nx + x;
                if (image[idx] >= 5.0f && weight[idx] > 0) {
                    mark[idx] = 1;
                } else {
                    mark[idx] = 0;
                }
            }
        }

        int nsb = 0;
        std::vector<double> xsb;
        std::vector<double> ysb;
        std::vector<float> flux_array;
        std::vector<int> order;

        constexpr int sizemax = 10000;
        constexpr int area_limit = 400;

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int start_idx = j * nx + i;
                if (mark[start_idx] == 1) {
                    int nbb = 0;
                    std::vector<std::pair<int, int>> buffer;
                    buffer.reserve(sizemax);
                    buffer.push_back({i, j});
                    mark[start_idx] = -1;

                    int xp = i;
                    int yp = j;
                    float peak = image[start_idx];
                    float flux = 0.0f;
                    int toobig = 0;

                    while (static_cast<int>(buffer.size()) > nbb) {
                        int k1 = nbb;
                        int k2 = buffer.size();
                        nbb = k2;
                        for (int k = k1; k < k2; ++k) {
                            int ix = buffer[k].first;
                            int iy = buffer[k].second;

                            int u_min = std::max(ix - 3, 0);
                            int u_max = std::min(ix + 3, nx - 1);
                            int v_min = std::max(iy - 3, 0);
                            int v_max = std::min(iy + 3, ny - 1);

                            for (int u = u_min; u <= u_max; ++u) {
                                for (int v = v_min; v <= v_max; ++v) {
                                    int nidx = v * nx + u;
                                    if (mark[nidx] == 1) {
                                        buffer.push_back({u, v});
                                        mark[nidx] = -1;
                                        flux += image[nidx];
                                        if (image[nidx] > peak) {
                                            peak = image[nidx];
                                            xp = u;
                                            yp = v;
                                        }
                                        if (buffer.size() == sizemax) {
                                            toobig = 1;
                                            goto label20;
                                        }
                                    } else if (mark[nidx] == 2) {
                                        toobig = 1;
                                        goto label20;
                                    }
                                }
                            }
                        }
                    }

                label20:
                    if (toobig == 1) {
                        for (const auto& pt : buffer) {
                            mark[pt.second * nx + pt.first] = 2;
                        }
                    } else {
                        int area = 0;
                        flux = 0.0f;
                        double xc = 0.0;
                        double yc = 0.0;
                        float temp = peak * 0.5f;
                        for (const auto& pt : buffer) {
                            int idx = pt.second * nx + pt.first;
                            if (image[idx] >= temp) {
                                area++;
                                flux += image[idx];
                                xc += image[idx] * (pt.first + 1.0); // 1-based coordinates
                                yc += image[idx] * (pt.second + 1.0); // 1-based coordinates
                            }
                        }
                        xc /= flux;
                        yc /= flux;

                        if (area <= area_limit) {
                            nsb++;
                            xsb.push_back(xc);
                            ysb.push_back(yc);
                            flux_array.push_back(flux);
                            order.push_back(nsb - 1);
                            if (nsb == ns_max) {
                                ns = 0;
                                return;
                            }
                        }
                    }
                }
            }
        }

        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return flux_array[a] > flux_array[b];
        });

        ns = 0;
        int limit = std::min(n_user_max, nsb);
        xs.resize(limit);
        ys.resize(limit);
        for (int i = 0; i < limit; ++i) {
            xs[i] = xsb[order[i]];
            ys[i] = ysb[order[i]];
            ns++;
        }
    }

    // ==========================================
    // Function: Match two astrometric point catalogs
    // Method: Exclude non-finite samples, preserve the F77 shift search, then run one checked affine refinement.
    // ==========================================
    void patternMatching(int np0, int n0, const std::vector<double>& x0, const std::vector<double>& y0,
                         int np1, int n1, const std::vector<double>& x1, const std::vector<double>& y1,
                         int shift_range, std::vector<int>& box_final) {
        std::vector<unsigned char> finite0(n0, 0);
        std::vector<unsigned char> finite1(n1, 0);
        std::size_t removed_samples = 0;
        for (int i = 0; i < n0; ++i) {
            finite0[i] = static_cast<unsigned char>(std::isfinite(x0[i]) && std::isfinite(y0[i]));
            if (finite0[i] == 0) {
                removed_samples++;
            }
        }
        for (int i = 0; i < n1; ++i) {
            finite1[i] = static_cast<unsigned char>(std::isfinite(x1[i]) && std::isfinite(y1[i]));
            if (finite1[i] == 0) {
                removed_samples++;
            }
        }
        int nn = shift_range;
        int dim = 2 * nn + 1;
        std::vector<float> mark1(dim * dim, 0.0f);
        std::vector<float> mark2(dim * dim, 0.0f);
        std::vector<int> mark3(dim * dim, 0);

        std::vector<std::vector<int>> box(n0);

        for (int i = 0; i < n0; ++i) {
            if (finite0[i] == 0) continue;
            for (int j = 0; j < n1; ++j) {
                if (finite1[j] == 0) continue;
                int dx = fortranShiftBin(x1[j] - x0[i]);
                int dy = fortranShiftBin(y1[j] - y0[i]);
                if (dx >= -nn && dx <= nn && dy >= -nn && dy <= nn) {
                    mark1[(dy + nn) * dim + (dx + nn)] += 1.0f;
                    box[i].push_back(j);
                }
            }
        }

        mark2 = mark1;

        constexpr int radius = 40;
        if (radius > 0) {
            double const_val = 1.0 / ((radius * 0.25) * (radius * 0.25));
            for (int dx_i = -nn; dx_i <= nn; ++dx_i) {
                for (int dy_i = -nn; dy_i <= nn; ++dy_i) {
                    float m1 = mark1[(dy_i + nn) * dim + (dx_i + nn)];
                    if (m1 > 0.0f) {
                        int u_min = std::max(dx_i - radius, -nn);
                        int u_max = std::min(dx_i + radius, nn);
                        int v_min = std::max(dy_i - radius, -nn);
                        int v_max = std::min(dy_i + radius, nn);
                        for (int u = u_min; u <= u_max; ++u) {
                            for (int v = v_min; v <= v_max; ++v) {
                                double dist_sq = (dx_i - u) * (dx_i - u) + (dy_i - v) * (dy_i - v);
                                mark2[(v + nn) * dim + (u + nn)] += m1 / (dist_sq * const_val + 1.0);
                            }
                        }
                    }
                }
            }
        }

        float peak = 0.0f;
        int xp = 0;
        int yp = 0;
        std::vector<float> arr;
        arr.reserve(dim * dim);
        for (int dx_i = -nn + radius; dx_i <= nn - radius; ++dx_i) {
            for (int dy_i = -nn + radius; dy_i <= nn - radius; ++dy_i) {
                float val = mark2[(dy_i + nn) * dim + (dx_i + nn)];
                arr.push_back(val);
                if (val > peak) {
                    peak = val;
                    xp = dx_i;
                    yp = dy_i;
                }
            }
        }

        std::sort(arr.begin(), arr.end());
        size_t median_idx = arr.size() / 2;
        if (median_idx > 0) {
            --median_idx;
        }
        float thresh = (arr[median_idx] + peak) * 0.5f;

        constexpr int r_fof = 1;
        int changed = 1;
        int xmin = xp;
        int xmax = xp;
        int ymin = yp;
        int ymax = yp;
        mark3[(yp + nn) * dim + (xp + nn)] = 1;

        while (changed == 1) {
            changed = 0;
            int current_xmin = xmin;
            int current_xmax = xmax;
            int current_ymin = ymin;
            int current_ymax = ymax;
            for (int i = current_xmin; i <= current_xmax; ++i) {
                for (int j = current_ymin; j <= current_ymax; ++j) {
                    if (mark3[(j + nn) * dim + (i + nn)] == 1) {
                        for (int u = i - r_fof; u <= i + r_fof; ++u) {
                            for (int v = j - r_fof; v <= j + r_fof; ++v) {
                                if (u >= -nn && u <= nn && v >= -nn && v <= nn) {
                                    int idx_uv = (v + nn) * dim + (u + nn);
                                    if (mark3[idx_uv] == 0 && mark2[idx_uv] > thresh) {
                                        mark3[idx_uv] = 1;
                                        xmin = std::min(xmin, u);
                                        xmax = std::max(xmax, u);
                                        ymin = std::min(ymin, v);
                                        ymax = std::max(ymax, v);
                                        changed = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        box_final.assign(n0, 0);
        for (int i = 0; i < n0; ++i) {
            if (finite0[i] == 0) continue;
            for (int k : box[i]) {
                int dx = fortranShiftBin(x1[k] - x0[i]);
                int dy = fortranShiftBin(y1[k] - y0[i]);
                if (dx >= -nn && dx <= nn && dy >= -nn && dy <= nn) {
                    if (mark3[(dy + nn) * dim + (dx + nn)] == 1) {
                        if (box_final[i] == 0) {
                            box_final[i] = k + 1; // 1-based index
                        } else {
                            box_final[i] = 0;
                            break;
                        }
                    }
                }
            }
        }

        int n_match = 0;
        std::vector<std::array<double, 2>> xx_fit;
        std::vector<std::array<double, 2>> xxt_fit;

        for (int i = 0; i < n0; ++i) {
            if (box_final[i] == 0) continue;
            int dpl = 0;
            for (int j = i + 1; j < n0; ++j) {
                if (box_final[j] == 0) continue;
                if (box_final[i] == box_final[j]) {
                    box_final[j] = 0;
                    dpl = 1;
                }
            }
            if (dpl == 1) {
                box_final[i] = 0;
                continue;
            }
            n_match++;
            xx_fit.push_back({x0[i], y0[i]});
            xxt_fit.push_back({x1[box_final[i] - 1], y1[box_final[i] - 1]});
        }

        if (n_match < 6) {
            std::fill(box_final.begin(), box_final.end(), 0);
            LinearSolve::reportFailure(
                "Astrometry::patternMatching", LinearSolve::SolveStatus::FailedRankDeficient,
                "valid_matches=" + std::to_string(n_match) +
                    " required=6 removed_samples=" + std::to_string(removed_samples) +
                    " action=CLEAR_MATCHES");
            return;
        }

        double coe[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        LinearSolve::SolveDiagnostics fit_diagnostics;
        LinearSolve::SolveStatus fit_status =
            UniversalUtils::fitLinear2D(xx_fit, xxt_fit, coe, &fit_diagnostics);
        if (fit_status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure(
                "Astrometry::fitLinear2D", fit_status,
                LinearSolve::diagnosticsContext(fit_diagnostics) +
                    " removed_samples=" + std::to_string(removed_samples) +
                    " action=CLEAR_MATCHES");
            std::fill(box_final.begin(), box_final.end(), 0);
            return;
        }

        std::fill(mark1.begin(), mark1.end(), 0.0f);
        std::fill(mark3.begin(), mark3.end(), 0);

        for (int i = 0; i < n0; ++i) {
            if (finite0[i] == 0) continue;
            for (int j = 0; j < n1; ++j) {
                if (finite1[j] == 0) continue;
                double predicted_x = x0[i] * coe[0][0] + y0[i] * coe[0][1] + coe[0][2];
                double predicted_y = x0[i] * coe[1][0] + y0[i] * coe[1][1] + coe[1][2];
                int dx = fortranShiftBin(x1[j] - predicted_x);
                int dy = fortranShiftBin(y1[j] - predicted_y);
                if (dx >= -nn && dx <= nn && dy >= -nn && dy <= nn) {
                    mark1[(dy + nn) * dim + (dx + nn)] += 1.0f;
                }
            }
        }

        float thresh_val = 0.5f;
        constexpr int crowdy = 1;
        if (crowdy == 1) {
            std::vector<float> arr_pass(dim * dim);
            int idx = 0;
            for (int dx_i = -nn; dx_i <= nn; ++dx_i) {
                for (int dy_i = -nn; dy_i <= nn; ++dy_i) {
                    arr_pass[idx++] = mark1[(dy_i + nn) * dim + (dx_i + nn)];
                }
            }
            std::sort(arr_pass.begin(), arr_pass.end());
            float peak_val = mark1[(0 + nn) * dim + (0 + nn)];
            size_t pass_median_idx = arr_pass.size() / 2;
            if (pass_median_idx > 0) {
                --pass_median_idx;
            }
            thresh_val = (arr_pass[pass_median_idx] + peak_val) * 0.5f;
        }

        changed = 1;
        xmin = 0;
        xmax = 0;
        ymin = 0;
        ymax = 0;
        mark3[(0 + nn) * dim + (0 + nn)] = 1;

        while (changed == 1) {
            changed = 0;
            int current_xmin = xmin;
            int current_xmax = xmax;
            int current_ymin = ymin;
            int current_ymax = ymax;
            for (int i = current_xmin; i <= current_xmax; ++i) {
                for (int j = current_ymin; j <= current_ymax; ++j) {
                    if (mark3[(j + nn) * dim + (i + nn)] == 1) {
                        for (int u = i - r_fof; u <= i + r_fof; ++u) {
                            for (int v = j - r_fof; v <= j + r_fof; ++v) {
                                if (u >= -nn && u <= nn && v >= -nn && v <= nn) {
                                    int idx_uv = (v + nn) * dim + (u + nn);
                                    if (mark3[idx_uv] == 0 && mark1[idx_uv] > thresh_val) {
                                        mark3[idx_uv] = 1;
                                        xmin = std::min(xmin, u);
                                        xmax = std::max(xmax, u);
                                        ymin = std::min(ymin, v);
                                        ymax = std::max(ymax, v);
                                        changed = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::fill(box_final.begin(), box_final.end(), 0);
        n_match = 0;
        for (int i = 0; i < n0; ++i) {
            if (finite0[i] == 0) continue;
            for (int j = 0; j < n1; ++j) {
                if (finite1[j] == 0) continue;
                double predicted_x = x0[i] * coe[0][0] + y0[i] * coe[0][1] + coe[0][2];
                double predicted_y = x0[i] * coe[1][0] + y0[i] * coe[1][1] + coe[1][2];
                int dx = fortranShiftBin(x1[j] - predicted_x);
                int dy = fortranShiftBin(y1[j] - predicted_y);
                if (dx >= -nn && dx <= nn && dy >= -nn && dy <= nn) {
                    if (mark3[(dy + nn) * dim + (dx + nn)] == 1) {
                        n_match++;
                        box_final[i] = j + 1; // 1-based index
                        break;
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Generate matched astrometry calibration data.
    // Method: Preserve F77 matching/output logic with 17-digit double serialization.
    // ==========================================
    void genAstrometryData(const std::string& catStandard, int nx, int ny,
                           const std::vector<float>& map, const std::vector<int>& weight,
                           WCSParams& wcs, const std::string& filename, int& procError) {
        if (procError == 1) {
            MainIO::OutputFile ofs(filename);
            if (ofs) {
                ofs << std::setprecision(17) << wcs.crpix[0] << " " << wcs.crpix[1] << " "
                    << wcs.crval[0] << " " << wcs.crval[1] << "\n";
                ofs << std::setprecision(17) << wcs.cd[0][0] << " " << wcs.cd[0][1] << " "
                    << wcs.cd[1][0] << " " << wcs.cd[1][1] << "\n";
                ofs << "0 0 0\n";
            }
            return;
        }

        constexpr double astrometry_shift_ratio = 0.2;
        double ra = 0.0, dra = 0.0;
        double dec[2] = {0.0, 0.0};

        getRaDecRange(nx, ny, ra, dec, dra, wcs.crpix, wcs.cd, wcs.crval, astrometry_shift_ratio);

        int n_ref = 0;
        constexpr int nss_max = 10000;
        constexpr int n_user_max = 200;

        std::vector<double> ra_r;
        std::vector<double> dec_r;
        std::vector<double> xr;
        std::vector<double> yr;

        std::ifstream ifs(catStandard);
        if (!ifs) {
            std::cerr << "Error / gen_astrometry_data catalog file error: " << catStandard << std::endl;
            procError = 1;
            MainIO::OutputFile ofs(filename);
            if (ofs) {
                ofs << std::setprecision(17) << wcs.crpix[0] << " " << wcs.crpix[1] << " "
                    << wcs.crval[0] << " " << wcs.crval[1] << "\n";
                ofs << std::setprecision(17) << wcs.cd[0][0] << " " << wcs.cd[0][1] << " "
                    << wcs.cd[1][0] << " " << wcs.cd[1][1] << "\n";
                ofs << "0 0 0\n";
            }
            return;
        }

        std::string header;
        std::getline(ifs, header); // skip first line
        
        std::string line;
        std::istringstream iss;
        double a = 0.0, d = 0.0;
        while (std::getline(ifs, line)) {
            std::replace(line.begin(), line.end(), ',', ' ');
            iss.clear();
            iss.str(line);
            if (!(iss >> a >> d)) {
                continue; 
            }
            
            if (std::abs(diffra(a, ra)) > dra * 0.5) continue;
            if (d < dec[0] || d > dec[1]) continue;

            n_ref++;
            if (n_ref > nss_max) {
                std::cerr << "n_ref is too large!! " << filename << std::endl;
                MainIO::OutputFile ofs(filename);
                if (ofs) {
                    ofs << std::setprecision(17) << wcs.crpix[0] << " " << wcs.crpix[1] << " "
                        << wcs.crval[0] << " " << wcs.crval[1] << "\n";
                    ofs << std::setprecision(17) << wcs.cd[0][0] << " " << wcs.cd[0][1] << " "
                        << wcs.cd[1][0] << " " << wcs.cd[1][1] << "\n";
                    ofs << "0 0 0\n";
                }
                return;
            }

            ra_r.push_back(a);
            dec_r.push_back(d);
            double x_pixel = 0.0, y_pixel = 0.0;
            coordinateTransferSimple(a, d, x_pixel, y_pixel, -1, wcs.crpix, wcs.cd, wcs.crval);
            xr.push_back(x_pixel);
            yr.push_back(y_pixel);
        }
        ifs.close();

        int n_user = 0;
        std::vector<double> xs, ys;
        getAstrometryCatalog(nx, ny, nx, ny, map, weight, n_user_max, nss_max, n_user, xs, ys);

        int astrometry_shift_range = static_cast<int>(std::max(nx, ny) * astrometry_shift_ratio);
        std::vector<int> box(n_ref, 0);

        patternMatching(nss_max, n_ref, xr, yr, nss_max, n_user, xs, ys, astrometry_shift_range, box);

        int nss = 0;
        std::vector<double> ra2, dec2, x2, y2;
        for (int i = 0; i < n_ref; ++i) {
            if (box[i] == 0) continue;
            nss++;
            ra2.push_back(ra_r[i]);
            dec2.push_back(dec_r[i]);
            int j = box[i] - 1; // back to 0-based
            x2.push_back(xs[j]);
            y2.push_back(ys[j]);
        }

        // std::cout << nss << " " << n_ref << " " << n_user << " " << filename << std::endl;

       MainIO::OutputFile ofs(filename);
        if (ofs) {
            ofs << std::setprecision(17) << wcs.crpix[0] << " " << wcs.crpix[1] << " "
                << wcs.crval[0] << " " << wcs.crval[1] << "\n";
            ofs << std::setprecision(17) << wcs.cd[0][0] << " " << wcs.cd[0][1] << " "
                << wcs.cd[1][0] << " " << wcs.cd[1][1] << "\n";
            ofs << nss << " " << n_user << " " << n_ref << "\n";
            for (int i = 0; i < nss; ++i) {
                ofs << std::setprecision(17) << ra2[i] << " " << dec2[i] << " " << x2[i] << " " << y2[i] << "\n";
            }
        }
    }

    // ==========================================
    // Function: Solve the global exposure astrometry with failure status
    // Method: Assemble the original least-squares design matrix and solve both coordinate RHS columns with one pivoted QR.
    // ==========================================
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
        LinearSolve::SolveDiagnostics* diagnostics) {
        (void)np;
        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};

        int tot_valid = 0;
        int total_rows = 0;
        for (int k = 0; k < nc; ++k) {
            if (valid[k] == 1) {
                tot_valid++;
                total_rows += n[k];
            }
        }

        int k_sys = npd + tot_valid * 3;
        if (total_rows < k_sys || k_sys <= 0) {
            diag.rows = total_rows;
            diag.cols = k_sys;
            diag.rank = std::min(total_rows, k_sys);
            diag.required_rank = k_sys;
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd design = Eigen::MatrixXd::Zero(total_rows, k_sys);
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(total_rows, 2);

        int ic = 0;
        int row = 0;
        for (int k = 0; k < nc; ++k) {
            if (valid[k] == 0) continue;
            int n_stars = n[k];
            for (int i = 0; i < n_stars; ++i) {
                double xi_val = 0.0, eta_val = 0.0;
                raDecToXiEta(ra[k][i], dec[k][i], xi_val, eta_val, cRVAL[0], cRVAL[1]);

                int px = 0;
                int py = 1;
                int order = 1;
                int nn = 0;
                while (nn < npd) {
                    if (py == order) {
                        order++;
                        px = order;
                        py = 0;
                    } else {
                        px--;
                        py++;
                    }
                    nn++;
                    design(row, nn - 1) = std::pow(xi_val, px) * std::pow(eta_val, py);
                }

                int j = npd + ic * 3;
                design(row, j) = x[k][i];
                design(row, j + 1) = y[k][i];
                design(row, j + 2) = 1.0;
                rhs(row, 0) = xi_val;
                rhs(row, 1) = eta_val;
                row++;
            }
            ic++;
        }

        if (!design.allFinite() || !rhs.allFinite()) {
            diag.rows = design.rows();
            diag.cols = design.cols();
            diag.required_rank = design.cols();
            return LinearSolve::SolveStatus::FailedSolver;
        }

        LinearSolve::LeastSquaresQR solver;
        LinearSolve::SolveStatus status = solver.factorize(design, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        Eigen::MatrixXd solution;
        status = solver.solve(rhs, solution);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }
        Eigen::VectorXd vec1 = solution.col(0);
        Eigen::VectorXd vec2 = solution.col(1);

        std::vector<std::array<double, 2>> fitted_pu(npd, {0.0, 0.0});
        std::vector<std::array<double, 2>> fitted_crpix = cRPIX;
        std::vector<std::array<std::array<double, 2>, 2>> fitted_cd = cD;

        // Unpack PU into temporary storage so a later chip failure cannot expose partial results.
        int px = 0;
        int py = 1;
        int order = 1;
        int nn = 0;
        while (nn < npd) {
            if (py == order) {
                order++;
                px = order;
                py = 0;
            } else {
                px--;
                py++;
            }
            nn++;
            const int eta_index = nn - 1 + order - py * 2;
            if (eta_index < 0 || eta_index >= npd) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
            fitted_pu[nn - 1][0] = vec1(nn - 1);
            fitted_pu[nn - 1][1] = vec2(eta_index);
        }

        // Unpack chip-specific cD and cRPIX into temporary storage.
        ic = 0;
        for (int k = 0; k < nc; ++k) {
            if (valid[k] == 0) continue;
            int j = npd + ic * 3;
            fitted_cd[k][0][0] = vec1(j);
            fitted_cd[k][0][1] = vec1(j + 1);
            fitted_cd[k][1][0] = vec2(j);
            fitted_cd[k][1][1] = vec2(j + 1);

            double V1 = vec1(j + 2);
            double V2 = vec2(j + 2);
            double det = fitted_cd[k][0][0] * fitted_cd[k][1][1] -
                         fitted_cd[k][1][0] * fitted_cd[k][0][1];
            if (!std::isfinite(det) || det == 0.0) {
                return LinearSolve::SolveStatus::FailedSolver;
            }

            fitted_crpix[k][0] =
                -(V1 * fitted_cd[k][1][1] - V2 * fitted_cd[k][0][1]) / det;
            fitted_crpix[k][1] =
                -(V2 * fitted_cd[k][0][0] - V1 * fitted_cd[k][1][0]) / det;
            if (!std::isfinite(fitted_crpix[k][0]) ||
                !std::isfinite(fitted_crpix[k][1])) {
                return LinearSolve::SolveStatus::FailedSolver;
            }

            ic++;
        }

        for (int term = 0; term < npd; ++term) {
            PU[0][term] = fitted_pu[term][0];
            PU[1][term] = fitted_pu[term][1];
        }
        cRPIX = std::move(fitted_crpix);
        cD = std::move(fitted_cd);
        return LinearSolve::SolveStatus::Normal;
    }

    void checkAstrometryGlobal(int np, const std::vector<int>& n, int nc,
                               const std::vector<std::vector<double>>& ra,
                               const std::vector<std::vector<double>>& dec,
                               const std::vector<std::vector<double>>& x,
                               const std::vector<std::vector<double>>& y,
                               std::vector<std::array<double, 2>>& cRPIX,
                               std::vector<std::array<std::array<double, 2>, 2>>& cD,
                               const double cRVAL[2], double PU[2][LensingConfig::npd],
                               int npd, std::vector<int>& valid) {
        constexpr double tolerate_shift = 0.1;
        for (int ic = 0; ic < nc; ++ic) {
            if (valid[ic] == 0) continue;

            int n_stars = n[ic];
            std::vector<double> a_stars(n_stars);
            std::vector<double> d_stars(n_stars);
            for (int i = 0; i < n_stars; ++i) {
                a_stars[i] = ra[ic][i];
                d_stars[i] = dec[ic][i];
            }

            double ra_c = 0.0, dra = 0.0, dec_c = 0.0, ddec = 0.0;
            getRaDecBound(np, n_stars, a_stars, d_stars, ra_c, dra, dec_c, ddec);

            double crp[2] = {cRPIX[ic][0], cRPIX[ic][1]};
            double cdd[2][2] = {{cD[ic][0][0], cD[ic][0][1]}, {cD[ic][1][0], cD[ic][1][1]}};

            for (int i = 0; i < n_stars; ++i) {
                double aa = 0.0, dd = 0.0;
                double x_val = x[ic][i];
                double y_val = y[ic][i];
                coordinateTransferPU(aa, dd, x_val, y_val, 1, crp, cdd, cRVAL, PU, npd);

                if (std::abs(diffra(aa, ra_c)) > dra * (0.5 + tolerate_shift) ||
                    std::abs(dd - dec_c) > ddec * (0.5 + tolerate_shift) ||
                    std::isnan(aa) || std::isnan(dd)) {
                    valid[ic] = 0;
                    break;
                }
            }
        }
    }

    void readAstrometryPara(const std::string& filename, int ichip,
                            double cRPIX[2], double cD[2][2], double cRVAL[2],
                            double PU[2][LensingConfig::npd], int npd, int& procError) {
        if (procError == 1) return;

        std::ifstream ifs(filename);
        if (!ifs) {
            procError = 1;
            return;
        }

        if (!(ifs >> cRVAL[0] >> cRVAL[1])) {
            procError = 1;
            return;
        }

        for (int i = 0; i < npd; ++i) {
            if (!(ifs >> PU[0][i] >> PU[1][i])) {
                procError = 1;
                return;
            }
        }

        std::string dummyLine;
        std::getline(ifs, dummyLine); // clear newline from last PU
        for (int k = 1; k < ichip; ++k) {
            if (!std::getline(ifs, dummyLine)) {
                procError = 1;
                return;
            }
        }

        int j = 0, valid = 0;
        if (!(ifs >> j >> valid >> cRPIX[0] >> cRPIX[1] >> cD[0][0] >> cD[0][1] >> cD[1][0] >> cD[1][1])) {
            procError = 1;
            return;
        }

        if (valid == 0) {
            procError = 1;
        }
    }

    // ==========================================
    // Function: Fit and write exposure astrometry parameters.
    // Method: Preserve F77 header layout with 17-digit double serialization.
    // ==========================================
    void getAstrometry(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput) {
        constexpr int nss_max = 10000;
        std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::vector<int> nss(nchip, 0);
        std::vector<std::vector<double>> ra2(nchip, std::vector<double>(nss_max, 0.0));
        std::vector<std::vector<double>> dec2(nchip, std::vector<double>(nss_max, 0.0));
        std::vector<std::vector<double>> x2(nchip, std::vector<double>(nss_max, 0.0));
        std::vector<std::vector<double>> y2(nchip, std::vector<double>(nss_max, 0.0));

        std::vector<std::array<double, 2>> cRPIX2(nchip);
        std::vector<std::array<std::array<double, 2>, 2>> cD2(nchip);
        double cRVAL2[2] = {0.0, 0.0};
        double PU[2][LensingConfig::npd];
        std::fill(&PU[0][0], &PU[0][0] + 2 * LensingConfig::npd, 0.0);

        std::vector<int> valid(nchip, 0);
        int tot_valid = 0;
        int tot_source = 0;

        for (int ichip = 1; ichip <= nchip; ++ichip) {
            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip - 1]);
            std::string filename = OutputLayout::chipPath(
                dirOutput, "astrometry/dat_Astro", prefix, "_astro.dat");

            std::ifstream ifs(filename);
            if (!ifs) {
                std::cout << "Error： astro.dat error" << std::endl;
                valid[ichip - 1] = 0;
                continue;
            }

            double crpix0 = 0.0;
            double crpix1 = 0.0;
            double crval0 = 0.0;
            double crval1 = 0.0;
            if (!(ifs >> crpix0 >> crpix1 >> crval0 >> crval1)) {
                std::cout << "Error： astro.dat error" << std::endl;
                valid[ichip - 1] = 0;
                continue;
            }
            if (!std::isfinite(crpix0) || !std::isfinite(crpix1) ||
                !std::isfinite(crval0) || !std::isfinite(crval1)) {
                LinearSolve::reportFailure(
                    "Astrometry::getAstrometry", LinearSolve::SolveStatus::FailedSolver,
                    "exposure=" + prefix_expo + " chip=" + std::to_string(ichip) +
                        " reason=NON_FINITE_INPUT action=MARK_CHIP_INVALID");
                valid[ichip - 1] = 0;
                continue;
            }
            cRPIX2[ichip - 1][0] = crpix0;
            cRPIX2[ichip - 1][1] = crpix1;
            cRVAL2[0] = crval0;
            cRVAL2[1] = crval1;
            double cd00 = 0.0;
            double cd01 = 0.0;
            double cd10 = 0.0;
            double cd11 = 0.0;
            if (!(ifs >> cd00 >> cd01 >> cd10 >> cd11)) {
                std::cout << "Error： astro.dat error" << std::endl;
                valid[ichip - 1] = 0;
                continue;
            }
            if (!std::isfinite(cd00) || !std::isfinite(cd01) ||
                !std::isfinite(cd10) || !std::isfinite(cd11)) {
                LinearSolve::reportFailure(
                    "Astrometry::getAstrometry", LinearSolve::SolveStatus::FailedSolver,
                    "exposure=" + prefix_expo + " chip=" + std::to_string(ichip) +
                        " reason=NON_FINITE_INPUT action=MARK_CHIP_INVALID");
                valid[ichip - 1] = 0;
                continue;
            }
            cD2[ichip - 1][0][0] = cd00;
            cD2[ichip - 1][0][1] = cd01;
            cD2[ichip - 1][1][0] = cd10;
            cD2[ichip - 1][1][1] = cd11;

            int n_user = 0, n_ref = 0;
            if (!(ifs >> nss[ichip - 1] >> n_user >> n_ref)) {
                std::cout << "Error： astro.dat error" << std::endl;
                valid[ichip - 1] = 0;
                continue;
            }

            int actual_nss = std::min(nss[ichip - 1], nss_max);
            int kept_nss = 0;
            int removed_non_finite = 0;
            bool astro_row_error = false;
            for (int i = 0; i < actual_nss; ++i) {
                double ra_value = 0.0;
                double dec_value = 0.0;
                double x_value = 0.0;
                double y_value = 0.0;
                if (!(ifs >> ra_value >> dec_value >> x_value >> y_value)) {
                    std::cout << "Error： astro.dat error" << std::endl;
                    valid[ichip - 1] = 0;
                    astro_row_error = true;
                    break;
                }
                if (!std::isfinite(ra_value) || !std::isfinite(dec_value) ||
                    !std::isfinite(x_value) || !std::isfinite(y_value)) {
                    removed_non_finite++;
                    continue;
                }
                ra2[ichip - 1][kept_nss] = ra_value;
                dec2[ichip - 1][kept_nss] = dec_value;
                x2[ichip - 1][kept_nss] = x_value;
                y2[ichip - 1][kept_nss] = y_value;
                kept_nss++;
            }
            if (astro_row_error) {
                continue;
            }
            nss[ichip - 1] = kept_nss;
            if (kept_nss >= 10) {
                valid[ichip - 1] = 1;
                tot_valid++;
                tot_source += kept_nss;
            } else {
                valid[ichip - 1] = 0;
                if (removed_non_finite > 0) {
                    LinearSolve::reportFailure(
                        "Astrometry::getAstrometry", LinearSolve::SolveStatus::FailedRankDeficient,
                        "exposure=" + prefix_expo + " chip=" + std::to_string(ichip) +
                            " valid_samples=" + std::to_string(kept_nss) +
                            " required=10 removed_samples=" +
                            std::to_string(removed_non_finite) +
                            " action=MARK_CHIP_INVALID");
                }
            }
        }

        tot_valid = 0;
        tot_source = 0;
        if (std::isfinite(cRVAL2[0]) && std::isfinite(cRVAL2[1])) {
            for (int ichip = 0; ichip < nchip; ++ichip) {
                if (valid[ichip] == 0) continue;
                int kept_nss = 0;
                int removed_projection = 0;
                for (int i = 0; i < nss[ichip]; ++i) {
                    double xi = 0.0;
                    double eta = 0.0;
                    raDecToXiEta(ra2[ichip][i], dec2[ichip][i], xi, eta, cRVAL2[0], cRVAL2[1]);
                    if (!std::isfinite(xi) || !std::isfinite(eta)) {
                        removed_projection++;
                        continue;
                    }
                    ra2[ichip][kept_nss] = ra2[ichip][i];
                    dec2[ichip][kept_nss] = dec2[ichip][i];
                    x2[ichip][kept_nss] = x2[ichip][i];
                    y2[ichip][kept_nss] = y2[ichip][i];
                    kept_nss++;
                }
                nss[ichip] = kept_nss;
                if (kept_nss >= 10) {
                    tot_valid++;
                    tot_source += kept_nss;
                } else {
                    valid[ichip] = 0;
                    if (removed_projection > 0) {
                        LinearSolve::reportFailure(
                            "Astrometry::getAstrometry", LinearSolve::SolveStatus::FailedRankDeficient,
                            "exposure=" + prefix_expo + " chip=" + std::to_string(ichip + 1) +
                                " valid_samples=" + std::to_string(kept_nss) +
                                " required=10 removed_samples=" +
                                std::to_string(removed_projection) +
                                " action=MARK_CHIP_INVALID");
                    }
                }
            }
        } else {
            LinearSolve::reportFailure(
                "Astrometry::getAstrometry", LinearSolve::SolveStatus::FailedSolver,
                "exposure=" + prefix_expo +
                    " reason=NON_FINITE_INPUT action=MARK_EXPOSURE_INVALID");
            std::fill(valid.begin(), valid.end(), 0);
        }

        int required_system_rows = (LensingConfig::npd + tot_valid * 3) * 3;
        if (tot_valid > 0 && tot_source >= required_system_rows) {
            LinearSolve::SolveDiagnostics fit_diagnostics;
            LinearSolve::SolveStatus fit_status = measureAstrometryGlobal(
                nss_max, nss, nchip, ra2, dec2, x2, y2, cRPIX2, cD2, cRVAL2,
                PU, LensingConfig::npd, valid, &fit_diagnostics);
            if (fit_status == LinearSolve::SolveStatus::Normal) {
                checkAstrometryGlobal(nss_max, nss, nchip, ra2, dec2, x2, y2,
                                      cRPIX2, cD2, cRVAL2, PU, LensingConfig::npd, valid);
            } else {
                LinearSolve::reportFailure(
                    "Astrometry::measureAstrometryGlobal", fit_status,
                    "exposure=" + prefix_expo + " " +
                        LinearSolve::diagnosticsContext(fit_diagnostics) +
                        " action=MARK_EXPOSURE_INVALID");
                std::fill(valid.begin(), valid.end(), 0);
            }
        } else {
            LinearSolve::reportFailure(
                "Astrometry::measureAstrometryGlobal",
                LinearSolve::SolveStatus::FailedRankDeficient,
                "exposure=" + prefix_expo + " valid_samples=" + std::to_string(tot_source) +
                    " required=" + std::to_string(required_system_rows) +
                    " action=MARK_EXPOSURE_INVALID");
            std::fill(valid.begin(), valid.end(), 0);
        }

        std::string out_head = dirOutput + "/astrometry/Head/" + prefix_expo + ".head";
        MainIO::OutputFile ofs(out_head);
        if (ofs) {
            ofs << std::setprecision(17) << cRVAL2[0] << " " << cRVAL2[1] << "\n";
            for (int i = 0; i < LensingConfig::npd; ++i) {
                ofs << std::setprecision(17) << PU[0][i] << " " << PU[1][i] << "\n";
            }
            for (int k = 1; k <= nchip; ++k) {
                ofs << k << " " << valid[k - 1] << " "
                    << std::setprecision(17)
                    << cRPIX2[k - 1][0] << " " << cRPIX2[k - 1][1] << " "
                    << cD2[k - 1][0][0] << " " << cD2[k - 1][0][1] << " "
                    << cD2[k - 1][1][0] << " " << cD2[k - 1][1][1] << "\n";
            }
        }
    }

    // ==========================================
    // Function: Apply exposure astrometry to each chip and write diagnostics.
    // Method: Preserve F77 check-file layout with 17-digit double serialization.
    // ==========================================
    void chipProcessAstrometry(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput) {
        getAstrometry(imageFiles, nchip, dirOutput);

        std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string check_filename = dirOutput + "/astrometry/dat_Chk/" + prefix_expo + "_check.dat";
        MainIO::OutputFile check_ofs(check_filename);
        if (!check_ofs) {
            std::cerr << "Error writing check file: " << check_filename << std::endl;
            return;
        }

        double cRPIX[2], cD[2][2], cRVAL[2];
        double PU[2][LensingConfig::npd];

        for (int ichip = 1; ichip <= nchip; ++ichip) {
            int proc_error = 0;
            std::string head_filename = dirOutput + "/astrometry/Head/" + prefix_expo + ".head";
            readAstrometryPara(head_filename, ichip, cRPIX, cD, cRVAL, PU, LensingConfig::npd, proc_error);

            if (proc_error == 0) {
                std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip - 1]);
                std::string norm_filename = OutputLayout::chipPath(
                    dirOutput, "stamps/Norm", prefix, "_norm.fits");
                
                WCSParams wcs;
                wcs.crpix[0] = cRPIX[0];
                wcs.crpix[1] = cRPIX[1];
                wcs.cd[0][0] = cD[0][0];
                wcs.cd[0][1] = cD[0][1];
                wcs.cd[1][0] = cD[1][0];
                wcs.cd[1][1] = cD[1][1];
                wcs.crval[0] = cRVAL[0];
                wcs.crval[1] = cRVAL[1];

                FitsIO::updatePara(norm_filename, wcs);

                std::string astro_filename = OutputLayout::chipPath(
                    dirOutput, "astrometry/dat_Astro", prefix, "_astro.dat");
                std::ifstream astro_ifs(astro_filename);
                if (astro_ifs) {
                    std::string dummyLine1, dummyLine2;
                    if (!std::getline(astro_ifs, dummyLine1) || !std::getline(astro_ifs, dummyLine2)) {
                        std::cout << "Error： astro.dat error" << std::endl;
                        continue;
                    }
                    int n = 0, j = 0, k = 0;
                    if (astro_ifs >> n >> j >> k) {
                        for (int i = 0; i < n; ++i) {
                            double ra_val = 0.0, dec_val = 0.0, x_val = 0.0, y_val = 0.0;
                            if (astro_ifs >> ra_val >> dec_val >> x_val >> y_val) {
                                double ra2_val = 0.0, dec2_val = 0.0;
                                coordinateTransferPU(ra2_val, dec2_val, x_val, y_val, 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd);
                                check_ofs << std::setprecision(17) << ra_val << " " << dec_val << " " << ra2_val << " " << dec2_val << "\n";
                            } else {
                                std::cout << "Error： astro.dat error" << std::endl;
                                break;
                            }
                        }
                    } else {
                        std::cout << "Error： astro.dat error" << std::endl;
                    }
                } else {
                    std::cout << "Error： astro.dat error" << std::endl;
                }
            }
        }
    }

    void procAstrometry(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = EXPO_FILE[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

        chipProcessAstrometry(image_files, image_files.size(), dir_output);
    }
}
