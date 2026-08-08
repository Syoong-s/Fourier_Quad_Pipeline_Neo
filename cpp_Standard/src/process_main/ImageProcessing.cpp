#include "ImageProcessing.hpp"
#include "NumericalRecipes.hpp"
#include "UniversalUtils.hpp"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <fftw3.h>

namespace ImageProcessing {

    // ==========================================
    // Function: Convert F77 stamp(x,y) coordinates to C++ row-major storage
    // Method: Preserve F77 first-dimension-is-x semantics in vectors stored as y*n+x.
    // ==========================================
    namespace {
        inline std::size_t stampIndex(int x, int y, int n) {
            return static_cast<std::size_t>(y) * n + x;
        }
    }

    // ==========================================
    // Function: Detect and mask the primary source on a stamp
    // Method: Mirror F77 mark_source with x as first stamp dimension and y as second dimension.
    // ==========================================
    void markSource(int n, const std::vector<float>& stamp, std::vector<int>& weight, double sig, 
                    double source_thresh, double core_thresh, int boundx[2], int boundy[2], 
                    double& total_flux, int& total_area, double& peak, double& half_light_flux, 
                    int& half_light_area, int& flag, double& radius, int& xp, int& yp) {
        
        int cc = n / 2;
        double r2min = n * n * 2.0;
        int ix = -1;
        int iy = -1;
        double r2_thresh = std::pow(n / 8.0, 2);

        double thresh1 = core_thresh * sig;
        double thresh2 = source_thresh * sig;

        std::vector<int> mark(n * n, 0);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (stamp[idx] >= thresh1 && weight[idx] > 0) {
                    mark[idx] = 1;
                    double r2 = std::pow(i - cc, 2) + std::pow(j - cc, 2);
                    if (r2 < r2min) {
                        r2min = r2;
                        ix = i;
                        iy = j;
                    }
                } else if (weight[idx] == 0) {
                    mark[idx] = 1;
                } else {
                    mark[idx] = 0;
                }
            }
        }

        if (r2min > r2_thresh) {
            flag = -2;
            return;
        }

        mark[stampIndex(ix, iy, n)] = 2;
        int changed = 1;
        int xmin = ix;
        int xmax = ix;
        int ymin = iy;
        int ymax = iy;

        total_flux = stamp[stampIndex(ix, iy, n)];
        total_area = 1;
        peak = stamp[stampIndex(ix, iy, n)];
        xp = ix;
        yp = iy;

        while (changed == 1) {
            changed = 0;
            for (int i = xmin; i <= xmax; ++i) {
                for (int j = ymin; j <= ymax; ++j) {
                    if (mark[stampIndex(i, j, n)] == 2) {
                        for (int u = std::max(i - 1, 0); u <= std::min(i + 1, n - 1); ++u) {
                            for (int v = std::max(j - 1, 0); v <= std::min(j + 1, n - 1); ++v) {
                                std::size_t uidx = stampIndex(u, v, n);
                                if (mark[uidx] < 2 && stamp[uidx] >= thresh2 && weight[uidx] > 0) {
                                    mark[uidx] = 2;
                                    changed = 1;
                                    xmin = std::min(xmin, u);
                                    ymin = std::min(ymin, v);
                                    xmax = std::max(xmax, u);
                                    ymax = std::max(ymax, v);
                                    if (stamp[uidx] > peak) {
                                        peak = stamp[uidx];
                                        xp = u;
                                        yp = v;
                                    }
                                    total_flux += stamp[uidx];
                                    total_area += 1;
                                } else if (weight[uidx] == 0) {
                                    flag = -1;
                                    return;
                                }
                            }
                        }
                    }
                }
            }
        }

        half_light_flux = 0.0;
        half_light_area = 0;
        double thresh = peak * 0.5;
        double r2max = 0.0;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (mark[idx] == 2) {
                    double r_2 = std::pow(i - xp, 2) + std::pow(j - yp, 2);
                    r2max = std::max(r2max, r_2);
                    if (stamp[idx] >= thresh) {
                        half_light_flux += stamp[idx];
                        half_light_area += 1;
                    }
                }
            }
        }

        radius = std::sqrt(r2max);
        total_flux /= sig;
        half_light_flux /= sig;
        peak /= sig;

        boundx[0] = xmin;
        boundx[1] = xmax;
        boundy[0] = ymin;
        boundy[1] = ymax;

        changed = 1;
        while (changed == 1) {
            changed = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::size_t idx = stampIndex(i, j, n);
                    if (mark[idx] != 1) continue;
                    mark[idx] = -1;
                    for (int u = std::max(i - 1, 0); u <= std::min(i + 1, n - 1); ++u) {
                        for (int v = std::max(j - 1, 0); v <= std::min(j + 1, n - 1); ++v) {
                            std::size_t uidx = stampIndex(u, v, n);
                            if (mark[uidx] == 0 && stamp[uidx] >= thresh2) {
                                mark[uidx] = 1;
                            }
                        }
                    }
                    changed = 1;
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (mark[idx] != -1) continue;
                mark[idx] = 1;
                for (int u = std::max(i - 2, 0); u <= std::min(i + 2, n - 1); ++u) {
                    for (int v = std::max(j - 2, 0); v <= std::min(j + 2, n - 1); ++v) {
                        std::size_t uidx = stampIndex(u, v, n);
                        if (mark[uidx] == 0) {
                            mark[uidx] = 1;
                        } else if (mark[uidx] == 2) {
                            flag = -1;
                            return;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (mark[idx] == 1) {
                    mark[idx] = 0;
                    weight[idx] = 0;
                }
            }
        }

        flag = 10;
        int s = 2;
        while (s <= 10) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::size_t idx = stampIndex(i, j, n);
                    if (mark[idx] == s) {
                        for (int u = std::max(i - 1, 0); u <= std::min(i + 1, n - 1); ++u) {
                            for (int v = std::max(j - 1, 0); v <= std::min(j + 1, n - 1); ++v) {
                                std::size_t uidx = stampIndex(u, v, n);
                                if (mark[uidx] == 0) {
                                    mark[uidx] = s + 1;
                                    if (weight[uidx] == 0) {
                                        flag = mark[uidx] - 2;
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            s++;
        }
    }

    // ==========================================
    // Function: Mask detected sources in a noise stamp
    // Method: Mirror F77 mark_noise with x-major loop order and row-major vector storage.
    // ==========================================
    void markNoise(int n, const std::vector<float>& stamp, std::vector<int>& weight, double sig, 
                   double source_thresh, double core_thresh) {
        double thresh1 = core_thresh * sig;
        double thresh2 = source_thresh * sig;

        std::vector<int> mark(n * n, 0);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (stamp[idx] >= thresh1 && weight[idx] > 0) {
                    mark[idx] = 1;
                } else if (weight[idx] == 0) {
                    mark[idx] = 1;
                } else {
                    mark[idx] = 0;
                }
            }
        }

        int changed = 1;
        while (changed == 1) {
            changed = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::size_t idx = stampIndex(i, j, n);
                    if (mark[idx] != 1) continue;
                    mark[idx] = -1;
                    for (int u = std::max(i - 1, 0); u <= std::min(i + 1, n - 1); ++u) {
                        for (int v = std::max(j - 1, 0); v <= std::min(j + 1, n - 1); ++v) {
                            std::size_t uidx = stampIndex(u, v, n);
                            if (mark[uidx] == 0 && stamp[uidx] >= thresh2) {
                                mark[uidx] = 1;
                            }
                        }
                    }
                    changed = 1;
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (mark[idx] != -1) continue;
                mark[idx] = 1;
                for (int u = std::max(i - 2, 0); u <= std::min(i + 2, n - 1); ++u) {
                    for (int v = std::max(j - 2, 0); v <= std::min(j + 2, n - 1); ++v) {
                        std::size_t uidx = stampIndex(u, v, n);
                        if (mark[uidx] == 0) {
                            mark[uidx] = 1;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::size_t idx = stampIndex(i, j, n);
                if (mark[idx] == 1) {
                    weight[idx] = 0;
                }
            }
        }
    }

    // ==========================================
    // Function: Subtract a fitted 2D background plane from a stamp
    // Method: Mirror F77 flatten_stamp_2D coordinates while using row-major storage.
    // ==========================================
    void flattenStamp2D(int ns, int nl, std::vector<float>& stamp, const std::vector<int>& weight, int& ierror) {
        ierror = 0;
        int d1 = (nl - ns) / 2;
        std::vector<Point3D> points;
        points.reserve(nl * nl);

        for (int i = 0; i < nl; ++i) {
            for (int j = 0; j < nl; ++j) {
                std::size_t idx = stampIndex(i, j, nl);
                if ((i < d1 || i >= nl - d1 || j < d1 || j >= nl - d1) && weight[idx] == 1) {
                    points.push_back({static_cast<double>(i + 1), static_cast<double>(j + 1), static_cast<double>(stamp[idx])});
                }
            }
        }

        int max_np = nl * nl - ns * ns;
        if (points.size() <= max_np * 0.3) {
            ierror = -1;
            return;
        }

        double aa = 0.0, bb = 0.0, cc = 0.0;
        UniversalUtils::findSlope2D(points, aa, bb, cc);

        for (int i = 0; i < nl; ++i) {
            for (int j = 0; j < nl; ++j) {
                std::size_t idx = stampIndex(i, j, nl);
                stamp[idx] -= static_cast<float>(aa + bb * (i + 1) + cc * (j + 1));
            }
        }
    }

    void flattenStampNew(int ns, int nl, std::vector<float>& stamp, const std::vector<int>& weight, int& ierror) {
        ierror = 0;
        int d1 = (nl - ns) / 2;
        std::vector<float> border_vals;
        border_vals.reserve(nl * nl - ns * ns);

        for (int i = 0; i < nl; ++i) {
            for (int j = 0; j < nl; ++j) {
                int idx = i * nl + j;
                if ((i < d1 || i >= nl - d1 || j < d1 || j >= nl - d1) && weight[idx] == 1) {
                    border_vals.push_back(stamp[idx]);
                }
            }
        }

        int max_np = nl * nl - ns * ns;
        if (border_vals.size() <= max_np * 0.3) {
            ierror = -1;
            return;
        }

        std::sort(border_vals.begin(), border_vals.end());
        int np = border_vals.size();
        double bg_median = 0.0;
        if (np % 2 == 0) {
            bg_median = (border_vals[np / 2 - 1] + border_vals[np / 2]) / 2.0;
        } else {
            bg_median = border_vals[np / 2];
        }

        for (int i = 0; i < nl; ++i) {
            for (int j = 0; j < nl; ++j) {
                stamp[i * nl + j] -= static_cast<float>(bg_median);
            }
        }
    }

    // ==========================================
    // Function: Fill masked pixels with Gaussian noise
    // Method: Preserve F77 decorate_stamp x-major traversal and row-major storage.
    // ==========================================
    void decorateStamp(int ns, double sig, const std::vector<int>& weights, std::vector<float>& stamp) {
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                std::size_t idx = stampIndex(i, j, ns);
                if (weights[idx] == 0) {
                    stamp[idx] = static_cast<float>(NumericalRecipes::gasdev() * sig);
                }
            }
        }
    }

    void smoothGrid33(std::vector<float>& f) {
        static const Eigen::Matrix<double, 6, 6> matx_inv = []() {
            Eigen::Matrix<double, 6, 6> mat = Eigen::Matrix<double, 6, 6>::Zero();
            std::array<double, 3> coords = {-1.0, 0.0, 1.0};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    double x = coords[i];
                    double y = coords[j];
                    Eigen::Matrix<double, 6, 1> vec;
                    vec << 1.0, x, y, x*x, x*y, y*y;
                    mat += vec * vec.transpose();
                }
            }
            return mat.inverse().eval();
        }();

        Eigen::Matrix<double, 6, 1> vec = Eigen::Matrix<double, 6, 1>::Zero();
        std::array<double, 3> coords = {-1.0, 0.0, 1.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double x = coords[i];
                double y = coords[j];
                double val = f[i * 3 + j];
                vec[0] += val;
                vec[1] += val * x;
                vec[2] += val * y;
                vec[3] += val * x * x;
                vec[4] += val * x * y;
                vec[5] += val * y * y;
            }
        }

        Eigen::Matrix<double, 6, 1> a = matx_inv * vec;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double x = coords[i];
                double y = coords[j];
                f[i * 3 + j] = static_cast<float>(a[0] + a[1]*x + a[2]*y + a[3]*x*x + a[4]*x*y + a[5]*y*y);
            }
        }
    }

    void smoothImage33(int nx, int ny, std::vector<float>& map) {
        std::vector<float> temp(nx * ny, 0.0f);
        std::vector<float> f(9);

        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                for (int u = 0; u < 3; ++u) {
                    for (int v = 0; v < 3; ++v) {
                        f[u * 3 + v] = map[(j + v - 1) * nx + (i + u - 1)];
                    }
                }
                smoothGrid33(f);

                temp[j * nx + i] = f[1 * 3 + 1];
                if (i == 1) {
                    temp[j * nx + 0] = f[0 * 3 + 1];
                }
                if (i == nx - 2) {
                    temp[j * nx + (nx - 1)] = f[2 * 3 + 1];
                }
                if (j == 1) {
                    temp[0 * nx + i] = f[1 * 3 + 0];
                }
                if (j == ny - 2) {
                    temp[(ny - 1) * nx + i] = f[1 * 3 + 2];
                }
                if (i == 1 && j == 1) {
                    temp[0 * nx + 0] = f[0 * 3 + 0];
                }
                if (i == 1 && j == ny - 2) {
                    temp[(ny - 1) * nx + 0] = f[0 * 3 + 2];
                }
                if (i == nx - 2 && j == 1) {
                    temp[0 * nx + (nx - 1)] = f[2 * 3 + 0];
                }
                if (i == nx - 2 && j == ny - 2) {
                    temp[(ny - 1) * nx + (nx - 1)] = f[2 * 3 + 2];
                }
            }
        }

        map = std::move(temp);
    }

    void smoothGrid55(std::vector<float>& f) {
        static const Eigen::Matrix<double, 6, 6> matx_inv = []() {
            Eigen::Matrix<double, 6, 6> mat = Eigen::Matrix<double, 6, 6>::Zero();
            std::array<double, 5> coords = {-2.0, -1.0, 0.0, 1.0, 2.0};
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    double x = coords[i];
                    double y = coords[j];
                    Eigen::Matrix<double, 6, 1> vec;
                    vec << 1.0, x, y, x*x, x*y, y*y;
                    mat += vec * vec.transpose();
                }
            }
            return mat.inverse().eval();
        }();

        Eigen::Matrix<double, 6, 1> vec = Eigen::Matrix<double, 6, 1>::Zero();
        std::array<double, 5> coords = {-2.0, -1.0, 0.0, 1.0, 2.0};
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double x = coords[i];
                double y = coords[j];
                double val = f[i * 5 + j];
                vec[0] += val;
                vec[1] += val * x;
                vec[2] += val * y;
                vec[3] += val * x * x;
                vec[4] += val * x * y;
                vec[5] += val * y * y;
            }
        }

        Eigen::Matrix<double, 6, 1> a = matx_inv * vec;

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double x = coords[i];
                double y = coords[j];
                f[i * 5 + j] = static_cast<float>(a[0] + a[1]*x + a[2]*y + a[3]*x*x + a[4]*x*y + a[5]*y*y);
            }
        }
    }

    void smoothGrid55_3rd_order(std::vector<float>& f) {
        static const Eigen::Matrix<double, 10, 10> matx_inv = []() {
            Eigen::Matrix<double, 10, 10> mat = Eigen::Matrix<double, 10, 10>::Zero();
            std::array<double, 5> coords = {-2.0, -1.0, 0.0, 1.0, 2.0};
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    double x = coords[i];
                    double y = coords[j];
                    Eigen::Matrix<double, 10, 1> v;
                    v << 1.0, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y;
                    mat += v * v.transpose();
                }
            }
            return mat.inverse().eval();
        }();

        Eigen::Matrix<double, 10, 1> vec = Eigen::Matrix<double, 10, 1>::Zero();
        std::array<double, 5> coords = {-2.0, -1.0, 0.0, 1.0, 2.0};
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double x = coords[i];
                double y = coords[j];
                double val = f[i * 5 + j];
                vec[0] += val;
                vec[1] += val * x;
                vec[2] += val * y;
                vec[3] += val * x * x;
                vec[4] += val * x * y;
                vec[5] += val * y * y;
                vec[6] += val * x * x * x;
                vec[7] += val * x * x * y;
                vec[8] += val * x * y * y;
                vec[9] += val * y * y * y;
            }
        }

        Eigen::Matrix<double, 10, 1> a = matx_inv * vec;

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double x = coords[i];
                double y = coords[j];
                f[i * 5 + j] = static_cast<float>(
                    a[0] + a[1]*x + a[2]*y + a[3]*x*x + a[4]*x*y + a[5]*y*y +
                    a[6]*x*x*x + a[7]*x*x*y + a[8]*x*y*y + a[9]*y*y*y
                );
            }
        }
    }

    void smoothImage55(int nx, int ny, std::vector<float>& map, int ord) {
        std::vector<float> temp(nx * ny, 0.0f);
        std::vector<float> f(25);

        for (int i = 2; i < nx - 2; ++i) {
            for (int j = 2; j < ny - 2; ++j) {
                for (int u = 0; u < 5; ++u) {
                    for (int v = 0; v < 5; ++v) {
                        f[u * 5 + v] = map[(j + v - 2) * nx + (i + u - 2)];
                    }
                }

                if (ord == 1) {
                    smoothGrid55(f);
                } else {
                    smoothGrid55_3rd_order(f);
                }

                temp[j * nx + i] = f[2 * 5 + 2];
                if (i == 2) {
                    temp[j * nx + 0] = f[0 * 5 + 2];
                    temp[j * nx + 1] = f[1 * 5 + 2];
                }
                if (i == nx - 3) {
                    temp[j * nx + (nx - 1)] = f[4 * 5 + 2];
                    temp[j * nx + (nx - 2)] = f[3 * 5 + 2];
                }
                if (j == 2) {
                    temp[0 * nx + i] = f[2 * 5 + 0];
                    temp[1 * nx + i] = f[2 * 5 + 1];
                }
                if (j == ny - 3) {
                    temp[(ny - 1) * nx + i] = f[2 * 5 + 4];
                    temp[(ny - 2) * nx + i] = f[2 * 5 + 3];
                }
                if (i == 2 && j == 2) {
                    temp[0 * nx + 0] = f[0 * 5 + 0];
                    temp[1 * nx + 0] = f[0 * 5 + 1];
                    temp[0 * nx + 1] = f[1 * 5 + 0];
                    temp[1 * nx + 1] = f[1 * 5 + 1];
                }
                if (i == 2 && j == ny - 3) {
                    temp[(ny - 1) * nx + 0] = f[0 * 5 + 4];
                    temp[(ny - 2) * nx + 0] = f[0 * 5 + 3];
                    temp[(ny - 1) * nx + 1] = f[1 * 5 + 4];
                    temp[(ny - 2) * nx + 1] = f[1 * 5 + 3];
                }
                if (i == nx - 3 && j == 2) {
                    temp[0 * nx + (nx - 1)] = f[4 * 5 + 0];
                    temp[0 * nx + (nx - 2)] = f[3 * 5 + 0];
                    temp[1 * nx + (nx - 1)] = f[4 * 5 + 1];
                    temp[1 * nx + (nx - 2)] = f[3 * 5 + 1];
                }
                if (i == nx - 3 && j == ny - 3) {
                    temp[(ny - 1) * nx + (nx - 1)] = f[4 * 5 + 4];
                    temp[(ny - 1) * nx + (nx - 2)] = f[3 * 5 + 4];
                    temp[(ny - 2) * nx + (nx - 1)] = f[4 * 5 + 3];
                    temp[(ny - 2) * nx + (nx - 2)] = f[3 * 5 + 3];
                }
            }
        }

        map = std::move(temp);
    }

    void smoothGrid55WithHole(const std::vector<float>& f, int xh, int yh, float& fc) {
        std::vector<std::pair<int, int>> valid_points;
        for (int i = -2; i <= 2; ++i) {
            for (int j = -2; j <= 2; ++j) {
                if (std::abs(i) == 2 && std::abs(j) == 2) continue;
                if (i == xh && j == yh) continue;
                valid_points.push_back({i, j});
            }
        }

        int n_pts = valid_points.size();
        Eigen::MatrixXd X(n_pts, 6);
        Eigen::VectorXd Z(n_pts);

        for (int k = 0; k < n_pts; ++k) {
            double x = valid_points[k].first;
            double y = valid_points[k].second;
            int grid_i = valid_points[k].first + 2;
            int grid_j = valid_points[k].second + 2;
            Z[k] = f[grid_i * 5 + grid_j];

            X(k, 0) = 1.0;
            X(k, 1) = x;
            X(k, 2) = y;
            X(k, 3) = x * x;
            X(k, 4) = x * y;
            X(k, 5) = y * y;
        }

        Eigen::VectorXd a = X.colPivHouseholderQr().solve(Z);
        fc = static_cast<float>(a[0]);
    }

    // ==========================================
    // Function: Smooth a periodic image while excluding a local hole.
    // Method: Preserve F77 (x,y) semantics in row-major y*nx+x storage.
    // ==========================================
    void smoothImage55Hole(int nx, int ny, std::vector<float>& map) {
        std::vector<float> temp(nx * ny, 0.0f);
        int cx = nx / 2;
        int cy = ny / 2;
        std::vector<float> f(25);

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int u = -2; u <= 2; ++u) {
                    for (int v = -2; v <= 2; ++v) {
                        int map_i = (i + u + nx) % nx;
                        int map_j = (j + v + ny) % ny;
                        f[(u + 2) * 5 + (v + 2)] = map[map_j * nx + map_i];
                    }
                }

                int xh = cx - i;
                int yh = cy - j;
                if (std::abs(xh) > 2 || std::abs(yh) > 2) {
                    xh = 2;
                    yh = 2;
                }

                float fc = 0.0f;
                smoothGrid55WithHole(f, xh, yh, fc);
                temp[j * nx + i] = fc;
            }
        }

        map = std::move(temp);
    }

    void smoothImage55HoleLn(int nx, int ny, std::vector<float>& map) {
        if (map.empty()) return;
        float mapmin = map[0];
        float mapmax = map[0];
        for (float val : map) {
            if (val < mapmin) mapmin = val;
            if (val > mapmax) mapmax = val;
        }

        float delta = 0.0001f * (mapmax - mapmin);
        for (float& val : map) {
            val = std::log(val + delta);
        }

        smoothImage55Hole(nx, ny, map);

        for (float& val : map) {
            val = std::exp(val) - delta;
        }
    }

    void processPowers(int n, std::vector<float>& sourcep, const std::vector<float>& noisep) {
        for (int i = 0; i < n * n; ++i) {
            sourcep[i] -= noisep[i];
        }

        double temp = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            temp += sourcep[i * n + 0] + sourcep[i * n + (n - 1)] + sourcep[0 * n + i] + sourcep[(n - 1) * n + i];
        }

        temp /= (4.0 * (n - 2));

        for (int i = 0; i < n * n; ++i) {
            sourcep[i] -= static_cast<float>(temp);
        }
    }

    // ==========================================
    // Function: Normalize a two-dimensional power spectrum.
    // Method: Preserve F77 center/neighbour coordinates in row-major storage.
    // ==========================================
    void regularizePower(int nx, int ny, std::vector<float>& power, int star_smooth) {
        int cx = nx / 2;
        int cy = ny / 2;

        if (star_smooth >= 1) {
            double temp = 1.0 / power[cy * nx + cx];
            for (float& val : power) {
                val = static_cast<float>(val * temp);
            }
        } else {
            double temp = 4.0 / (power[cy * nx + (cx + 1)] +
                                 power[cy * nx + (cx - 1)] +
                                 power[(cy + 1) * nx + cx] +
                                 power[(cy - 1) * nx + cx]);
            for (float& val : power) {
                val = static_cast<float>(val * temp);
            }
            power[cy * nx + cx] = 1.0f;
        }
    }

    // ==========================================
    // Function: Measure thresholded power-spectrum ellipticity.
    // Method: Preserve F77 x/y flood-fill coordinates in row-major storage.
    // ==========================================
    void getPowerShape(int nx, int ny, const std::vector<float>& power, double& e, double thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        double thresh = power[cy * nx + cx] * thresh_ratio;

        double e1 = 0.0;
        double e2 = 0.0;
        double norm = 0.0;

        std::vector<int> mark(nx * ny, 0);
        mark[cy * nx + cx] = 1;

        std::vector<std::pair<int, int>> stack;
        stack.reserve(nx * ny);
        stack.push_back({cx, cy});

        size_t area0 = 0;
        while (stack.size() > area0) {
            size_t tempi = stack.size();
            for (size_t k = area0; k < tempi; ++k) {
                int x = stack[k].first;
                int y = stack[k].second;

                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack.push_back({u, v});

                            double kx = u - cx;
                            double ky = v - cy;
                            e1 += (kx * kx - ky * ky);
                            e2 += 2.0 * kx * ky;
                            norm += (kx * kx + ky * ky);
                        }
                    }
                }
            }
            area0 = tempi;
        }

        if (norm > 0.0) {
            e1 /= norm;
            e2 /= norm;
            e = e1 * e1 + e2 * e2;
        } else {
            e = 0.0;
        }
    }

    // ==========================================
    // Function: Estimate robust sigma and median from random samples
    // Method: Match F77 get_sig_med with exactly 1000 sampled values and 1-based percentile indices.
    // ==========================================
    void getSigMed(int nx, int ny, const std::vector<float>& image, float& sig, float& med) {
        const int npp = 1000;
        const int margin = 2;
        std::vector<float> pix(npp);

        for (int i = 0; i < npp; ++i) {
            int ix = static_cast<int>(NumericalRecipes::ran1() * (nx - 2 * margin)) + margin - 1;
            int iy = static_cast<int>(NumericalRecipes::ran1() * (ny - 2 * margin)) + margin - 1;
            pix[i] = image[iy * nx + ix];
        }

        std::sort(pix.begin(), pix.end());

        sig = 0.5f * (pix[5 * npp / 6 - 1] - pix[npp / 6 - 1]);
        med = pix[npp / 2 - 1];
    }

    void getEntropy(int nx, int ny, const std::vector<float>& image, double sig, double med, int r, std::vector<float>& entropy) {
        entropy.assign(nx * ny, 0.0f);
        std::array<double, 5> a = {0.25 * sig, 0.52 * sig, 0.85 * sig, 1.28 * sig, 100.0 * sig};

        auto wrap = [](int val, int limit) {
            int res = val % limit;
            return res < 0 ? res + limit : res;
        };

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                std::array<double, 5> p = {0.0, 0.0, 0.0, 0.0, 0.0};

                for (int ii = i - r; ii <= i + r; ++ii) {
                    int ix = wrap(ii, nx);
                    for (int jj = j - r; jj <= j + r; ++jj) {
                        int iy = wrap(jj, ny);
                        double val = std::abs(image[iy * nx + ix] - med);
                        int k = 0;
                        while (k < 4 && val > a[k]) {
                            k++;
                        }
                        p[k] += 1.0;
                    }
                }

                double ent = 0.0;
                for (int k = 0; k < 5; ++k) {
                    if (p[k] > 0.0) {
                        ent += p[k] * std::log(p[k]);
                    }
                }
                entropy[j * nx + i] = static_cast<float>(ent);
            }
        }
    }

    void removeContinuous(int nx, int ny, int npx, int npy, std::vector<float>& map, 
                          const std::function<double(double, int)>& func, int ord) {
        std::vector<float> maps(nx * ny, 0.0f);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                maps[j * nx + i] = static_cast<float>(func(map[j * npx + i], 1));
            }
        }

        if (ord == 5) {
            smoothImage55(nx, ny, maps, 2);
        } else if (ord == 4) {
            smoothImage55(nx, ny, maps, 1);
        } else {
            smoothImage33(nx, ny, maps);
        }

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                float val = static_cast<float>(func(maps[j * nx + i], -1));
                map[j * npx + i] -= val;
            }
        }
    }

    // ==========================================
    // Function: Perform F77-equivalent 2D complex FFT
    // Method: FFTPACK CFFT2F normalizes forward transforms by n1*n2; FFTW does not, so apply the scale here.
    // ==========================================
    void FFT2D(int n1, int n2, std::vector<std::complex<float>>& arr, int direction) {
        int fftw_dir = (direction == 1) ? FFTW_FORWARD : FFTW_BACKWARD;
        fftwf_complex* data = reinterpret_cast<fftwf_complex*>(arr.data());
        fftwf_plan plan = fftwf_plan_dft_2d(n2, n1, data, data, fftw_dir, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        if (direction == 1) {
            float norm = 1.0f / static_cast<float>(n1 * n2);
            for (auto& val : arr) {
                val *= norm;
            }
        }
    }

    // ==========================================
    // Function: Compute F77-equivalent shifted FFT power spectrum
    // Method: Use FFT2D's F77-equivalent forward normalization, then compute |FFT|^2 and fftshift.
    // ==========================================
    void getPower(int n1, int n2, const std::vector<float>& map, std::vector<float>& power, int smooth, double& pc) {
        int n1_2 = n1 / 2;
        int n2_2 = n2 / 2;

        std::vector<std::complex<float>> arr(n1 * n2);
        for (int i = 0; i < n1 * n2; ++i) {
            arr[i] = std::complex<float>(map[i], 0.0f);
        }

        FFT2D(n1, n2, arr, 1);

        power.resize(n1 * n2);
        for (int i = 0; i < n1; ++i) {
            int ii = (i + n1_2) % n1;
            for (int j = 0; j < n2; ++j) {
                int jj = (j + n2_2) % n2;
                float val = std::abs(arr[j * n1 + i]);
                power[jj * n1 + ii] = static_cast<float>(val * val);
            }
        }

        pc = power[n2_2 * n1 + n1_2];

        if (smooth == 1) {
            smoothImage55Hole(n1, n2, power);
        } else if (smooth == 2) {
            smoothImage55HoleLn(n1, n2, power);
        }
    }

    // ==========================================
    // Function: Draw a filled dot at a one-based image coordinate.
    // Method: Map F77 map(x,y) to row-major map[(y-1)*nx+(x-1)].
    // ==========================================
    void drawDot(int nx, int ny, std::vector<float>& map, double x, double y, double intensity, double thickness) {
        int ix = static_cast<int>(x + 0.5);
        int iy = static_cast<int>(y + 0.5);

        int i_start = static_cast<int>(ix - thickness + 0.5);
        int i_end = static_cast<int>(ix + thickness + 0.5);
        int j_start = static_cast<int>(iy - thickness + 0.5);
        int j_end = static_cast<int>(iy + thickness + 0.5);

        for (int i = i_start; i <= i_end; ++i) {
            for (int j = j_start; j <= j_end; ++j) {
                if (i >= 1 && i <= nx && j >= 1 && j <= ny) {
                    map[(j - 1) * nx + (i - 1)] = static_cast<float>(intensity);
                }
            }
        }
    }

    void drawLine(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity, double thickness) {
        double r = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        if (r > 0.0) {
            double cosx = (x2 - x1) / r;
            double cosy = (y2 - y1) / r;
            int ndots = static_cast<int>(r);
            if (ndots > 0) {
                double dr = r / ndots;
                double x = x1;
                double y = y1;
                drawDot(nx, ny, map, x, y, intensity, thickness);
                for (int i = 0; i < ndots; ++i) {
                    x += dr * cosx;
                    y += dr * cosy;
                    drawDot(nx, ny, map, x, y, intensity, thickness);
                }
            } else {
                drawDot(nx, ny, map, x1, y1, intensity, thickness);
            }
        }
    }

    void drawRectangle(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity, double thickness) {
        drawLine(nx, ny, map, x1, y1, x1, y2, intensity, thickness);
        drawLine(nx, ny, map, x2, y1, x2, y2, intensity, thickness);
        drawLine(nx, ny, map, x1, y1, x2, y1, intensity, thickness);
        drawLine(nx, ny, map, x1, y2, x2, y2, intensity, thickness);
    }

    // ==========================================
    // Function: Draw a filled axis-aligned box.
    // Method: Map F77 map(x,y) to row-major map[(y-1)*nx+(x-1)].
    // ==========================================
    void drawBoxFill(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity) {
        int ix1 = std::max(1, static_cast<int>(x1 + 0.5));
        int ix2 = std::min(nx, static_cast<int>(x2 + 0.5));
        int iy1 = std::max(1, static_cast<int>(y1 + 0.5));
        int iy2 = std::min(ny, static_cast<int>(y2 + 0.5));

        for (int i = ix1; i <= ix2; ++i) {
            for (int j = iy1; j <= iy2; ++j) {
                map[(j - 1) * nx + (i - 1)] = static_cast<float>(intensity);
            }
        }
    }

    void reverseColor(int nx, int ny, std::vector<float>& map) {
        for (float& val : map) {
            val = 255.0f - val;
        }
    }

    void drawShearExpo(int n, std::vector<float>& map, 
                       const std::vector<std::array<double, 4>>& pc, 
                       const std::vector<std::array<double, 5>>& sk, 
                       double intensity, double thickness) {
        map.assign(n * n, 0.0f);
        double xmin = 1e10;
        double xmax = -1e10;
        double ymin = 1e10;
        double ymax = -1e10;
        std::vector<std::array<double, 4>> valid_pc;
        valid_pc.reserve(pc.size());
        for (const auto& chip : pc) {
            if (chip[0] == 0.0 && chip[1] == 0.0 && chip[2] == 0.0 && chip[3] == 0.0) {
                continue;
            }
            valid_pc.push_back(chip);
            xmin = std::min({chip[0], chip[2], xmin});
            xmax = std::max({chip[0], chip[2], xmax});
            ymin = std::min({chip[1], chip[3], ymin});
            ymax = std::max({chip[1], chip[3], ymax});
        }

        int nc = static_cast<int>(valid_pc.size());
        if (nc == 0) {
            reverseColor(n, n, map);
            return;
        }

        double margin_val = 0.05;
        double dx = (xmax - xmin) * margin_val;
        double dy = (ymax - ymin) * margin_val;
        xmin -= dx;
        xmax += dx;
        ymin -= dy;
        ymax += dy;
        double ratiox = (n - 1.0) / (xmax - xmin);
        double ratioy = (n - 1.0) / (ymax - ymin);

        for (int i = 0; i < nc; ++i) {
            double x1 = (valid_pc[i][0] - xmin) * ratiox + 1.0;
            double x2 = (valid_pc[i][2] - xmin) * ratiox + 1.0;
            double y1 = (valid_pc[i][1] - ymin) * ratioy + 1.0;
            double y2 = (valid_pc[i][3] - ymin) * ratioy + 1.0;
            drawRectangle(n, n, map, x1, y1, x2, y2, intensity, thickness);
        }

        int ns = sk.size();
        if (ns == 0) {
            reverseColor(n, n, map);
            return;
        }

        auto sk_local = sk;
        for (int i = 0; i < ns; ++i) {
            sk_local[i][0] = (sk[i][0] - xmin) * ratiox + 1.0;
            sk_local[i][1] = (sk[i][1] - ymin) * ratioy + 1.0;
        }

        double tmp = ns;
        double dd = n / std::sqrt(tmp) * 0.25;

        double kmin = 10000.0;
        double kmax = -kmin;
        std::vector<double> e(ns);
        std::vector<double> tt(ns);

        for (int i = 0; i < ns; ++i) {
            kmin = std::min(kmin, sk_local[i][2]);
            kmax = std::max(kmax, sk_local[i][2]);
            e[i] = std::sqrt(sk_local[i][3] * sk_local[i][3] + sk_local[i][4] * sk_local[i][4]);
            tt[i] = e[i];
        }

        std::sort(tt.begin(), tt.end());
        tmp = std::max(std::abs(kmax), std::abs(kmin)) * 0.01;
        double ratio = intensity * 0.7 / (kmax - kmin + tmp);

        for (int i = 0; i < ns; ++i) {
            double tmp_val = (sk_local[i][2] - kmin) * ratio + 1.0;
            drawBoxFill(n, n, map, sk_local[i][0] - dd, sk_local[i][1] - dd, sk_local[i][0] + dd, sk_local[i][1] + dd, tmp_val);
        }

        ratio = dd / tt[std::max(0, ns * 3 / 4 - 1)];
        for (int i = 0; i < ns; ++i) {
            if (e[i] <= 0.0) continue;
            double cos2t = sk_local[i][3] / e[i];
            double cost = std::sqrt((1.0 + cos2t) * 0.5);
            double sint = std::sqrt((1.0 - cos2t) * 0.5);
            if (sk_local[i][4] > 0.0) {
                sint = -sint;
            }
            double x1 = sk_local[i][0] - e[i] * ratio * cost;
            double x2 = sk_local[i][0] + e[i] * ratio * cost;
            double y1 = sk_local[i][1] - e[i] * ratio * sint;
            double y2 = sk_local[i][1] + e[i] * ratio * sint;
            drawLine(n, n, map, x1, y1, x2, y2, intensity, 0.0);
        }

        reverseColor(n, n, map);
    }
}
