#include "general/NumericalRecipes.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace NumericalRecipes {

    namespace {

    struct RngState {
        int idum = -123;
        int iy = 0;
        std::array<int, 32> iv{};
        int gasdev_iset = 0;
        float gasdev_gset = 0.0f;
    };

    thread_local RngState rng_state;

    }  // namespace

    // ==========================================
    // Function: Initialize a distinct ran1 stream for one MPI rank
    // Method: Pack a monotonic clock count/rank and apply the F77 bijective fifth-power map.
    // ==========================================
    unsigned int initializeRan1Seed(int rank, int numProcs) {
        constexpr std::int64_t prime = 2147483647LL;
        constexpr std::int64_t limit = prime - 1LL;

        if (numProcs <= 0) {
            throw std::invalid_argument("Invalid MPI process count for ran1");
        }
        if (rank < 0 || rank >= numProcs) {
            throw std::invalid_argument("Invalid MPI rank for ran1");
        }

        const std::int64_t nproc = static_cast<std::int64_t>(numProcs);
        if (nproc > limit) {
            throw std::invalid_argument("Too many MPI ranks for ran1");
        }

        const std::int64_t slots = limit / nproc;
        const std::int64_t clockCount =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now().time_since_epoch()).count();
        std::int64_t clockSlot = clockCount % slots;
        if (clockSlot < 0) {
            clockSlot += slots;
        }

        // Encode rank in the residue modulo numProcs, so ranks cannot collide.
        const std::int64_t packed =
            clockSlot * nproc + static_cast<std::int64_t>(rank) + 1LL;

        // The fifth-power map is one-to-one because gcd(5, prime - 1) is one.
        const std::int64_t power2 = (packed * packed) % prime;
        const std::int64_t power4 = (power2 * power2) % prime;
        const std::int64_t seedValue = (power4 * packed) % prime;

        if (seedValue <= 0LL || seedValue >= prime) {
            throw std::runtime_error("Generated ran1 seed is outside the valid range");
        }

        const unsigned int seed = static_cast<unsigned int>(seedValue);
        seedRandom(seed);
        return seed;
    }

    // ==========================================
    // Function: Seed F77-equivalent random generator
    // Method: Normalize to ran1's valid Park-Miller state range and reset all saved RNG state.
    // ==========================================
    void seedRandom(unsigned int seed) {
        constexpr unsigned int max_valid_seed = 2147483646u;
        const unsigned int normalized_seed =
            (seed == 0u) ? 1u : ((seed - 1u) % max_valid_seed) + 1u;
        rng_state.idum = -static_cast<int>(normalized_seed);
        rng_state.iy = 0;
        rng_state.iv.fill(0);
        rng_state.gasdev_iset = 0;
        rng_state.gasdev_gset = 0.0f;
    }

    // ==========================================
    // Function: Uniform random number
    // Method: Exact C++ translation of F77 press.f ran1, including REAL return precision.
    // ==========================================
    double ran1() {
        constexpr int IA = 16807;
        constexpr int IM = 2147483647;
        constexpr float AM = 1.0f / static_cast<float>(IM);
        constexpr int IQ = 127773;
        constexpr int IR = 2836;
        constexpr int NTAB = 32;
        constexpr int NDIV = 1 + (IM - 1) / NTAB;
        constexpr float EPS = 1.2e-7f;
        constexpr float RNMX = 1.0f - EPS;

        int k = 0;
        if (rng_state.idum <= 0 || rng_state.iy == 0) {
            rng_state.idum = std::max(-rng_state.idum, 1);
            for (int j = NTAB + 7; j >= 0; --j) {
                k = rng_state.idum / IQ;
                rng_state.idum = IA * (rng_state.idum - k * IQ) - IR * k;
                if (rng_state.idum < 0) rng_state.idum += IM;
                if (j < NTAB) rng_state.iv[static_cast<size_t>(j)] = rng_state.idum;
            }
            rng_state.iy = rng_state.iv[0];
        }

        k = rng_state.idum / IQ;
        rng_state.idum = IA * (rng_state.idum - k * IQ) - IR * k;
        if (rng_state.idum < 0) rng_state.idum += IM;
        int j = rng_state.iy / NDIV;
        rng_state.iy = rng_state.iv[static_cast<size_t>(j)];
        rng_state.iv[static_cast<size_t>(j)] = rng_state.idum;
        float value = std::min(AM * static_cast<float>(rng_state.iy), RNMX);
        return static_cast<double>(value);
    }

    // ==========================================
    // Function: Gaussian random number
    // Method: Exact Box-Muller cache logic from F77 press.f gasdev using REAL intermediates.
    // ==========================================
    double gasdev() {
        if (rng_state.gasdev_iset == 0) {
            float v1 = 0.0f;
            float v2 = 0.0f;
            float rsq = 0.0f;
            do {
                v1 = 2.0f * static_cast<float>(ran1()) - 1.0f;
                v2 = 2.0f * static_cast<float>(ran1()) - 1.0f;
                rsq = v1 * v1 + v2 * v2;
            } while (rsq >= 1.0f || rsq == 0.0f);

            float fac = std::sqrt(-2.0f * std::log(rsq) / rsq);
            rng_state.gasdev_gset = v1 * fac;
            rng_state.gasdev_iset = 1;
            return static_cast<double>(v2 * fac);
        }

        rng_state.gasdev_iset = 0;
        return static_cast<double>(rng_state.gasdev_gset);
    }

    // ==========================================
    // Function: Pair of Gaussian random numbers
    // Method: Draw two sequential F77-equivalent gasdev() values.
    // ==========================================
    void gasdev2(double& x, double& y) {
        x = gasdev();
        y = gasdev();
    }

    void sort(std::vector<float>& arr) {
        std::sort(arr.begin(), arr.end());
    }

    void sort(int n, std::vector<float>& arr) {
        if (n > 0 && n <= static_cast<int>(arr.size())) {
            std::sort(arr.begin(), arr.begin() + n);
        }
    }

    void sortDoub(std::vector<double>& arr) {
        std::sort(arr.begin(), arr.end());
    }

    void sortDoub(int n, std::vector<double>& arr) {
        if (n > 0 && n <= static_cast<int>(arr.size())) {
            std::sort(arr.begin(), arr.begin() + n);
        }
    }

    void indexx(const std::vector<float>& arr, std::vector<int>& indx) {
        indexx(static_cast<int>(arr.size()), arr, indx);
    }

    void indexx(int n, const std::vector<float>& arr, std::vector<int>& indx) {
        indx.resize(n);
        std::iota(indx.begin(), indx.end(), 0);
        std::sort(indx.begin(), indx.end(), [&arr](int i1, int i2) {
            return arr[i1] < arr[i2];
        });
    }

    void sort2i(std::vector<float>& arr, std::vector<int>& brr) {
        sort2i(static_cast<int>(arr.size()), arr, brr);
    }

    void sort2i(int n, std::vector<float>& arr, std::vector<int>& brr) {
        if (n <= 0) return;
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&arr](int i1, int i2) {
            return arr[i1] < arr[i2];
        });
        std::vector<float> temp_arr(n);
        std::vector<int> temp_brr(n);
        for (int i = 0; i < n; ++i) {
            temp_arr[i] = arr[idx[i]];
            temp_brr[i] = brr[idx[i]];
        }
        for (int i = 0; i < n; ++i) {
            arr[i] = temp_arr[i];
            brr[i] = temp_brr[i];
        }
    }

    void getPeakWidth(const std::vector<float>& arr, float& p, float& sig, int& status, int direc) {
        int n = static_cast<int>(arr.size());
        if (n <= 1) {
            p = 0;
            sig = 0;
            status = 0;
            return;
        }

        std::vector<float> sorted_arr = arr;
        std::sort(sorted_arr.begin(), sorted_arr.end());

        if (direc > 0) {
            n = (n * 9) / 10;
        }

        constexpr int nb = 200;
        constexpr int smooth = 3;

        float arr_min = sorted_arr[0];
        float arr_max = sorted_arr[n - 1];
        float d = (arr_max - arr_min) / (nb - 1.0f);

        std::vector<float> posi(nb);
        std::vector<float> den(nb + 2, 0.0f); // 0 to nb+1
        std::vector<int> mark(nb + 2, 0);

        for (int i = 0; i < nb; ++i) {
            posi[i] = arr_min + i * d;
        }

        for (int i = 0; i < n; ++i) {
            float tmp = (sorted_arr[i] - arr_min) / d + 1.0f; // 1-based tmp
            int ip = static_cast<int>(tmp + 0.5f);
            int start_j = std::max(ip - 4 * smooth, 1);
            int end_j = std::min(ip + 4 * smooth, nb);
            for (int j = start_j; j <= end_j; ++j) {
                float val = (tmp - j) / smooth;
                den[j] += std::exp(-0.5f * val * val);
            }
        }

        float vol1 = 0;
        float vol2 = 0;
        int bb[2][2] = {{0, 0}, {0, 0}};

        for (int i = 1; i <= nb; ++i) {
            if (mark[i] > 0) continue;
            mark[i] = 1;
            if (den[i] > den[i - 1] && den[i] > den[i + 1]) {
                int bound1 = i;
                int bound2 = i;
                int ip = i;
                float thresh = den[i] * 0.5f;

                while (bound2 < nb && (den[bound2 + 1] > thresh || den[bound2 + 1] < den[bound2])) {
                    bound2++;
                    mark[bound2] = 2;
                    if (den[bound2] * 0.5f > thresh) {
                        thresh = den[bound2] * 0.5f;
                        ip = bound2;
                    }
                }
                while (bound1 > 1 && mark[bound1 - 1] != 2 && (den[bound1 - 1] > thresh || den[bound1 - 1] < den[bound1])) {
                    bound1--;
                    if (den[bound1] * 0.5f > thresh) {
                        thresh = den[bound1] * 0.5f;
                        ip = bound1;
                    }
                }

                thresh *= 0.5f;
                int b2 = ip;
                while (b2 < bound2 && den[b2] > thresh) {
                    b2++;
                }
                int b1 = ip;
                while (b1 > bound1 && den[b1] > thresh) {
                    b1--;
                }
                bound1 = b1;
                bound2 = b2;

                float vol = 0.0f;
                for (int j = bound1; j <= bound2; ++j) {
                    vol += den[j];
                }

                if (vol > vol1) {
                    vol2 = vol1;
                    bb[1][0] = bb[0][0];
                    bb[1][1] = bb[0][1];
                    vol1 = vol;
                    bb[0][0] = bound1;
                    bb[0][1] = bound2;
                } else if (vol > vol2) {
                    vol2 = vol;
                    bb[1][0] = bound1;
                    bb[1][1] = bound2;
                }
            }
        }

        int bound1 = 0, bound2 = 0;
        // Map 1-based Fortran arrays to 0-based indices for bb
        // Fortran: bb(1,1) -> C++: bb[0][0]
        if (direc > 0) {
            if (bb[0][0] > bb[1][0]) {
                bound1 = bb[0][0];
                bound2 = bb[0][1];
                status = 1;
            } else {
                bound1 = bb[1][0];
                bound2 = bb[1][1];
                status = 2;
            }
        } else {
            if (bb[0][0] > bb[1][0]) {
                bound1 = bb[1][0];
                bound2 = bb[1][1];
                status = 2;
            } else {
                bound1 = bb[0][0];
                bound2 = bb[0][1];
                status = 1;
            }
        }

        // Adjust bounds to 0-based index for posi
        int b1_idx = bound1 - 1;
        int b2_idx = bound2 - 1;
        if (b1_idx < 0) b1_idx = 0;
        if (b2_idx >= nb) b2_idx = nb - 1;

        double sum_p = 0;
        double sum_sig = 0;
        int num = 0;
        for (int i = 0; i < n; ++i) {
            if (sorted_arr[i] < posi[b1_idx] || sorted_arr[i] > posi[b2_idx]) continue;
            sum_p += sorted_arr[i];
            sum_sig += sorted_arr[i] * sorted_arr[i];
            num++;
        }

        if (num > 0) {
            p = static_cast<float>(sum_p / num);
            sig = static_cast<float>(std::max(std::sqrt(std::max(sum_sig / num - p * p, 0.0)), 0.02 * (posi[b2_idx] - posi[b1_idx])));
        } else {
            p = 0;
            sig = 0;
        }
    }

    void getPeakWidthLowSide(const std::vector<float>& arr, float& p, float& sig) {
        int n = static_cast<int>(arr.size());
        if (n <= 4) {
            p = 0;
            sig = 0;
            return;
        }

        std::vector<float> sorted_arr = arr;
        std::sort(sorted_arr.begin(), sorted_arr.end());

        std::vector<float> den(n);
        for (int i = 1; i < n - 1; ++i) {
            float val = (sorted_arr[i + 1] - sorted_arr[i - 1]) / 2.0f;
            den[i] = val * val;
        }
        float val_first = sorted_arr[1] - sorted_arr[0];
        den[0] = val_first * val_first;
        float val_last = sorted_arr[n - 1] - sorted_arr[n - 2];
        den[n - 1] = val_last * val_last;

        std::vector<float> den2(n, 0.0f);
        for (int i = 2; i < n - 2; ++i) {
            float sum_val = 0.0f;
            for (int j = i - 2; j <= i + 2; ++j) {
                sum_val += den[j];
            }
            den2[i] = 1.0f / std::sqrt(sum_val / 5.0f);
        }

        den2[0] = 1.0f / std::sqrt((den[0] + den[1] + den[2]) / 3.0f);
        den2[1] = 1.0f / std::sqrt((den[0] + den[1] + den[2] + den[3]) / 4.0f);
        den2[n - 1] = 1.0f / std::sqrt((den[n - 1] + den[n - 2] + den[n - 3]) / 3.0f);
        den2[n - 2] = 1.0f / std::sqrt((den[n - 1] + den[n - 2] + den[n - 3] + den[n - 4]) / 4.0f);

        int ip = 0;
        p = den2[0];
        for (int i = 1; i < n; ++i) {
            if (den2[i] > p) {
                p = den2[i];
                ip = i;
            }
        }

        float thresh = p * 0.5f;
        int i_idx = ip + 1;
        for (; i_idx < n; ++i_idx) {
            if (den2[i_idx] < thresh) break;
        }
        if (i_idx >= n) i_idx = n - 1;

        p = sorted_arr[ip];
        sig = sorted_arr[i_idx] - sorted_arr[ip];
    }

    double gammln(double xx) {
        return std::lgamma(xx);
    }

    // Helper functions for gammq (incomplete gamma)
    static void gser(double& gamser, double a, double x, double& gln) {
        gln = std::lgamma(a);
        if (x <= 0.0) {
            gamser = 0.0;
            return;
        }
        double sum = 1.0 / a;
        double del = sum;
        double ap = a;
        for (int n = 1; n <= 100; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 3e-7) {
                gamser = sum * std::exp(-x + a * std::log(x) - gln);
                return;
            }
        }
        gamser = sum * std::exp(-x + a * std::log(x) - gln); // fallback
    }

    static void gcf(double& gammcf, double a, double x, double& gln) {
        gln = std::lgamma(a);
        double gold = 0.0;
        double a0 = 1.0;
        double a1 = x;
        double b0 = 0.0;
        double b1 = 1.0;
        double fac = 1.0;
        for (int n = 1; n <= 100; ++n) {
            double an = n;
            double ana = an - a;
            a0 = (a1 + a0 * ana) * fac;
            b0 = (b1 + b0 * ana) * fac;
            double anf = an * fac;
            a1 = x * a0 + anf * a1;
            b1 = x * b0 + anf * b1;
            if (b1 != 0.0) {
                fac = 1.0 / b1;
                double g = a1 * fac;
                if (std::abs((g - gold) / g) < 3e-7) {
                    gammcf = std::exp(-x + a * std::log(x) - gln) * g;
                    return;
                }
                gold = g;
            }
        }
        gammcf = std::exp(-x + a * std::log(x) - gln) * (a1 * fac); // fallback
    }

    double gammq(double a, double x) {
        if (x < 0.0 || a <= 0.0) return 0.0;
        double gln;
        if (x < a + 1.0) {
            double gamser;
            gser(gamser, a, x, gln);
            return 1.0 - gamser;
        } else {
            double gammcf;
            gcf(gammcf, a, x, gln);
            return gammcf;
        }
    }
}
