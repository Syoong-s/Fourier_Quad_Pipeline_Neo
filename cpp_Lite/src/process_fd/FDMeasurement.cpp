#include "process_fd/FDMeasurement.hpp"
#include "FDConfig.hpp"
#include "process_fd/QuadraticFitting.hpp"
#include "process_main/NumericalRecipes.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace fc = FDConfig;

FDMeasurement::FDMeasurement(int rank, int num_procs)
    : rank_(rank), num_procs_(num_procs), xbin_(fc::NMAX, 0.0) {}

// ==========================================
// Function: statis
// Method:
//   SWSE mode (3): c = 2*sum(y*ww)/sum(ww), dc = 0  (ratio estimator)
//   PDF  mode: equal-probability binning -> chi2 sign test scan ->
//              quadratic fitting -> c_best, sigma = 1/sqrt(2*|a1|)
// ==========================================
void FDMeasurement::statis(int n, int nt, const float* x, const float* de,
                            int nbin, float& c, float& dc) {
    if constexpr (fc::FD_USE_PDF_STATIS) {
        // ---- PDF statis (Mode 1 & 2): chi2 sign test + quadratic fitting ----
        std::fill(xbin_.begin(), xbin_.end(), 0.0);

        // Equal-probability boundary initialization
        int effnode = (rank_ == 0) ? 1 : 0;
        std::vector<float> xx(n);
        if (rank_ != 0) {
            if (n < (nbin + 1)) {
                effnode = 1;
            } else {
                int dn = n / (nbin + 1);
                for (int i = 0; i < n; ++i) xx[i] = std::fabs(x[i]);
                std::sort(xx.begin(), xx.begin() + n);
                for (int i = 0; i < nbin; ++i)
                    xbin_[i] = xx[dn * (i + 1) - 1];  // 0-based: xx[dn*(i+1)-1]
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<float> xbin_tmp(fc::NMAX, 0.0);
        MPI_Reduce(xbin_.data(), xbin_tmp.data(), fc::NMAX, MPI_FLOAT,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        int effnodet = 0;
        MPI_Reduce(&effnode, &effnodet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank_ == 0) {
            int eff = num_procs_ - effnodet;
            if (eff > 0)
                for (int i = 0; i < fc::NMAX; ++i)
                    xbin_[i] = xbin_tmp[i] / float(eff);
        }
        MPI_Bcast(xbin_.data(), fc::NMAX, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Boundary iteration for balance
        for (int i = 0; i < nbin; ++i) {
            float vc = sourceAccumulate(n, nt, xx.data(), i, nbin);
            if (std::fabs(vc) < 0.03) continue;
            float xb1, xb2;
            if (vc > 0.0) {
                if (i - 1 < 0) xb1 = xbin_[i] * 1.5 - xbin_[i + 1] * 0.5;
                else xb1 = (xbin_[i] + xbin_[i - 1]) * 0.5;
                xb2 = xbin_[i];
            } else {
                xb1 = xbin_[i];
                if (i + 1 >= nbin) xb2 = xbin_[i] * 1.5 - xbin_[i - 1] * 0.5;
                else xb2 = (xbin_[i] + xbin_[i + 1]) * 0.5;
            }
            int change = 1, iter_count = 0;
            while (change == 1) {
                iter_count++;
                if (iter_count > 20) break;
                change = 0;
                xbin_[i] = xb1;
                float v1 = sourceAccumulate(n, nt, xx.data(), i, nbin);
                if (std::fabs(v1) < 0.03) break;
                xbin_[i] = xb2;
                float v2 = sourceAccumulate(n, nt, xx.data(), i, nbin);
                if (std::fabs(v2) < 0.03) break;
                change = 1;
                if (v1 * v2 < 0.0) {
                    if (std::fabs(v1) > std::fabs(v2)) xb1 = xb1 * 0.75 + xb2 * 0.25;
                    else xb2 = xb2 * 0.75 + xb1 * 0.25;
                } else {
                    if (v1 > 0.0) {
                        float tmp = xb1;
                        xb1 = xb1 * 2.0 - xb2;
                        xb2 = tmp * 1.3 - xb2 * 0.3;
                    } else {
                        float tmp = xb2;
                        xb2 = xb2 * 2.0 - xb1;
                        xb1 = tmp * 1.3 - xb1 * 0.3;
                    }
                }
            }
        }

        // Coarse chi2 scan to narrow [c1, c2]
        c = 0.0; dc = 0.2;
        float c1 = c - dc, c2 = c + dc;
        int change = 1;
        float thresh = 20.0;
        while (change == 1) {
            change = 0;
            float cc = (c1 + c2) * 0.5;
            float v1 = chi2SignTest(n, x, de, nbin, (c1 + cc) * 0.5);
            float v2 = chi2SignTest(n, x, de, nbin, (c2 + cc) * 0.5);
            float vc = chi2SignTest(n, x, de, nbin, cc);
            if (v1 > vc + thresh) { c1 = (c1 + cc) * 0.5; change = 1; }
            if (v2 > vc + thresh) { c2 = (c2 + cc) * 0.5; change = 1; }
        }

        // Fine grid + quadratic fitting
        dc = (c2 - c1) / float(fc::NMAX - 1);
        std::vector<float> cdata(2 * fc::NMAX);
        for (int i = 0; i < fc::NMAX; ++i) {
            float ci = c1 + i * dc;
            cdata[2 * i] = ci;
            cdata[2 * i + 1] = chi2SignTest(n, x, de, nbin, ci);
        }
        std::vector<float> a(3, 0.0);
        QuadraticFitting::fit(fc::NMAX, cdata, a);
        if (std::fabs(a[0]) > 1e-30) {
            c = -a[1] / (2.0f * a[0]);
            dc = 1.0f / std::sqrt(2.0f * std::fabs(a[0]));
        } else {
            c = 0.0;
            dc = 0.0;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        // ---- SWSE statis (Mode 3): simple ratio estimator ----
        float nume_loc = 0.0, deno_loc = 0.0;
        for (int i = 0; i < n; ++i) {
            nume_loc += x[i];
            deno_loc += de[i];
        }
        float nume_all = 0.0, deno_all = 0.0;
        MPI_Allreduce(&nume_loc, &nume_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&deno_loc, &deno_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        c = (deno_all != 0.0) ? 2.0f * nume_all / deno_all : 0.0f;
        dc = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// ==========================================
// Function: chi2SignTest
// Method: Non-parametric sign-test chi2: sum over bins of D^2/(2N),
//         where D = N+ - N- per bin.
// ==========================================
float FDMeasurement::chi2SignTest(int n, const float* x, const float* de,
                                   int nbin, float c) {
    std::vector<float> diff(fc::NMAX, 0.0), summ(fc::NMAX, 0.0);
    std::vector<float> diff_t(fc::NMAX, 0.0), summ_t(fc::NMAX, 0.0);

    if (rank_ != 0) {
        for (int i = 0; i < n; ++i) {
            float y = x[i] - c * de[i];
            float abs_y = std::fabs(y);
            int ibin;
            if (abs_y > xbin_[nbin - 1]) {
                ibin = nbin;  // overflow bin (0-based)
            } else {
                ibin = 0;
                while (ibin < nbin - 1 && abs_y > xbin_[ibin]) ibin++;
            }
            summ[ibin] += 1.0;
            if (y > 0) diff[ibin] += 1.0;
            if (y < 0) diff[ibin] -= 1.0;
        }
    }

    MPI_Reduce(diff.data(), diff_t.data(), fc::NMAX, MPI_FLOAT, MPI_SUM,
               0, MPI_COMM_WORLD);
    MPI_Reduce(summ.data(), summ_t.data(), fc::NMAX, MPI_FLOAT, MPI_SUM,
               0, MPI_COMM_WORLD);

    float temp = 0.0;
    if (rank_ == 0) {
        for (int ibin = 0; ibin <= nbin; ++ibin) {
            if (summ_t[ibin] > 0.0)
                temp += diff_t[ibin] * diff_t[ibin] / (2.0f * summ_t[ibin]);
        }
    }
    MPI_Bcast(&temp, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    return temp;
}

// ==========================================
// Function: sourceAccumulate
// Method: Bin-balance metric v = N_bin/N_tot - 1/(nbin+1)
// ==========================================
float FDMeasurement::sourceAccumulate(int n, int ntot, const float* x,
                                       int sbin, int nbin) {
    int summ = 0;
    if (rank_ != 0) {
        for (int i = 0; i < n; ++i) {
            int ibin;
            if (x[i] > xbin_[nbin - 1]) ibin = nbin;
            else {
                ibin = 0;
                while (ibin < nbin - 1 && x[i] > xbin_[ibin]) ibin++;
            }
            if (ibin == sbin) summ++;
        }
    }
    int summt = 0;
    MPI_Allreduce(&summ, &summt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return float(summt) / float(ntot) - 1.0f / float(nbin + 1);
}

// ==========================================
// Function: plotComparison
// Method: Spatial bin by field distortion, then per bin:
//   Jack: jackknife resampling -> c from ratio, sigma from jackknife variance
//   PDF:  single statis call -> c from chi2 fitting, sigma = 1/sqrt(2|a1|)
// ==========================================
void FDMeasurement::plotComparison(int n, int nbin,
                                    const std::vector<float>& x,
                                    const std::vector<float>& y,
                                    const std::vector<float>& de,
                                    const std::vector<int>& labels,
                                    int nb, int njack,
                                    std::vector<float>& arr) {
    float xmin = -float(fc::gf_lim);
    float xmax = float(fc::gf_lim);
    float dx = (xmax - xmin) / float(nb);

    // Temporaries for bin sources
    std::vector<float> yy(n), dd(n);
    std::vector<int> lbl(n);

    for (int i = 0; i < nb; ++i) {
        arr[i * 4 + 0] = xmin + dx * (i + 0.5);  // gf_center
        int is = 0;
        for (int j = 0; j < n; ++j) {
            if (x[j] >= xmin + dx * i && x[j] < xmin + dx * (i + 1)) {
                yy[is] = y[j];
                dd[is] = de[j];
                lbl[is] = labels[j];
                is++;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int ntot = 0;
        MPI_Allreduce(&is, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if constexpr (!fc::FD_USE_JACKKNIFE) {
            // ---- Mode 1 (PDF_SIGMA): single statis call, c and sigma from PDF ----
            float c = 0.0, dc = 0.0;
            statis(is, ntot, yy.data(), dd.data(), nbin - 1, c, dc);
            arr[i * 4 + 1] = c;    // c_best
            arr[i * 4 + 2] = dc;   // sigma (PDF)
        } else {
            // ---- Mode 2 & 3: jackknife resampling for sigma ----
            std::vector<float> yy_jk(n), dd_jk(n);
            std::vector<float> gg_jk(njack, 0.0);
            std::vector<float> w_jk(njack, 0.0);

            for (int ijack = 0; ijack < njack; ++ijack) {
                int igal_jk = 0;
                for (int igal = 0; igal < is; ++igal) {
                    if (lbl[igal] == ijack) continue;  // leave-one-out
                    yy_jk[igal_jk] = yy[igal];
                    dd_jk[igal_jk] = dd[igal];
                    igal_jk++;
                }
                float ggpt = 0.0, sspt = 0.0;
                statis(igal_jk, ntot, yy_jk.data(), dd_jk.data(),
                       fc::FD_USE_PDF_STATIS ? (nbin - 1) : 0, ggpt, sspt);
                gg_jk[ijack] = ggpt;
                int w_p = igal_jk;
                int w_tot = 0;
                MPI_Allreduce(&w_p, &w_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                w_jk[ijack] = float(ntot - w_tot);
            }

            // Full sample
            float gn = 0.0, sspt = 0.0;
            statis(is, ntot, yy.data(), dd.data(),
                   fc::FD_USE_PDF_STATIS ? (nbin - 1) : 0, gn, sspt);

            // Count non-empty regions
            int g_prime = 0;
            for (int ijack = 0; ijack < njack; ++ijack)
                if (w_jk[ijack] > 0.0) g_prime++;

            // Jackknife variance
            float exterm = 0.0;
            if (g_prime > 0) {
                for (int ijack = 0; ijack < njack; ++ijack)
                    if (w_jk[ijack] > 0.0)
                        exterm += (1.0f - w_jk[ijack] / float(ntot)) * gg_jk[ijack];
            }
            float summ = 0.0;
            if (g_prime >= 2) {
                for (int ijack = 0; ijack < njack; ++ijack) {
                    if (w_jk[ijack] > 0.0) {
                        float h_j = float(ntot) / w_jk[ijack];
                        float t2 = ((h_j - float(g_prime)) * gn
                                    - (h_j - 1.0f) * gg_jk[ijack] + exterm);
                        t2 = t2 * t2 / (h_j - 1.0f);
                        summ += t2;
                    }
                }
                sspt = std::sqrt(summ / float(g_prime));
            } else {
                sspt = 0.0;
            }
            arr[i * 4 + 1] = gn;     // c_best
            arr[i * 4 + 2] = sspt;   // sigma (jackknife)
        }

        if (rank_ == 0) {
            arr[i * 4 + 3] = float(ntot);
            std::cout << i << "  " << arr[i * 4 + 0] << "  "
                      << arr[i * 4 + 1] << "  " << arr[i * 4 + 2] << "  "
                      << arr[i * 4 + 3] << std::endl;
        }
    }
}
