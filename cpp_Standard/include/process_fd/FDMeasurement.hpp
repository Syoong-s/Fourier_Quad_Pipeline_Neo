#ifndef FD_MEASUREMENT_HPP
#define FD_MEASUREMENT_HPP

#include "process_fd/FDData.hpp"

#include <vector>

// ==========================================
// FDMeasurement - spatial binning + shear recovery
// Method: For each spatial bin (by field distortion), recover the mean
//         shear and its uncertainty using either:
//   - PDF mode: chi2 sign test + quadratic fitting
//   - Jackknife mode: ratio estimator + jackknife variance
// ==========================================
class FDMeasurement {
public:
    FDMeasurement();

    // plot_comparison_MPI: bin sources by gf, measure per bin
    // arr[nb][4]: {gf_center, c_best, sigma, ntot}
    void plotComparison(int n, int nbin,
                        const std::vector<float>& x,
                        const std::vector<float>& y,
                        const std::vector<float>& de,
                        const std::vector<int>& labels,
                        int nb, int njack,
                        std::vector<float>& arr);

private:
    std::vector<float> xbin_;  // equal-probability bin boundaries (chi2_mpi_pass)

    // statis_MPI: compute best-fit c and uncertainty dc
    void statis(int n, int nt, const float* x, const float* de,
                int nbin, float& c, float& dc);

    // chi2_MPI: sign-test chi2 for a trial c
    float chi2SignTest(int n, const float* x, const float* de,
                       int nbin, float c);

    // source_accumulate: bin balance metric for boundary iteration
    float sourceAccumulate(int n, int ntot, const float* x,
                           int sbin, int nbin);
};

#endif  // FD_MEASUREMENT_HPP
