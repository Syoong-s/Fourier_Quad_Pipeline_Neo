#ifndef PRE_PROCESS_LEGACY_HPP
#define PRE_PROCESS_LEGACY_HPP

#include <vector>

namespace PreProcessLegacy {

// ==========================================
// Function: Apply the historical F77 two-level background estimator
// Method: Use random middle-third rough fitting, random rank-median blocks,
//         legacy percentile rejection, and serialize the total model in the
//         normalized coordinate basis consumed by the current Stage 3.
// ==========================================
bool setBackground(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    int nx,
    int ny,
    std::vector<float>& image,
    int blocksize,
    int nct,
    int ncx,
    std::vector<double>& normalized_coefficients);

// ==========================================
// Function: Apply the historical F77 random-triple sigma estimator
// Method: Fit the middle third of 2000 squared-difference samples on absolute
//         one-based coordinates and apply the unscaled positive linear plane.
// ==========================================
bool setSig(
    int x_start,
    int x_end,
    int y_start,
    int y_end,
    int nx,
    int ny,
    std::vector<float>& image,
    double& aa,
    double& bb,
    double& cc);

}  // namespace PreProcessLegacy

#endif  // PRE_PROCESS_LEGACY_HPP
