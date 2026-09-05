#ifndef PSF_RECONS_HPP
#define PSF_RECONS_HPP

#include <string>
#include <vector>
#include <array>
#include "process_main/LinearSolve.hpp"

namespace PSFRecons {
    namespace Internal {
        // ==========================================
        // Function: Build the physical-CCD to Science-list index for one exposure
        // Method: Read every Science FITS CCDNUM, reject invalid or duplicate
        //         identities, and leave genuinely absent CCD slots at -1.
        // ==========================================
        std::vector<int> buildChipImageIndex(
            const std::vector<std::string>& image_files,
            int max_chip_id);

        // ==========================================
        // Function: Resolve one physical CCD through a flattened exposure index
        // Method: Return its Science path, or null when that CCD is absent.
        // ==========================================
        const std::string* indexedChipImage(
            const std::vector<std::string>& image_files,
            const std::vector<int>& flattened_indices,
            int exposure_index, int chip_id, int max_chip_id);
    }

    // Stage 6 main entry: coordinates PSF fitting and reconstruction across chips and exposures
    void chipPSFRecons(int nexpo);

    // ==========================================
    // Function: Fit residual PCA for one physical CCD
    // Method: Reuse the broadcast CCDNUM-to-Science index for both residual passes.
    // ==========================================
    void chipResPCAFit(int ichip, int nexpo,
                       const std::vector<int>& chip_image_indices,
                       int max_chip_id);

    // Plot residuals and map modified residuals for a specific exposure
    void plotResidualsV2(int iexpo);

    // Hierarchical PSF model reconstruction at a specific position (x, y) on a CCD chip
    void getPSFModelHierarchical(int i_ccd, double x, double y, float refactor, 
                                 const std::vector<double>& local_coe, std::vector<float>& psf_model);

    // PSF interpolation with covariance matrix (future use: error propagation)
    LinearSolve::SolveStatus itpNormPSFCov(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        std::vector<float>& sigmarr, std::vector<float>& comat,
        LinearSolve::SolveDiagnostics* diagnostics = nullptr);

    // Polynomial interpolation of 6th degree
    LinearSolve::SolveStatus interpolate_6th(
        int nsam, const std::vector<float>& x, const std::vector<float>& y,
        const std::vector<float>& z, int nc, std::vector<float>& coef,
        LinearSolve::SolveDiagnostics* diagnostics = nullptr);
}

#endif // PSF_RECONS_HPP
