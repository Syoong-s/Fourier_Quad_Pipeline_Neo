#ifndef PSF_MODEL_HPP
#define PSF_MODEL_HPP

#include <string>
#include <vector>
#include <array>
#include "LensingConfig.hpp"
#include "LinearSolve.hpp"

namespace PSFModel {
    // Global PSF PCA components and coefficients loaded from disk (replacing psf_storage_mod)
    extern std::vector<double> global_components; // [NMAX_CHIP][nsns][n_pcs]
    extern std::vector<double> global_mean_psf;   // [NMAX_CHIP][nsns]
    extern std::vector<float> global_poly_coefs;  // [NMAX_CHIP][2][2][n_pcs][npp6th]
    extern bool is_data_loaded;

    // Load and broadcast PSF PCA components
    void initAndLoadAllPSF(const std::string& dirOutput, int myRank);
    void freePSFMemory();

    // Stage 5 main entry
    void procPSF(int iexpo);

    // Coordinate scaling helpers (used by both PSFModel and PSFRecons)
    void PSF_rescale(std::vector<float>& resi, float refactor);
    void PSF_unscale(std::vector<float>& model_in, float scale_factor);

    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics = nullptr);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);
    void getPSFModelVeryLocal(const std::vector<float>& psfmap, double x, double y, std::vector<float>& model, double& dmax, int stride = LensingConfig::npx);
    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio);

    // Helpers to access global components flat layout
    inline size_t getCompIndex(int ccd_idx, int k_idx, int j_idx) {
        return (ccd_idx * LensingConfig::nsns + k_idx) * LensingConfig::n_pcs + j_idx;
    }

    inline size_t getMeanIndex(int ccd_idx, int k_idx) {
        return ccd_idx * LensingConfig::nsns + k_idx;
    }

    inline size_t getPolyIndex(int ccd_idx, int bx_idx, int by_idx, int u_idx, int j_idx) {
        return (((ccd_idx * 2 + bx_idx) * 2 + by_idx) * LensingConfig::n_pcs + u_idx) * LensingConfig::npp6th + j_idx;
    }
}

#endif // PSF_MODEL_HPP
