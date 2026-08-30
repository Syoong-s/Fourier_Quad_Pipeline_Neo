#ifndef PSF_MODEL_HPP
#define PSF_MODEL_HPP

#include <string>
#include <vector>
#include <array>
#include "LensingConfig.hpp"
#include "process_main/LinearSolve.hpp"

namespace PSFModel {
    struct PcaCacheState {
        std::vector<double> components;
        std::vector<double> mean_psf;
        std::vector<float> poly_coefs;
        bool data_loaded = false;

        // ==========================================
        // Function: Release the loaded PCA model cache
        // Method: Clear every aligned array and reset the load marker together.
        // ==========================================
        void clear() {
            components.clear();
            mean_psf.clear();
            poly_coefs.clear();
            data_loaded = false;
        }
    };

    extern PcaCacheState pca_cache;

    // Load and broadcast PSF PCA components
    void initAndLoadAllPSF(const std::string& dirOutput);
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
        LinearSolve::SolveDiagnostics* diagnostics = nullptr,
        std::vector<double>* leverage = nullptr);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);
    void getPSFModelVeryLocal(const std::vector<float>& psfmap, double x,
                              double y, std::vector<float>& model,
                              double& dmax, int stride);
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
