#ifndef PSF_MODEL_HPP
#define PSF_MODEL_HPP

#include <string>
#include <vector>
#include <array>
#include "LensingConfig.hpp"
#include "LinearSolve.hpp"

namespace PSFModel {
    // cpp_lite: the global PSF PCA storage (global_components / global_mean_psf /
    // global_poly_coefs / is_data_loaded, initAndLoadAllPSF, freePSFMemory, their flat-layout
    // index helpers, PSF_rescale, PSF_unscale and getPSFModelVeryLocal) belonged exclusively to
    // the PSF_Ms=1 / PSF_type=2 paths and has been removed together with PSFRecons.

    // Stage 5 main entry
    void procPSF(int iexpo);

    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics = nullptr,
        std::vector<double>* leverage = nullptr);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);
    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio);
}

#endif // PSF_MODEL_HPP
