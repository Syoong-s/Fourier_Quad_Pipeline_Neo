#include "LensingConfig.hpp"
#include "RuntimeConfig.hpp"
#include "general/OutputLayout.hpp"
#include "process_main/PSFModel.hpp"
#include "process_main/PSFRecons.hpp"
#include "process_main/ProcessMainState.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

// ==========================================
// Function: Stop the direct-path regression on a failed requirement
// Method: Print one focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFRecons direct-path test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated synthetic Stage-6 output tree
// Method: Create a unique temporary root and remove it on scope exit.
// ==========================================
class TemporaryTree {
public:
    TemporaryTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_psfrecons_direct_"
                   + std::to_string(
                       std::chrono::steady_clock::now().time_since_epoch().count()))) {
        std::filesystem::create_directories(root_);
    }
    ~TemporaryTree() {
        PSFModel::pca_cache.clear();
        ProcessMain::state.clear();
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }
    const std::filesystem::path& root() const { return root_; }
private:
    std::filesystem::path root_;
};

// ==========================================
// Function: Write one finite local-polynomial PSF coefficient product
// Method: Match the four-field producer header and place a Gaussian in each constant term.
// ==========================================
void writeLocalCoefficients(const std::filesystem::path& root,
                            const std::string& prefix) {
    const std::filesystem::path path = OutputLayout::chipPath(
        root.string(), "stamps/dat_PsfFit", prefix, "_PSF_coe_local.dat");
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create local PSF coefficients");
    output << "1 1 0 0\n";
    const int ns = LensingConfig::ns;
    const int npl = LensingConfig::npl;
    const double center = 0.5 * static_cast<double>(ns - 1);
    for (int x = 0; x < ns; ++x) {
        for (int y = 0; y < ns; ++y) {
            const double radius2 = (x - center) * (x - center)
                                 + (y - center) * (y - center);
            output << std::exp(-0.2 * radius2);
            for (int term = 1; term < npl + 1; ++term) output << " 0";
            output << '\n';
        }
    }
}

// ==========================================
// Function: Initialize a PCA cache that requests polynomial-only fallback
// Method: Allocate the physical CCD range and mark every first component invalid.
// ==========================================
void initializeFallbackPcaCache(int nmax_chip) {
    PSFModel::pca_cache.components.assign(
        static_cast<std::size_t>(nmax_chip)
            * static_cast<std::size_t>(LensingConfig::nsns)
            * static_cast<std::size_t>(LensingConfig::n_pcs),
        -1.0e30);
    PSFModel::pca_cache.mean_psf.clear();
    PSFModel::pca_cache.poly_coefs.clear();
    PSFModel::pca_cache.data_loaded = true;
}

// ==========================================
// Function: Build a reordered, gapped exposure and its direct CCDNUM-7 products
// Method: Omit Science FITS creation so any positional/header lookup is guaranteed to fail.
// ==========================================
void prepareExposure(const TemporaryTree& tree) {
    const std::filesystem::path image1 =
        tree.root() / "science" / "exposure" / "exposure_1.fits";
    const std::filesystem::path image7 =
        tree.root() / "science" / "exposure" / "exposure_7.fits";
    const std::filesystem::path list = tree.root() / "exposure.list";
    {
        std::ofstream output(list);
        require(static_cast<bool>(output), "cannot create exposure list");
        output << image7.string() << '\n' << image1.string() << '\n';
    }
    ProcessMain::state.exposure_files = {list.string()};

    std::filesystem::create_directories(tree.root() / "stamps" / "dat_Rescale");
    {
        std::ofstream output(
            tree.root() / "stamps" / "dat_Rescale" / "exposure_factor.dat");
        require(static_cast<bool>(output), "cannot create rescale factor");
        output << "1\n";
    }

    std::filesystem::create_directories(tree.root() / "stamps" / "dat_StarComp");
    {
        std::ofstream output(
            tree.root() / "stamps" / "dat_StarComp"
            / "exposure_star_comp_expo.dat");
        require(static_cast<bool>(output), "cannot create StarComp product");
        output << "7 1 1\n";
        output << "100 100 2 0.1 0.2 2 0.1 0.2\n";
    }
    std::filesystem::create_directories(tree.root() / "stamps" / "dat_StarCompV2");
    writeLocalCoefficients(tree.root(), "exposure_7");
}

// ==========================================
// Function: Verify residual reconstruction uses StarComp CCDNUM directly
// Method: Process CCDNUM 7 with only two reordered list entries and inspect the output key.
// ==========================================
void testDirectPhysicalPath(const TemporaryTree& tree) {
    PSFRecons::plotResidualsV2(1);
    const std::filesystem::path output_path =
        tree.root() / "stamps" / "dat_StarCompV2"
        / "exposure_star_comp_expo_v2.dat";
    std::ifstream input(output_path);
    require(static_cast<bool>(input), "StarCompV2 output is missing");
    int ccdnum = 0;
    int nstar = 0;
    int valid = 0;
    require(static_cast<bool>(input >> ccdnum >> nstar >> valid),
            "StarCompV2 header is incomplete");
    require(ccdnum == 7 && nstar == 1 && valid == 1,
            "StarCompV2 did not preserve the physical CCDNUM key");
    double value = 0.0;
    for (int field = 0; field < 8; ++field) {
        require(static_cast<bool>(input >> value) && std::isfinite(value),
                "StarCompV2 row is incomplete or nonfinite");
    }
}

}  // namespace

// ==========================================
// Function: Run the Standard direct-CCDNUM PSFRecons regression
// Method: Initialize runtime geometry, force PCA fallback, and reconstruct a gapped identity.
// ==========================================
int main() {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    const int nmax_chip = config.lensing.nmax_chip;
    std::string error;
    require(RuntimeConfigStore::initialize(std::move(config), error),
            "cannot initialize runtime configuration: " + error);

    TemporaryTree tree;
    prepareExposure(tree);
    initializeFallbackPcaCache(nmax_chip);
    testDirectPhysicalPath(tree);
    std::cout << "PSFRecons direct-path tests passed\n";
    return EXIT_SUCCESS;
}
