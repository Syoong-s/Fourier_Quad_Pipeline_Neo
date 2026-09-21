#include "RuntimeConfig.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/PSFModelState.hpp"

#include <array>
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

namespace PSFModel {
void plotStars(
    int nchip, const std::vector<std::string>& imageFiles,
    const std::string& dirOutput, int nc,
    const std::vector<std::array<double, 4>>& p_chip,
    Internal::ExposurePSFState& state);
void makePSFLocalFit(
    int nchip, const std::vector<std::string>& imageFiles,
    const std::string& dirOutput, Internal::ExposurePSFState& state);
void getPSFModel(
    int ns, int npp, const std::vector<double>& coefficients,
    double x, double y, std::vector<float>& model,
    std::vector<float>& constant_model);
void getPowerAll(
    int nx, int ny, const std::vector<float>& power,
    std::array<double, 2>& ellipticity, double& size,
    float threshold_ratio);
}

namespace {

// ==========================================
// Function: Stop the exposure-table regression on a failed requirement
// Method: Print one focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF exposure-table CCDNUM test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated Stage-5 table-output tree
// Method: Create a unique native-Linux root and remove it on scope exit.
// ==========================================
class TemporaryTree {
public:
    TemporaryTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_psf_table_ccdnum_"
                   + std::to_string(
                       std::chrono::steady_clock::now().time_since_epoch().count()))) {
        std::filesystem::create_directories(root_ / "stamps" / "dat_StarInfo");
        std::filesystem::create_directories(root_ / "stamps" / "fits_PsfSrc");
        std::filesystem::create_directories(root_ / "stamps" / "dat_StarComp");
        std::filesystem::create_directories(
            root_ / "stamps" / "dat_PsfFit" / "exposure");
        std::filesystem::create_directories(
            root_ / "stamps" / "fits_StarCanP" / "exposure");
    }

    ~TemporaryTree() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    const std::filesystem::path& root() const { return root_; }

private:
    std::filesystem::path root_;
};

// ==========================================
// Function: Read every first-column chip key from one Stage-5 exposure table
// Method: Optionally discard the StarInfo header, then consume only zero-star rows.
// ==========================================
std::vector<int> readZeroStarKeys(
    const std::filesystem::path& path, bool has_header) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "expected exposure table is missing");
    if (has_header) {
        std::string header;
        require(static_cast<bool>(std::getline(input, header)),
                "StarInfo header is missing");
    }

    std::vector<int> keys;
    int ccdnum = 0;
    int nstar = 0;
    while (input >> ccdnum >> nstar) {
        keys.push_back(ccdnum);
        require(nstar == 0, "empty synthetic state produced a nonzero star count");
        std::string remainder;
        std::getline(input, remainder);
    }
    return keys;
}

// ==========================================
// Function: Verify Stage-5 exposure tables preserve shuffled physical CCDNUM
// Method: Publish empty-chip StarInfo and StarComp rows from dense state slots and inspect keys.
// ==========================================
void testExposureTables(const TemporaryTree& tree) {
    const std::vector<std::string> image_files = {
        (tree.root() / "science" / "exposure" / "exposure_7.fits").string(),
        (tree.root() / "science" / "exposure" / "exposure_1.fits").string(),
        (tree.root() / "science" / "exposure" / "exposure_3.fits").string(),
    };
    PSFModel::Internal::ExposurePSFState state(3);
    const std::vector<std::array<double, 4>> chip_bounds(3);

    PSFModel::plotStars(
        3, image_files, tree.root().string(), 0, chip_bounds, state);
    PSFModel::makePSFLocalFit(
        3, image_files, tree.root().string(), state);

    const std::vector<int> expected = {7, 1, 3};
    require(readZeroStarKeys(
                tree.root() / "stamps" / "dat_StarInfo"
                    / "exposure_star_info_expo.dat",
                true) == expected,
            "StarInfo first column follows dense/list position");
    require(readZeroStarKeys(
                tree.root() / "stamps" / "dat_StarComp"
                    / "exposure_star_comp_expo.dat",
                false) == expected,
            "StarComp first column follows dense/list position");
}

// ==========================================
// Function: Verify local-fit StarComp rows use the all-star fitted model
// Method: Publish a cached constant polynomial whose ordinary and analytic-LOO
//         shapes differ, then compare the live writer with both predictions.
// ==========================================
void testFullFitModelDiagnostics(const TemporaryTree& tree) {
    constexpr int star_count = LensingConfig::nstar_min_local;
    constexpr int ns = LensingConfig::ns;
    constexpr int npl = LensingConfig::npl;
    const std::vector<std::string> image_files = {
        (tree.root() / "science" / "exposure" / "exposure_7.fits").string()};
    PSFModel::Internal::ExposurePSFState state(1);
    auto& chip = state.chips[0];

    std::vector<float> fitted_power(static_cast<std::size_t>(ns) * ns, 0.0f);
    const int center = ns / 2;
    for (int y = 0; y < ns; ++y) {
        for (int x = 0; x < ns; ++x) {
            const double dx = static_cast<double>(x - center);
            const double dy = static_cast<double>(y - center);
            fitted_power[static_cast<std::size_t>(y) * ns + x] =
                static_cast<float>(std::exp(-dx * dx / 40.0 - dy * dy / 24.0));
        }
    }
    std::vector<float> observed_power = fitted_power;
    observed_power[static_cast<std::size_t>(center) * ns + center + 4] += 1.5f;

    chip.stars.reserve(star_count);
    chip.selection.resize(star_count);
    chip.fit.valid = true;
    chip.fit.initial_star_count = star_count;
    chip.fit.coefficients.assign(
        static_cast<std::size_t>(ns) * ns * (npl + 1), 0.0);
    for (int pixel = 0; pixel < ns * ns; ++pixel) {
        chip.fit.coefficients[static_cast<std::size_t>(pixel) * (npl + 1)] =
            fitted_power[pixel];
        chip.fit.coefficients[static_cast<std::size_t>(pixel) * (npl + 1) + npl] =
            fitted_power[pixel];
    }

    std::vector<float> stamp_cube;
    stamp_cube.reserve(static_cast<std::size_t>(star_count) * ns * ns);
    for (int star_index = 0; star_index < star_count; ++star_index) {
        PSFModel::Internal::ChipPSFState::StarRow row{};
        row[1] = 500.0 + star_index;
        row[2] = 900.0 + star_index;
        row[4] = 1.0;
        row[7] = 1.0;
        row[8] = 0.0;
        row[9] = 0.0;
        chip.stars.push_back(row);
        chip.fit.star_indices.push_back(star_index);
        chip.fit.leverage.push_back(0.5);
        stamp_cube.insert(
            stamp_cube.end(), observed_power.begin(), observed_power.end());
    }

    const std::filesystem::path cube_path =
        tree.root() / "stamps" / "fits_StarCanP" / "exposure"
        / "exposure_7_star_can_power.fits";
    require(FitsIO::writeStampCube(
                cube_path.string(), ns, ns, star_count, stamp_cube),
            "cannot write the synthetic fitted-star cube");
    PSFModel::makePSFLocalFit(
        1, image_files, tree.root().string(), state);

    std::ifstream input(
        tree.root() / "stamps" / "dat_StarComp"
            / "exposure_star_comp_expo.dat");
    int ccdnum = 0;
    int serialized_count = 0;
    int status = 0;
    require(static_cast<bool>(input >> ccdnum >> serialized_count >> status),
            "cannot read the synthetic StarComp chip header");
    require(ccdnum == 7 && serialized_count == star_count && status == 1,
            "synthetic StarComp header lost CCDNUM or fit status");

    double px = 0.0;
    double py = 0.0;
    double observed_size = 0.0;
    double observed_e1 = 0.0;
    double observed_e2 = 0.0;
    double serialized_size = 0.0;
    double serialized_e1 = 0.0;
    double serialized_e2 = 0.0;
    require(static_cast<bool>(
                input >> px >> py >> observed_size >> observed_e1 >> observed_e2
                      >> serialized_size >> serialized_e1 >> serialized_e2),
            "cannot read the first synthetic StarComp model row");

    const LensingRuntimeConfig& lensing = RuntimeConfigStore::get().lensing;
    const double normalized_x =
        2.0 * (px / static_cast<double>(lensing.chipnx)) - 1.0;
    const double normalized_y =
        2.0 * (py / static_cast<double>(lensing.chipny)) - 1.0;
    std::vector<float> full_model;
    std::vector<float> constant_model;
    PSFModel::getPSFModel(
        ns, npl, chip.fit.coefficients, normalized_x, normalized_y,
        full_model, constant_model);
    std::array<double, 2> expected_shape = {0.0, 0.0};
    double expected_size = 0.0;
    PSFModel::getPowerAll(
        ns, ns, full_model, expected_shape, expected_size, 0.02f);

    std::vector<float> old_loo_model(full_model.size(), 0.0f);
    for (std::size_t pixel = 0; pixel < full_model.size(); ++pixel) {
        old_loo_model[pixel] =
            2.0f * full_model[pixel] - observed_power[pixel];
    }
    std::array<double, 2> old_loo_shape = {0.0, 0.0};
    double old_loo_size = 0.0;
    PSFModel::getPowerAll(
        ns, ns, old_loo_model, old_loo_shape, old_loo_size, 0.02f);
    const bool fixture_distinguishes_loo =
        expected_size != old_loo_size
        || std::abs(expected_shape[0] - old_loo_shape[0]) > 1.0e-8
        || std::abs(expected_shape[1] - old_loo_shape[1]) > 1.0e-8;
    require(fixture_distinguishes_loo,
            "synthetic fixture does not distinguish full-fit from LOO output");
    require(std::abs(serialized_size - expected_size) < 1.0e-12
                && std::abs(serialized_e1 - expected_shape[0]) < 1.0e-12
                && std::abs(serialized_e2 - expected_shape[1]) < 1.0e-12,
            "StarComp row does not serialize the ordinary full-fit model");
}

}  // namespace

// ==========================================
// Function: Run Stage-5 persistent CCDNUM table regressions
// Method: Initialize runtime defaults and publish a shuffled zero-star exposure.
// ==========================================
int main() {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    std::string error;
    require(RuntimeConfigStore::initialize(std::move(config), error),
            "cannot initialize runtime configuration: " + error);
    TemporaryTree tree;
    testExposureTables(tree);
    testFullFitModelDiagnostics(tree);
    std::cout << "PSF exposure-table CCDNUM tests passed\n";
    return EXIT_SUCCESS;
}
