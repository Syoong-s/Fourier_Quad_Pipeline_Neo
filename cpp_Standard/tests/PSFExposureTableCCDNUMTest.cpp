#include "RuntimeConfig.hpp"
#include "process_main/PSFModelState.hpp"

#include <array>
#include <chrono>
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
    std::cout << "PSF exposure-table CCDNUM tests passed\n";
    return EXIT_SUCCESS;
}
