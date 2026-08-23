#include "CatalogCombiner.hpp"
#include "FitsIO.hpp"
#include "LensingConfig.hpp"
#include "OutputLayout.hpp"
#include "UniversalUtils.hpp"
#include "Universalblock.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

std::vector<std::string> EXPO_FILE;

namespace ExposureInfo {
std::vector<float> expo_para;
}

namespace {

// ==========================================
// Function: Stop the catalog-lifecycle test on a failed requirement
// Method: Print one focused diagnostic and terminate with failure status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "CatalogCombiner lifecycle test failed: " << message
                  << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated Stage-9 input and output tree
// Method: Create production-layout chip paths and remove them after the test.
// ==========================================
class TemporaryCatalogTree {
public:
    // ==========================================
    // Function: Initialize the synthetic exposure tree
    // Method: Create a minimal chip FITS and every parent directory used by Stage 9.
    // ==========================================
    TemporaryCatalogTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_catalog_lifecycle_" + std::to_string(::getpid()))),
          image_file_((root_ / "input" / "exposure_1.fits").string()),
          prefix_(UniversalUtils::getPrefix(image_file_)),
          shear_file_(OutputLayout::chipPath(
              root_.string(), "stamps/dat_Shear", prefix_, "_shear.dat")),
          orig_file_(OutputLayout::chipPath(
              root_.string(), "stamps/cat_Orig", prefix_, "_orig.cat")),
          output_file_((root_ / "result" / "exposure_all.cat").string()) {
        std::filesystem::create_directories(
            std::filesystem::path(image_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(shear_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(orig_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(output_file_).parent_path());
        FitsIO::writeImage(image_file_, 1, 1, std::vector<float>{0.0f});
        writeNorm(-1.0f);
    }

    // ==========================================
    // Function: Remove the synthetic exposure tree
    // Method: Clean every temporary product without propagating teardown errors.
    // ==========================================
    ~TemporaryCatalogTree() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    // ==========================================
    // Function: Replace the Stage-1 norm sentinel
    // Method: Write one valid 2D norm product at the production-derived path.
    // ==========================================
    void writeNorm(float sentinel) const {
        const std::string filename =
            Universalblock::normFilename(image_file_, root_.string());
        std::filesystem::create_directories(
            std::filesystem::path(filename).parent_path());
        FitsIO::writeImage(filename, 1, 1, std::vector<float>{sentinel});
    }

    // ==========================================
    // Function: Write independently populated shear and original catalogs
    // Method: Always write nonblank headers and optionally one data row per input.
    // ==========================================
    void writeCatalogs(bool with_shear_data, bool with_original_data) const {
        std::ofstream shear(shear_file_, std::ios::trunc);
        std::ofstream orig(orig_file_, std::ios::trunc);
        require(static_cast<bool>(shear) && static_cast<bool>(orig),
                "synthetic Stage-9 inputs must open");

        for (int column = 0; column < LensingConfig::shear_cat_ncols; ++column) {
            shear << (column == 0 ? "s0" : " s" + std::to_string(column));
        }
        shear << '\n';
        orig << "ra dec\n";
        if (with_shear_data) {
            for (int column = 0; column < LensingConfig::shear_cat_ncols; ++column) {
                shear << (column == 0 ? "0" : " 0");
            }
            shear << '\n';
        }
        if (with_original_data) {
            orig << "180 0\n";
        }
    }

    // ==========================================
    // Function: Seed a stale combined catalog
    // Method: Replace the exposure output with recognizable prior-run content.
    // ==========================================
    void writeStaleOutput() const {
        std::ofstream output(output_file_, std::ios::trunc);
        require(static_cast<bool>(output), "stale combined output must open");
        output << "stale\n";
    }

    // ==========================================
    // Function: Invoke Stage 9 for the synthetic exposure
    // Method: Pass the single live chip through the production catalog combiner.
    // ==========================================
    void combine(float chi2) const {
        CatalogCombiner::combineExpoCatalog(
            1, std::vector<std::string>{image_file_}, root_.string(), chi2);
    }

    const std::string& outputFile() const noexcept { return output_file_; }

private:
    std::filesystem::path root_;
    std::string image_file_;
    std::string prefix_;
    std::string shear_file_;
    std::string orig_file_;
    std::string output_file_;
};

// ==========================================
// Function: Verify stale output removal without replacement data
// Method: Exercise invalid, header-only-shear, and high-Chi2 inputs independently.
// ==========================================
void testNoOutputCases(TemporaryCatalogTree& tree) {
    tree.writeNorm(1.0f);
    tree.writeStaleOutput();
    tree.combine(0.0f);
    require(!std::filesystem::exists(tree.outputFile()),
            "all-invalid exposure must remove stale output");

    tree.writeNorm(-1.0f);
    tree.writeCatalogs(false, true);
    tree.writeStaleOutput();
    tree.combine(0.0f);
    require(!std::filesystem::exists(tree.outputFile()),
            "header-only shear with populated original catalog must be skipped");

    tree.writeCatalogs(true, true);
    tree.writeStaleOutput();
    tree.combine(static_cast<float>(LensingConfig::chi2_thresh + 1.0));
    require(!std::filesystem::exists(tree.outputFile()),
            "high-Chi2 exposure must not retain output");
}

// ==========================================
// Function: Verify lazy creation and terminal Chi2 serialization
// Method: Combine one aligned data row and inspect both output lines.
// ==========================================
void testLiveOutput(TemporaryCatalogTree& tree) {
    constexpr float chi2 = 0.005f;
    tree.writeCatalogs(true, true);
    tree.combine(chi2);
    require(std::filesystem::exists(tree.outputFile()),
            "live exposure must create a combined catalog");

    std::ifstream output(tree.outputFile());
    std::string header;
    std::string row;
    std::string extra;
    require(std::getline(output, header) && std::getline(output, row)
                && !std::getline(output, extra),
            "combined catalog must contain one header and one data row");
    require(header.size() >= 4
                && header.substr(header.size() - 4) == "Chi2",
            "combined catalog header must end in Chi2");

    std::istringstream values(row);
    double value = 0.0;
    double last_value = 0.0;
    int columns = 0;
    while (values >> value) {
        last_value = value;
        ++columns;
    }
    require(columns > 0 && std::abs(last_value - chi2) < 1.0e-7,
            "combined catalog data must end in the exposure Chi2");
}

}  // namespace

// ==========================================
// Function: Run Stage-9 output lifecycle regressions
// Method: Cover stale removal, asymmetric sentinel skip, lazy output, and Chi2 output.
// ==========================================
int main() {
    TemporaryCatalogTree tree;
    testNoOutputCases(tree);
    testLiveOutput(tree);
    std::cout << "CatalogCombiner lifecycle tests passed\n";
    return EXIT_SUCCESS;
}
