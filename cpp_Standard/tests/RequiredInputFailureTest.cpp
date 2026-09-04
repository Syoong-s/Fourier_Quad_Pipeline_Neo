#include "RuntimeConfig.hpp"
#include "general/OutputLayout.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/FourierTransformSt1.hpp"
#include "process_main/FourierTransformSt2.hpp"
#include "process_main/SourceExtractor.hpp"
#include "process_main/Universalblock.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

namespace {

enum class Stage { StarFft, GalaxyFft };

// ==========================================
// Function: Stop the test program when a requirement is not met
// Method: Print one focused failure and return a non-zero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Required-input failure test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own a synthetic chip and its sharded Stage-3 products
// Method: Create isolated norm/catalog paths and remove the complete tree on exit.
// ==========================================
class TemporaryInputTree {
public:
    TemporaryInputTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_required_input_" + std::to_string(::getpid()))),
          image_file_((root_ / "input" / "exposure_3.fits").string()),
          norm_template_file_((root_ / "input" / "norm_template.fits").string()),
          norm_file_(Universalblock::normFilename(image_file_, root_.string())),
          star_file_(OutputLayout::chipPath(
              root_.string(), "stamps/dat_StarCanInfo", "exposure_3",
              "_star_can_info.dat")),
          source_file_(OutputLayout::chipPath(
              root_.string(), "stamps/dat_SrcInfo", "exposure_3",
              "_source_info.dat")) {
        std::filesystem::create_directories(
            std::filesystem::path(image_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(norm_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(star_file_).parent_path());
        std::filesystem::create_directories(
            std::filesystem::path(source_file_).parent_path());
    }

    ~TemporaryInputTree() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    // ==========================================
    // Function: Replace the synthetic norm with one requested sentinel
    // Method: Write a valid 2D FITS image whose first pixel controls NormStatus.
    // ==========================================
    void writeNorm(float sentinel) const {
        require(FitsIO::writeImage(
                    norm_file_, 2, 2, {sentinel, -1.0f, -1.0f, -1.0f}),
                "cannot write synthetic norm FITS");
    }

    // ==========================================
    // Function: Replace the synthetic science image with one valid FITS image
    // Method: Write finite pixels at the requested dimensions for Stage-3 input tests.
    // ==========================================
    void writeScience(int nx, int ny) const {
        require(FitsIO::writeImage(
                    image_file_, nx, ny,
                    std::vector<float>(static_cast<std::size_t>(nx) * ny, 0.0f)),
                "cannot write synthetic science FITS");
    }

    // ==========================================
    // Function: Remove the synthetic science image
    // Method: Preserve the derived product prefix while forcing the science read to fail.
    // ==========================================
    void removeScience() const {
        std::error_code error;
        std::filesystem::remove(image_file_, error);
    }

    // ==========================================
    // Function: Write a complete valid-metadata norm at requested dimensions
    // Method: Build exact background/sigma coefficient arrays and copy a matching template HDU.
    // ==========================================
    void writeCompleteNorm(int nx, int ny) const {
        require(FitsIO::writeImage(
                    norm_template_file_, nx, ny,
                    std::vector<float>(static_cast<std::size_t>(nx) * ny, -1.0f)),
                "cannot write norm template FITS");
        const int split = RuntimeConfigStore::get().lensing.ccd_split;
        std::vector<double> background(
            static_cast<std::size_t>(split) * LensingConfig::nct, 0.0);
        std::vector<double> sigma(static_cast<std::size_t>(split) * 3U, 0.0);
        require(FitsIO::writeNormHDU(
                    norm_template_file_, norm_file_, nx, ny,
                    std::vector<float>(static_cast<std::size_t>(nx) * ny, -1.0f),
                    background, sigma, split, LensingConfig::nct),
                "cannot write complete synthetic norm FITS");
    }

    // ==========================================
    // Function: Replace one required Stage-3 catalog
    // Method: Remove it for Missing, truncate it for empty, or write one header line.
    // ==========================================
    void prepareCatalog(Stage stage, int state) const {
        const std::string& filename = catalog(stage);
        std::error_code error;
        std::filesystem::remove(filename, error);
        if (state < 0) {
            return;
        }
        std::ofstream output(filename, std::ios::trunc);
        require(output.is_open(), "cannot create synthetic catalog");
        if (state > 0) {
            output << "synthetic header\n";
        }
    }

    // ==========================================
    // Function: Return the synthetic science-image path
    // Method: Supply only the stable chip prefix needed by the tested stages.
    // ==========================================
    const std::string& imageFile() const noexcept { return image_file_; }

    // ==========================================
    // Function: Return the synthetic process-main output root
    // Method: Expose the root shared by all generated test products.
    // ==========================================
    std::string outputRoot() const { return root_.string(); }

private:
    // ==========================================
    // Function: Select the required catalog for one FFT stage
    // Method: Map Stage 4 to star candidates and Stage 6 to source information.
    // ==========================================
    const std::string& catalog(Stage stage) const {
        return stage == Stage::StarFft ? star_file_ : source_file_;
    }

    std::filesystem::path root_;
    std::string image_file_;
    std::string norm_template_file_;
    std::string norm_file_;
    std::string star_file_;
    std::string source_file_;
};

// ==========================================
// Function: Invoke one chip-level FFT stage
// Method: Dispatch directly to the production Stage-4 or Stage-6 entry point.
// ==========================================
void invokeStage(Stage stage, const TemporaryInputTree& tree) {
    if (stage == Stage::StarFft) {
        FourierTransformSt1::chipProcessFourierTSt1(
            tree.imageFile(), tree.outputRoot());
    } else {
        FourierTransformSt2::chipProcessFourierTSt2(
            tree.imageFile(), tree.outputRoot());
    }
}

// ==========================================
// Function: Run one production stage in an isolated child process
// Method: Observe fail-fast exit status without terminating the test controller.
// ==========================================
bool stageFails(Stage stage, const TemporaryInputTree& tree) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        invokeStage(stage, tree);
        std::_Exit(EXIT_SUCCESS);
    }
    int status = 0;
    require(::waitpid(child, &status, 0) == child, "waitpid failed");
    return !WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS;
}

// ==========================================
// Function: Run the production Stage-3 chip entry in an isolated child
// Method: Observe fail-fast status without allowing an expected abort to end the suite.
// ==========================================
bool sourceStageFails(const TemporaryInputTree& tree) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        SourceExtractor::chipProcessSource(
            {tree.imageFile()}, 1, tree.outputRoot());
        std::_Exit(EXIT_SUCCESS);
    }
    int status = 0;
    require(::waitpid(child, &status, 0) == child, "waitpid failed");
    return !WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS;
}

// ==========================================
// Function: Verify Stage-3 science and full-norm fail-fast boundaries
// Method: Preserve Invalid skip, then exercise missing science, metadata, and geometry failures.
// ==========================================
void testSourceInputStateMachine(TemporaryInputTree& tree) {
    tree.removeScience();
    tree.writeNorm(1.0f);
    require(!sourceStageFails(tree), "Invalid norm must skip missing science");

    tree.writeNorm(-1.0f);
    require(sourceStageFails(tree), "valid norm plus missing science must fail");

    tree.writeScience(4, 4);
    tree.writeNorm(-1.0f);
    require(sourceStageFails(tree), "valid norm without metadata must fail full read");

    tree.writeCompleteNorm(3, 4);
    require(sourceStageFails(tree), "science/norm geometry mismatch must fail");
}

// ==========================================
// Function: Verify one FFT stage's required-input state machine
// Method: Preserve Invalid skip, fail on missing/empty input, and accept a header-only catalog.
// ==========================================
void testStage(Stage stage, TemporaryInputTree& tree) {
    tree.writeNorm(1.0f);
    tree.prepareCatalog(stage, -1);
    require(!stageFails(stage, tree), "Invalid norm must skip missing catalog");

    tree.writeNorm(-1.0f);
    tree.prepareCatalog(stage, -1);
    require(stageFails(stage, tree), "valid norm plus missing catalog must fail");

    tree.prepareCatalog(stage, 0);
    require(stageFails(stage, tree), "zero-byte catalog must fail");

    tree.prepareCatalog(stage, 1);
    require(!stageFails(stage, tree), "header-only catalog must remain valid");
}

}  // namespace

// ==========================================
// Function: Run Stage-4 and Stage-6 required-input failure regressions
// Method: Initialize compiled defaults once and exercise both chip entry points.
// ==========================================
int main() {
    std::string error;
    require(RuntimeConfigStore::initialize(makeDefaultRuntimeConfig(), error),
            "cannot initialize runtime config: " + error);
    TemporaryInputTree tree;
    testSourceInputStateMachine(tree);
    testStage(Stage::StarFft, tree);
    testStage(Stage::GalaxyFft, tree);
    std::cout << "Required-input failure tests passed\n";
    return EXIT_SUCCESS;
}
