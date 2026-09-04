#include "process_main/FitsIO.hpp"
#include "process_main/Universalblock.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

namespace {

// ==========================================
// Class: Own one temporary sharded norm-product tree
// Method: Create an isolated process-main layout and remove it after every test run.
// ==========================================
class TemporaryNormTree {
public:
    TemporaryNormTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_universalblock_" + std::to_string(::getpid()))),
          imageFile_((root_ / "input" / "exposure_1.fits").string()),
          normFile_(Universalblock::normFilename(imageFile_, root_.string())) {
        std::filesystem::create_directories(
            std::filesystem::path(normFile_).parent_path());
    }

    ~TemporaryNormTree() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    // ==========================================
    // Function: Return the synthetic chip image path
    // Method: Expose a stable prefix source without requiring a real raw image.
    // ==========================================
    const std::string& imageFile() const noexcept {
        return imageFile_;
    }

    // ==========================================
    // Function: Return the synthetic process-main output root
    // Method: Expose the owned temporary root for Universalblock checks.
    // ==========================================
    const std::string outputRoot() const {
        return root_.string();
    }

    // ==========================================
    // Function: Return the derived norm-product path
    // Method: Expose the module-generated sharded pathname to test setup.
    // ==========================================
    const std::string& normFile() const noexcept {
        return normFile_;
    }

    // ==========================================
    // Function: Replace the norm product with one valid 2D FITS image
    // Method: Put the requested sentinel in the first pixel and benign values elsewhere.
    // ==========================================
    void writeNorm(float sentinel) const {
        const std::vector<float> pixels = {sentinel, 2.0f, 3.0f, 4.0f};
        FitsIO::writeImage(normFile_, 2, 2, pixels);
    }

private:
    std::filesystem::path root_;
    std::string imageFile_;
    std::string normFile_;
};

// ==========================================
// Function: Stop the test program when a requirement is not met
// Method: Print a focused failure and return a non-zero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Universalblock test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Verify valid and invalid norm sentinel classification
// Method: Rewrite only the first FITS pixel and require invalid checks to remain silent.
// ==========================================
void testSentinelClassification(TemporaryNormTree& tree) {
    tree.writeNorm(-1.0f);
    float sentinel = 0.0f;
    require(Universalblock::checkNorm(
                tree.imageFile(), tree.outputRoot(), &sentinel)
                == Universalblock::NormStatus::Valid,
            "-1 sentinel must be valid");
    require(sentinel == -1.0f, "valid sentinel output mismatch");

    const std::vector<float> invalidSentinels = {
        1.0f,
        -99999.0f,
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity()
    };
    for (float value : invalidSentinels) {
        tree.writeNorm(value);
        std::ostringstream capturedError;
        std::streambuf* originalError = std::cerr.rdbuf(capturedError.rdbuf());
        const Universalblock::NormStatus status =
            Universalblock::checkNorm(tree.imageFile(), tree.outputRoot());
        std::cerr.rdbuf(originalError);
        require(status == Universalblock::NormStatus::Invalid,
                "failed/non-finite sentinel must be invalid");
        require(capturedError.str().empty(),
                "invalid sentinel classification must be silent");
    }
}

// ==========================================
// Function: Verify missing and malformed norm input classification
// Method: Remove the file for Missing, then replace it with non-FITS bytes for ReadError.
// ==========================================
void testInputFailures(TemporaryNormTree& tree) {
    std::filesystem::remove(tree.normFile());
    const Universalblock::NormStatus missing_status =
        Universalblock::checkNorm(tree.imageFile(), tree.outputRoot());
    require(missing_status == Universalblock::NormStatus::Missing,
            "absent norm path must be Missing");
    require(Universalblock::normErrorDetail(
                missing_status, tree.imageFile(), tree.outputRoot())
                == "missing norm FITS: " + tree.normFile(),
            "Missing diagnostic must identify the derived norm path");

    std::ofstream malformed(tree.normFile(), std::ios::binary | std::ios::trunc);
    malformed << "not a FITS image";
    malformed.close();
    const Universalblock::NormStatus read_error_status =
        Universalblock::checkNorm(tree.imageFile(), tree.outputRoot());
    require(read_error_status == Universalblock::NormStatus::ReadError,
            "malformed existing norm file must be ReadError");
    require(Universalblock::normErrorDetail(
                read_error_status, tree.imageFile(), tree.outputRoot())
                == "unreadable norm FITS: " + tree.normFile(),
            "ReadError diagnostic must identify the derived norm path");
}

// ==========================================
// Function: Verify compact FITS failure state and a true 1x1 Invalid norm
// Method: Read one missing image, then serialize and classify the Stage-1
//         unreadable-input product without a source HDU template.
// ==========================================
void testCompactReadFailureAndInvalidNorm(TemporaryNormTree& tree) {
    int nx = 7;
    int ny = 9;
    std::vector<float> pixels(8, 0.0f);
    require(!FitsIO::readImage(tree.imageFile(), nx, ny, pixels),
            "missing FITS input must fail");
    require(pixels.size() == 1 && pixels[0] == -99999.0f,
            "failed FITS read must return one sentinel pixel");

    require(FitsIO::writeImage(tree.normFile(), 1, 1, {1.0f}),
            "1x1 Invalid norm write failed");
    nx = 0;
    ny = 0;
    pixels.clear();
    require(FitsIO::readImage(tree.normFile(), nx, ny, pixels)
                && nx == 1 && ny == 1
                && pixels.size() == 1 && pixels[0] == 1.0f,
            "1x1 Invalid norm did not round-trip");
    require(Universalblock::checkNorm(tree.imageFile(), tree.outputRoot())
                == Universalblock::NormStatus::Invalid,
            "1x1 positive sentinel must classify as Invalid");
}

}  // namespace

// ==========================================
// Function: Run the focused Universalblock norm-validity regression suite
// Method: Exercise path derivation, valid/invalid sentinels, silence, missing input, and corruption.
// ==========================================
int main() {
    TemporaryNormTree tree;
    testCompactReadFailureAndInvalidNorm(tree);
    testSentinelClassification(tree);
    testInputFailures(tree);
    std::cout << "Universalblock tests passed\n";
    return EXIT_SUCCESS;
}
