#include "process_init/FitsExtractor.hpp"

#include <fitsio.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

enum class KeyKind {
    Missing,
    Integer,
    Malformed,
};

// ==========================================
// Function: Stop the initializer regression on a failed requirement
// Method: Print one focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Initializer CCDNUM test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated initializer test tree
// Method: Create a unique native-Linux temporary root and remove it on scope exit.
// ==========================================
class TemporaryTree {
public:
    explicit TemporaryTree(const std::string& label)
        : root_(std::filesystem::temp_directory_path()
                / ("fq_init_ccdnum_" + label + "_"
                   + std::to_string(
                       std::chrono::steady_clock::now().time_since_epoch().count()))) {
        std::filesystem::create_directories(root_);
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
// Function: Close a synthetic FITS handle with checked status
// Method: Use a fresh close status and fail the regression on any CFITSIO error.
// ==========================================
void closeFits(fitsfile*& file, int status, const std::string& context) {
    int close_status = 0;
    if (file != nullptr) fits_close_file(file, &close_status);
    file = nullptr;
    require(status == 0 && close_status == 0, context);
}

// ==========================================
// Function: Append one two-dimensional synthetic image HDU
// Method: Write either an integer CCDNUM, no key, or an intentionally malformed key.
// ==========================================
void appendImage(fitsfile* file, KeyKind key_kind, int ccdnum) {
    int status = 0;
    long axes[2] = {2, 2};
    fits_create_img(file, FLOAT_IMG, 2, axes, &status);
    if (key_kind == KeyKind::Integer) {
        fits_update_key(file, TINT, "CCDNUM", &ccdnum, nullptr, &status);
    } else if (key_kind == KeyKind::Malformed) {
        char value[] = "not-an-integer";
        fits_update_key(file, TSTRING, "CCDNUM", value, nullptr, &status);
    }
    float pixels[4] = {
        static_cast<float>(ccdnum), 1.0f, 2.0f, 3.0f};
    fits_write_img(file, TFLOAT, 1, 4, pixels, &status);
    require(status == 0, "cannot append synthetic image HDU");
}

// ==========================================
// Function: Create one multi-HDU initializer input archive
// Method: Write a primary HDU followed by caller-specified two-dimensional extensions.
// ==========================================
void createArchive(
    const std::filesystem::path& path,
    const std::vector<std::pair<KeyKind, int>>& extensions) {
    fitsfile* file = nullptr;
    int status = 0;
    const std::string create_name = "!" + path.string();
    fits_create_file(&file, create_name.c_str(), &status);
    fits_create_img(file, BYTE_IMG, 0, nullptr, &status);
    require(status == 0, "cannot create synthetic archive primary HDU");
    for (const auto& extension : extensions) {
        appendImage(file, extension.first, extension.second);
    }
    closeFits(file, status, "cannot close synthetic archive");
}

// ==========================================
// Function: Read the CCDNUM stored in one extracted chip FITS
// Method: Open the primary image and require an integer identity keyword.
// ==========================================
int readCCDNUM(const std::filesystem::path& path) {
    fitsfile* file = nullptr;
    int status = 0;
    fits_open_file(&file, path.string().c_str(), READONLY, &status);
    int ccdnum = 0;
    if (status == 0) {
        fits_read_key(file, TINT, "CCDNUM", &ccdnum, nullptr, &status);
    }
    closeFits(file, status, "cannot read extracted CCDNUM");
    return ccdnum;
}

// ==========================================
// Function: Replace one extracted chip's CCDNUM
// Method: Corrupt only the identity keyword to exercise resume validation and re-extraction.
// ==========================================
void replaceCCDNUM(const std::filesystem::path& path, int ccdnum) {
    fitsfile* file = nullptr;
    int status = 0;
    fits_open_file(&file, path.string().c_str(), READWRITE, &status);
    if (status == 0) {
        fits_update_key(file, TINT, "CCDNUM", &ccdnum, nullptr, &status);
    }
    closeFits(file, status, "cannot corrupt extracted CCDNUM");
}

// ==========================================
// Function: Verify CCDNUM-native Science and DQ extraction
// Method: Use shuffled, gapped identities plus a missing-key HDU and inspect names and headers.
// ==========================================
void testCanonicalExtractionAndResume() {
    TemporaryTree tree("canonical");
    const std::filesystem::path science = tree.root() / "exposure_ood.fits.fz";
    createArchive(science, {
        {KeyKind::Integer, 7},
        {KeyKind::Missing, 0},
        {KeyKind::Integer, 1},
        {KeyKind::Integer, 3},
    });

    const std::filesystem::path final_dir = tree.root() / "science";
    const std::filesystem::path staging_dir = tree.root() / "stage_science";
    fqinit::ExtractionResult science_result = fqinit::extractArchive(
        science, fqinit::ProductKind::Science, final_dir, staging_dir,
        fqinit::ExistingPolicy::Overwrite);
    require(science_result.success && science_result.output_paths.size() == 3,
            "Science extraction must retain all keyed HDUs");

    const std::vector<int> expected_order = {7, 1, 3};
    for (std::size_t index = 0; index < expected_order.size(); ++index) {
        const int ccdnum = expected_order[index];
        require(science_result.output_paths[index].filename()
                    == "exposure_ood_" + std::to_string(ccdnum) + ".fits",
                "Science filename must use physical CCDNUM in archive order");
        require(readCCDNUM(science_result.output_paths[index]) == ccdnum,
                "Science output header must preserve physical CCDNUM");
    }

    const std::filesystem::path dq_final = tree.root() / "dq";
    fqinit::ExtractionResult dq_result = fqinit::extractArchive(
        science, fqinit::ProductKind::DqMask, dq_final,
        tree.root() / "stage_dq", fqinit::ExistingPolicy::Overwrite);
    require(dq_result.success && dq_result.output_paths.size() == 3,
            "DQ extraction must retain all keyed HDUs");
    const std::string dq_stem = fqinit::dqOutputStem(science);
    require(dq_result.output_paths.front().filename()
                == dq_stem + "_7.fits",
            "DQ filename must use mapped exposure stem and physical CCDNUM");

    replaceCCDNUM(final_dir / "exposure_ood_7.fits", 99);
    fqinit::ExtractionResult resumed = fqinit::extractArchive(
        science, fqinit::ProductKind::Science, final_dir,
        tree.root() / "stage_resume", fqinit::ExistingPolicy::Resume);
    require(resumed.success && resumed.output_paths.size() == 3,
            "partial resume must publish a complete output set");
    require(readCCDNUM(final_dir / "exposure_ood_7.fits") == 7,
            "resume must re-extract an output whose header identity is wrong");
}

// ==========================================
// Function: Verify malformed persistent identities fail before publication
// Method: Exercise duplicate, non-positive, and unconvertible CCDNUM values independently.
// ==========================================
void testMalformedIdentities() {
    const auto extractionFails = [](
        const std::string& label,
        const std::vector<std::pair<KeyKind, int>>& extensions) {
        TemporaryTree tree(label);
        const std::filesystem::path archive = tree.root() / "bad.fits.fz";
        createArchive(archive, extensions);
        const fqinit::ExtractionResult result = fqinit::extractArchive(
            archive, fqinit::ProductKind::Science, tree.root() / "out",
            tree.root() / "stage", fqinit::ExistingPolicy::Overwrite);
        return !result.success && !result.error.empty();
    };

    require(extractionFails("duplicate", {
                {KeyKind::Integer, 3}, {KeyKind::Integer, 3}}),
            "duplicate CCDNUM must fail the archive");
    require(extractionFails("zero", {{KeyKind::Integer, 0}}),
            "non-positive CCDNUM must fail the archive");
    require(extractionFails("malformed", {{KeyKind::Malformed, 0}}),
            "malformed CCDNUM must fail the archive");
}

}  // namespace

// ==========================================
// Function: Run focused process_init CCDNUM regressions
// Method: Exercise canonical extraction, DQ naming, resume repair, and fatal identities.
// ==========================================
int main() {
    testCanonicalExtractionAndResume();
    testMalformedIdentities();
    std::cout << "Initializer CCDNUM tests passed\n";
    return EXIT_SUCCESS;
}
