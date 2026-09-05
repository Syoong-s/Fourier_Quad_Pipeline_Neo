#include "general/MPIScheduler.hpp"
#include "process_init/Initializer.hpp"

#include <fitsio.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace {

// ==========================================
// Function: Stop the publication regression on a failed requirement
// Method: Print one focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Initializer publication CCDNUM test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated clean-output initializer tree
// Method: Create unique Science, DQ, and output roots and remove them on scope exit.
// ==========================================
class TemporaryTree {
public:
    TemporaryTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_init_publish_ccdnum_"
                   + std::to_string(
                       std::chrono::steady_clock::now().time_since_epoch().count()))) {
        std::filesystem::create_directories(scienceRoot());
        std::filesystem::create_directories(dqRoot());
        std::filesystem::create_directories(outputRoot());
    }

    ~TemporaryTree() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    std::filesystem::path scienceRoot() const { return root_ / "science_input"; }
    std::filesystem::path dqRoot() const { return root_ / "dq_input"; }
    std::filesystem::path outputRoot() const { return root_ / "output"; }

private:
    std::filesystem::path root_;
};

// ==========================================
// Function: Create one synthetic multi-HDU archive with shuffled physical chips
// Method: Write CCDNUM 7, 1, and 3 in archive order to exercise unsorted publication.
// ==========================================
void createArchive(const std::filesystem::path& path) {
    fitsfile* file = nullptr;
    int status = 0;
    const std::string create_name = "!" + path.string();
    fits_create_file(&file, create_name.c_str(), &status);
    fits_create_img(file, BYTE_IMG, 0, nullptr, &status);
    for (int ccdnum : {7, 1, 3}) {
        long axes[2] = {2, 2};
        fits_create_img(file, FLOAT_IMG, 2, axes, &status);
        fits_update_key(file, TINT, "CCDNUM", &ccdnum, nullptr, &status);
        float pixels[4] = {
            static_cast<float>(ccdnum), 1.0f, 2.0f, 3.0f};
        fits_write_img(file, TFLOAT, 1, 4, pixels, &status);
    }
    int close_status = 0;
    if (file != nullptr) fits_close_file(file, &close_status);
    require(status == 0 && close_status == 0,
            "cannot create synthetic initializer archive");
}

// ==========================================
// Function: Read nonempty text lines from one published initializer artifact
// Method: Preserve row order while discarding only empty lines.
// ==========================================
std::vector<std::string> readLines(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "published text artifact is missing");
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}

// ==========================================
// Function: Verify manifest schema and unsorted per-exposure publication
// Method: Run the production initializer on a clean root and inspect all public identity files.
// ==========================================
void testPublication(const TemporaryTree& tree) {
    const std::string exposure = "demo_v1_ooi";
    createArchive(tree.scienceRoot() / (exposure + ".fits.fz"));
    createArchive(tree.dqRoot() / "demo_v1_ood.fits.fz");

    fqinit::Config config;
    config.science_root = tree.scienceRoot();
    config.dq_root = tree.dqRoot();
    config.output_root = tree.outputRoot();
    config.target = "target";
    config.filename_prefix = "demo";
    config.filename_tokens = {"v1"};
    config.existing_policy = fqinit::ExistingPolicy::Fail;
    fqinit::normalizeAndValidateConfig(config);
    require(fqinit::runInitializer(config) == 0,
            "production initializer returned failure");

    const std::filesystem::path target_root = tree.outputRoot() / "target";
    const std::vector<std::string> rows = readLines(
        target_root / "stamps" / (exposure + ".list"));
    require(rows.size() == 3, "per-exposure list has the wrong chip count");
    const std::vector<int> expected = {7, 1, 3};
    for (std::size_t index = 0; index < expected.size(); ++index) {
        require(std::filesystem::path(rows[index]).filename()
                    == exposure + "_" + std::to_string(expected[index]) + ".fits",
                "per-exposure publication imposed positional/lexical chip order");
    }
    require(!std::filesystem::exists(
                target_root / "science" / exposure / (exposure + "_2.fits")),
            "missing CCDNUM 2 was synthesized from list position");

    const std::vector<std::string> top_rows = readLines(
        tree.outputRoot() / "expo_target.list");
    require(top_rows.size() == 1
                && top_rows.front().find("     3") != std::string::npos,
            "top-level exposure list lost the physical chip count");

    std::ifstream manifest_input(
        tree.outputRoot() / "init_target_manifest.json");
    require(static_cast<bool>(manifest_input), "initializer manifest is missing");
    std::ostringstream manifest;
    manifest << manifest_input.rdbuf();
    require(manifest.str().find("\"schema_version\": 3") != std::string::npos,
            "initializer manifest did not publish schema version 3");
    require(manifest.str().find("\"science_numbering\": \"CCDNUM\"")
                != std::string::npos,
            "Science numbering metadata is not CCDNUM");
    require(manifest.str().find("\"dq_numbering\": \"CCDNUM\"")
                != std::string::npos,
            "DQ numbering metadata is not CCDNUM");
}

}  // namespace

// ==========================================
// Function: Run clean-output initializer publication regressions
// Method: Initialize one-rank MPI, exercise the production path, then finalize MPI.
// ==========================================
int main(int argc, char** argv) {
    MPIScheduler::init(argc, argv);
    TemporaryTree tree;
    testPublication(tree);
    MPIScheduler::finalize();
    std::cout << "Initializer publication CCDNUM tests passed\n";
    return EXIT_SUCCESS;
}
