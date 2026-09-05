#include "general/OutputLayout.hpp"
#include "process_main/Astrometry.hpp"
#include "process_main/UniversalUtils.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

namespace {

// ==========================================
// Function: Stop the identity regression on a failed requirement
// Method: Print one focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "CCDNUM identity test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own one isolated cross-stage identity test tree
// Method: Create a unique temporary root and remove it on scope exit.
// ==========================================
class TemporaryTree {
public:
    TemporaryTree()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_ccdnum_identity_"
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
// Function: Observe one malformed canonical filename in a child process
// Method: Let the production fail-fast parser exit only the forked child.
// ==========================================
bool chipIdFails(const std::string& path) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        (void)UniversalUtils::getChipId(path);
        std::_Exit(EXIT_SUCCESS);
    }
    int status = 0;
    require(::waitpid(child, &status, 0) == child, "waitpid failed");
    return !WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS;
}

// ==========================================
// Function: Verify filename-only physical identity parsing
// Method: Accept positive decimal suffixes and reject missing, malformed, or zero suffixes.
// ==========================================
void testFilenameIdentity() {
    require(UniversalUtils::getChipId("/tmp/exposure_name_37.fits") == 37,
            "positive CCDNUM suffix was not parsed");
    require(UniversalUtils::getChipId("/tmp/exposure_name_001.fits") == 1,
            "leading-zero CCDNUM suffix was not parsed");
    require(chipIdFails("/tmp/exposure.fits"),
            "missing CCDNUM suffix must fail");
    require(chipIdFails("/tmp/exposure_bad.fits"),
            "non-decimal CCDNUM suffix must fail");
    require(chipIdFails("/tmp/exposure_0.fits"),
            "zero CCDNUM suffix must fail");
}

// ==========================================
// Function: Write one shuffled keyed astrometry header
// Method: Preserve the production shared projection block and emit non-positional chip rows.
// ==========================================
void writeKeyedHead(const std::filesystem::path& path) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create keyed astrometry header");
    output << "10 20\n";
    for (int term = 0; term < LensingConfig::npd; ++term) output << "0 0\n";
    output << "7 1 70 71 7 0 0 7\n";
    output << "1 1 10 11 1 0 0 1\n";
    output << "3 0 30 31 3 0 0 3\n";
}

// ==========================================
// Function: Verify `.head` consumers scan by physical identity
// Method: Query shuffled valid, invalid, and absent CCDNUM rows independently.
// ==========================================
void testKeyedHeadLookup(const std::filesystem::path& root) {
    const std::filesystem::path head = root / "keyed.head";
    writeKeyedHead(head);
    double crpix[2] = {0.0, 0.0};
    double cd[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double crval[2] = {0.0, 0.0};
    double pu[2][LensingConfig::npd] = {{0.0}, {0.0}};

    int error = 0;
    Astrometry::readAstrometryPara(
        head.string(), 1, crpix, cd, crval, pu, LensingConfig::npd, error);
    require(error == 0 && crpix[0] == 10.0 && cd[0][0] == 1.0,
            "CCDNUM 1 resolved a positional row");

    error = 0;
    Astrometry::readAstrometryPara(
        head.string(), 7, crpix, cd, crval, pu, LensingConfig::npd, error);
    require(error == 0 && crpix[0] == 70.0 && cd[0][0] == 7.0,
            "CCDNUM 7 did not resolve its keyed row");

    error = 0;
    Astrometry::readAstrometryPara(
        head.string(), 3, crpix, cd, crval, pu, LensingConfig::npd, error);
    require(error == 1, "invalid keyed row must remain invalid");
    error = 0;
    Astrometry::readAstrometryPara(
        head.string(), 9, crpix, cd, crval, pu, LensingConfig::npd, error);
    require(error == 1, "absent keyed row must fail lookup");
}

// ==========================================
// Function: Write one chip-level astrometry match file
// Method: Emit ten finite rows so the chip is initially valid before the exposure fit gate.
// ==========================================
void writeAstroData(const std::filesystem::path& root,
                    const std::string& prefix, double crval) {
    const std::filesystem::path path = OutputLayout::chipPath(
        root.string(), "astrometry/dat_Astro", prefix, "_astro.dat");
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create chip astrometry data");
    output << "100 200 " << crval << " 0\n";
    output << "0.001 0 0 0.001\n";
    output << "10 10 10\n";
    for (int row = 0; row < 10; ++row) {
        output << crval + 1.0e-5 * row << " " << 1.0e-5 * row
               << " " << 100 + row << " " << 200 + row << "\n";
    }
}

// ==========================================
// Function: Read the shared CRVAL and keyed row identities from one generated header
// Method: Skip the projection block and return the first column of every chip row.
// ==========================================
std::pair<double, std::vector<int>> readGeneratedHead(
    const std::filesystem::path& path, int row_count) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "generated astrometry header is missing");
    double crval = 0.0;
    double dec = 0.0;
    require(static_cast<bool>(input >> crval >> dec), "generated CRVAL is missing");
    double discarded = 0.0;
    for (int term = 0; term < 2 * LensingConfig::npd; ++term) {
        require(static_cast<bool>(input >> discarded), "generated PU block is incomplete");
    }
    std::vector<int> ccdnums;
    for (int row = 0; row < row_count; ++row) {
        int ccdnum = 0;
        int valid = 0;
        require(static_cast<bool>(input >> ccdnum >> valid), "generated keyed row is missing");
        ccdnums.push_back(ccdnum);
        for (int field = 0; field < 6; ++field) {
            require(static_cast<bool>(input >> discarded), "generated keyed row is incomplete");
        }
    }
    return {crval, ccdnums};
}

// ==========================================
// Function: Verify normal Astrometry is invariant to exposure-list order
// Method: Use gapped CCDNUM values with distinct CRVALs and require CCDNUM 1 as reference.
// ==========================================
void testNormalAstrometryOrder(const std::filesystem::path& root) {
    const std::vector<std::string> first_order = {
        (root / "science" / "expo" / "exposure_7.fits").string(),
        (root / "science" / "expo" / "exposure_1.fits").string(),
        (root / "science" / "expo" / "exposure_3.fits").string(),
    };
    writeAstroData(root, "exposure_7", 107.0);
    writeAstroData(root, "exposure_1", 101.0);
    writeAstroData(root, "exposure_3", 103.0);
    std::filesystem::create_directories(root / "astrometry" / "Head");
    Astrometry::getAstrometry(first_order, 3, root.string());
    const auto first = readGeneratedHead(
        root / "astrometry" / "Head" / "exposure.head", 3);
    require(std::fabs(first.first - 101.0) < 1.0e-12,
            "normal Astrometry did not choose the smallest valid CCDNUM CRVAL");
    require(first.second == std::vector<int>({7, 1, 3}),
            "normal Astrometry wrote dense row numbers instead of CCDNUM");

    const std::vector<std::string> second_order = {
        first_order[2], first_order[0], first_order[1]};
    Astrometry::getAstrometry(second_order, 3, root.string());
    const auto second = readGeneratedHead(
        root / "astrometry" / "Head" / "exposure.head", 3);
    require(std::fabs(second.first - 101.0) < 1.0e-12,
            "reordered list changed the shared CRVAL");
    require(second.second == std::vector<int>({3, 7, 1}),
            "reordered output rows lost their physical CCDNUM identity");
}

#ifdef FQ_STANDARD_VARIANT
// ==========================================
// Function: Verify trivial Astrometry excludes unreadable rows from identity selection
// Method: Omit one chip file, shuffle the list, and require the smallest readable CCDNUM CRVAL.
// ==========================================
void testTrivialAstrometryValidity(const std::filesystem::path& root) {
    const std::vector<std::string> images = {
        (root / "science" / "trivial" / "trivial_7.fits").string(),
        (root / "science" / "trivial" / "trivial_3.fits").string(),
        (root / "science" / "trivial" / "trivial_1.fits").string(),
    };
    writeAstroData(root, "trivial_7", 207.0);
    writeAstroData(root, "trivial_1", 201.0);
    Astrometry::getAstrometryTrivial(images, 3, root.string());
    const auto generated = readGeneratedHead(
        root / "astrometry" / "Head" / "trivial.head", 3);
    require(std::fabs(generated.first - 201.0) < 1.0e-12,
            "trivial Astrometry selected an unreadable or non-minimum CCDNUM");
    require(generated.second == std::vector<int>({7, 3, 1}),
            "trivial Astrometry did not write keyed CCDNUM rows");
}
#endif

}  // namespace

// ==========================================
// Function: Run cross-stage CCDNUM identity regressions
// Method: Exercise filename parsing, keyed lookup, deterministic writing, and list reordering.
// ==========================================
int main() {
    TemporaryTree tree;
    testFilenameIdentity();
    testKeyedHeadLookup(tree.root());
    testNormalAstrometryOrder(tree.root());
#ifdef FQ_STANDARD_VARIANT
    testTrivialAstrometryValidity(tree.root());
#endif
    std::cout << "CCDNUM identity tests passed\n";
    return EXIT_SUCCESS;
}
