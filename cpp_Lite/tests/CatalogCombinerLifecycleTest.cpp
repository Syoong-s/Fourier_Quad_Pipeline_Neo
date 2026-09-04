#include "process_main/CatalogCombiner.hpp"
#include "process_main/ExposureInfo.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/FitsIO.hpp"
#include "LensingConfig.hpp"
#include "general/OutputLayout.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"

#include <cmath>
#include <cerrno>
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

ProcessMain::State ProcessMain::state;

namespace ExposureInfo {
State state;
}

namespace {

constexpr int kExposureIndex = 7;

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
// Function: Require one production Stage-9 failure and return its diagnostic text
// Method: Redirect child stderr through a pipe, then accept only nonzero exit or signal.
// ==========================================
template <typename Operation>
std::string requireAbort(Operation operation, const std::string& message) {
    int diagnostic_pipe[2] = {-1, -1};
    require(::pipe(diagnostic_pipe) == 0,
            "pipe failed for fatal Stage-9 case");

    const pid_t child = ::fork();
    require(child >= 0, "fork failed for fatal Stage-9 case");
    if (child == 0) {
        ::close(diagnostic_pipe[0]);
        if (::dup2(diagnostic_pipe[1], STDERR_FILENO) < 0) {
            std::_Exit(EXIT_FAILURE);
        }
        ::close(diagnostic_pipe[1]);
        operation();
        std::_Exit(EXIT_SUCCESS);
    }

    ::close(diagnostic_pipe[1]);
    std::string diagnostics;
    char buffer[4096];
    while (true) {
        const ssize_t count = ::read(diagnostic_pipe[0], buffer, sizeof(buffer));
        if (count > 0) {
            diagnostics.append(buffer, static_cast<std::size_t>(count));
            continue;
        }
        if (count < 0 && errno == EINTR) {
            continue;
        }
        require(count == 0, "read failed for fatal Stage-9 diagnostic");
        break;
    }
    ::close(diagnostic_pipe[0]);

    int status = 0;
    require(::waitpid(child, &status, 0) == child,
            "waitpid failed for fatal Stage-9 case");
    require((WIFEXITED(status) && WEXITSTATUS(status) != EXIT_SUCCESS)
                || WIFSIGNALED(status),
            message);
    return diagnostics;
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
    // Function: Remove the Stage-1 norm product
    // Method: Delete the production-derived path to synthesize a Missing state.
    // ==========================================
    void removeNorm() const {
        std::filesystem::remove(
            Universalblock::normFilename(image_file_, root_.string()));
    }

    // ==========================================
    // Function: Replace the Stage-1 norm product with malformed bytes
    // Method: Truncate the production-derived path to synthesize a ReadError state.
    // ==========================================
    void writeMalformedNorm() const {
        std::ofstream malformed(
            Universalblock::normFilename(image_file_, root_.string()),
            std::ios::binary | std::ios::trunc);
        require(static_cast<bool>(malformed),
                "malformed norm fixture must open");
        malformed << "not a FITS image";
    }

    // ==========================================
    // Function: Format one complete Stage-7 shear data row
    // Method: Fill the live schema width while allowing the first value to vary.
    // ==========================================
    static std::string shearRow(float first_shear_value = 0.0f) {
        std::ostringstream row;
        for (int column = 0; column < LensingConfig::shear_cat_ncols; ++column) {
            row << (column == 0 ? "" : " ")
                << (column == 0 ? first_shear_value : 0.0f);
        }
        return row.str();
    }

    // ==========================================
    // Function: Write explicit shear and original catalog data rows
    // Method: Keep live headers while allowing valid, invalid, and malformed fixtures.
    // ==========================================
    void writeCatalogRows(const std::vector<std::string>& shear_rows,
                          const std::vector<std::string>& original_rows) const {
        std::ofstream shear(shear_file_, std::ios::trunc);
        std::ofstream orig(orig_file_, std::ios::trunc);
        require(static_cast<bool>(shear) && static_cast<bool>(orig),
                "synthetic Stage-9 inputs must open");

        for (int column = 0; column < LensingConfig::shear_cat_ncols; ++column) {
            shear << (column == 0 ? "s0" : " s" + std::to_string(column));
        }
        shear << '\n';
        orig << "ra dec\n";
        for (const std::string& row : shear_rows) {
            shear << row << '\n';
        }
        for (const std::string& row : original_rows) {
            orig << row << '\n';
        }
    }

    // ==========================================
    // Function: Replace the shear catalog with a zero-line invalid file
    // Method: Truncate only the shear path so preflight must reject the missing header.
    // ==========================================
    void writeEmptyShearCatalog() const {
        std::ofstream shear(shear_file_, std::ios::trunc);
        require(static_cast<bool>(shear),
                "empty synthetic shear catalog must open");
    }

    // ==========================================
    // Function: Write the common zero-or-one-row lifecycle fixtures
    // Method: Adapt the established cases to the explicit paired-row writer.
    // ==========================================
    void writeCatalogs(bool with_shear_data, bool with_original_data,
                       float first_shear_value = 0.0f) const {
        std::vector<std::string> shear_rows;
        std::vector<std::string> original_rows;
        if (with_shear_data) {
            shear_rows.push_back(shearRow(first_shear_value));
        }
        if (with_original_data) {
            original_rows.emplace_back("180 0");
        }
        writeCatalogRows(shear_rows, original_rows);
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
            1, std::vector<std::string>{image_file_}, root_.string(),
            kExposureIndex, chi2);
    }

    const std::string& outputFile() const noexcept { return output_file_; }
    int chipIndex() const { return UniversalUtils::getChipId(image_file_); }

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
// Function: Verify missing and unreadable norm products are fatal
// Method: Exercise both integrity states through Stage 9, then restore a valid
//         sentinel so the ordinary output lifecycle cases can continue.
// ==========================================
void testNormIntegrityFailures(TemporaryCatalogTree& tree) {
    tree.writeCatalogs(true, true);
    tree.removeNorm();
    requireAbort(
        [&]() { tree.combine(0.0f); },
        "Missing norm must terminate Stage 9");

    tree.writeMalformedNorm();
    requireAbort(
        [&]() { tree.combine(0.0f); },
        "ReadError norm must terminate Stage 9");
    tree.writeNorm(-1.0f);
}

// ==========================================
// Function: Verify lazy creation and the complete Stage-9 schema
// Method: Require exposure/CCD header and data order, exact width, and terminal Chi2.
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
    std::istringstream header_stream(header);
    std::vector<std::string> header_fields;
    std::string header_field;
    while (header_stream >> header_field) {
        header_fields.push_back(header_field);
    }

    std::istringstream values_stream(row);
    std::vector<double> values;
    double value = 0.0;
    while (values_stream >> value) {
        values.push_back(value);
    }
    const std::size_t expected_columns =
        2U + 2U
        + static_cast<std::size_t>(LensingConfig::shear_cat_ncols) + 1U;
    require(header_fields.size() == expected_columns
                && values.size() == expected_columns,
            "combined header and row must have the complete schema width");
    require(header_fields[2] == "EXPO_NUM"
                && header_fields[3] == "ccD_NUM",
            "EXPO_NUM must immediately precede ccD_NUM in the header");
    require(std::abs(values[2] - kExposureIndex) < 1.0e-7
                && std::abs(values[3] - tree.chipIndex()) < 1.0e-7,
            "data exposure and CCD identities must match the header order");
    require(header_fields.back() == "Chi2"
                && std::abs(values.back() - chi2) < 1.0e-7,
            "combined catalog data must end in the exposure Chi2");
}

// ==========================================
// Function: Verify scientific rejection advances both paired input streams
// Method: Reject pair one, retain pair two, and inspect its external fields.
// ==========================================
void testScientificInvalidPairConsumption(TemporaryCatalogTree& tree) {
    tree.writeCatalogRows(
        {TemporaryCatalogTree::shearRow(-99999.0f),
         TemporaryCatalogTree::shearRow()},
        {"101 1", "202 2"});
    tree.combine(0.0f);

    std::ifstream output(tree.outputFile());
    std::string header;
    std::string row;
    std::string extra;
    require(std::getline(output, header) && std::getline(output, row)
                && !std::getline(output, extra),
            "one rejected and one valid pair must produce exactly one data row");

    std::istringstream values(row);
    double external_ra = 0.0;
    double external_dec = 0.0;
    require(static_cast<bool>(values >> external_ra >> external_dec),
            "retained pair must contain its external catalog fields");
    require(std::abs(external_ra - 202.0) < 1.0e-12
                && std::abs(external_dec - 2.0) < 1.0e-12,
            "valid shear row must remain paired with the second original row");
}

// ==========================================
// Function: Verify Stage-9 sentinel rejection preserves row pairing
// Method: Combine one full-width sentinel shear row with its original row and
//         require a header-only science catalog.
// ==========================================
void testSentinelOutput(TemporaryCatalogTree& tree) {
    tree.writeCatalogs(true, true, -99999.0f);
    tree.combine(0.0f);

    std::ifstream output(tree.outputFile());
    std::string header;
    std::string row;
    require(std::getline(output, header) && !std::getline(output, row),
            "sentinel shear row must be consumed but omitted from science output");
}

// ==========================================
// Function: Verify structural pairing failures cannot silently desynchronize rows
// Method: Verify two preflight attempts plus incomplete, blank, and empty row failures.
// ==========================================
void testFatalPairingCases(TemporaryCatalogTree& tree) {
    const std::string valid = TemporaryCatalogTree::shearRow();

    const std::string short_orig_diagnostics = requireAbort(
        [&]() {
            tree.writeCatalogRows({valid, valid}, {"101 1"});
            tree.combine(0.0f);
        },
        "external catalog ending before shear must fail fast");
    for (const char* expected : {
             "attempt1_shear_lines=3", "attempt1_orig_lines=2",
             "attempt2_shear_lines=3", "attempt2_orig_lines=2"}) {
        require(short_orig_diagnostics.find(expected) != std::string::npos,
                std::string("short-orig preflight diagnostic must contain ")
                    + expected);
    }

    const std::string long_orig_diagnostics = requireAbort(
        [&]() {
            tree.writeCatalogRows({valid}, {"101 1", "202 2"});
            tree.combine(0.0f);
        },
        "unpaired trailing external catalog row must fail fast");
    for (const char* expected : {
             "attempt1_shear_lines=2", "attempt1_orig_lines=3",
             "attempt2_shear_lines=2", "attempt2_orig_lines=3"}) {
        require(long_orig_diagnostics.find(expected) != std::string::npos,
                std::string("long-orig preflight diagnostic must contain ")
                    + expected);
    }

    requireAbort(
        [&]() {
            tree.writeCatalogRows({valid, "0 0", valid},
                                  {"101 1", "202 2", "303 3"});
            tree.combine(0.0f);
        },
        "incomplete paired shear row must fail fast");

    requireAbort(
        [&]() {
            tree.writeCatalogRows({valid, "", valid},
                                  {"101 1", "202 2", "303 3"});
            tree.combine(0.0f);
        },
        "blank paired shear row must fail fast");

    requireAbort(
        [&]() {
            tree.writeCatalogRows({valid, valid}, {"101 1", ""});
            tree.combine(0.0f);
        },
        "empty paired external catalog row must fail fast");

    requireAbort(
        [&]() {
            tree.writeEmptyShearCatalog();
            tree.combine(0.0f);
        },
        "zero-line shear catalog must fail before row subtraction");
}

}  // namespace

// ==========================================
// Function: Run Stage-9 output lifecycle regressions
// Method: Cover lifecycle, lockstep consumption, schema, and fatal mismatches.
// ==========================================
int main() {
    TemporaryCatalogTree tree;
    testNoOutputCases(tree);
    testNormIntegrityFailures(tree);
    testLiveOutput(tree);
    testScientificInvalidPairConsumption(tree);
    testSentinelOutput(tree);
    testFatalPairingCases(tree);
    std::cout << "CatalogCombiner lifecycle tests passed\n";
    return EXIT_SUCCESS;
}
