#include "CatalogRowCount.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>

#include <sys/wait.h>
#include <unistd.h>

namespace {

// ==========================================
// Class: Own one temporary physical-row catalog
// Method: Write a required header plus an exact number of data lines and remove
//         the file after each test case.
// ==========================================
class TemporaryCatalog {
public:
    TemporaryCatalog(const std::string& label, std::size_t data_rows)
        : path_(std::filesystem::temp_directory_path()
                / ("fq_catalog_rows_" + label + "_"
                   + std::to_string(::getpid()) + ".dat")) {
        std::ofstream output(path_);
        output << "header\n";
        for (std::size_t row = 0; row < data_rows; ++row) {
            output << row << '\n';
        }
    }

    ~TemporaryCatalog() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    // ==========================================
    // Function: Return the temporary catalog path
    // Method: Expose the immutable path to the production row-count helper.
    // ==========================================
    std::string path() const {
        return path_.string();
    }

private:
    std::filesystem::path path_;
};

// ==========================================
// Function: Stop the row-count test when one invariant fails
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "CatalogRowCount test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Require the production mismatch guard to terminate a child process
// Method: Fork one isolated caller, invoke the fail-fast guard, and accept only
//         a nonzero exit or terminating signal as evidence of abort behavior.
// ==========================================
void requireMismatchAbort(const TemporaryCatalog& shear,
                          const TemporaryCatalog& orig,
                          const std::string& message) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed for mismatch case");
    if (child == 0) {
        CatalogCombiner::Internal::requireMatchingCatalogDataRows(
            shear.path(), orig.path());
        std::_Exit(EXIT_SUCCESS);
    }

    int status = 0;
    require(::waitpid(child, &status, 0) == child,
            "waitpid failed for mismatch case");
    require((WIFEXITED(status) && WEXITSTATUS(status) != EXIT_SUCCESS)
                || WIFSIGNALED(status),
            message);
}

// ==========================================
// Function: Exercise equal, shorter, longer, and zero-row catalog pairs
// Method: Count with the production helper and verify the exact Stage-9
//         equality decisions required by the preview plan.
// ==========================================
void testRowCountCases() {
    const TemporaryCatalog shear_equal("shear_equal", 100);
    const TemporaryCatalog orig_equal("orig_equal", 100);
    const TemporaryCatalog orig_short("orig_short", 99);
    const TemporaryCatalog orig_long("orig_long", 101);
    const TemporaryCatalog shear_zero("shear_zero", 0);
    const TemporaryCatalog orig_zero("orig_zero", 0);

    const std::size_t shear_rows =
        CatalogCombiner::Internal::countCatalogDataRows(shear_equal.path());
    require(shear_rows == 100, "shear catalog must contain 100 data rows");
    CatalogCombiner::Internal::requireMatchingCatalogDataRows(
        shear_equal.path(), orig_equal.path());
    requireMismatchAbort(
        shear_equal, orig_short,
        "100/99 rows must terminate through the production mismatch guard");
    requireMismatchAbort(
        shear_equal, orig_long,
        "100/101 rows must terminate through the production mismatch guard");
    require(CatalogCombiner::Internal::countCatalogDataRows(shear_zero.path())
                == CatalogCombiner::Internal::countCatalogDataRows(
                    orig_zero.path()),
            "0/0 rows must pass equality");
    CatalogCombiner::Internal::requireMatchingCatalogDataRows(
        shear_zero.path(), orig_zero.path());
}

}  // namespace

// ==========================================
// Function: Run focused Stage-9 physical row-count regression cases
// Method: Execute all synthetic catalog pairs and report one success line.
// ==========================================
int main() {
    testRowCountCases();
    std::cout << "CatalogRowCount tests passed\n";
    return EXIT_SUCCESS;
}
