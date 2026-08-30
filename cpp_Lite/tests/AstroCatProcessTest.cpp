#include "process_astrocat/process_astrocat.hpp"

#include "general/MPIUtils.hpp"
#include "general/MPIScheduler.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

namespace fs = std::filesystem;

// ==========================================
// Function: Record one process_astrocat test assertion
// Method: Print a focused rank-zero diagnostic and preserve cumulative status.
// ==========================================
bool require(bool condition, const std::string& message) {
    if (!condition && MPIScheduler::state.rank == 0) {
        std::cerr << "AstroCatProcessTest: " << message << std::endl;
    }
    return condition;
}

// ==========================================
// Function: Write one synthetic raw Gaia fixture
// Method: Truncate the exact test path and serialize caller-supplied text.
// ==========================================
bool writeFixture(const fs::path& path, const std::string& content) {
    std::ofstream output(path, std::ios::trunc);
    output << content;
    return output.good();
}

// ==========================================
// Function: Format a double with round-trip precision
// Method: Use max_digits10 so nextafter fixtures survive text ingestion.
// ==========================================
std::string formatDouble(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10)
           << value;
    return stream.str();
}

// ==========================================
// Function: Count and load every generated Gaia output row
// Method: Scan direct .dat tiles, require the normalized header, and parse all
//         coordinate pairs into one deterministic vector.
// ==========================================
bool loadGeneratedRows(const fs::path& output_directory,
                       std::vector<std::pair<double, double>>& rows,
                       std::vector<std::string>& filenames,
                       std::string& error) {
    rows.clear();
    filenames.clear();
    for (const fs::directory_entry& entry :
         fs::directory_iterator(output_directory)) {
        if (!entry.is_regular_file()
            || entry.path().extension() != ".dat") {
            continue;
        }
        filenames.push_back(entry.path().filename().string());
        std::ifstream input(entry.path());
        std::string header;
        if (!std::getline(input, header) || header != "RA    DEC") {
            error = "generated tile has the wrong header: "
                    + entry.path().string();
            return false;
        }
        double ra = 0.0;
        double dec = 0.0;
        while (input >> ra >> dec) {
            rows.emplace_back(ra, dec);
        }
        if (!input.eof()) {
            error = "generated tile contains a malformed row: "
                    + entry.path().string();
            return false;
        }
    }
    std::sort(rows.begin(), rows.end());
    std::sort(filenames.begin(), filenames.end());
    error.clear();
    return true;
}

// ==========================================
// Function: Test one coordinate's exact presence
// Method: Compare parsed doubles exactly because fixtures use round-trip text.
// ==========================================
bool containsRow(const std::vector<std::pair<double, double>>& rows,
                 double ra,
                 double dec) {
    return std::find(rows.begin(), rows.end(), std::make_pair(ra, dec))
           != rows.end();
}

// ==========================================
// Function: Test one generated filename's exact presence
// Method: Search the sorted basename vector for the Type-2 lookup contract.
// ==========================================
bool containsFilename(const std::vector<std::string>& filenames,
                      const std::string& filename) {
    return std::binary_search(filenames.begin(), filenames.end(), filename);
}

}  // namespace

// ==========================================
// Function: Exercise process_astrocat MPI tiling and lifecycle contracts
// Method: Build synthetic whole-file jobs covering headers, ULP de-duplication,
//         tile edges, overwrite cleanup, and independent output paths.
// ==========================================
int main(int argc, char** argv) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::state.rank;
    bool local_ok = true;
    std::string error;

    std::string root_text;
    if (rank == 0) {
        root_text = (fs::temp_directory_path()
                     / ("fourier_quad_astrocat_test_"
                        + std::to_string(getpid())))
                        .string();
    }
    local_ok = require(
                   MPIUtils::broadcastString(root_text, 0, error),
                   "temporary-root broadcast failed: " + error)
               && local_ok;
    const fs::path root(root_text);
    const fs::path input = root / "raw";
    const fs::path output = root / "tiles";

    const double same_ra = 10.25;
    const double grouped_ra = 20.25;
    const double grouped_ra_next = std::nextafter(
        grouped_ra, std::numeric_limits<double>::infinity());
    const double two_ulp_ra = 30.25;
    const double two_ulp_ra_next = std::nextafter(
        std::nextafter(two_ulp_ra,
                       std::numeric_limits<double>::infinity()),
        std::numeric_limits<double>::infinity());
    const double ra_boundary_left = std::nextafter(
        180.0, -std::numeric_limits<double>::infinity());
    const double dec_boundary_down = std::nextafter(
        20.0, -std::numeric_limits<double>::infinity());

    if (rank == 0) {
        fs::create_directories(input);
        std::ostringstream first;
        first << formatDouble(same_ra) << " 0.25\n"
              << formatDouble(grouped_ra) << " 5.25\n"
              << formatDouble(grouped_ra) << " 5.75\n"
              << formatDouble(two_ulp_ra) << " 6.25\n"
              << formatDouble(ra_boundary_left) << " 10.25\n"
              << "40.25 " << formatDouble(dec_boundary_down) << "\n"
              << "50 90\n"
              << "360 1.25\n"
              << "malformed row\n";
        std::ostringstream second;
        second << formatDouble(same_ra) << ",0.25\n"
               << formatDouble(grouped_ra_next) << ",5.25\n"
               << formatDouble(two_ulp_ra_next) << ",6.25\n"
               << "180,10.25\n"
               << "40.25,20\n";
        local_ok = writeFixture(input / "a.dat", first.str())
                   && writeFixture(input / "b.csv", second.str());
    }
    MPIScheduler::barrier();

    ProcessAstrocat::Config config;
    config.input_directory = input;
    config.output_directory = output;
    config.add_header = true;
    config.existing_policy = ProcessAstrocat::ExistingPolicy::Fail;
    local_ok = require(process_astrocat(config) == 0,
                       "primary no-header run failed")
               && local_ok;

    if (rank == 0) {
        std::vector<std::pair<double, double>> rows;
        std::vector<std::string> filenames;
        local_ok = require(loadGeneratedRows(
                               output, rows, filenames, error), error)
                   && local_ok;
        local_ok = require(rows.size() == 9,
                           "exact, one-ULP, and boundary duplicates were not removed")
                   && local_ok;
        local_ok = require(containsRow(rows, same_ra, 0.25),
                           "exact duplicate canonical row is missing")
                   && local_ok;
        local_ok = require(containsRow(rows, grouped_ra, 5.25)
                               && containsRow(rows, grouped_ra, 5.75)
                               && !containsRow(rows, grouped_ra_next, 5.25),
                           "non-adjacent same-tile one-ULP handling is incorrect")
                   && local_ok;
        local_ok = require(containsRow(rows, two_ulp_ra, 6.25)
                               && containsRow(rows, two_ulp_ra_next, 6.25),
                           "two-ULP coordinates were incorrectly merged")
                   && local_ok;
        local_ok = require(containsRow(rows, ra_boundary_left, 10.25)
                               && !containsRow(rows, 180.0, 10.25),
                           "cross-RA-tile canonical de-duplication failed")
                   && local_ok;
        local_ok = require(containsRow(rows, 40.25, dec_boundary_down)
                               && !containsRow(rows, 40.25, 20.0),
                           "cross-Dec-tile canonical de-duplication failed")
                   && local_ok;
        local_ok = require(containsRow(rows, 0.0, 1.25),
                           "RA=360 was not normalized to the zero-degree tile")
                   && local_ok;
        local_ok = require(containsFilename(
                               filenames,
                               "des_y6_RA_000_001_Dec_p01_p02.dat")
                               && containsFilename(
                                   filenames,
                                   "des_y6_RA_050_051_Dec_p89_p90.dat"),
                           "Type-2 filename boundary contract is incorrect")
                   && local_ok;
    }
    MPIScheduler::barrier();

    local_ok = require(process_astrocat(config) != 0,
                       "existing-policy fail accepted generated tiles")
               && local_ok;
    if (rank == 0) {
        local_ok = writeFixture(output / "keep.txt", "unrelated\n")
                   && local_ok;
    }
    MPIScheduler::barrier();
    config.existing_policy = ProcessAstrocat::ExistingPolicy::Overwrite;
    local_ok = require(process_astrocat(config) == 0,
                       "existing-policy overwrite failed")
               && local_ok;
    if (rank == 0) {
        local_ok = require(fs::exists(output / "keep.txt"),
                           "overwrite removed an unrelated output file")
                   && local_ok;
        for (const fs::directory_entry& entry : fs::directory_iterator(output)) {
            const std::string name = entry.path().filename().string();
            local_ok = require(name.find("staging") == std::string::npos
                                   && entry.path().extension() != ".part",
                               "an intermediate output artifact was created")
                       && local_ok;
        }
    }
    MPIScheduler::barrier();

    const fs::path header_input = root / "header_raw";
    const fs::path header_output = root / "header_tiles";
    if (rank == 0) {
        fs::create_directories(header_input);
        local_ok = writeFixture(header_input / "with_header.dat",
                                "RA DEC\n70.25 -5.25\n")
                   && local_ok;
    }
    MPIScheduler::barrier();
    ProcessAstrocat::Config header_config;
    header_config.input_directory = header_input;
    header_config.output_directory = header_output;
    header_config.add_header = false;
    local_ok = require(process_astrocat(header_config) == 0,
                       "input-header skip run failed")
               && local_ok;
    if (rank == 0) {
        std::vector<std::pair<double, double>> rows;
        std::vector<std::string> filenames;
        local_ok = require(loadGeneratedRows(
                               header_output, rows, filenames, error), error)
                   && local_ok;
        local_ok = require(rows.size() == 1
                               && containsRow(rows, 70.25, -5.25),
                           "input header was not skipped exactly once")
                   && local_ok;
    }

    ProcessAstrocat::Config overlapping;
    overlapping.input_directory = input;
    overlapping.output_directory = input / "nested";
    bool rejected_overlap = false;
    try {
        ProcessAstrocat::normalizeAndValidateConfig(overlapping);
    } catch (const std::invalid_argument&) {
        rejected_overlap = true;
    }
    local_ok = require(rejected_overlap,
                       "nested input/output paths were accepted")
               && local_ok;

    bool global_ok = false;
    const bool reduction_ok = MPIUtils::allRanksSucceeded(
        local_ok, global_ok, error);
    MPIScheduler::barrier();
    if (rank == 0) {
        fs::remove_all(root);
    }
    MPIScheduler::finalize();
    return reduction_ok && global_ok ? 0 : 1;
}
