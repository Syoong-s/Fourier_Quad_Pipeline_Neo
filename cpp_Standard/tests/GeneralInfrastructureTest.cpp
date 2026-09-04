#include "general/ExposureList.hpp"
#include "general/MPIScheduler.hpp"
#include "general/MPIUtils.hpp"
#include "general/PathUtils.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

// ==========================================
// Function: Record one test assertion
// Method: Print the failed contract and return its Boolean result to the caller.
// ==========================================
bool require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "GeneralInfrastructureTest: " << message << std::endl;
    }
    return condition;
}

// ==========================================
// Function: Write one synthetic text fixture
// Method: Truncate the exact temporary path and publish the supplied content.
// ==========================================
bool writeFixture(const std::filesystem::path& path,
                  const std::string& content) {
    std::ofstream output(path);
    output << content;
    return output.good();
}

}  // namespace

// ==========================================
// Function: Exercise shared exposure-list, path, and MPI helpers
// Method: Use deterministic temporary fixtures and singleton-safe collectives
//         to validate normal and rejected protocol cases.
// ==========================================
int main(int argc, char** argv) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::state.rank;

    std::string root_text;
    std::string error;
    if (rank == 0) {
        root_text = (std::filesystem::temp_directory_path()
                     / ("fourier_quad_general_test_" + std::to_string(getpid())))
                        .string();
    }
    bool local_ok = require(
        MPIUtils::broadcastString(root_text, 0, error),
        "temporary-root broadcast failed: " + error);
    const std::filesystem::path root =
        std::filesystem::path(root_text);
    if (rank == 0) {
        std::filesystem::create_directories(root / "dataset" / "science" / "expo");
        local_ok = writeFixture(root / "pipeline.list",
                                "\"/tmp/exposure_a.list\" 62\n"
                                "/tmp/exposure_b.list 61\n")
                   && writeFixture(root / "paths.list",
                                   "\"/tmp/catalog_a.cat\"\n"
                                   "/tmp/catalog_b.cat\n")
                   && writeFixture(root / "malformed.list",
                                   "/tmp/exposure_a.list bad\n")
                   && writeFixture(root / "truncated.list",
                                   "/tmp/exposure_a.list\n")
                   && writeFixture(root / "empty.list", "");
    }
    MPIScheduler::barrier();

    std::vector<ExposureList::Entry> entries;
    local_ok = require(
                   ExposureList::loadPipelineList(
                       (root / "pipeline.list").string(), entries, 2, error),
                   error)
               && local_ok;
    local_ok = require(entries.size() == 2
                           && entries[0].path == "/tmp/exposure_a.list"
                           && entries[0].chip_count == 62,
                       "pipeline records were not preserved")
               && local_ok;
    local_ok = require(
                   ExposureList::loadPipelineList(
                       (root / "pipeline.list").string(), entries, 0, error)
                       && entries.size() == 2,
                   "zero must select unlimited pipeline-list loading: " + error)
               && local_ok;
    local_ok = require(
                   !ExposureList::loadPipelineList(
                       (root / "truncated.list").string(), entries, 2, error),
                   "truncated pipeline record was accepted")
               && local_ok;
    local_ok = require(
                   !ExposureList::loadPipelineList(
                       (root / "pipeline.list").string(), entries, 1, error),
                   "maximum-entry guard accepted oversized input")
               && local_ok;
    local_ok = require(
                   !ExposureList::loadPipelineList(
                       (root / "malformed.list").string(), entries, 2, error),
                   "malformed pipeline record was accepted")
               && local_ok;
    local_ok = require(
                   !ExposureList::loadPipelineList(
                       (root / "empty.list").string(), entries, 2, error),
                   "empty pipeline list was accepted")
               && local_ok;

    std::vector<std::string> paths;
    local_ok = require(
                   ExposureList::loadPathList(
                       (root / "paths.list").string(), paths, 2, error),
                   error)
               && local_ok;
    local_ok = require(paths.size() == 2 && paths[0] == "/tmp/catalog_a.cat",
                       "path-only records were not preserved")
               && local_ok;

    const std::filesystem::path normalized_root =
        PathUtils::normalizedAbsolute(root);
    const std::filesystem::path image =
        normalized_root / "dataset" / "science" / "expo" / "chip.fits";
    std::filesystem::path parent;
    local_ok = require(PathUtils::parentAtLevel(image, 3, parent, error)
                           && parent == normalized_root / "dataset",
                       "fixed-level parent resolution failed")
               && local_ok;
    local_ok = require(!PathUtils::parentAtLevel(image, 100, parent, error),
                       "insufficient parent depth was accepted")
               && local_ok;
    local_ok = require(PathUtils::isPathWithin(image, normalized_root)
                           && !PathUtils::isPathWithin(normalized_root,
                                                       image.parent_path()),
                       "component-wise path containment failed")
               && local_ok;

    std::string broadcast_value = rank == 0 ? "shared-value" : std::string();
    local_ok = require(MPIUtils::broadcastString(
                           broadcast_value, 0, error)
                           && broadcast_value == "shared-value",
                       "string broadcast failed: " + error)
               && local_ok;
    std::string empty_broadcast = rank == 0 ? std::string() : "stale";
    local_ok = require(MPIUtils::broadcastString(
                           empty_broadcast, 0, error)
                           && empty_broadcast.empty(),
                       "empty-string broadcast failed: " + error)
               && local_ok;
    std::vector<std::string> broadcast_values =
        rank == 0 ? std::vector<std::string>{"", "alpha", "beta"}
                  : std::vector<std::string>{};
    local_ok = require(MPIUtils::broadcastStrings(
                           broadcast_values, 0, error)
                           && broadcast_values.size() == 3
                           && broadcast_values[1] == "alpha",
                       "string-vector broadcast failed: " + error)
               && local_ok;
    std::vector<std::string> empty_broadcast_values =
        rank == 0 ? std::vector<std::string>{}
                  : std::vector<std::string>{"stale"};
    local_ok = require(MPIUtils::broadcastStrings(
                           empty_broadcast_values, 0, error)
                           && empty_broadcast_values.empty(),
                       "empty string-vector broadcast failed: " + error)
               && local_ok;

    const int size = MPIScheduler::state.size;
    bool expected_global_success = true;
    const bool expected_reduction_ok = MPIUtils::allRanksSucceeded(
        rank == 0, expected_global_success, error);
    local_ok = require(expected_reduction_ok
                           && expected_global_success == (size == 1),
                       "collective failure propagation was incorrect: " + error)
               && local_ok;

    bool global_ok = false;
    const bool reduction_ok = MPIUtils::allRanksSucceeded(
        local_ok, global_ok, error);
    if (rank == 0) {
        std::filesystem::remove_all(root);
    }
    MPIScheduler::finalize();
    return reduction_ok && global_ok ? 0 : 1;
}
