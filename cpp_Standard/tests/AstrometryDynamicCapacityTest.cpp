#include "Astrometry.hpp"
#include "LensingConfig.hpp"
#include "OutputLayout.hpp"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> EXPO_FILE;

namespace {

// ==========================================
// Function: Stop the dynamic-capacity test when a requirement is not met.
// Method: Print a focused failure message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Astrometry dynamic-capacity test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Create an isolated temporary directory for one file-level case.
// Method: Combine the system temporary root with a monotonic-clock token.
// ==========================================
std::filesystem::path makeTemporaryDirectory(const std::string& case_name) {
    const auto token = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("fq_astrometry_" + case_name + "_" + std::to_string(token));
    std::filesystem::create_directories(directory);
    return directory;
}

// ==========================================
// Function: Verify source rows grow beyond every former count reservation.
// Method: Place 10,100 isolated detections farther apart than the segmentation
//         linking radius and require every source to reach the output catalog.
// ==========================================
void testSourceRowGrowth() {
    constexpr int columns = 101;
    constexpr int rows = 100;
    constexpr int spacing = 8;
    constexpr int margin = 3;
    constexpr int expected_sources = columns * rows;
    constexpr int nx = 2 * margin + (columns - 1) * spacing + 1;
    constexpr int ny = 2 * margin + (rows - 1) * spacing + 1;

    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 0.0f);
    std::vector<int> weight(image.size(), 1);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int x = margin + column * spacing;
            const int y = margin + row * spacing;
            image[static_cast<std::size_t>(y) * nx + x] = 10.0f;
        }
    }

    int source_count = 0;
    std::vector<double> x;
    std::vector<double> y;
    Astrometry::getAstrometryCatalog(
        nx, ny, image, weight, source_count, x, y);

    require(source_count == expected_sources,
            "all 10,100 isolated sources must survive dynamic growth");
    require(x.size() == expected_sources && y.size() == expected_sources,
            "coordinate vectors must contain every detected source");
}

// ==========================================
// Function: Verify one connected component can exceed 10,000 pixels.
// Method: Fill a 101 by 101 component but keep only its central peak above the
//         half-peak morphology threshold so the component remains selectable.
// ==========================================
void testConnectedComponentGrowth() {
    constexpr int nx = 101;
    constexpr int ny = 101;
    constexpr int center = 50;

    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 5.0f);
    std::vector<int> weight(image.size(), 1);
    image[static_cast<std::size_t>(center) * nx + center] = 20.0f;

    int source_count = 0;
    std::vector<double> x;
    std::vector<double> y;
    Astrometry::getAstrometryCatalog(
        nx, ny, image, weight, source_count, x, y);

    require(source_count == 1 && x.size() == 1 && y.size() == 1,
            "the 10,201-pixel component must remain selectable");
    require(std::fabs(x[0] - (center + 1.0)) < 1.0e-12 &&
                std::fabs(y[0] - (center + 1.0)) < 1.0e-12,
            "the large component must preserve its one-based peak centroid");
}

// ==========================================
// Function: Verify reference-catalog ingestion exceeds 10,000 rows.
// Method: Generate 10,001 in-bounds sky rows with no image detections and read
//         the reference count serialized by the production Stage-1 function.
// ==========================================
void testReferenceCatalogGrowth() {
    constexpr int reference_count = 10001;
    constexpr int nx = 512;
    constexpr int ny = 512;
    const std::filesystem::path directory =
        makeTemporaryDirectory("reference_rows");
    const std::filesystem::path catalog_path = directory / "reference.cat";
    const std::filesystem::path output_path = directory / "chip_astro.dat";

    {
        std::ofstream catalog(catalog_path);
        require(static_cast<bool>(catalog), "temporary reference catalog must open");
        catalog << "ra dec\n";
        for (int row = 0; row < reference_count; ++row) {
            catalog << "180 0\n";
        }
    }

    WCSParams wcs;
    wcs.crpix[0] = 256.0;
    wcs.crpix[1] = 256.0;
    wcs.crval[0] = 180.0;
    wcs.crval[1] = 0.0;
    wcs.cd[0][0] = 0.001;
    wcs.cd[1][1] = 0.001;

    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 0.0f);
    std::vector<int> weight(image.size(), 1);
    int process_error = 0;
    Astrometry::genAstrometryData(
        catalog_path.string(), nx, ny, image, weight, wcs,
        output_path.string(), process_error);

    std::ifstream output(output_path);
    require(process_error == 0 && static_cast<bool>(output),
            "production reference ingestion must write an output catalog");
    double discarded_value = 0.0;
    for (int value = 0; value < 8; ++value) {
        require(static_cast<bool>(output >> discarded_value),
                "astrometry output must contain both WCS rows");
    }
    int matched_count = 0;
    int user_count = 0;
    int serialized_reference_count = 0;
    require(static_cast<bool>(output >> matched_count >> user_count >>
                              serialized_reference_count),
            "astrometry output must contain catalog counts");
    require(serialized_reference_count == reference_count,
            "all 10,001 reference rows must be serialized");

    output.close();
    std::filesystem::remove_all(directory);
}

// ==========================================
// Function: Verify the post-detection astrometry user-source selection
// Method: Keep the dynamic catalog complete, then require genAstrometryData to serialize
//         min(n_detected, n_user_max) as the count passed to pattern matching.
// ==========================================
void testUserSelectionLimit() {
    constexpr int columns = 20;
    constexpr int rows = 15;
    constexpr int spacing = 8;
    constexpr int margin = 3;
    constexpr int nx = 2 * margin + (columns - 1) * spacing + 1;
    constexpr int ny = 2 * margin + (rows - 1) * spacing + 1;
    const std::filesystem::path directory = makeTemporaryDirectory("user_limit");
    const std::filesystem::path catalog_path = directory / "reference.cat";

    {
        std::ofstream catalog(catalog_path);
        require(static_cast<bool>(catalog), "selection reference catalog must open");
        catalog << "ra dec\n";
    }

    WCSParams wcs{};
    wcs.crpix[0] = 0.5 * nx;
    wcs.crpix[1] = 0.5 * ny;
    wcs.crval[0] = 180.0;
    wcs.crval[1] = 0.0;
    wcs.cd[0][0] = 0.001;
    wcs.cd[1][1] = 0.001;

    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 0.0f);
    std::vector<int> weight(image.size(), 1);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int x = margin + column * spacing;
            const int y = margin + row * spacing;
            image[static_cast<std::size_t>(y) * nx + x] =
                10.0f + static_cast<float>(row * columns + column);
        }
    }

    int detected = 0;
    std::vector<double> xs;
    std::vector<double> ys;
    Astrometry::getAstrometryCatalog(nx, ny, image, weight, detected, xs, ys);
    require(detected == columns * rows && xs.size() == static_cast<std::size_t>(detected)
                && ys.size() == static_cast<std::size_t>(detected),
            "selection regression must detect all 300 sources before truncation");

    const auto readUserCount = [&](const std::filesystem::path& output_path,
                                   int expected_count) {
        int process_error = 0;
        Astrometry::genAstrometryData(
            catalog_path.string(), nx, ny, image, weight, wcs,
            output_path.string(), process_error);
        std::ifstream output(output_path);
        require(process_error == 0 && static_cast<bool>(output),
                "selection regression must write an astrometry output");
        double discarded_value = 0.0;
        for (int value = 0; value < 8; ++value) {
            require(static_cast<bool>(output >> discarded_value),
                    "selection output must contain both WCS rows");
        }
        int matched_count = 0;
        int user_count = 0;
        int reference_count = 0;
        require(static_cast<bool>(output >> matched_count >> user_count >> reference_count),
                "selection output must contain catalog counts");
        require(user_count == expected_count && reference_count == 0,
                "_astro.dat must record the number of selected image detections");
    };

    readUserCount(directory / "many_astro.dat", LensingConfig::n_user_max);

    std::fill(image.begin(), image.end(), 0.0f);
    for (int index = 0; index < 37; ++index) {
        const int row = index / columns;
        const int column = index % columns;
        const int x = margin + column * spacing;
        const int y = margin + row * spacing;
        image[static_cast<std::size_t>(y) * nx + x] = 10.0f + static_cast<float>(index);
    }
    readUserCount(directory / "few_astro.dat", 37);
    std::filesystem::remove_all(directory);
}

// ==========================================
// Function: Verify the Stage-2 reader consumes rows beyond 10,000.
// Method: Put a malformed sentinel at row 10,001 and require the production
//         reader to report that row instead of silently stopping at row 10,000.
// ==========================================
void testMatchedCatalogReaderGrowth() {
    constexpr int declared_rows = 10001;
    const std::filesystem::path directory =
        makeTemporaryDirectory("matched_rows");
    const std::filesystem::path output_root = directory / "output";
    const std::filesystem::path head_directory =
        output_root / "astrometry" / "Head";
    std::filesystem::create_directories(head_directory);

    const std::string image_path = (directory / "exposure_01.fits").string();
    const std::filesystem::path catalog_path =
        OutputLayout::chipPath(
            output_root.string(), "astrometry/dat_Astro", "exposure_01",
            "_astro.dat");
    std::filesystem::create_directories(catalog_path.parent_path());
    {
        std::ofstream catalog(catalog_path);
        require(static_cast<bool>(catalog), "temporary matched catalog must open");
        catalog << "256 256 180 0\n";
        catalog << "0.001 0 0 0.001\n";
        catalog << declared_rows << " " << declared_rows << " "
                << declared_rows << '\n';
        for (int row = 0; row < declared_rows - 1; ++row) {
            catalog << "180 0 256 256\n";
        }
        catalog << "malformed sentinel row\n";
    }
    require(std::filesystem::is_regular_file(catalog_path),
            "matched catalog must exist at the production chip path");

    std::ostringstream captured_output;
    std::streambuf* original_output = std::cout.rdbuf(captured_output.rdbuf());
    Astrometry::getAstrometry({image_path}, 1, output_root.string());
    std::cout.rdbuf(original_output);

    require(captured_output.str().find("astro.dat error") != std::string::npos,
            "the Stage-2 reader must attempt to parse row 10,001");
    std::filesystem::remove_all(directory);
}

}  // namespace

// ==========================================
// Function: Run Astrometry dynamic-capacity regression cases.
// Method: Exercise source, component, reference, and matched-row paths past old limits.
// ==========================================
int main() {
    testSourceRowGrowth();
    testConnectedComponentGrowth();
    testReferenceCatalogGrowth();
    testUserSelectionLimit();
    testMatchedCatalogReaderGrowth();
    std::cout << "Astrometry dynamic-capacity tests passed\n";
    return EXIT_SUCCESS;
}
