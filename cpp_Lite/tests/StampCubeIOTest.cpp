#include "FitsIO.hpp"

#include <fitsio.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

// ==========================================
// Class: Remove one temporary FITS artifact after a test case
// Method: Delete the owned path during normal return or exception unwinding.
// ==========================================
class TemporaryFitsFile {
public:
    explicit TemporaryFitsFile(const std::string& label)
        : path_((std::filesystem::temp_directory_path()
                 / ("fq_" + label + "_" + std::to_string(::getpid()) + ".fits"))
                    .string()) {}

    ~TemporaryFitsFile() {
        std::remove(path_.c_str());
    }

    // ==========================================
    // Function: Return the temporary FITS pathname
    // Method: Expose the immutable generated path to the test operations.
    // ==========================================
    const std::string& path() const noexcept {
        return path_;
    }

private:
    std::string path_;
};

// ==========================================
// Function: Stop the test program when a requirement is not met
// Method: Print a focused failure message and return a non-zero process code.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "StampCubeIO test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Construct one orientation-sensitive stamp collection
// Method: Encode stamp, row, and column independently in every float value.
// ==========================================
std::vector<float> makeCube(int nx, int ny, int count) {
    std::vector<float> cube(static_cast<std::size_t>(count) * ny * nx);
    for (int stamp = 0; stamp < count; ++stamp) {
        for (int row = 0; row < ny; ++row) {
            for (int col = 0; col < nx; ++col) {
                const std::size_t index =
                    (static_cast<std::size_t>(stamp) * ny + row) * nx + col;
                cube[index] = static_cast<float>(10000 * stamp + 100 * row + col);
            }
        }
    }
    return cube;
}

// ==========================================
// Function: Verify the FITS cube axes and self-description keywords
// Method: Read the raw primary-header values through CFITSIO and compare them
//         with the public stamp-cube contract.
// ==========================================
void verifyHeader(const std::string& filename, int nx, int ny, int count) {
    fitsfile* fptr = nullptr;
    int status = 0;
    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    require(status == 0, "cannot open generated FITS header");

    int naxis = 0;
    long naxes[3] = {0, 0, 0};
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 3, naxes, &status);

    char format[FLEN_VALUE] = {};
    char order[FLEN_VALUE] = {};
    fits_read_key(fptr, TSTRING, "FQFMT", format, nullptr, &status);
    fits_read_key(fptr, TSTRING, "FQORDER", order, nullptr, &status);

    int closeStatus = 0;
    fits_close_file(fptr, &closeStatus);
    require(status == 0 && closeStatus == 0, "cannot read or close FITS header");
    require(naxis == 3, "NAXIS must equal 3");
    require(naxes[0] == nx && naxes[1] == ny && naxes[2] == count,
            "FITS axes must be (nx, ny, count)");
    require(std::string(format) == "STAMP_CUBE", "FQFMT keyword mismatch");
    require(std::string(order) == "X,Y,STAMP", "FQORDER keyword mismatch");
}

// ==========================================
// Function: Exercise non-square cube round-trip and header orientation
// Method: Write and read a 5x3x4 asymmetric cube and compare every value.
// ==========================================
void testNonSquareRoundTrip() {
    constexpr int nx = 5;
    constexpr int ny = 3;
    constexpr int count = 4;
    const std::vector<float> expected = makeCube(nx, ny, count);
    TemporaryFitsFile file("stamp_cube_round_trip");

    require(FitsIO::writeStampCube(file.path(), nx, ny, count, expected),
            "cube write returned false");
    verifyHeader(file.path(), nx, ny, count);

    FitsIO::StampCubeShape shape;
    std::vector<float> actual;
    require(FitsIO::readStampCube(file.path(), shape, actual),
            "cube read returned false");
    require(shape.matches(nx, ny, count), "decoded cube shape mismatch");
    require(actual == expected, "non-square cube changed pixel orientation");
    require(actual[(static_cast<std::size_t>(2) * ny + 1) * nx + 4] == 20104.0f,
            "orientation sentinel stamp=2,row=1,col=4 mismatch");
}

// ==========================================
// Function: Exercise selected-plane compaction
// Method: Select input stamps 1, 3, and 4 and require a dense three-plane cube
//         in the original input order.
// ==========================================
void testSelectedPlanes() {
    constexpr int nx = 5;
    constexpr int ny = 3;
    constexpr int inputCount = 5;
    const std::vector<float> input = makeCube(nx, ny, inputCount);
    const std::vector<int> selection = {0, 1, 0, 1, 1};
    TemporaryFitsFile file("selected_stamp_cube");

    require(FitsIO::writeSelectedStampCube(
                file.path(), nx, ny, input, inputCount, selection, 1),
            "selected cube write returned false");

    FitsIO::StampCubeShape shape;
    std::vector<float> actual;
    require(FitsIO::readStampCube(file.path(), shape, actual),
            "selected cube read returned false");
    require(shape.matches(nx, ny, 3), "selected cube must contain three planes");

    const std::size_t planeSize = static_cast<std::size_t>(nx) * ny;
    const int expectedStamp[3] = {1, 3, 4};
    for (int outputStamp = 0; outputStamp < 3; ++outputStamp) {
        const std::size_t inputOffset =
            static_cast<std::size_t>(expectedStamp[outputStamp]) * planeSize;
        const std::size_t outputOffset =
            static_cast<std::size_t>(outputStamp) * planeSize;
        require(std::equal(input.begin() + inputOffset,
                           input.begin() + inputOffset + planeSize,
                           actual.begin() + outputOffset),
                "selected plane content or order mismatch");
    }
}

}  // namespace

// ==========================================
// Function: Run the focused three-dimensional stamp-cube regression suite
// Method: Validate orientation, header metadata, exact shape recovery, and
//         selected-plane compaction.
// ==========================================
int main() {
    testNonSquareRoundTrip();
    testSelectedPlanes();
    std::cout << "StampCubeIO tests passed\n";
    return EXIT_SUCCESS;
}
