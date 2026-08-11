#ifndef FITSIO_HPP
#define FITSIO_HPP

#include <string>
#include <vector>

struct WCSParams {
    double crpix[2] = {0.0, 0.0};
    double crval[2] = {0.0, 0.0};
    double cd[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
};

namespace FitsIO {
    // ==========================================
    // Enum: Classify a focused FITS pixel read
    // Method: Separate an absent path from an existing file that cannot supply a pixel.
    // ==========================================
    enum class PixelReadStatus {
        Ok,
        Missing,
        ReadError
    };

    // ==========================================
    // Structure: Describe one contiguous FITS stamp cube
    // Method: Map FITS axes (x, y, stamp) to the in-memory
    //         [stamp][row][col] collection contract.
    // ==========================================
    struct StampCubeShape {
        int nx = 0;
        int ny = 0;
        int count = 0;

        // ==========================================
        // Function: Compare a decoded cube shape with one pipeline contract
        // Method: Require exact x, y, and plane-count equality.
        // ==========================================
        bool matches(int expectedNx, int expectedNy,
                     int expectedCount) const noexcept {
            return nx == expectedNx && ny == expectedNy
                && count == expectedCount;
        }
    };

    // Utility to print cfitsio errors
    void printError(int status);

    // Read CCD number from primary header
    bool readCCDNUM(const std::string& filename, int& ccdNum);

    // Read image dimensions
    bool readPara(const std::string& filename, int& nx, int& ny);

    // Read 2D image data
    bool readImage(const std::string& filename, int& nx, int& ny, std::vector<float>& data);

    // ==========================================
    // Function: Read the first pixel of one two-dimensional FITS image
    // Method: Validate the image shape and transfer a single float without loading the CCD.
    // ==========================================
    PixelReadStatus readFirstPixel(const std::string& filename, float& value);

    // Read 2D image data and WCS WCSParams
    bool readImagePara(const std::string& filename, int& nx, int& ny, std::vector<float>& data, WCSParams& wcs);

    // Update WCS keywords in a FITS file
    bool updatePara(const std::string& filename, const WCSParams& wcs);

    // Create a new FITS file, copy headers from file1, and write array
    bool writeImageCopyHDU(const std::string& templateFile, const std::string& filename, int nx, int ny, const std::vector<float>& data);

    // Write a standard 2D float image
    bool writeImage(const std::string& filename, int nx, int ny, const std::vector<float>& data);

    // ==========================================
    // Function: Read one contiguous three-dimensional stamp cube
    // Method: Decode (x, y, stamp) axes from the FITS header and return pixels
    //         in [stamp][row][col] order without caller-supplied mosaic geometry.
    // ==========================================
    bool readStampCube(const std::string& filename, StampCubeShape& shape,
                       std::vector<float>& stamps);

    // ==========================================
    // Function: Write one contiguous three-dimensional stamp cube
    // Method: Require exactly count * ny * nx pixels and store axes
    //         (x, y, stamp) with self-describing Fourier_Quad keywords.
    // ==========================================
    bool writeStampCube(const std::string& filename, int nx, int ny, int count,
                        const std::vector<float>& stamps);

    // ==========================================
    // Function: Compact and write selected planes from one stamp collection
    // Method: Preserve input order, retain planes matching selectedValue, and
    //         publish an exact-size three-dimensional cube.
    // ==========================================
    bool writeSelectedStampCube(const std::string& filename, int nx, int ny,
                                const std::vector<float>& stamps,
                                int inputCount,
                                const std::vector<int>& selection,
                                int selectedValue);

    // Stateful class for writing serial stamps to a FITS file
    class FitsSerialWriter {
    public:
        FitsSerialWriter();
        ~FitsSerialWriter();

        bool init(const std::string& filename);
        bool writeStamp(int nx, int ny, const std::vector<float>& data, bool newHdu);
        bool writeKey(const std::string& keyName, int val, const std::string& comment);
        void close();

    private:
        void* fptr; // fitsfile* typecast to void* to avoid exposing cfitsio headers in FitsIO.hpp
        int status;
        std::string outputFilename;
    };
}

#endif // FITSIO_HPP
