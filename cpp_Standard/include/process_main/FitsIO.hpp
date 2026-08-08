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
    // Utility to print cfitsio errors
    void printError(int status);

    // Read CCD number from primary header
    bool readCCDNUM(const std::string& filename, int& ccdNum);

    // Read image dimensions
    bool readPara(const std::string& filename, int& nx, int& ny);

    // Read 2D image data
    bool readImage(const std::string& filename, int& nx, int& ny, std::vector<float>& data);

    // Read 2D image data and WCS WCSParams
    bool readImagePara(const std::string& filename, int& nx, int& ny, std::vector<float>& data, WCSParams& wcs);

    // Update WCS keywords in a FITS file
    bool updatePara(const std::string& filename, const WCSParams& wcs);

    // Create a new FITS file, copy headers from file1, and write array
    bool writeImageCopyHDU(const std::string& templateFile, const std::string& filename, int nx, int ny, const std::vector<float>& data);

    // Write a standard 2D float image
    bool writeImage(const std::string& filename, int nx, int ny, const std::vector<float>& data);

    // Read multiple stamp sub-images from a large 2D FITS image
    bool readStamps(int np, int nstart, int n, int nsx, int nsy, std::vector<float>& stamps, int n1, int n2, const std::string& filename);

    // Write multiple stamp sub-images to a large 2D FITS image
    bool writeStamps(int np, int nstart, int n, int nsx, int nsy, const std::vector<float>& stamps, int n1, int n2, const std::string& filename);

    // Write stamps with options
    bool writeStamps2(int np, int n, int nsx, int nsy, const std::vector<float>& stamps, const std::vector<int>& opt, int val, int n1, int n2, const std::string& filename);

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
