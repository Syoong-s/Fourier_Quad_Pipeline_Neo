#include "FitsIO.hpp"
#include "OutputFile.hpp"
#include <fitsio.h>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <filesystem>

namespace FitsIO {

    // ==========================================
    // Function: Format one CFITSIO status and diagnostic stack
    // Method: Capture the primary status text and drain all queued messages into
    //         one reason string suitable for read diagnostics or fatal output errors.
    // ==========================================
    static std::string fitsErrorMessage(int status) {
        char errtext[FLEN_STATUS] = {};
        fits_get_errstatus(status, errtext);
        std::ostringstream message;
        message << "CFITSIO status " << status << ": " << errtext;
        char errmessage[FLEN_ERRMSG] = {};
        while (fits_read_errmsg(errmessage) != 0) {
            message << "; " << errmessage;
        }
        return message.str();
    }

    // ==========================================
    // Function: Print one non-fatal CFITSIO diagnostic
    // Method: Reuse the complete formatted status stack for input-side failures.
    // ==========================================
    void printError(int status) {
        if (status != 0) {
            std::cerr << fitsErrorMessage(status) << std::endl;
        }
    }

    // ==========================================
    // Function: Terminate after a CFITSIO output failure
    // Method: Convert the complete CFITSIO status stack into the shared MPI-wide
    //         fail-fast output diagnostic.
    // ==========================================
    [[noreturn]] static void failFitsOutput(const std::string& operation,
                                            const std::string& filename,
                                            int status) {
        MainIO::failOutput(operation, filename, fitsErrorMessage(status));
    }

    // ==========================================
    // Function: Mark image read failure
    // Method: Match F77 readimage/readimage_para by placing -99999 in the first pixel.
    // ==========================================
    static void markReadFailure(int nx, int ny, std::vector<float>& data) {
        size_t n = (nx > 0 && ny > 0) ? static_cast<size_t>(nx) * static_cast<size_t>(ny) : 1u;
        data.assign(n, 0.0f);
        data[0] = -99999.0f;
    }

    // ==========================================
    // Function: Close FITS file after a failed read
    // Method: Use a fresh status so CFITSIO does not ignore close after an earlier error.
    // ==========================================
    static void closeAfterFailure(fitsfile* fptr) {
        if (fptr == nullptr) return;
        int closeStatus = 0;
        fits_close_file(fptr, &closeStatus);
    }

    // ==========================================
    // Function: Compute the exact storage size of one stamp cube
    // Method: Reject non-positive axes and guard both size_t and CFITSIO
    //         LONGLONG multiplication before allocating or performing I/O.
    // ==========================================
    static bool stampCubeElementCount(int nx, int ny, int count,
                                      std::size_t& elementCount,
                                      LONGLONG& fitsElementCount) {
        elementCount = 0;
        fitsElementCount = 0;
        if (nx <= 0 || ny <= 0 || count <= 0) {
            return false;
        }

        const std::size_t sx = static_cast<std::size_t>(nx);
        const std::size_t sy = static_cast<std::size_t>(ny);
        const std::size_t sc = static_cast<std::size_t>(count);
        if (sx > std::numeric_limits<std::size_t>::max() / sy) {
            return false;
        }
        const std::size_t planeSize = sx * sy;
        if (planeSize > std::numeric_limits<std::size_t>::max() / sc) {
            return false;
        }
        elementCount = planeSize * sc;
        if (elementCount > static_cast<std::size_t>(
                std::numeric_limits<LONGLONG>::max())) {
            elementCount = 0;
            return false;
        }
        fitsElementCount = static_cast<LONGLONG>(elementCount);
        return true;
    }

    bool readCCDNUM(const std::string& filename, int& ccdNum) {
        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            printError(status);
            return false;
        }
        fits_read_key(fptr, TINT, "CCDNUM", &ccdNum, nullptr, &status);
        if (status != 0) {
            printError(status);
            closeAfterFailure(fptr);
            return false;
        }
        int closeStatus = 0;
        fits_close_file(fptr, &closeStatus);
        if (closeStatus != 0) {
            printError(closeStatus);
            return false;
        }
        return true;
    }

    bool readPara(const std::string& filename, int& nx, int& ny) {
        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            printError(status);
            return false;
        }
        long naxes[2] = {0, 0};
        int nfound = 0;
        fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read NAXIS keywords of: " << filename << std::endl;
            fits_close_file(fptr, &status);
            return false;
        }
        nx = static_cast<int>(naxes[0]);
        ny = static_cast<int>(naxes[1]);
        fits_close_file(fptr, &status);
        return (status == 0);
    }

    bool readImage(const std::string& filename, int& nx, int& ny, std::vector<float>& data) {
        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            std::cerr << "Error opening file: " << filename << std::endl;
            printError(status);
            return false;
        }
        long naxes[2] = {0, 0};
        int nfound = 0;
        fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read NAXIS keywords of: " << filename << std::endl;
            fits_close_file(fptr, &status);
            return false;
        }
        nx = static_cast<int>(naxes[0]);
        ny = static_cast<int>(naxes[1]);

        data.resize(nx * ny);
        long fpixel[2] = {1, 1};
        float nullval = 0.0f;
        int anynull = 0;
        fits_read_pix(fptr, TFLOAT, fpixel, nx * ny, &nullval, data.data(), &anynull, &status);
        fits_close_file(fptr, &status);
        if (status != 0) {
            printError(status);
            return false;
        }
        return true;
    }

    // ==========================================
    // Function: Read the first pixel of one two-dimensional FITS image
    // Method: Distinguish a missing path before opening, validate positive 2D axes, and read
    //         exactly one float while mapping every existing-file failure to ReadError.
    // ==========================================
    PixelReadStatus readFirstPixel(const std::string& filename, float& value) {
        std::error_code filesystemError;
        const bool exists = std::filesystem::exists(filename, filesystemError);
        if (filesystemError) {
            return PixelReadStatus::ReadError;
        }
        if (!exists) {
            return PixelReadStatus::Missing;
        }

        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            fits_clear_errmsg();
            return PixelReadStatus::ReadError;
        }

        int naxis = 0;
        long naxes[2] = {0, 0};
        fits_get_img_dim(fptr, &naxis, &status);
        if (status == 0 && naxis == 2) {
            fits_get_img_size(fptr, 2, naxes, &status);
        }
        if (status != 0 || naxis != 2 || naxes[0] <= 0 || naxes[1] <= 0) {
            closeAfterFailure(fptr);
            fits_clear_errmsg();
            return PixelReadStatus::ReadError;
        }

        long firstPixel[2] = {1, 1};
        float nullValue = std::numeric_limits<float>::quiet_NaN();
        int anyNull = 0;
        fits_read_pix(fptr, TFLOAT, firstPixel, 1, &nullValue, &value,
                      &anyNull, &status);
        if (status != 0) {
            closeAfterFailure(fptr);
            fits_clear_errmsg();
            return PixelReadStatus::ReadError;
        }

        int closeStatus = 0;
        fits_close_file(fptr, &closeStatus);
        if (closeStatus != 0) {
            fits_clear_errmsg();
            return PixelReadStatus::ReadError;
        }
        return PixelReadStatus::Ok;
    }

    // ==========================================
    // Function: Read image and WCS parameters
    // Method: Preserve F77 readimage_para sentinel semantics on read failure.
    // ==========================================
    bool readImagePara(const std::string& filename, int& nx, int& ny, std::vector<float>& data, WCSParams& wcs) {
        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            std::cerr << "Error opening file: " << filename << std::endl;
            printError(status);
            markReadFailure(nx, ny, data);
            return false;
        }

        int nfound = 0;
        long naxes[2] = {0, 0};
        fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read NAXIS keywords of: " << filename << std::endl;
            markReadFailure(nx, ny, data);
            closeAfterFailure(fptr);
            return false;
        }
        nx = static_cast<int>(naxes[0]);
        ny = static_cast<int>(naxes[1]);

        fits_read_keys_dbl(fptr, "CRPIX", 1, 2, wcs.crpix, &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read CRPIX keywords of: " << filename << std::endl;
            markReadFailure(nx, ny, data);
            closeAfterFailure(fptr);
            return false;
        }

        fits_read_keys_dbl(fptr, "CRVAL", 1, 2, wcs.crval, &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read CRVAL keywords of: " << filename << std::endl;
            markReadFailure(nx, ny, data);
            closeAfterFailure(fptr);
            return false;
        }

        fits_read_keys_dbl(fptr, "CD1_", 1, 2, wcs.cd[0], &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read CD1_ keywords of: " << filename << std::endl;
            markReadFailure(nx, ny, data);
            closeAfterFailure(fptr);
            return false;
        }

        fits_read_keys_dbl(fptr, "CD2_", 1, 2, wcs.cd[1], &nfound, &status);
        if (status != 0 || nfound != 2) {
            std::cerr << "Failed to read CD2_ keywords of: " << filename << std::endl;
            markReadFailure(nx, ny, data);
            closeAfterFailure(fptr);
            return false;
        }

        data.resize(nx * ny);
        long fpixel[2] = {1, 1};
        float nullval = 0.0f;
        int anynull = 0;
        fits_read_pix(fptr, TFLOAT, fpixel, nx * ny, &nullval, data.data(), &anynull, &status);
        fits_close_file(fptr, &status);
        if (status != 0) {
            printError(status);
            markReadFailure(nx, ny, data);
            return false;
        }
        return true;
    }

    bool updatePara(const std::string& filename, const WCSParams& wcs) {
        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READWRITE, &status);
        if (status != 0) {
            printError(status);
            return false;
        }

        fits_update_key(fptr, TDOUBLE, "CRPIX1", const_cast<double*>(&wcs.crpix[0]), "replaced", &status);
        fits_update_key(fptr, TDOUBLE, "CRPIX2", const_cast<double*>(&wcs.crpix[1]), "replaced", &status);
        fits_update_key(fptr, TDOUBLE, "CD1_1", const_cast<double*>(&wcs.cd[0][0]), "replaced", &status);
        fits_update_key(fptr, TDOUBLE, "CD1_2", const_cast<double*>(&wcs.cd[0][1]), "replaced", &status);
        fits_update_key(fptr, TDOUBLE, "CD2_1", const_cast<double*>(&wcs.cd[1][0]), "replaced", &status);
        fits_update_key(fptr, TDOUBLE, "CD2_2", const_cast<double*>(&wcs.cd[1][1]), "replaced", &status);

        fits_close_file(fptr, &status);
        if (status != 0) {
            printError(status);
            return false;
        }
        return true;
    }

    // ==========================================
    // Function: Write an image while preserving the template primary HDU
    // Method: Treat template-open failure as an input error, but terminate the
    //         complete MPI program for output creation, write, or close failure.
    // ==========================================
    bool writeImageCopyHDU(const std::string& templateFile, const std::string& filename, int nx, int ny, const std::vector<float>& data) {
        fitsfile* infptr = nullptr;
        fitsfile* outfptr = nullptr;
        int status = 0;

        fits_open_file(&infptr, templateFile.c_str(), READONLY, &status);
        if (status != 0) {
            printError(status);
            return false;
        }

        const std::string create_name = "!" + filename;
        fits_create_file(&outfptr, create_name.c_str(), &status);
        if (status != 0) {
            int close_status = 0;
            fits_close_file(infptr, &close_status);
            failFitsOutput("create FITS output", filename, status);
        }

        // Copy primary HDU
        fits_copy_hdu(infptr, outfptr, 0, &status);
        
        // Modify bitpix to float
        int bitpix = -32;
        fits_update_key(outfptr, TINT, "BITPIX", &bitpix, nullptr, &status);

        long fpixel[2] = {1, 1};
        fits_write_pix(outfptr, TFLOAT, fpixel, nx * ny, const_cast<float*>(data.data()), &status);
        if (status != 0) {
            failFitsOutput("write FITS output", filename, status);
        }

        int input_close_status = 0;
        fits_close_file(infptr, &input_close_status);
        if (input_close_status != 0) {
            printError(input_close_status);
            return false;
        }

        int output_close_status = 0;
        fits_close_file(outfptr, &output_close_status);
        if (output_close_status != 0) {
            failFitsOutput("close FITS output", filename, output_close_status);
        }
        return true;
    }

    // ==========================================
    // Function: Write one standard two-dimensional float FITS image
    // Method: Use CFITSIO overwrite syntax and terminate the complete MPI
    //         program for creation, header, pixel, or close failure.
    // ==========================================
    bool writeImage(const std::string& filename, int nx, int ny, const std::vector<float>& data) {
        fitsfile* fptr = nullptr;
        int status = 0;
        const std::string create_name = "!" + filename;
        fits_create_file(&fptr, create_name.c_str(), &status);
        if (status != 0) {
            failFitsOutput("create FITS output", filename, status);
        }

        int bitpix = -32;
        int naxis = 2;
        long naxes[2] = {nx, ny};
        fits_write_imghdr(fptr, bitpix, naxis, naxes, &status);

        long fpixel[2] = {1, 1};
        fits_write_pix(fptr, TFLOAT, fpixel, nx * ny, const_cast<float*>(data.data()), &status);
        if (status != 0) {
            failFitsOutput("write FITS output", filename, status);
        }

        int close_status = 0;
        fits_close_file(fptr, &close_status);
        if (close_status != 0) {
            failFitsOutput("close FITS output", filename, close_status);
        }
        return true;
    }

    // ==========================================
    // Function: Read one contiguous three-dimensional stamp cube
    // Method: Derive all axes from the FITS header, require positive
    //         (x, y, stamp) dimensions, and read the row-major planes once.
    // ==========================================
    bool readStampCube(const std::string& filename, StampCubeShape& shape,
                       std::vector<float>& stamps) {
        shape = {};
        stamps.clear();

        fitsfile* fptr = nullptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        if (status != 0) {
            std::cerr << "Error opening stamp cube: " << filename << std::endl;
            printError(status);
            return false;
        }

        int naxis = 0;
        fits_get_img_dim(fptr, &naxis, &status);
        if (status != 0 || naxis != 3) {
            if (status != 0) {
                printError(status);
            } else {
                std::cerr << "Stamp cube must have NAXIS=3: " << filename
                          << " (actual NAXIS=" << naxis << ")" << std::endl;
            }
            closeAfterFailure(fptr);
            return false;
        }

        long naxes[3] = {0, 0, 0};
        fits_get_img_size(fptr, 3, naxes, &status);
        if (status != 0 || naxes[0] <= 0 || naxes[1] <= 0 || naxes[2] <= 0
            || naxes[0] > std::numeric_limits<int>::max()
            || naxes[1] > std::numeric_limits<int>::max()
            || naxes[2] > std::numeric_limits<int>::max()) {
            if (status != 0) {
                printError(status);
            } else {
                std::cerr << "Stamp cube axes must be positive int values: "
                          << filename << std::endl;
            }
            closeAfterFailure(fptr);
            return false;
        }

        shape.nx = static_cast<int>(naxes[0]);
        shape.ny = static_cast<int>(naxes[1]);
        shape.count = static_cast<int>(naxes[2]);
        std::size_t elementCount = 0;
        LONGLONG fitsElementCount = 0;
        if (!stampCubeElementCount(shape.nx, shape.ny, shape.count,
                                   elementCount, fitsElementCount)) {
            std::cerr << "Stamp cube element count overflows supported storage: "
                      << filename << std::endl;
            closeAfterFailure(fptr);
            shape = {};
            return false;
        }

        stamps.resize(elementCount);
        long fpixel[3] = {1, 1, 1};
        float nullval = 0.0f;
        int anynull = 0;
        fits_read_pix(fptr, TFLOAT, fpixel, fitsElementCount, &nullval,
                      stamps.data(), &anynull, &status);
        if (status != 0) {
            printError(status);
            closeAfterFailure(fptr);
            shape = {};
            stamps.clear();
            return false;
        }

        int closeStatus = 0;
        fits_close_file(fptr, &closeStatus);
        if (closeStatus != 0) {
            printError(closeStatus);
            shape = {};
            stamps.clear();
            return false;
        }
        return true;
    }

    // ==========================================
    // Function: Write one contiguous three-dimensional stamp cube
    // Method: Validate the exact collection size, create axes
    //         (x, y, stamp), add format metadata, and write all pixels once.
    // ==========================================
    bool writeStampCube(const std::string& filename, int nx, int ny, int count,
                        const std::vector<float>& stamps) {
        std::size_t elementCount = 0;
        LONGLONG fitsElementCount = 0;
        if (!stampCubeElementCount(nx, ny, count, elementCount,
                                   fitsElementCount)) {
            MainIO::failOutput(
                "validate stamp cube", filename,
                "nx, ny, and count must be positive and their product must fit supported storage");
        }
        if (stamps.size() != elementCount) {
            MainIO::failOutput(
                "validate stamp cube", filename,
                "stamp collection size does not equal count * ny * nx");
        }

        fitsfile* fptr = nullptr;
        int status = 0;
        const std::string createName = "!" + filename;
        fits_create_file(&fptr, createName.c_str(), &status);
        if (status != 0) {
            failFitsOutput("create stamp cube", filename, status);
        }

        int bitpix = FLOAT_IMG;
        int naxis = 3;
        long naxes[3] = {nx, ny, count};
        fits_write_imghdr(fptr, bitpix, naxis, naxes, &status);

        char format[] = "STAMP_CUBE";
        char order[] = "X,Y,STAMP";
        fits_update_key(fptr, TSTRING, "FQFMT", format,
                        "Fourier_Quad stamp collection format", &status);
        fits_update_key(fptr, TSTRING, "FQORDER", order,
                        "FITS axis semantic order", &status);

        long fpixel[3] = {1, 1, 1};
        fits_write_pix(fptr, TFLOAT, fpixel, fitsElementCount,
                       const_cast<float*>(stamps.data()), &status);
        if (status != 0) {
            failFitsOutput("write stamp cube", filename, status);
        }

        int closeStatus = 0;
        fits_close_file(fptr, &closeStatus);
        if (closeStatus != 0) {
            failFitsOutput("close stamp cube", filename, closeStatus);
        }
        return true;
    }

    // ==========================================
    // Function: Compact selected stamp planes into one output cube
    // Method: Preserve input order, copy only matching planes, and delegate
    //         the single contiguous FITS write to writeStampCube.
    // ==========================================
    bool writeSelectedStampCube(const std::string& filename, int nx, int ny,
                                const std::vector<float>& stamps,
                                int inputCount,
                                const std::vector<int>& selection,
                                int selectedValue) {
        std::size_t inputElementCount = 0;
        LONGLONG ignoredFitsCount = 0;
        if (!stampCubeElementCount(nx, ny, inputCount, inputElementCount,
                                   ignoredFitsCount)) {
            MainIO::failOutput(
                "validate selected stamp cube", filename,
                "nx, ny, and inputCount must be positive and fit supported storage");
        }
        if (stamps.size() != inputElementCount
            || selection.size() != static_cast<std::size_t>(inputCount)) {
            MainIO::failOutput(
                "validate selected stamp cube", filename,
                "input stamp or selection size does not match inputCount");
        }

        const int selectedCount = static_cast<int>(std::count(
            selection.begin(), selection.end(), selectedValue));
        if (selectedCount <= 0) {
            MainIO::failOutput(
                "validate selected stamp cube", filename,
                "selection contains no output stamp");
        }

        const std::size_t planeSize = static_cast<std::size_t>(nx) * ny;
        std::vector<float> selected;
        selected.reserve(static_cast<std::size_t>(selectedCount) * planeSize);
        for (int stamp = 0; stamp < inputCount; ++stamp) {
            if (selection[stamp] != selectedValue) {
                continue;
            }
            const std::size_t offset = static_cast<std::size_t>(stamp) * planeSize;
            selected.insert(selected.end(), stamps.begin() + offset,
                            stamps.begin() + offset + planeSize);
        }
        return writeStampCube(filename, nx, ny, selectedCount, selected);
    }

    // ==========================================
    // Function: Initialize one serial FITS writer
    // Method: Start with no open output and no retained CFITSIO status.
    // ==========================================
    FitsSerialWriter::FitsSerialWriter() : fptr(nullptr), status(0) {}

    // ==========================================
    // Function: Finalize one serial FITS writer
    // Method: Route scope-driven close through the checked output-close path.
    // ==========================================
    FitsSerialWriter::~FitsSerialWriter() {
        close();
    }

    // ==========================================
    // Function: Create one serial multi-HDU FITS output
    // Method: Use CFITSIO overwrite syntax and terminate all MPI ranks when
    //         output creation fails.
    // ==========================================
    bool FitsSerialWriter::init(const std::string& filename) {
        close();
        outputFilename = filename;
        status = 0;
        fitsfile* f = nullptr;
        const std::string create_name = "!" + filename;
        fits_create_file(&f, create_name.c_str(), &status);
        fptr = static_cast<void*>(f);
        if (status != 0) {
            failFitsOutput("create serial FITS output", outputFilename, status);
        }
        return true;
    }

    // ==========================================
    // Function: Append one image HDU to a serial FITS output
    // Method: Validate writer initialization and terminate all MPI ranks on
    //         any CFITSIO image-header or pixel-write failure.
    // ==========================================
    bool FitsSerialWriter::writeStamp(int nx, int ny, const std::vector<float>& data, bool newHdu) {
        fitsfile* f = static_cast<fitsfile*>(fptr);
        if (!f) {
            MainIO::failOutput(
                "write serial FITS output", outputFilename, "serial writer is not initialized");
        }

        if (newHdu) {
            fits_create_img(f, -32, 2, nullptr, &status); // create extension HDU
        }

        int bitpix = -32;
        int naxis = 2;
        long naxes[2] = {nx, ny};
        fits_write_imghdr(f, bitpix, naxis, naxes, &status);

        long fpixel[2] = {1, 1};
        fits_write_pix(f, TFLOAT, fpixel, nx * ny, const_cast<float*>(data.data()), &status);

        if (status != 0) {
            failFitsOutput("write serial FITS output", outputFilename, status);
        }
        return true;
    }

    // ==========================================
    // Function: Append one integer keyword to a serial FITS output
    // Method: Validate writer initialization and terminate all MPI ranks on
    //         any CFITSIO keyword-write failure.
    // ==========================================
    bool FitsSerialWriter::writeKey(const std::string& keyName, int val, const std::string& comment) {
        fitsfile* f = static_cast<fitsfile*>(fptr);
        if (!f) {
            MainIO::failOutput(
                "write serial FITS keyword", outputFilename, "serial writer is not initialized");
        }

        fits_write_key(f, TINT, keyName.c_str(), &val, comment.c_str(), &status);
        if (status != 0) {
            failFitsOutput("write serial FITS keyword", outputFilename, status);
        }
        return true;
    }

    // ==========================================
    // Function: Close one serial FITS output
    // Method: Use a fresh close status and terminate all MPI ranks when buffered
    //         output cannot be finalized.
    // ==========================================
    void FitsSerialWriter::close() {
        fitsfile* f = static_cast<fitsfile*>(fptr);
        if (f) {
            int close_status = 0;
            fits_close_file(f, &close_status);
            fptr = nullptr;
            if (close_status != 0) {
                failFitsOutput("close serial FITS output", outputFilename, close_status);
            }
            status = 0;
        }
    }
}
