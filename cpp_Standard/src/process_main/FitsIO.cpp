#include "FitsIO.hpp"
#include "OutputFile.hpp"
#include <fitsio.h>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <sstream>

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

    bool readStamps(int np, int nstart, int n, int nsx, int nsy, std::vector<float>& stamps, int n1, int n2, const std::string& filename) {
        int lnx = 0, lny = 0;
        std::vector<float> largeStamp;
        if (!readImage(filename, lnx, lny, largeStamp)) {
            std::cerr << "Error opening large file: " << filename << std::endl;
            return false;
        }

        if (lnx != n1 || lny != n2) {
            std::cerr << "Warning: FITS size (" << lnx << "x" << lny << ") does not match requested (" << n1 << "x" << n2 << ")." << std::endl;
        }

        stamps.resize(static_cast<size_t>(np) * nsx * nsy, 0.0f);

        int offx = 0;
        int offy = 0;

        // Fortran nstart is 1-based, we map to 0-based k
        for (int k = nstart - 1; k < n; ++k) {
            if (offy + nsy > lny) {
                std::cerr << "large_stamp is too small for reading!" << std::endl;
                return false;
            }
            for (int y = 0; y < nsy; ++y) {
                for (int x = 0; x < nsx; ++x) {
                    size_t stampIdx = static_cast<size_t>(k) * nsx * nsy + y * nsx + x;
                    size_t largeIdx = static_cast<size_t>(y + offy) * lnx + (x + offx);
                    stamps[stampIdx] = largeStamp[largeIdx];
                }
            }
            offx += nsx;
            if (offx + nsx > lnx) {
                offx = 0;
                offy += nsy;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Pack and write a collection of image stamps.
    // Method: Match F77 write_stamps packing and terminate on fixed-buffer overflow.
    // ==========================================
    bool writeStamps(int np, int nstart, int n, int nsx, int nsy, const std::vector<float>& stamps, int n1, int n2, const std::string& filename) {
        constexpr int f77Nmax = 7000;
        if (n2 > f77Nmax) {
            MainIO::failOutput(
                "pack FITS stamps", filename,
                "requested output height exceeds the retained 7000-row stamp buffer contract");
        }

        std::vector<float> largeStamp(static_cast<size_t>(n1) * n2, 0.0f);

        int offx = 0;
        int offy = 0;

        for (int k = nstart - 1; k < n; ++k) {
            if (offy + nsy > n2) {
                MainIO::failOutput(
                    "pack FITS stamps", filename,
                    "stamp collection exceeds the requested FITS image dimensions");
            }
            for (int y = 0; y < nsy; ++y) {
                for (int x = 0; x < nsx; ++x) {
                    size_t stampIdx = static_cast<size_t>(k) * nsx * nsy + y * nsx + x;
                    size_t largeIdx = static_cast<size_t>(y + offy) * n1 + (x + offx);
                    largeStamp[largeIdx] = stamps[stampIdx];
                }
            }
            offx += nsx;
            if (offx + nsx > n1) {
                offx = 0;
                offy += nsy;
            }
        }

        return writeImage(filename, n1, n2, largeStamp);
    }

    // ==========================================
    // Function: Pack selected image stamps into one checked FITS image
    // Method: Abort the complete MPI program on packing overflow and delegate
    //         creation/write failure handling to writeImage.
    // ==========================================
    bool writeStamps2(int np, int n, int nsx, int nsy, const std::vector<float>& stamps, const std::vector<int>& opt, int val, int n1, int n2, const std::string& filename) {
        std::vector<float> largeStamp(static_cast<size_t>(n1) * n2, 0.0f);

        int offx = 0;
        int offy = 0;

        for (int k = 0; k < n; ++k) {
            if (opt[k] != val) continue;
            if (offy + nsy > n2) {
                MainIO::failOutput(
                    "pack selected FITS stamps", filename,
                    "selected stamp collection exceeds the requested FITS image dimensions");
            }
            for (int y = 0; y < nsy; ++y) {
                for (int x = 0; x < nsx; ++x) {
                    size_t stampIdx = static_cast<size_t>(k) * nsx * nsy + y * nsx + x;
                    size_t largeIdx = static_cast<size_t>(y + offy) * n1 + (x + offx);
                    largeStamp[largeIdx] = stamps[stampIdx];
                }
            }
            offx += nsx;
            if (offx + nsx > n1) {
                offx = 0;
                offy += nsy;
            }
        }

        return writeImage(filename, n1, n2, largeStamp);
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
