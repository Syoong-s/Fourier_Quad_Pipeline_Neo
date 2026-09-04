#include "process_init/FitsExtractor.hpp"
#include "Initialize.hpp"

#include <fitsio.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fqinit {
namespace {

constexpr std::size_t kCopyBufferBytes = 16U * 1024U * 1024U;

struct PlannedImage {
    int hdu_number = 0;
    int expected_ccdnum = 0;
    bool check_ccdnum = false;
    std::filesystem::path final_path;
    std::filesystem::path staged_path;
};

// ==========================================
// Function: Convert a CFITSIO status stack into one diagnostic string
// Method: Read the primary status text and drain all queued CFITSIO messages.
// ==========================================
std::string fitsDiagnostic(int status) {
    char status_text[FLEN_STATUS] = {};
    fits_get_errstatus(status, status_text);
    std::ostringstream message;
    message << status_text;
    char detail[FLEN_ERRMSG] = {};
    while (fits_read_errmsg(detail) != 0) {
        message << "; " << detail;
    }
    return message.str();
}

// ==========================================
// Function: Close a FITS handle even after an earlier CFITSIO failure
// Method: Use a fresh close status because CFITSIO ignores calls after error.
// ==========================================
void closeFits(fitsfile*& file) {
    if (file == nullptr) {
        return;
    }
    int close_status = 0;
    fits_close_file(file, &close_status);
    file = nullptr;
}

// ==========================================
// Function: Replace all non-overlapping occurrences in a filename stem
// Method: Advance past each replacement so the configured DQ stem mapping is exact.
// ==========================================
void replaceAll(std::string& value, const std::string& from, const std::string& to) {
    if (from.empty()) {
        return;
    }
    std::size_t position = 0;
    while ((position = value.find(from, position)) != std::string::npos) {
        value.replace(position, from.size(), to);
        position += to.size();
    }
}

// ==========================================
// Function: Map FITS BITPIX codes to CFITSIO memory datatypes
// Method: Preserve signed, unsigned, floating, and 64-bit image representations.
// ==========================================
std::pair<int, std::size_t> imageDatatype(int bitpix) {
    switch (bitpix) {
        case BYTE_IMG:
            return {TBYTE, 1U};
        case SBYTE_IMG:
            return {TSBYTE, 1U};
        case SHORT_IMG:
            return {TSHORT, 2U};
        case USHORT_IMG:
            return {TUSHORT, 2U};
        case LONG_IMG:
            return {TINT, 4U};
        case ULONG_IMG:
            return {TUINT, 4U};
        case LONGLONG_IMG:
            return {TLONGLONG, 8U};
        case ULONGLONG_IMG:
            return {TULONGLONG, 8U};
        case FLOAT_IMG:
            return {TFLOAT, 4U};
        case DOUBLE_IMG:
            return {TDOUBLE, 8U};
        default:
            throw std::runtime_error("unsupported FITS BITPIX code: " + std::to_string(bitpix));
    }
}

// ==========================================
// Function: Copy the current logical FITS image into one uncompressed FITS file
// Method: Follow the CFITSIO imcopy algorithm: create a normal image, copy only
//         non-structural/non-compression cards, and stream raw pixels in chunks.
// ==========================================
void copyCurrentImage(fitsfile* input, const std::filesystem::path& output_path) {
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_get_img_param(input, 9, &bitpix, &naxis, naxes, &status);
    if (status != 0) {
        throw std::runtime_error("cannot read logical image parameters: " + fitsDiagnostic(status));
    }
    if (naxis != 2) {
        throw std::runtime_error("initializer only extracts two-dimensional image HDUs");
    }

    fitsfile* output = nullptr;
    const std::string create_name = "!" + output_path.string();
    fits_create_file(&output, create_name.c_str(), &status);
    fits_create_img(output, bitpix, naxis, naxes, &status);
    if (status != 0) {
        closeFits(output);
        throw std::runtime_error("cannot create uncompressed output: " + fitsDiagnostic(status));
    }

    int key_count = 0;
    fits_get_hdrspace(input, &key_count, nullptr, &status);
    for (int key_index = 1; key_index <= key_count && status == 0; ++key_index) {
        char card[FLEN_CARD] = {};
        fits_read_record(input, key_index, card, &status);
        if (status != 0) {
            break;
        }
        const int key_class = fits_get_keyclass(card);
        if (key_class > TYP_CMPRS_KEY && key_class != TYP_CKSUM_KEY) {
            fits_write_record(output, card, &status);
        }
    }
    if (status != 0) {
        closeFits(output);
        throw std::runtime_error("cannot copy FITS header cards: " + fitsDiagnostic(status));
    }

    const auto [datatype, bytes_per_pixel] = imageDatatype(bitpix);
    LONGLONG total_pixels = 1;
    for (int axis = 0; axis < naxis; ++axis) {
        if (naxes[axis] <= 0) {
            closeFits(output);
            throw std::runtime_error("image has a non-positive axis length");
        }
        total_pixels *= static_cast<LONGLONG>(naxes[axis]);
    }

    const LONGLONG chunk_pixels = std::max<LONGLONG>(
        1, static_cast<LONGLONG>(kCopyBufferBytes / bytes_per_pixel));
    const std::size_t buffer_bytes = static_cast<std::size_t>(chunk_pixels) * bytes_per_pixel;
    std::vector<double> aligned_buffer((buffer_bytes + sizeof(double) - 1U) / sizeof(double));
    double null_value = 0.0;
    int any_null = 0;
    fits_set_bscale(input, 1.0, 0.0, &status);
    fits_set_bscale(output, 1.0, 0.0, &status);

    LONGLONG first_pixel = 1;
    LONGLONG remaining = total_pixels;
    while (remaining > 0 && status == 0) {
        const LONGLONG count = std::min(remaining, chunk_pixels);
        fits_read_img(input, datatype, first_pixel, count, &null_value,
                      aligned_buffer.data(), &any_null, &status);
        fits_write_img(output, datatype, first_pixel, count, aligned_buffer.data(), &status);
        first_pixel += count;
        remaining -= count;
    }

    int output_close_status = 0;
    fits_close_file(output, &output_close_status);
    output = nullptr;
    if (status == 0 && output_close_status != 0) {
        status = output_close_status;
    }
    if (status != 0) {
        throw std::runtime_error("cannot stream FITS pixels: " + fitsDiagnostic(status));
    }
}

// ==========================================
// Function: Verify a committed or resumed chip image
// Method: Require a two-dimensional uncompressed primary image and, for DQ,
//         require the configured DQ chip identifier used in its output filename.
// ==========================================
bool validateOutput(const PlannedImage& plan, std::string& error) {
    fitsfile* file = nullptr;
    int status = 0;
    fits_open_file(&file, plan.final_path.string().c_str(), READONLY, &status);
    if (status != 0) {
        error = "cannot open existing output " + plan.final_path.string() + ": "
                + fitsDiagnostic(status);
        closeFits(file);
        return false;
    }

    int naxis = 0;
    fits_get_img_dim(file, &naxis, &status);
    const int compressed = fits_is_compressed_image(file, &status);
    if (status == 0 && (naxis != 2 || compressed != 0)) {
        error = "existing output is not an uncompressed two-dimensional FITS image: "
                + plan.final_path.string();
        closeFits(file);
        return false;
    }
    if (status == 0 && plan.check_ccdnum) {
        int ccdnum = 0;
        fits_read_key(file, TINT, Initialize::CCDNUM_KEYWORD,
                      &ccdnum, nullptr, &status);
        if (status == 0 && ccdnum != plan.expected_ccdnum) {
            error = "existing DQ output has the wrong "
                    + std::string(Initialize::CCDNUM_KEYWORD) + ": "
                    + plan.final_path.string();
            closeFits(file);
            return false;
        }
    }
    if (status != 0) {
        error = "cannot validate existing output " + plan.final_path.string() + ": "
                + fitsDiagnostic(status);
        closeFits(file);
        return false;
    }
    closeFits(file);
    return true;
}

// ==========================================
// Function: Discover the output name of every extractable HDU in one archive
// Method: Number science images by two-dimensional HDU occurrence and DQ images
//         by the configured chip-keyword header without changing either convention.
// ==========================================
std::vector<PlannedImage> planImages(fitsfile* input,
                                     const std::filesystem::path& source,
                                     ProductKind kind,
                                     const std::filesystem::path& final_directory,
                                     const std::filesystem::path& staging_directory) {
    int status = 0;
    int hdu_count = 0;
    fits_get_num_hdus(input, &hdu_count, &status);
    if (status != 0) {
        throw std::runtime_error("cannot count FITS HDUs: " + fitsDiagnostic(status));
    }

    const std::string output_stem = kind == ProductKind::Science
                                        ? archiveStem(source)
                                        : dqOutputStem(source);
    int science_index = 0;
    std::set<std::string> unique_outputs;
    std::vector<PlannedImage> plans;

    for (int hdu_number = 1; hdu_number <= hdu_count; ++hdu_number) {
        int hdu_type = 0;
        int hdu_status = 0;
        fits_movabs_hdu(input, hdu_number, &hdu_type, &hdu_status);
        if (hdu_status != 0) {
            fits_clear_errmsg();
            continue;
        }
        if (hdu_type != IMAGE_HDU) {
            continue;
        }
        int naxis = 0;
        int dim_status = 0;
        fits_get_img_dim(input, &naxis, &dim_status);
        if (dim_status != 0) {
            fits_clear_errmsg();
            continue;
        }
        if (naxis != 2) {
            continue;
        }

        int output_number = 0;
        bool check_ccdnum = false;
        if (kind == ProductKind::Science) {
            output_number = ++science_index;
        } else {
            int key_status = 0;
            fits_read_key(input, TINT, Initialize::CCDNUM_KEYWORD,
                          &output_number, nullptr, &key_status);
            if (key_status != 0) {
                fits_clear_errmsg();
                continue;
            }
            check_ccdnum = true;
        }

        const std::string filename = output_stem + "_" + std::to_string(output_number) + ".fits";
        if (!unique_outputs.insert(filename).second) {
            continue;
        }
        plans.push_back({hdu_number, output_number, check_ccdnum,
                         final_directory / filename, staging_directory / filename});
    }
    if (plans.empty()) {
        throw std::runtime_error("archive contains no extractable two-dimensional image HDUs");
    }
    return plans;
}

// ==========================================
// Function: Commit all staged images for one source archive
// Method: Use same-filesystem POSIX rename so each published chip file appears
//         atomically.
// ==========================================
void commitImages(const std::vector<PlannedImage>& plans) {
    for (const PlannedImage& plan : plans) {
        if (std::rename(plan.staged_path.c_str(), plan.final_path.c_str()) != 0) {
            throw std::runtime_error("cannot commit " + plan.final_path.string() + ": "
                                     + std::strerror(errno));
        }
    }
}

}  // namespace

// ==========================================
// Function: Return the exposure stem used by science chip output names
// Method: Remove the exact configured archive suffix from the source basename.
// ==========================================
std::string archiveStem(const std::filesystem::path& source) {
    const std::string filename = source.filename().string();
    const std::size_t suffix_length = std::strlen(Initialize::ARCHIVE_SUFFIX);
    if (filename.size() <= suffix_length
        || filename.compare(filename.size() - suffix_length,
                            suffix_length, Initialize::ARCHIVE_SUFFIX) != 0) {
        throw std::runtime_error(
            "archive does not end in " + std::string(Initialize::ARCHIVE_SUFFIX)
            + ": " + source.string());
    }
    return filename.substr(0, filename.size() - suffix_length);
}

// ==========================================
// Function: Return the DQ exposure stem expected by the pipeline
// Method: Remove the archive suffix and apply the configured DQ stem replacement.
// ==========================================
std::string dqOutputStem(const std::filesystem::path& source) {
    std::string stem = archiveStem(source);
    replaceAll(stem, Initialize::DQ_STEM_REPLACE_FROM,
               Initialize::DQ_STEM_REPLACE_TO);
    return stem;
}

// ==========================================
// Function: Extract one multi-HDU FITS/FZ archive into pipeline chip images
// Method: Read the source archive in place with CFITSIO, stage uncompressed
//         two-dimensional HDUs, then commit the completed output set.
//         Skip individual HDUs that fail extraction and continue with the rest.
// ==========================================
ExtractionResult extractArchive(const std::filesystem::path& source,
                                ProductKind kind,
                                const std::filesystem::path& final_directory,
                                const std::filesystem::path& staging_directory,
                                ExistingPolicy policy) {
    ExtractionResult result;
    fitsfile* input = nullptr;
    try {
        int status = 0;
        fits_open_file(&input, source.string().c_str(), READONLY, &status);
        if (status != 0) {
            throw std::runtime_error("cannot open archive: " + fitsDiagnostic(status));
        }

        const std::vector<PlannedImage> plans = planImages(
            input, source, kind, final_directory, staging_directory);

        // Partition plans into existing (resumable) and to-extract.
        std::vector<PlannedImage> to_extract;
        bool any_existing = false;
        for (const PlannedImage& plan : plans) {
            if (std::filesystem::exists(plan.final_path)) {
                any_existing = true;
            }
        }
        if (any_existing && policy == ExistingPolicy::Fail) {
            result.error = "output already exists: "
                           + plans.front().final_path.parent_path().string();
            closeFits(input);
            return result;
        }
        for (const PlannedImage& plan : plans) {
            if (std::filesystem::exists(plan.final_path) && policy == ExistingPolicy::Resume) {
                std::string validation_error;
                if (validateOutput(plan, validation_error)) {
                    result.output_paths.push_back(plan.final_path);
                    continue;
                }
            }
            to_extract.push_back(plan);
        }
        if (to_extract.empty()) {
            result.success = true;
            result.resumed = true;
            closeFits(input);
            return result;
        }

        std::filesystem::create_directories(staging_directory);
        std::filesystem::create_directories(final_directory);
        std::vector<PlannedImage> successful_plans;
        int skipped = 0;
        for (const PlannedImage& plan : to_extract) {
            try {
                int hdu_type = 0;
                int hdu_status = 0;
                fits_movabs_hdu(input, plan.hdu_number, &hdu_type, &hdu_status);
                if (hdu_status != 0 || hdu_type != IMAGE_HDU) {
                    throw std::runtime_error("cannot revisit planned FITS image HDU: "
                                             + fitsDiagnostic(hdu_status));
                }
                copyCurrentImage(input, plan.staged_path);
                successful_plans.push_back(plan);
            } catch (const std::exception& exception) {
                ++skipped;
                fits_clear_errmsg();
                std::error_code remove_error;
                std::filesystem::remove(plan.staged_path, remove_error);
            }
        }
        closeFits(input);

        if (successful_plans.empty()) {
            result.error = "all planned HDUs failed extraction";
            std::error_code cleanup_error;
            std::filesystem::remove_all(staging_directory, cleanup_error);
            return result;
        }

        commitImages(successful_plans);
        for (const PlannedImage& plan : successful_plans) {
            result.output_paths.push_back(plan.final_path);
        }
        result.success = true;
        result.skipped_hdus = skipped;
        std::error_code cleanup_error;
        std::filesystem::remove_all(staging_directory, cleanup_error);
        return result;
    } catch (const std::exception& exception) {
        closeFits(input);
        result.error = exception.what();
        std::error_code cleanup_error;
        std::filesystem::remove_all(staging_directory, cleanup_error);
        return result;
    }
}

}  // namespace fqinit
