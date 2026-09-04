#include "process_main/CatalogCombiner.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/OutputFile.hpp"
#include "process_main/MPIFailure.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/ExposureInfo.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <system_error>
#include <cstddef>


namespace CatalogCombiner {

namespace {

// ==========================================
// Function: Trim trailing catalog whitespace
// Method: Remove spaces and line-ending characters without changing prefixes.
// ==========================================
std::string trimRight(std::string str) {
    while (!str.empty() && (str.back() == ' ' || str.back() == '\r'
                            || str.back() == '\n' || str.back() == '\t')) {
        str.pop_back();
    }
    return str;
}

// ==========================================
// Function: Count every physical line in one Stage-9 input catalog
// Method: Use a fresh getline stream and distinguish clean EOF from an I/O failure.
// ==========================================
std::size_t countCatalogLines(const std::string& filename,
                              const std::string& role) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        MPIFailure::abortWorld(
            "open catalog for row-count preflight", role + "=" + filename);
    }

    std::size_t line_count = 0;
    std::string line;
    while (std::getline(input, line)) {
        ++line_count;
    }
    if (input.bad()) {
        MPIFailure::abortWorld(
            "read catalog for row-count preflight", role + "=" + filename);
    }
    return line_count;
}

// ==========================================
// Function: Determine the fixed number of paired Stage-9 data rows
// Method: Count shear then orig, retry both with fresh streams after one mismatch,
//         and preserve a one-line shear catalog as the header-only sentinel.
// ==========================================
std::size_t determinePairedDataRows(const std::string& filename_shear,
                                    const std::string& filename_orig,
                                    const std::string& prefix) {
    const std::size_t shear_1 = countCatalogLines(filename_shear, "shear");
    if (shear_1 == 0) {
        MPIFailure::abortWorld(
            "preflight Stage 7 shear catalog",
            "shear catalog contains no header prefix=" + prefix
                + " shear=" + filename_shear);
    }
    if (shear_1 == 1) {
        return 0;
    }

    const std::size_t orig_1 = countCatalogLines(filename_orig, "orig");
    if (shear_1 == orig_1) {
        return shear_1 - 1;
    }

    const std::size_t shear_2 = countCatalogLines(filename_shear, "shear");
    if (shear_2 == 0) {
        MPIFailure::abortWorld(
            "preflight Stage 7 shear catalog",
            "shear catalog contains no header prefix=" + prefix
                + " shear=" + filename_shear);
    }
    if (shear_2 == 1) {
        return 0;
    }

    const std::size_t orig_2 = countCatalogLines(filename_orig, "orig");
    if (shear_2 != orig_2) {
        std::ostringstream detail;
        detail << "prefix=" << prefix
               << " attempt1_shear_lines=" << shear_1
               << " attempt1_orig_lines=" << orig_1
               << " attempt2_shear_lines=" << shear_2
               << " attempt2_orig_lines=" << orig_2
               << " shear=" << filename_shear
               << " orig=" << filename_orig;
        MPIFailure::abortWorld(
            "combine catalog row-count preflight", detail.str());
    }

    std::cout << "CATALOG_ROWCOUNT_RECOVERED"
              << " prefix=" << prefix
              << " attempt1_shear_lines=" << shear_1
              << " attempt1_orig_lines=" << orig_1
              << " attempt2_shear_lines=" << shear_2
              << " attempt2_orig_lines=" << orig_2 << std::endl;
    return shear_2 - 1;
}

// ==========================================
// Function: Parse one Stage-7 shear row
// Method: Reuse the caller-owned live-schema buffer and reject incomplete rows.
// ==========================================
bool parseShearRow(const std::string& line, int num_cols,
                   std::vector<float>& cat) {
    std::stringstream stream(line);
    for (int column = 0; column < num_cols; ++column) {
        if (!(stream >> cat[column])) return false;
    }
    return true;
}

// ==========================================
// Function: Apply the Lite combined-catalog source cuts
// Method: Preserve the existing peak-boundary and invalid-PSF rejection order.
// ==========================================
bool passesCombinedCatalogCuts(const std::vector<float>& cat) {
    if (cat[LensingConfig::i_imax] >= LensingConfig::ns
        || cat[LensingConfig::i_jmax] >= LensingConfig::ns) {
        return false;
    }
    return cat[0] >= -900.0f;
}

// ==========================================
// Function: Apply the Lite external-catalog shear calibration
// Method: Execute the frozen zero-correction formula so finite and exceptional
//         floating-point propagation stays identical to the Standard branch.
// ==========================================
void applyLiteCatalogCalibration(std::vector<float>& cat) {
    constexpr double g1c = 0.0;
    constexpr double g2c = 0.0;
    cat[LensingConfig::ig1] = static_cast<float>(
        cat[LensingConfig::ig1]
        - g1c * cat[LensingConfig::ide]
        + g1c * cat[LensingConfig::ih1]
        + g2c * cat[LensingConfig::ih2]);
    cat[LensingConfig::ig2] = static_cast<float>(
        cat[LensingConfig::ig2]
        - g2c * cat[LensingConfig::ide]
        + g1c * cat[LensingConfig::ih2]
        - g2c * cat[LensingConfig::ih1]);
}

}  // namespace

// ==========================================
// Function: Combine one exposure's chip catalogs into the final result catalog
// Method: Preflight paired physical row counts with one fresh retry, preserve
//         header-only shear sentinels, and consume matched pairs in a fixed loop.
// ==========================================
void combineExpoCatalog(int nchip, const std::vector<std::string>& imageFiles,
                        const std::string& dirOutput, int expo_index,
                        float chi2) {
    if (nchip <= 0 || imageFiles.empty()) {
        MPIFailure::abortWorld(
            "combine exposure catalog", "exposure contains no chip paths");
    }

    const std::string prefix_expo =
        UniversalUtils::getPrefixExpo(imageFiles[0]);
    const std::string out_filename =
        dirOutput + "/result/" + prefix_expo + "_all.cat";

    std::error_code filesystem_error;
    std::filesystem::remove(out_filename, filesystem_error);
    if (filesystem_error) {
        MPIFailure::abortWorld(
            "remove stale combined catalog",
            out_filename + ": " + filesystem_error.message());
    }

    MainIO::OutputFile fout20;

    int accepted_count = 0;
    int rejected_count = 0;
    const int num_cols = LensingConfig::shear_cat_ncols;
    bool output_opened = false;
    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        const Universalblock::NormStatus normStatus =
            Universalblock::checkNorm(imageFiles[ichip], dirOutput);
        if (normStatus == Universalblock::NormStatus::Invalid) {
            continue;
        }
        if (normStatus != Universalblock::NormStatus::Valid) {
            MPIFailure::abortWorld(
                "check Stage 1 norm before catalog combination",
                Universalblock::normErrorDetail(
                    normStatus, imageFiles[ichip], dirOutput));
        }

        const int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        const std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        const std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        const std::string filename_orig = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
        const std::size_t paired_data_rows = determinePairedDataRows(
            filename_shear, filename_orig, prefix);
        if (paired_data_rows == 0) {
            continue;
        }

        std::ifstream shear_input(filename_shear);
        if (!shear_input.is_open()) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }

        std::string shear_header;
        if (!std::getline(shear_input, shear_header)) {
            MPIFailure::abortWorld(
                "read Stage 7 shear catalog header", filename_shear);
        }
        shear_header = trimRight(shear_header);
        if (shear_header.empty()) {
            MPIFailure::abortWorld(
                "parse Stage 7 shear catalog header", filename_shear);
        }

        std::ifstream original_input(filename_orig);
        if (!original_input.is_open()) {
            MPIFailure::abortWorld("read external source catalog", filename_orig);
        }

        std::string original_header;
        if (!std::getline(original_input, original_header)) {
            MPIFailure::abortWorld(
                "read external source catalog header", filename_orig);
        }
        original_header = trimRight(original_header);
        if (original_header.empty()) {
            MPIFailure::abortWorld(
                "parse external source catalog header", filename_orig);
        }

        if (chi2 > LensingConfig::chi2_thresh) {
            std::cout << prefix << " contains no valid sources!" << std::endl;
            return;
        }

        if (!output_opened) {
            fout20.open(out_filename);
            fout20 << std::setprecision(10);
            fout20 << original_header << " EXPO_NUM ccD_NUM "
                   << shear_header << " Chi2\n";
            output_opened = true;
        }

        std::vector<float> cat(num_cols);
        for (std::size_t pair_index = 0;
             pair_index < paired_data_rows; ++pair_index) {
            std::string shear_line;
            std::string original_line;
            const bool shear_ok = static_cast<bool>(
                std::getline(shear_input, shear_line));
            const bool orig_ok = static_cast<bool>(
                std::getline(original_input, original_line));
            if (!shear_ok || !orig_ok) {
                std::ostringstream detail;
                detail << "prefix=" << prefix
                       << " row_index_zero_based=" << pair_index
                       << " row_index_one_based=" << (pair_index + 1)
                       << " expected_data_rows=" << paired_data_rows
                       << " shear_read_ok=" << shear_ok
                       << " orig_read_ok=" << orig_ok
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "read fixed paired catalog row", detail.str());
            }

            original_line = trimRight(original_line);
            if (original_line.empty()) {
                std::ostringstream detail;
                detail << "empty external catalog data row prefix=" << prefix
                       << " pair_index=" << (pair_index + 1)
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired external catalog row", detail.str());
            }

            if (shear_line.empty()) {
                std::ostringstream detail;
                detail << "empty shear data row prefix=" << prefix
                       << " pair_index=" << (pair_index + 1)
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired Stage 7 shear row", detail.str());
            }
            if (!parseShearRow(shear_line, num_cols, cat)) {
                std::ostringstream detail;
                detail << "incomplete shear data row prefix=" << prefix
                       << " pair_index=" << (pair_index + 1)
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired Stage 7 shear row", detail.str());
            }

            // A parsed shear row and its original catalog row form one pair.
            // Consume both before scientific rejection so omitted pairs cannot
            // shift the remaining external fields out of alignment.
            if (!passesCombinedCatalogCuts(cat)) {
                ++rejected_count;
                continue;
            }

            ++accepted_count;
            applyLiteCatalogCalibration(cat);
            fout20 << original_line << " " << expo_index << " "
                   << chip_index;
            for (int column = 0; column < num_cols; ++column) {
                fout20 << " " << cat[column];
            }
            fout20 << " " << chi2 << "\n";
        }

    }

    std::cout << (last_prefix.empty() ? prefix_expo : last_prefix)
              << " " << accepted_count << " "
              << rejected_count << std::endl;
    if (output_opened) {
        fout20.close();
    }
}

// ==========================================
// Function: Run Stage-9 catalog combination for one exposure
// Method: Resolve chip paths, obtain the reduced Stage-8 chi2, and invoke the
//         Lite external-catalog combiner.
// ==========================================
void procComb(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    const std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    float chi2 = 0.0f;
    if (ExposureInfo::state.parameters.size() >= static_cast<std::size_t>(iexpo) * 6) {
        chi2 = ExposureInfo::state.parameters[(iexpo - 1) * 6 + 2];
    }

    combineExpoCatalog(static_cast<int>(image_files.size()), image_files,
                       dir_output, iexpo, chi2);
}

}  // namespace CatalogCombiner
