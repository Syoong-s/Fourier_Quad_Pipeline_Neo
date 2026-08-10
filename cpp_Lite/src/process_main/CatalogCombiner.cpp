#include "CatalogCombiner.hpp"
#include "OutputFile.hpp"
#include "MPIFailure.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "ExposureInfo.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

extern std::vector<std::string> EXPO_FILE;

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
    return !std::isnan(cat[0]) && cat[0] >= -900.0f;
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
// Method: Parse and filter Stage-7 rows once, preserve external-row alignment,
//         and publish the Lite external-catalog layout through checked output.
// ==========================================
void combineExpoCatalog(int nchip, const std::vector<std::string>& imageFiles,
                        const std::string& dirOutput, float chi2) {
    const std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    const std::string out_filename =
        dirOutput + "/result/" + prefix_expo + "_all.cat";

    MainIO::OutputFile fout20(out_filename);
    fout20 << std::setprecision(10);

    std::string original_header;
    for (int ichip = 0; ichip < nchip; ++ichip) {
        const std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        const std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
        std::ifstream input(filename);
        if (input.is_open() && std::getline(input, original_header)) {
            std::string first_row;
            if (std::getline(input, first_row)) {
                original_header = trimRight(original_header);
                break;
            }
        }
    }

    int accepted_count = 0;
    int rejected_count = 0;
    const int num_cols = LensingConfig::shear_cat_ncols;
    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        const int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        const std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        const std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        std::ifstream shear_input(filename_shear);
        if (!shear_input.is_open()) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }

        std::string shear_header;
        std::getline(shear_input, shear_header);
        shear_header = trimRight(shear_header);

        const std::string filename_orig = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
        std::ifstream original_input(filename_orig);
        if (!original_input.is_open()) {
            MPIFailure::abortWorld("read external source catalog", filename_orig);
        }

        std::string ignored_original_header;
        std::getline(original_input, ignored_original_header);

        if (ichip == 0) {
            fout20 << original_header << " ccD_NUM " << shear_header
                   << " Chi2\n";
            if (chi2 > LensingConfig::chi2_thresh) {
                shear_input.close();
                original_input.close();
                fout20.close();
                std::cout << prefix << " contains no valid sources!" << std::endl;
                return;
            }
        }

        std::vector<float> cat(num_cols);
        std::string shear_line;
        while (std::getline(shear_input, shear_line)) {
            if (shear_line.empty()) continue;
            if (!parseShearRow(shear_line, num_cols, cat)) continue;

            std::string original_line;
            if (!std::getline(original_input, original_line)) break;
            original_line = trimRight(original_line);

            if (!passesCombinedCatalogCuts(cat)) {
                ++rejected_count;
                continue;
            }

            ++accepted_count;
            applyLiteCatalogCalibration(cat);
            fout20 << original_line << " " << chip_index;
            for (int column = 0; column < num_cols; ++column) {
                fout20 << " " << cat[column];
            }
            fout20 << " " << chi2 << "\n";
        }
        original_input.close();
        shear_input.close();
    }

    std::cout << last_prefix << " " << accepted_count << " "
              << rejected_count << std::endl;
    fout20.close();
}

// ==========================================
// Function: Run Stage-9 catalog combination for one exposure
// Method: Resolve chip paths, obtain the reduced Stage-8 chi2, and invoke the
//         Lite external-catalog combiner.
// ==========================================
void procComb(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    const std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    float chi2 = 0.0f;
    if (ExposureInfo::expo_para.size() >= static_cast<std::size_t>(iexpo) * 6) {
        chi2 = ExposureInfo::expo_para[(iexpo - 1) * 6 + 2];
    }

    combineExpoCatalog(static_cast<int>(image_files.size()), image_files,
                       dir_output, chi2);
}

}  // namespace CatalogCombiner
