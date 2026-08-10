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
    while (!str.empty() && (str.back() == ' ' || str.back() == '\r' || str.back() == '\n' || str.back() == '\t')) {
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
// Function: Apply the common combined-catalog source cuts
// Method: Preserve the existing peak-boundary and invalid-PSF rejection order.
// ==========================================
bool passesCombinedCatalogCuts(const std::vector<float>& cat) {
    if (cat[LensingConfig::i_imax] >= LensingConfig::ns ||
        cat[LensingConfig::i_jmax] >= LensingConfig::ns) {
        return false;
    }
    return !std::isnan(cat[0]) && cat[0] >= -900.0f;
}

// ==========================================
// Function: Apply the active combined-catalog shear calibration
// Method: Preserve zero correction for external catalogs and gf plus compiled
//         additive correction for the non-external branch.
// ==========================================
void applyCombinedCatalogCalibration(std::vector<float>& cat,
                                     bool use_external_catalog) {
    double g1c = 0.0;
    double g2c = 0.0;
    if (!use_external_catalog) {
        g1c = static_cast<double>(cat[LensingConfig::igf1])
            + LensingConfig::g1_c;
        g2c = static_cast<double>(cat[LensingConfig::igf2])
            + LensingConfig::g2_c;
    }

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
// Method: Share Stage-7 parsing, cuts, calibration, counters, and shear-column
//         output while conditionally aligning external-catalog rows.
// ==========================================
void combineExpoCatalog(int nchip, const std::vector<std::string>& imageFiles,
                        const std::string& dirOutput, float chi2) {
    constexpr bool use_external_catalog = LensingConfig::ext_cat == 1;
    std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    std::string out_filename = dirOutput + "/result/" + prefix_expo + "_all.cat";

    MainIO::OutputFile fout20(out_filename);
    fout20 << std::setprecision(10);

    std::string original_header;
    if constexpr (use_external_catalog) {
        for (int ichip = 0; ichip < nchip; ++ichip) {
            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
            std::string filename = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
            std::ifstream fin15(filename);
            if (fin15.is_open()) {
                if (std::getline(fin15, original_header)) {
                    std::string cat_content;
                    if (std::getline(fin15, cat_content)) {
                        original_header = trimRight(original_header);
                        break;
                    }
                }
            }
        }
    }

    int n = 0;
    int m = 0;
    int num_cols = LensingConfig::shear_cat_ncols;

    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        const std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        std::ifstream fin10(filename_shear);
        if (!fin10.is_open()) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }

        std::string cat_list1;
        std::getline(fin10, cat_list1);
        cat_list1 = trimRight(cat_list1);

        std::ifstream fin15;
        if constexpr (use_external_catalog) {
            const std::string filename_orig = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
            fin15.open(filename_orig);
            if (!fin15.is_open()) {
                MPIFailure::abortWorld("read external source catalog", filename_orig);
            }

            std::string dummy_orig_header;
            std::getline(fin15, dummy_orig_header);
        }

        if (ichip == 0) {
            if constexpr (use_external_catalog) {
                fout20 << original_header << " ccD_NUM " << cat_list1 << " Chi2\n";
            } else {
                fout20 << " ccD_NUM " << cat_list1 << "\n";
            }
            if (chi2 > LensingConfig::chi2_thresh) {
                fin10.close();
                if constexpr (use_external_catalog) {
                    fin15.close();
                }
                fout20.close();
                std::cout << prefix << " contains no valid sources!" << std::endl;
                return;
            }
        }

        std::vector<float> cat(num_cols);
        std::string line10;
        while (std::getline(fin10, line10)) {
            if (line10.empty()) continue;
            if (!parseShearRow(line10, num_cols, cat)) continue;

            std::string cat_content;
            if constexpr (use_external_catalog) {
                if (!std::getline(fin15, cat_content)) {
                    break;
                }
                cat_content = trimRight(cat_content);
            }

            if (!passesCombinedCatalogCuts(cat)) {
                ++m;
                continue;
            }

            ++n;
            applyCombinedCatalogCalibration(cat, use_external_catalog);

            if constexpr (use_external_catalog) {
                fout20 << cat_content << " " << chip_index;
            } else {
                fout20 << chip_index;
            }
            for (int u = 0; u < num_cols; ++u) {
                fout20 << " " << cat[u];
            }
            if constexpr (use_external_catalog) {
                fout20 << " " << chi2;
            }
            fout20 << "\n";
        }

        if constexpr (use_external_catalog) {
            fin15.close();
        }
        fin10.close();
    }

    std::cout << last_prefix << " " << n << " " << m << std::endl;
    fout20.close();
}

// ==========================================
// Function: Run Stage-9 catalog combination for one exposure
// Method: Resolve chip paths, obtain the reduced Stage-8 chi2, and invoke the
//         shared external/non-external catalog combiner.
// ==========================================
void procComb(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    float chi2 = 0.0f;
    if (ExposureInfo::expo_para.size() >= static_cast<size_t>(iexpo) * 6) {
        chi2 = ExposureInfo::expo_para[(iexpo - 1) * 6 + 2]; // 3rd element in Fortran, index 2
    }

    combineExpoCatalog(static_cast<int>(image_files.size()), image_files, dir_output, chi2);
}

} // namespace CatalogCombiner
