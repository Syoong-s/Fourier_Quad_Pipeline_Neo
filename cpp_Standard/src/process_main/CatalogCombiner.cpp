#include "CatalogCombiner.hpp"
#include "CatalogRowCount.hpp"
#include "OutputFile.hpp"
#include "MPIFailure.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "Universalblock.hpp"
#include "ExposureInfo.hpp"
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

enum class ShearCatalogStatus {
    HasSources,
    Empty,
    Missing,
    ReadError
};

struct ShearCatalogProbe {
    ShearCatalogStatus status = ShearCatalogStatus::ReadError;
    std::string header;
};

// ==========================================
// Function: Classify one Stage-7 shear catalog without opening its paired original catalog
// Method: Require a nonempty header and scan only until the first nonblank data row,
//         distinguishing a valid zero-source file from missing or unreadable input.
// ==========================================
ShearCatalogProbe probeShearCatalog(const std::string& filename) {
    ShearCatalogProbe result;
    std::ifstream input(filename);
    if (!input.is_open()) {
        result.status = ShearCatalogStatus::Missing;
        return result;
    }
    if (!std::getline(input, result.header)) {
        result.status = ShearCatalogStatus::ReadError;
        return result;
    }
    result.header = trimRight(result.header);
    if (result.header.empty()) {
        result.status = ShearCatalogStatus::ReadError;
        return result;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (!trimRight(line).empty()) {
            result.status = ShearCatalogStatus::HasSources;
            return result;
        }
    }
    result.status = input.bad()
        ? ShearCatalogStatus::ReadError
        : ShearCatalogStatus::Empty;
    return result;
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
// Method: Remove stale output, gate every chip by norm and shear data presence, then lazily
//         create the exposure catalog from the first contributing chip's live headers.
// ==========================================
void combineExpoCatalog(int nchip, const std::vector<std::string>& imageFiles,
                        const std::string& dirOutput, float chi2) {
    constexpr bool use_external_catalog = LensingConfig::ext_cat == 1;
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

    int n = 0;
    int m = 0;
    int num_cols = LensingConfig::shear_cat_ncols;
    bool output_opened = false;

    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        const Universalblock::NormStatus normStatus =
            Universalblock::checkNorm(imageFiles[ichip], dirOutput);
        if (normStatus == Universalblock::NormStatus::Invalid) {
            continue;
        }
        if (normStatus != Universalblock::NormStatus::Valid) {
            Universalblock::reportNormError(
                normStatus, imageFiles[ichip], dirOutput);
            continue;
        }

        int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        const std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        const ShearCatalogProbe shear_probe =
            probeShearCatalog(filename_shear);
        if (shear_probe.status == ShearCatalogStatus::Missing) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }
        if (shear_probe.status == ShearCatalogStatus::ReadError) {
            MPIFailure::abortWorld("parse Stage 7 shear catalog", filename_shear);
        }

        std::string filename_orig;
        if constexpr (use_external_catalog) {
            filename_orig = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
            Internal::requireMatchingCatalogDataRows(
                filename_shear, filename_orig);
        }
        if (shear_probe.status == ShearCatalogStatus::Empty) {
            continue;
        }
        if (chi2 > LensingConfig::chi2_thresh) {
            std::cout << prefix << " contains no valid sources!" << std::endl;
            return;
        }

        std::ifstream fin10(filename_shear);
        std::string ignored_shear_header;
        if (!fin10.is_open() || !std::getline(fin10, ignored_shear_header)) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }

        std::ifstream fin15;
        std::string original_header;
        if constexpr (use_external_catalog) {
            fin15.open(filename_orig);
            if (!fin15.is_open()) {
                MPIFailure::abortWorld("read external source catalog", filename_orig);
            }

            if (!std::getline(fin15, original_header)) {
                MPIFailure::abortWorld(
                    "read external source catalog header", filename_orig);
            }
            original_header = trimRight(original_header);
            if (original_header.empty()) {
                MPIFailure::abortWorld(
                    "parse external source catalog header", filename_orig);
            }
        }

        if (!output_opened) {
            fout20.open(out_filename);
            fout20 << std::setprecision(10);
            if constexpr (use_external_catalog) {
                fout20 << original_header << " ccD_NUM "
                       << shear_probe.header << " Chi2\n";
            } else {
                fout20 << " ccD_NUM " << shear_probe.header << " Chi2\n";
            }
            output_opened = true;
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
                fout20 << cat_content << " ";
            }

            fout20 << chip_index;

            for (int u = 0; u < num_cols; ++u) {
                fout20 << " " << cat[u];
            }
            fout20 << " " << chi2 << "\n";
        }

        if constexpr (use_external_catalog) {
            fin15.close();
        }
        fin10.close();
    }

    std::cout << (last_prefix.empty() ? prefix_expo : last_prefix)
              << " " << n << " " << m << std::endl;
    if (output_opened) {
        fout20.close();
    }
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
