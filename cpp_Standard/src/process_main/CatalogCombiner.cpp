#include "CatalogCombiner.hpp"
#include "OutputFile.hpp"
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
#include <algorithm>

extern std::vector<std::string> EXPO_FILE;

namespace CatalogCombiner {

static inline std::string trimRight(std::string str) {
    while (!str.empty() && (str.back() == ' ' || str.back() == '\r' || str.back() == '\n' || str.back() == '\t')) {
        str.pop_back();
    }
    return str;
}

// ==========================================
// Function: Combine one exposure's chip catalogs into the final result catalog
// Method: Read aligned original/shear rows and publish the merged text output
//         through the checked main-process writer.
// ==========================================
void combineExpoCatalog(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, float chi2) {
    std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    std::string out_filename = dirOutput + "/result/" + prefix_expo + "_all.cat";

    MainIO::OutputFile fout20(out_filename);
    if (!fout20.is_open()) {
        std::cerr << "Error: cannot open combined catalog output: " << out_filename << std::endl;
        std::exit(1);
    }
    fout20 << std::setprecision(10);

    std::string cat_list2;
    if (LensingConfig::ext_cat == 1) {
        for (int ichip = 0; ichip < nchip; ++ichip) {
            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
            std::string filename = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
            std::ifstream fin15(filename);
            if (fin15.is_open()) {
                if (std::getline(fin15, cat_list2)) {
                    std::string cat_content;
                    if (std::getline(fin15, cat_content)) {
                        cat_list2 = trimRight(cat_list2);
                        fin15.close();
                        break;
                    }
                }
                fin15.close();
            }
        }
    }

    int n = 0;
    int m = 0;
    int num_cols = LensingConfig::iparity + 1; // 24 columns

    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        std::ifstream fin10(filename_shear);
        if (!fin10.is_open()) {
            std::cerr << filename_shear << " is missing!" << std::endl;
            std::exit(1);
        }

        std::string cat_list1;
        std::getline(fin10, cat_list1);
        cat_list1 = trimRight(cat_list1);

        if (LensingConfig::ext_cat == 1) {
            std::string filename_orig = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
            std::ifstream fin15(filename_orig);
            if (!fin15.is_open()) {
                std::cerr << filename_orig << " is missing!" << std::endl;
                std::exit(1);
            }

            std::string dummy_orig_header;
            std::getline(fin15, dummy_orig_header);

            if (ichip == 0) {
                fout20 << cat_list2 << " ccD_NUM " << cat_list1 << " Chi2\n";
                if (chi2 > LensingConfig::chi2_thresh) {
                    fin10.close();
                    fin15.close();
                    fout20.close();
                    std::cout << prefix << " contains no valid sources!" << std::endl;
                    return;
                }
            }

            std::string line10;
            while (std::getline(fin10, line10)) {
                if (line10.empty()) continue;
                std::stringstream ss(line10);
                std::vector<float> cat(num_cols);
                bool read_ok = true;
                for (int u = 0; u < num_cols; ++u) {
                    if (!(ss >> cat[u])) {
                        read_ok = false;
                        break;
                    }
                }
                if (!read_ok) continue;

                std::string cat_content;
                if (!std::getline(fin15, cat_content)) {
                    break;
                }
                cat_content = trimRight(cat_content);

                if (cat[LensingConfig::i_imax] >= LensingConfig::ns || cat[LensingConfig::i_jmax] >= LensingConfig::ns) {
                    m++;
                    continue;
                }
                if (std::isnan(cat[0]) || cat[0] < -900.0f) {
                    m++;
                    continue;
                }

                n++;
                double g1c = 0.0;
                double g2c = 0.0;
                // double g1c = static_cast<double>(cat[LensingConfig::igf1]) + LensingConfig::g1_c;
                // double g2c = static_cast<double>(cat[LensingConfig::igf2]) + LensingConfig::g2_c;
                cat[LensingConfig::ig1] = static_cast<float>(cat[LensingConfig::ig1] - g1c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih1] + g2c * cat[LensingConfig::ih2]);
                cat[LensingConfig::ig2] = static_cast<float>(cat[LensingConfig::ig2] - g2c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih2] - g2c * cat[LensingConfig::ih1]);

                fout20 << cat_content << " " << chip_index;
                for (int u = 0; u < num_cols; ++u) {
                    fout20 << " " << cat[u];
                }
                fout20 << " " << chi2 << "\n";
            }
            fin15.close();
        } else {
            if (ichip == 0) {
                fout20 << " ccD_NUM " << cat_list1 << "\n";
                if (chi2 > LensingConfig::chi2_thresh) {
                    fin10.close();
                    fout20.close();
                    std::cout << prefix << " contains no valid sources!" << std::endl;
                    return;
                }
            }

            std::string line10;
            while (std::getline(fin10, line10)) {
                if (line10.empty()) continue;
                std::stringstream ss(line10);
                std::vector<float> cat(num_cols);
                bool read_ok = true;
                for (int u = 0; u < num_cols; ++u) {
                    if (!(ss >> cat[u])) {
                        read_ok = false;
                        break;
                    }
                }
                if (!read_ok) continue;

                if (cat[LensingConfig::i_imax] >= LensingConfig::ns || cat[LensingConfig::i_jmax] >= LensingConfig::ns) {
                    m++;
                    continue;
                }
                if (std::isnan(cat[0]) || cat[0] < -900.0f) {
                    m++;
                    continue;
                }

                n++;
                double g1c = cat[LensingConfig::igf1] + LensingConfig::g1_c;
                double g2c = cat[LensingConfig::igf2] + LensingConfig::g2_c;
                cat[LensingConfig::ig1] = static_cast<float>(cat[LensingConfig::ig1] - g1c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih1] + g2c * cat[LensingConfig::ih2]);
                cat[LensingConfig::ig2] = static_cast<float>(cat[LensingConfig::ig2] - g2c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih2] - g2c * cat[LensingConfig::ih1]);

                fout20 << chip_index;
                for (int u = 0; u < num_cols; ++u) {
                    fout20 << " " << cat[u];
                }
                fout20 << "\n";
            }
        }
        fin10.close();
    }

    std::cout << last_prefix << " " << n << " " << m << std::endl;
    fout20.close();
}

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
