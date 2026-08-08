#include "ExposureInfo.hpp"
#include "OutputFile.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "Astrometry.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

extern std::vector<std::string> EXPO_FILE;

namespace ExposureInfo {

std::vector<float> expo_para;

// ==========================================
// Function: Aggregate exposure-level chip diagnostics
// Method: Match F77 get_expo_info and stop immediately if astrometry head reading fails.
// ==========================================
void getExpoInfo(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput, float para[6]) {
    std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    std::string fstar = dirOutput + "/stamps/dat_StarInfo/" + prefix_expo + "_star_info_expo.dat";
    std::string fastro = dirOutput + "/astrometry/Head/" + prefix_expo + ".head";
    std::string fexpo = dirOutput + "/stamps/dat_ExpoInfo/" + prefix_expo + "_expo_info.dat";

    MainIO::OutputFile fout10(fexpo);
    if (!fout10.is_open()) {
        std::cerr << "Error: cannot open output expo info file: " << fexpo << std::endl;
        std::exit(1);
    }
    fout10 << "# ichip nstar FWHM e1 e2 chi_d cRPIX_1 cRPIX_2 cD_11 cD_12 cD_21 cD_22\n";

    std::ifstream fin20(fstar);
    if (!fin20.is_open()) {
        std::cerr << "Error: cannot open star info file: " << fstar << std::endl;
        std::exit(1);
    }
    
    std::string header;
    std::getline(fin20, header); // skip header line

    int nvalid = 0;
    float FWHM_AVE = 0.0f;
    float chi_d_AVE = 0.0f;
    float nstar_AVE = 0.0f;
    double cRVAL1 = 0.0;
    double cRVAL2 = 0.0;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        int i = 0;
        int nstar = 0;
        float FWHM = 0.0f, e1 = 0.0f, e2 = 0.0f, chi_d = 0.0f;
        
        if (!(fin20 >> i >> nstar >> FWHM >> e1 >> e2 >> chi_d)) {
            std::cerr << "Error reading chip info from star info file at chip: " << ichip + 1 << std::endl;
            std::exit(1);
        }

        if (nstar == 0) {
            fout10 << ichip + 1 << " 0 -99.0 -99.0 -99.0 -99.0 -99.0 -99.0 -99.0 -99.0 -99.0 -99.0\n";
            continue;
        }

        nvalid++;
        FWHM_AVE += FWHM;
        chi_d_AVE += chi_d;
        nstar_AVE += static_cast<float>(nstar);

        double cRPIX[2] = {0.0, 0.0};
        double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double cRVAL[2] = {0.0, 0.0};
        double PU[2][LensingConfig::npd] = {{0.0}, {0.0}};
        int ierror = 0;

        Astrometry::readAstrometryPara(fastro, ichip + 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd, ierror);
        if (ierror != 0) {
            std::cerr << "Error reading astrometry head file: " << fastro
                      << " at chip: " << ichip + 1 << std::endl;
            std::exit(1);
        }
        
        cRVAL1 = cRVAL[0];
        cRVAL2 = cRVAL[1];

        fout10 << std::setprecision(10)
               << ichip + 1 << " " << nstar << " " << FWHM << " " << e1 << " " << e2 << " " << chi_d << " "
               << std::setprecision(17)
               << cRPIX[0] << " " << cRPIX[1] << " "
               << cD[0][0] << " " << cD[0][1] << " "
               << cD[1][0] << " " << cD[1][1] << "\n";
    }
    fin20.close();
    fout10.close();

    if (nvalid > 0) {
        FWHM_AVE /= nvalid;
        chi_d_AVE /= nvalid;
        nstar_AVE /= nvalid;
    }

    para[0] = static_cast<float>(nvalid);
    para[1] = FWHM_AVE;
    para[2] = chi_d_AVE;
    para[3] = nstar_AVE;
    para[4] = static_cast<float>(cRVAL1);
    para[5] = static_cast<float>(cRVAL2);
}

void procInfo(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    float para[6] = {0.0f};
    getExpoInfo(image_files, image_files.size(), dir_output, para);

    // Initialize expo_para if it is not allocated or empty (it is defined in main.cpp, but here we can write directly to it)
    if (expo_para.size() < static_cast<size_t>(LensingConfig::NMAX_EXPO) * 6) {
        expo_para.resize(static_cast<size_t>(LensingConfig::NMAX_EXPO) * 6, 0.0f);
    }

    for (int i = 0; i < 6; ++i) {
        // F77: expo_para(i,iexpo)=para(i)
        // Memory index: (iexpo - 1) * 6 + i
        expo_para[(iexpo - 1) * 6 + i] = para[i];
    }
}

} // namespace ExposureInfo
