#include "FourierTransformSt2.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

extern std::vector<std::string> EXPO_FILE;

namespace FourierTransformSt2 {

// ==========================================
// Function: Transform one chip's source stamps into Fourier-space products
// Method: Read Stage-3 source data, update source diagnostics, and publish all
//         text/FITS outputs through checked main-process writers.
// ==========================================
void chipProcessFourierTSt2(const std::string& imageFile, const std::string& dirOutput) {
    int ns = LensingConfig::ns;

    std::string raw_prefix = UniversalUtils::getPrefix(imageFile);
    // PREFIX inlined: per-type stamps/ subdirs (reorganized layout)

    int nsource = 0;
    std::string info_filename = OutputLayout::chipPath(
        dirOutput, "stamps/dat_SrcInfo", raw_prefix, "_source_info.dat");
    std::vector<std::vector<float>> source_para;

    std::ifstream fin(info_filename);
    if (!fin.is_open()) {
        std::cerr << "Error / FFT2 source_info catalog file error!! " << info_filename << std::endl;
        return;
    }

    std::string header;
    std::getline(fin, header); // skip header line

    std::string line;
    int iflag_col = LensingConfig::iflag + 1; // 10 columns
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<float> row(LensingConfig::src_npara, 0.0f);
        bool success = true;
        for (int i = 0; i < iflag_col; ++i) {
            if (!(ss >> row[i])) {
                success = false;
                break;
            }
        }
        if (success) {
            source_para.push_back(row);
        }
    }
    fin.close();

    nsource = source_para.size();
    if (nsource == 0) {
        // Write header and return
        MainIO::OutputFile fout(info_filename);
        fout << "ig xp yp sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F\n";
        fout.close();
        return;
    }

    int len_g = LensingConfig::len_g;
    int ngal_max = LensingConfig::ngal_max;
    int nn1 = ns * len_g;
    int nn2 = ns * (nsource / len_g + 1);

    std::vector<float> source_coll;
    std::string source_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Src", raw_prefix, "_source.fits");
    if (!FitsIO::readStamps(ngal_max, 1, nsource, ns, ns, source_coll, nn1, nn2, source_fits)) {
        std::cerr << "Error reading source stamps: " << source_fits << std::endl;
        return;
    }

    std::vector<float> noise_coll;
    std::string noise_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Noise", raw_prefix, "_noise.fits");
    if (!FitsIO::readStamps(ngal_max, 1, nsource, ns, ns, noise_coll, nn1, nn2, noise_fits)) {
        std::cerr << "Error reading noise stamps: " << noise_fits << std::endl;
        return;
    }

    std::vector<float> power_coll(static_cast<size_t>(ngal_max) * ns * ns, 0.0f);
    const std::size_t stamp_size =
        static_cast<std::size_t>(ns) * static_cast<std::size_t>(ns);
    std::vector<float> source(stamp_size);
    std::vector<float> noise(stamp_size);
    std::vector<float> source_p(stamp_size);
    std::vector<float> noise_p(stamp_size);

    for (int i = 0; i < nsource; ++i) {
        const std::size_t offset = static_cast<std::size_t>(i) * stamp_size;
        std::copy_n(source_coll.data() + offset, stamp_size, source.data());
        std::copy_n(noise_coll.data() + offset, stamp_size, noise.data());
        double pc = 0.0;

        // SNR calculation uses star_smooth (2)
        ImageProcessing::getPower(ns, ns, source, source_p,
                                  LensingConfig::star_smooth, pc);

        int ns_2 = LensingConfig::ns_2;
        float cen_val = source_p[ns_2 * ns + ns_2];
        source_para[i][10] = std::sqrt(std::max(static_cast<float>(pc), cen_val));
        source_para[i][11] = source_para[i][10] / source_para[i][3] * ns;

        // Main power spectrum computation uses gal_smooth
        ImageProcessing::getPower(ns, ns, source, source_p, LensingConfig::gal_smooth, pc);
        ImageProcessing::getPower(ns, ns, noise, noise_p, LensingConfig::gal_smooth, pc);
        ImageProcessing::processPowers(ns, source_p, noise_p);

        std::copy_n(source_p.data(), stamp_size, power_coll.data() + offset);
    }

    MainIO::OutputFile fout(info_filename);
    fout << "ig xp yp sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F\n";
    for (int i = 0; i < nsource; ++i) {
        for (int j = 0; j <= LensingConfig::iSNR_F; ++j) {
            fout << source_para[i][j] << (j == LensingConfig::iSNR_F ? "" : " ");
        }
        fout << "\n";
    }
    fout.close();

    std::string power_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_SrcP", raw_prefix, "_source_p.fits");
    if (!FitsIO::writeStamps(ngal_max, 1, nsource, ns, ns, power_coll, nn1, nn2, power_fits)) {
        MainIO::failOutput("write FFT2 source power FITS", power_fits,
                           "FitsIO::writeStamps returned false");
    }
}

// ==========================================
// Function: Run the Stage-6 galaxy FFT for one exposure
// Method: Resolve the exposure's chip list and process each chip serially.
// ==========================================
void procFourierTSt2(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    for (const auto& image_file : image_files) {
        chipProcessFourierTSt2(image_file, dir_output);
    }
}

} // namespace FourierTransformSt2
