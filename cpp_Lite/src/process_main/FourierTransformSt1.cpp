#include "FourierTransformSt1.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

extern std::vector<std::string> EXPO_FILE;

namespace FourierTransformSt1 {

    void procFourierTSt1(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = EXPO_FILE[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

        for (const auto& image_file : image_files) {
            chipProcessFourierTSt1(image_file, dir_output);
        }
    }

    void chipProcessFourierTSt1(const std::string& imageFile, const std::string& dirOutput) {
        std::string raw_prefix = UniversalUtils::getPrefix(imageFile);
        // PREFIX inlined: per-type stamps/ subdirs (reorganized layout)

        int nsource = 0;
        std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/dat_StarCanInfo", raw_prefix, "_star_can_info.dat");

        std::ifstream fin(filename);
        if (!fin.is_open()) {
            std::cerr << filename << "\n";
            std::cerr << "Error / FFT1 star_can_info catalog file error!!\n";
            return;
        }

        std::string header;
        std::getline(fin, header); // skip header line

        double val1, val2, val3, val4;
        while (fin >> val1 >> val2 >> val3 >> val4) {
            nsource++;
        }
        fin.close();

        if (nsource > 0) {
            int ns = LensingConfig::ns;
            int len_s = LensingConfig::len_s;
            int ngal_max = LensingConfig::ngal_max;

            int nn1 = ns * len_s;
            int nn2 = ns * (nsource / len_s + 1);

            std::vector<float> source_coll;
            std::vector<float> noise_coll;

            std::string filename_star_can = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCan", raw_prefix, "_star_can.fits");
            if (!FitsIO::readStamps(ngal_max, 1, nsource, ns, ns, source_coll, nn1, nn2, filename_star_can)) {
                std::cerr << "Error reading stamps: " << filename_star_can << "\n";
                return;
            }

            std::string filename_star_can_noise = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanN", raw_prefix, "_star_can_noise.fits");
            if (!FitsIO::readStamps(ngal_max, 1, nsource, ns, ns, noise_coll, nn1, nn2, filename_star_can_noise)) {
                std::cerr << "Error reading stamps: " << filename_star_can_noise << "\n";
                return;
            }

            std::vector<float> power_coll(static_cast<size_t>(ngal_max) * ns * ns, 0.0f);

            for (int i = 0; i < nsource; ++i) {
                std::vector<float> source(ns * ns);
                std::vector<float> noise(ns * ns);

                std::copy(source_coll.begin() + i * ns * ns, source_coll.begin() + (i + 1) * ns * ns, source.begin());
                std::copy(noise_coll.begin() + i * ns * ns, noise_coll.begin() + (i + 1) * ns * ns, noise.begin());

                std::vector<float> source_p(ns * ns);
                std::vector<float> noise_p(ns * ns);
                double source_pc = 0.0;
                double noise_pc = 0.0;

                ImageProcessing::getPower(ns, ns, source, source_p, LensingConfig::star_smooth, source_pc);
                ImageProcessing::getPower(ns, ns, noise, noise_p, LensingConfig::star_smooth, noise_pc);
                ImageProcessing::processPowers(ns, source_p, noise_p);
                ImageProcessing::regularizePower(ns, ns, source_p, LensingConfig::star_smooth);

                std::copy(source_p.begin(), source_p.end(), power_coll.begin() + i * ns * ns);
            }

            std::string filename_star_can_power = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", raw_prefix, "_star_can_power.fits");
            FitsIO::writeStamps(ngal_max, 1, nsource, ns, ns, power_coll, nn1, nn2, filename_star_can_power);
        }
    }
}
