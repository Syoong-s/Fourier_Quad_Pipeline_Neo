#include "FourierTransformSt1.hpp"
#include "MPIFailure.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "Universalblock.hpp"
#include "FitsIO.hpp"
#include "ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

extern std::vector<std::string> EXPO_FILE;

namespace FourierTransformSt1 {

    // ==========================================
    // Function: Run the Stage-4 star FFT for one exposure
    // Method: Resolve the exposure's chip list and process each chip serially.
    // ==========================================
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

    // ==========================================
    // Function: Transform one chip's star-candidate stamps to Fourier power
    // Method: Apply the shared norm gate before chip-product reads, then reuse fixed-size scratch
    //         vectors while preserving the existing FFT, subtraction, and regularization.
    // ==========================================
    void chipProcessFourierTSt1(const std::string& imageFile,
                                const std::string& dirOutput) {
        const Universalblock::NormStatus normStatus =
            Universalblock::checkNorm(imageFile, dirOutput);
        if (normStatus == Universalblock::NormStatus::Invalid) {
            return;
        }
        if (normStatus != Universalblock::NormStatus::Valid) {
            Universalblock::reportNormError(normStatus, imageFile, dirOutput);
            return;
        }

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

            std::vector<float> source_coll;
            std::vector<float> noise_coll;

            std::string filename_star_can = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCan", raw_prefix, "_star_can.fits");
            FitsIO::StampCubeShape sourceShape;
            if (!FitsIO::readStampCube(
                    filename_star_can, sourceShape, source_coll)
                || !sourceShape.matches(ns, ns, nsource)) {
                MPIFailure::abortWorld(
                    "read star-candidate cube with expected shape "
                        + std::to_string(ns) + "x" + std::to_string(ns) + "x"
                        + std::to_string(nsource),
                    filename_star_can);
            }

            std::string filename_star_can_noise = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanN", raw_prefix, "_star_can_noise.fits");
            FitsIO::StampCubeShape noiseShape;
            if (!FitsIO::readStampCube(
                    filename_star_can_noise, noiseShape, noise_coll)
                || !noiseShape.matches(ns, ns, nsource)) {
                MPIFailure::abortWorld(
                    "read star-candidate noise cube with expected shape "
                        + std::to_string(ns) + "x" + std::to_string(ns) + "x"
                        + std::to_string(nsource),
                    filename_star_can_noise);
            }

            const std::size_t stamp_size =
                static_cast<std::size_t>(ns) * static_cast<std::size_t>(ns);
            std::vector<float> power_coll(
                static_cast<std::size_t>(nsource) * stamp_size, 0.0f);
            std::vector<float> source(stamp_size);
            std::vector<float> noise(stamp_size);
            std::vector<float> source_p(stamp_size);
            std::vector<float> noise_p(stamp_size);

            for (int i = 0; i < nsource; ++i) {
                const std::size_t offset = static_cast<std::size_t>(i) * stamp_size;
                std::copy_n(source_coll.data() + offset, stamp_size, source.data());
                std::copy_n(noise_coll.data() + offset, stamp_size, noise.data());
                double source_pc = 0.0;
                double noise_pc = 0.0;

                ImageProcessing::getPower(ns, ns, source, source_p, LensingConfig::star_smooth, source_pc);
                ImageProcessing::getPower(ns, ns, noise, noise_p, LensingConfig::star_smooth, noise_pc);
                ImageProcessing::processPowers(ns, source_p, noise_p);
                ImageProcessing::regularizePower(ns, ns, source_p, LensingConfig::star_smooth);

                std::copy_n(source_p.data(), stamp_size,
                            power_coll.data() + offset);
            }

            std::string filename_star_can_power = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", raw_prefix, "_star_can_power.fits");
            FitsIO::writeStampCube(filename_star_can_power, ns, ns, nsource,
                                   power_coll);
        }
    }
}
