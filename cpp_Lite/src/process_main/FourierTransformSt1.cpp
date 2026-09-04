#include "process_main/FourierTransformSt1.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/MPIFailure.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "RuntimeConfig.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>


namespace FourierTransformSt1 {

    // ==========================================
    // Function: Run the Stage-4 star FFT for one exposure
    // Method: Resolve the exposure's chip list and process each chip serially.
    // ==========================================
    void procFourierTSt1(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

        for (const auto& image_file : image_files) {
            chipProcessFourierTSt1(image_file, dir_output);
        }
    }

    // ==========================================
    // Function: Transform one chip's star-candidate stamps to Fourier power
    // Method: Prepare the configured noise product, build one shared corrected-power path,
    //         regularize it, and reuse fixed-size scratch vectors.
    // ==========================================
    void chipProcessFourierTSt1(const std::string& imageFile,
                                const std::string& dirOutput) {
        const int star_smooth =
            RuntimeConfigStore::get().lensing.star_smooth;
        const Universalblock::NormStatus normStatus =
            Universalblock::checkNorm(imageFile, dirOutput);
        if (normStatus == Universalblock::NormStatus::Invalid) {
            return;
        }
        if (normStatus != Universalblock::NormStatus::Valid) {
            MPIFailure::abortWorld(
                "check Stage 1 norm before star FFT",
                Universalblock::normErrorDetail(
                    normStatus, imageFile, dirOutput));
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
            std::vector<float> noise_product(stamp_size);
            std::vector<float> source_p(stamp_size);
            std::vector<float> noise_p(stamp_size);

            for (int i = 0; i < nsource; ++i) {
                const std::size_t offset = static_cast<std::size_t>(i) * stamp_size;
                std::copy_n(source_coll.data() + offset, stamp_size, source.data());
                std::copy_n(
                    noise_coll.data() + offset, stamp_size, noise_product.data());
                double source_pc = 0.0;

                if (!ImageProcessing::prepareNoisePower(
                        ns, noise_product, noise_p)) {
                    MPIFailure::abortWorld(
                        "prepare star-candidate noise power",
                        filename_star_can_noise);
                }
                if (!ImageProcessing::buildCorrectedPower(
                        ns, ns, source, noise_p, star_smooth,
                        source_p, source_pc)) {
                    MPIFailure::abortWorld(
                        "build corrected star-candidate power",
                        filename_star_can_noise);
                }
                ImageProcessing::regularizePower(ns, ns, source_p, star_smooth);

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
