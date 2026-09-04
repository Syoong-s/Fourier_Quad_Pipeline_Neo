#include "process_main/FourierTransformSt2.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/MPIFailure.hpp"
#include "process_main/OutputFile.hpp"
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
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>


namespace FourierTransformSt2 {

// ==========================================
// Function: Transform one chip's source stamps into Fourier-space products
// Method: Gate on the Stage-1 norm, require the Stage-3 catalog and header, keep
//         source-only smooth-2 diagnostics, and publish through checked writers.
// ==========================================
void chipProcessFourierTSt2(const std::string& imageFile, const std::string& dirOutput) {
    const LensingRuntimeConfig& lensing =
        RuntimeConfigStore::get().lensing;
    const Universalblock::NormStatus normStatus =
        Universalblock::checkNorm(imageFile, dirOutput);
    if (normStatus == Universalblock::NormStatus::Invalid) {
        return;
    }
    if (normStatus != Universalblock::NormStatus::Valid) {
        MPIFailure::abortWorld(
            "check Stage 1 norm before galaxy FFT",
            Universalblock::normErrorDetail(normStatus, imageFile, dirOutput));
    }

    int ns = LensingConfig::ns;

    std::string raw_prefix = UniversalUtils::getPrefix(imageFile);
    // PREFIX inlined: per-type stamps/ subdirs (reorganized layout)

    int nsource = 0;
    std::string info_filename = OutputLayout::chipPath(
        dirOutput, "stamps/dat_SrcInfo", raw_prefix, "_source_info.dat");
    std::vector<std::vector<float>> source_para;
    std::vector<bool> bad_source_para;

    std::ifstream fin(info_filename);
    if (!fin.is_open()) {
        MPIFailure::abortWorld(
            "read source info for galaxy FFT", info_filename);
    }

    std::string header;
    if (!std::getline(fin, header)) {
        MPIFailure::abortWorld(
            "read source-info header for galaxy FFT", info_filename);
    }

    std::string line;
    int iflag_col = LensingConfig::iflag + 1; // 10 columns
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<float> row(LensingConfig::src_npara, 0.0f);
        bool bad_source = false;
        for (int i = 0; i < iflag_col; ++i) {
            std::string token;
            if (!(ss >> token)) {
                MPIFailure::abortWorld(
                    "parse Stage-3 source-info row", info_filename);
            }

            char* token_end = nullptr;
            const float value = std::strtof(token.c_str(), &token_end);
            if (token_end != token.c_str() + token.size()) {
                MPIFailure::abortWorld(
                    "parse Stage-3 source-info token", info_filename);
            }
            row[i] = value;
            bad_source = bad_source || !std::isfinite(value);
        }
        if (bad_source) {
            std::fill(row.begin(), row.end(), -99999.0f);
        }
        source_para.push_back(std::move(row));
        bad_source_para.push_back(bad_source);
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

    std::vector<float> source_coll;
    std::string source_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Src", raw_prefix, "_source.fits");
    FitsIO::StampCubeShape sourceShape;
    if (!FitsIO::readStampCube(source_fits, sourceShape, source_coll)
        || !sourceShape.matches(ns, ns, nsource)) {
        MPIFailure::abortWorld(
            "read source cube with expected shape " + std::to_string(ns) + "x"
                + std::to_string(ns) + "x" + std::to_string(nsource),
            source_fits);
    }

    std::vector<float> noise_coll;
    std::string noise_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Noise", raw_prefix, "_noise.fits");
    FitsIO::StampCubeShape noiseShape;
    if (!FitsIO::readStampCube(noise_fits, noiseShape, noise_coll)
        || !noiseShape.matches(ns, ns, nsource)) {
        MPIFailure::abortWorld(
            "read noise cube with expected shape " + std::to_string(ns) + "x"
                + std::to_string(ns) + "x" + std::to_string(nsource),
            noise_fits);
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
        if (bad_source_para[i]) {
            continue;
        }

        const std::size_t offset = static_cast<std::size_t>(i) * stamp_size;
        std::copy_n(source_coll.data() + offset, stamp_size, source.data());
        std::copy_n(
            noise_coll.data() + offset, stamp_size, noise_product.data());
        double pc = 0.0;

        // SNR retains the main-branch source-only smooth-2 definition.
        ImageProcessing::getPower(ns, ns, source, source_p, 2, pc);

        int ns_2 = LensingConfig::ns_2;
        float cen_val = source_p[ns_2 * ns + ns_2];
        source_para[i][10] = std::sqrt(std::max(static_cast<float>(pc), cen_val));
        source_para[i][11] = source_para[i][10] / source_para[i][3] * ns;
        for (int j = 0; j <= LensingConfig::iSNR_F; ++j) {
            if (!std::isfinite(source_para[i][j])) {
                std::fill(source_para[i].begin(), source_para[i].end(), -99999.0f);
                bad_source_para[i] = true;
                break;
            }
        }
        if (bad_source_para[i]) {
            continue;
        }

        if (!ImageProcessing::prepareNoisePower(
                ns, noise_product, LensingConfig::NstampType, noise_p)) {
            MPIFailure::abortWorld("prepare source noise power", noise_fits);
        }
        if (!ImageProcessing::buildCorrectedPower(
                ns, ns, source, noise_p, lensing.gal_smooth, source_p, pc)) {
            MPIFailure::abortWorld("build corrected source power", noise_fits);
        }

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
    if (!FitsIO::writeStampCube(power_fits, ns, ns, nsource, power_coll)) {
        MainIO::failOutput("write FFT2 source power FITS", power_fits,
                           "FitsIO::writeStampCube returned false");
    }
}

// ==========================================
// Function: Run the Stage-6 galaxy FFT for one exposure
// Method: Resolve the exposure's chip list and process each chip serially.
// ==========================================
void procFourierTSt2(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    for (const auto& image_file : image_files) {
        chipProcessFourierTSt2(image_file, dir_output);
    }
}

} // namespace FourierTransformSt2
