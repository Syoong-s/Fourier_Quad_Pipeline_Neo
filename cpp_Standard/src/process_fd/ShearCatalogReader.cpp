#include "process_fd/ShearCatalogReader.hpp"
#include "FDConfig.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace fc = FDConfig;

// ==========================================
// Function: readExposure
// Method: Read one exposure's _all.cat file, apply quality cuts, deduplicate
//         overlapping detections, and append valid sources to the shared
//         FDData arrays.  Faithful translation of Fortran read_shear_cat_v2.
// ==========================================
void ShearCatalogReader::readExposure(int iexpo, FDData& data,
                                      const std::vector<std::string>& expo_files,
                                      int rank) {
    if (iexpo < 1 || iexpo > static_cast<int>(expo_files.size())) {
        if (rank == 0)
            std::cerr << "Invalid exposure index: " << iexpo << std::endl;
        return;
    }

    const std::string& filename = expo_files[iexpo - 1];
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << filename << " does not exist!!" << std::endl;
        return;
    }

    // Skip header line
    std::string header;
    std::getline(file, header);

    using CatalogRow = std::array<float, fc::ICHI2>;
    using DuplicateBuffer = std::array<CatalogRow, fc::MAX_DUP>;

    DuplicateBuffer dup_buf{};
    CatalogRow item{};
    int ndup = 0;
    float last_dec = -999.0, last_ra = -999.0;
    constexpr float multisam_thrsh = 1e-7f;

    // ==========================================
    // Function: Flush one duplicate group into the FD arrays
    // Method: Preserve first-MAX_DUP ordering and all derived shear, weight,
    //         magnitude, size, diagnostic, coordinate, and exposure fields.
    // ==========================================
    auto flushDuplicates = [&]() -> bool {
        const int count = std::min(ndup, fc::MAX_DUP);
        for (int i = 0; i < count; ++i) {
            const int idx = data.ng;
            if (idx >= fc::nmax_per_core) {
                std::cerr << "nmax_per_core is too small!" << std::endl;
                return false;
            }

            const CatalogRow& row = dup_buf[i];
            data.x1[idx] = row[fc::col_gf1];
            data.y1[idx] = row[fc::col_g1];
            data.de1[idx] = row[fc::col_de] - row[fc::col_h1];
            data.x2[idx] = row[fc::col_gf2];
            data.y2[idx] = row[fc::col_g2];
            data.de2[idx] = row[fc::col_de] + row[fc::col_h1];

            const float de_val = row[fc::col_de];
            const float y1j = row[fc::col_g1] / de_val;
            const float y2j = row[fc::col_g2] / de_val;
            const float gmag = std::sqrt(y1j * y1j + y2j * y2j);
            const float sign = de_val >= 0.0f ? 1.0f : -1.0f;
            data.ww[idx] = gmag > 0.0f ? sign / gmag : 0.0f;

            data.magr[idx] = row[fc::col_mag_r];
            data.magg[idx] = row[fc::col_mag_g];
            data.magi[idx] = row[fc::col_mag_i];
            data.sizerel[idx] =
                ((std::sqrt(row[fc::col_h_area] / LensingConfig::pi)
                  * LensingConfig::pixel_size * 2.0)
                 - row[fc::col_PSF])
                / row[fc::col_PSF];
            data.src_snr[idx] =
                row[fc::col_h_flux] / std::sqrt(row[fc::col_h_area]);
            data.delta_chi2[idx] = row[fc::col_delta_chi2];
            data.orth_ext[idx] = row[fc::col_orth_ext];
            data.rra[idx] = row[fc::col_ra];
            data.ddec[idx] = row[fc::col_dec];

            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                data.iexpo[idx] =
                    static_cast<int>(std::lround(row[fc::col_chi2]));
                data.snrf[idx] = row[fc::col_SNR_F];
            }
            ++data.ng;
        }

        ndup = 0;
        return true;
    };

    std::string line;
    std::istringstream iss;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        iss.clear();
        iss.str(line);
        bool parse_ok = true;
        for (int u = 0; u < fc::ICHI2; ++u) {
            if (!(iss >> item[u])) { parse_ok = false; break; }
        }
        if (!parse_ok) continue;

        // NaN / Inf check
        bool undefined = false;
        for (int u = 0; u < fc::ICHI2; ++u) {
            if (std::isnan(item[u])) { undefined = true; break; }
            if (std::fabs(item[u]) > 1.0e30) { undefined = true; break; }
        }
        if (undefined) continue;

        // Pixel coordinates for chip-edge masking
        int ix = static_cast<int>(item[fc::col_pixx]);  // pixel x
        int iy = static_cast<int>(item[fc::col_pixy]);  // pixel y

        // Bad CCD check
        int ccd_val = static_cast<int>(std::lround(item[fc::col_ccd]));
        bool bad_ccd = false;
        for (int i = 0; i < fc::n_bad_ccds; ++i) {
            if (ccd_val == fc::bad_ccds[i]) { bad_ccd = true; break; }
        }
        if (bad_ccd) continue;

        // Chip-edge masking (DES)
        if (ix < fc::chip_xmin || ix > fc::chip_xmax ||
            iy < fc::chip_ymin || iy > fc::chip_ymax)
            continue;

        // External-catalog flag cuts
        if (fc::ft_cut >= 0.0 && std::fabs(item[fc::col_flags_ft] - fc::ft_cut) > 1e-3) continue;
        if (fc::fg_cut >= 0.0 && std::fabs(item[fc::col_flags_fg] - fc::fg_cut) > 1e-3) continue;
        if (fc::gold_cut >= 0.0 && std::fabs(item[fc::col_flags_gold] - fc::gold_cut) > 1e-3) continue;
        if (fc::ext_cut >= 0.0 && std::fabs(item[fc::col_ext_mash] - fc::ext_cut) > 1e-3) continue;

        // SNR_F cut
        if (item[fc::col_SNR_F] < fc::snrfcut) continue;

        // SNR cut
        float snr = item[fc::col_h_flux] / std::sqrt(item[fc::col_h_area]);
        if (fc::snrlow > 0.0 && snr < fc::snrlow) continue;
        if (fc::snrhigh > 0.0 && snr > fc::snrhigh) continue;

        // Half-light radius cut
        if (fc::r_half_thresh > 0.0) {
            float r_half = std::sqrt(item[fc::col_h_area] / LensingConfig::pi) * LensingConfig::pixel_size * 2.0;
            if (r_half <= fc::r_half_thresh * item[fc::col_PSF]) continue;
        }

        // Magnitude cut
        if (item[fc::col_mag_i] < fc::mag_min_val ||
            item[fc::col_mag_i] > fc::mag_max_val) {
            continue;
        }

        // PSF polychi2 cut
        if (item[fc::col_polychi2] > fc::psf_chi2_mltp) continue;

        // Star classification cut
        if (item[fc::col_star] < fc::starcut) continue;

        // Chi2 cut (single mode only; per-exposure uses col_chi2 for exposure index)
        if (!fc::FD_PER_EXPOSURE_STAR_BAR) {
            if (item[fc::col_chi2] > fc::chi2_thresh) continue;
        }

        // Flag cut
        if (item[fc::col_flag] <= fc::flagcut) continue;

        // Imax / Jmax cut
        if (item[fc::col_imax] >= fc::imaxcut) continue;
        if (item[fc::col_jmax] >= fc::jmaxcut) continue;

        // Zero-point cut
        if (item[fc::col_zp] <= fc::zplow) continue;
        if (item[fc::col_zp] >= fc::zphigh) continue;

        // Field-distortion range cut
        if (std::fabs(item[fc::col_gf1]) > fc::gf_lim) continue;
        if (std::fabs(item[fc::col_gf2]) > fc::gf_lim) continue;

        // --- Duplicate detection ---
        bool new_galaxy = false;
        if (std::fabs(item[fc::col_dec] - last_dec) > multisam_thrsh ||
            std::fabs(item[fc::col_ra] - last_ra) > multisam_thrsh) {
            new_galaxy = true;
            last_dec = item[fc::col_dec];
            last_ra = item[fc::col_ra];
        }

        // Flush previous duplicates when a new galaxy is detected
        if (new_galaxy && ndup > 0) {
            if (!flushDuplicates()) return;
        }

        // Store current row in duplicate buffer
        if (ndup < fc::MAX_DUP) {
            dup_buf[ndup] = item;
        }
        ++ndup;
    }

    // Flush remaining duplicates
    if (ndup > 0) {
        if (!flushDuplicates()) return;
    }
    file.close();
}
