#include "process_fd/ShearCatalogReader.hpp"
#include "FDConfig.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace fc = FDConfig;
namespace lc = LensingConfig;

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

    // Duplicate-detection buffer (stores up to MAX_DUP rows per source)
    std::vector<std::vector<float>> dup_buf(fc::MAX_DUP,
                                            std::vector<float>(fc::ICHI2));
    int ndup = 0;
    float last_dec = -999.0, last_ra = -999.0;
    const float multisam_thrsh = 1e-7;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<float> item(fc::ICHI2);
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
        if (item[fc::col_mag_i] < 10.0 || item[fc::col_mag_i] > 30.0) continue;

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
        if (std::fabs(item[fc::col_gf1]) > 0.0015) continue;
        if (std::fabs(item[fc::col_gf2]) > 0.0015) continue;

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
            for (int i = 0; i < ndup && i < fc::MAX_DUP; ++i) {
                int idx = data.ng;
                if (idx >= fc::nmax_per_core) {
                    std::cerr << "nmax_per_core is too small!" << std::endl;
                    file.close();
                    return;
                }
                data.x1[idx]  = dup_buf[i][fc::col_gf1];
                data.y1[idx]  = dup_buf[i][fc::col_g1];
                data.de1[idx] = dup_buf[i][fc::col_de] - dup_buf[i][fc::col_h1];
                data.x2[idx]  = dup_buf[i][fc::col_gf2];
                data.y2[idx]  = dup_buf[i][fc::col_g2];
                data.de2[idx] = dup_buf[i][fc::col_de] + dup_buf[i][fc::col_h1];

                // Jackknife weight: ww = ±1/sqrt((g1/de)² + (g2/de)²)
                {
                    float de_val = dup_buf[i][fc::col_de];
                    float y1j = dup_buf[i][fc::col_g1] / de_val;
                    float y2j = dup_buf[i][fc::col_g2] / de_val;
                    float gmag = std::sqrt(y1j * y1j + y2j * y2j);
                    float sign = (de_val >= 0.0) ? 1.0 : -1.0;
                    data.ww[idx] = (gmag > 0.0) ? sign / gmag : 0.0;
                }

                data.magr[idx]    = dup_buf[i][fc::col_mag_r];
                data.magg[idx]    = dup_buf[i][fc::col_mag_g];
                data.magi[idx]    = dup_buf[i][fc::col_mag_i];
                data.sizerel[idx] =
                    ((std::sqrt(dup_buf[i][fc::col_h_area] / LensingConfig::pi) * LensingConfig::pixel_size * 2.0)
                     - dup_buf[i][fc::col_PSF]) / dup_buf[i][fc::col_PSF];
                data.src_snr[idx] = dup_buf[i][fc::col_h_flux] /
                                     std::sqrt(dup_buf[i][fc::col_h_area]);
                data.rra[idx]  = dup_buf[i][fc::col_ra];
                data.ddec[idx] = dup_buf[i][fc::col_dec];

                if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                    data.iexpo[idx] = static_cast<int>(std::lround(dup_buf[i][fc::col_chi2]));
                    data.snrf[idx]  = dup_buf[i][fc::col_SNR_F];
                }
                data.ng++;
            }
            ndup = 0;
        }

        // Store current row in duplicate buffer
        if (ndup < fc::MAX_DUP) {
            for (int u = 0; u < fc::ICHI2; ++u)
                dup_buf[ndup][u] = item[u];
        }
        ndup++;
    }

    // Flush remaining duplicates
    if (ndup > 0) {
        for (int i = 0; i < ndup && i < fc::MAX_DUP; ++i) {
            int idx = data.ng;
            if (idx >= fc::nmax_per_core) {
                std::cerr << "nmax_per_core is too small!" << std::endl;
                break;
            }
            data.x1[idx]  = dup_buf[i][fc::col_gf1];
            data.y1[idx]  = dup_buf[i][fc::col_g1];
            data.de1[idx] = dup_buf[i][fc::col_de] - dup_buf[i][fc::col_h1];
            data.x2[idx]  = dup_buf[i][fc::col_gf2];
            data.y2[idx]  = dup_buf[i][fc::col_g2];
            data.de2[idx] = dup_buf[i][fc::col_de] + dup_buf[i][fc::col_h1];

            // Jackknife weight: ww = ±1/sqrt((g1/de)² + (g2/de)²)
            {
                float de_val = dup_buf[i][fc::col_de];
                float y1j = dup_buf[i][fc::col_g1] / de_val;
                float y2j = dup_buf[i][fc::col_g2] / de_val;
                float gmag = std::sqrt(y1j * y1j + y2j * y2j);
                float sign = (de_val >= 0.0) ? 1.0 : -1.0;
                data.ww[idx] = (gmag > 0.0) ? sign / gmag : 0.0;
            }

            data.magr[idx]    = dup_buf[i][fc::col_mag_r];
            data.magg[idx]    = dup_buf[i][fc::col_mag_g];
            data.magi[idx]    = dup_buf[i][fc::col_mag_i];
            data.sizerel[idx] =
                ((std::sqrt(dup_buf[i][fc::col_h_area] / LensingConfig::pi) * LensingConfig::pixel_size * 2.0)
                 - dup_buf[i][fc::col_PSF]) / dup_buf[i][fc::col_PSF];
            data.src_snr[idx] = dup_buf[i][fc::col_h_flux] /
                                 std::sqrt(dup_buf[i][fc::col_h_area]);
            data.rra[idx]  = dup_buf[i][fc::col_ra];
            data.ddec[idx] = dup_buf[i][fc::col_dec];

            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                data.iexpo[idx] = static_cast<int>(std::lround(dup_buf[i][fc::col_chi2]));
                data.snrf[idx]  = dup_buf[i][fc::col_SNR_F];
            }
            data.ng++;
        }
    }
    file.close();
}
