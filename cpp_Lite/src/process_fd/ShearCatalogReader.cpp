#include "process_fd/ShearCatalogReader.hpp"
#include "FDConfig.hpp"
#include "LensingConfig.hpp"
#include "RuntimeConfig.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

namespace fc = FDConfig;

// ==========================================
// Function: Parse one exact-width numeric FD catalog row
// Method: Require the runtime layout width, finite values, and no trailing columns.
// ==========================================
bool ShearCatalogReader::parseCatalogRow(const std::string& line,
                                         std::size_t column_count,
                                         std::vector<float>& values) {
    if (column_count == 0) {
        return false;
    }

    std::istringstream input(line);
    std::vector<float> parsed(column_count);
    for (std::size_t column = 0; column < column_count; ++column) {
        if (!(input >> parsed[column]) || !std::isfinite(parsed[column])
            || std::fabs(parsed[column]) > 1.0e30f) {
            return false;
        }
    }
    std::string extra;
    if (input >> extra) {
        return false;
    }
    values = std::move(parsed);
    return true;
}

// ==========================================
// Function: readExposure
// Method: Read one exposure's _all.cat file, apply quality cuts, deduplicate
//         overlapping detections, and append valid sources to the shared
//         FDData arrays.  Faithful translation of Fortran read_shear_cat_v2.
// ==========================================
void ShearCatalogReader::readExposure(int iexpo, FDData& data,
                                      const std::vector<std::string>& expo_files,
                                      const PipelineCatalog::CatalogLayout& layout,
                                      std::size_t magnitude_column,
                                      int rank) {
    const double pixel_size = RuntimeConfigStore::get().lensing.pixel_size;
    if (iexpo < 1 || iexpo > static_cast<int>(expo_files.size())) {
        if (rank == 0)
            std::cerr << "Invalid exposure index: " << iexpo << std::endl;
        return;
    }
    if (magnitude_column >= layout.external_columns) {
        if (rank == 0) {
            std::cerr << "Invalid FD magnitude column: "
                      << magnitude_column << std::endl;
        }
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

    using CatalogRow = std::vector<float>;
    std::vector<CatalogRow> dup_buf(
        fc::MAX_DUP, CatalogRow(layout.all_columns));
    CatalogRow item;
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
            data.x1[idx] = row[layout.source.gf1];
            data.y1[idx] = row[layout.source.g1];
            data.de1[idx] = row[layout.source.de] - row[layout.source.h1];
            data.x2[idx] = row[layout.source.gf2];
            data.y2[idx] = row[layout.source.g2];
            data.de2[idx] = row[layout.source.de] + row[layout.source.h1];

            const float de_val = row[layout.source.de];
            const float y1j = row[layout.source.g1] / de_val;
            const float y2j = row[layout.source.g2] / de_val;
            const float gmag = std::sqrt(y1j * y1j + y2j * y2j);
            const float sign = de_val >= 0.0f ? 1.0f : -1.0f;
            data.ww[idx] = gmag > 0.0f ? sign / gmag : 0.0f;

            data.star_mag[idx] = row[magnitude_column];
            data.sizerel[idx] =
                ((std::sqrt(row[layout.source.h_area] / LensingConfig::pi)
                  * pixel_size * 2.0)
                 - row[layout.source.psf])
                / row[layout.source.psf];
            data.src_snr[idx] =
                row[layout.source.h_flux]
                / std::sqrt(row[layout.source.h_area]);
            data.delta_chi2[idx] = row[layout.source.delta_chi2];
            data.orth_ext[idx] = row[layout.source.orth_ext];
            data.rra[idx] = row[layout.source.ra];
            data.ddec[idx] = row[layout.source.dec];

            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                data.iexpo[idx] =
                    static_cast<int>(std::lround(row[layout.source.chi2]));
                data.snrf[idx] = row[layout.source.snr_f];
            }
            ++data.ng;
        }

        ndup = 0;
        return true;
    };

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (!parseCatalogRow(line, layout.all_columns, item)) continue;

        // Pixel coordinates for chip-edge masking
        int ix = static_cast<int>(item[layout.source.pixx]);  // pixel x
        int iy = static_cast<int>(item[layout.source.pixy]);  // pixel y

        // Bad CCD check
        int ccd_val = static_cast<int>(std::lround(item[layout.ccd]));
        bool bad_ccd = false;
        for (int i = 0; i < fc::n_bad_ccds; ++i) {
            if (ccd_val == fc::bad_ccds[i]) { bad_ccd = true; break; }
        }
        if (bad_ccd) continue;

        // Chip-edge masking (DES)
        if (ix < fc::chip_xmin || ix > fc::chip_xmax ||
            iy < fc::chip_ymin || iy > fc::chip_ymax)
            continue;

        // SNR_F cut
        if (item[layout.source.snr_f] < fc::snrfcut) continue;

        // SNR cut
        float snr = item[layout.source.h_flux]
                    / std::sqrt(item[layout.source.h_area]);
        if (fc::snrlow > 0.0 && snr < fc::snrlow) continue;
        if (fc::snrhigh > 0.0 && snr > fc::snrhigh) continue;

        // Half-light radius cut
        if (fc::r_half_thresh > 0.0) {
            const float r_half = static_cast<float>(
                std::sqrt(item[layout.source.h_area] / LensingConfig::pi)
                * pixel_size * 2.0);
            if (r_half <= fc::r_half_thresh * item[layout.source.psf]) continue;
        }

        // Magnitude cut
        if (item[magnitude_column] < fc::mag_min_val ||
            item[magnitude_column] > fc::mag_max_val) {
            continue;
        }

        // PSF polychi2 cut
        if (item[layout.source.polychi2] > fc::psf_chi2_mltp) continue;

        // Star classification cut
        if (item[layout.source.star] < fc::starcut) continue;

        // Chi2 cut (single mode only; per-exposure uses source.chi2 as exposure index)
        if (!fc::FD_PER_EXPOSURE_STAR_BAR) {
            if (item[layout.source.chi2] > fc::chi2_thresh) continue;
        }

        // Flag cut
        if (item[layout.source.flag] <= fc::flagcut) continue;

        // Imax / Jmax cut
        if (item[layout.source.imax] >= fc::imaxcut) continue;
        if (item[layout.source.jmax] >= fc::jmaxcut) continue;

        // Zero-point cut
        if (item[layout.external.zp] <= fc::zplow) continue;
        if (item[layout.external.zp] >= fc::zphigh) continue;

        // Field-distortion range cut
        if (std::fabs(item[layout.source.gf1]) > fc::gf_lim) continue;
        if (std::fabs(item[layout.source.gf2]) > fc::gf_lim) continue;

        // --- Duplicate detection ---
        bool new_galaxy = false;
        if (std::fabs(item[layout.source.dec] - last_dec) > multisam_thrsh ||
            std::fabs(item[layout.source.ra] - last_ra) > multisam_thrsh) {
            new_galaxy = true;
            last_dec = item[layout.source.dec];
            last_ra = item[layout.source.ra];
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
