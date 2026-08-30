#include "process_fd/ShearCatalogReader.hpp"
#include "FDConfig.hpp"
#include "LensingConfig.hpp"
#include "RuntimeConfig.hpp"
#include "general/NumericalRecipes.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
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
// Method: Group raw valid rows by external-catalog coordinates, uniformly
//         retain at most MAX_DUP measurements, then apply the existing FD cuts
//         and append accepted measurements with their serialized exposure ID.
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
    bool have_group = false;
    constexpr float multisam_thrsh = 1e-7f;

    // ==========================================
    // Function: Flush one duplicate group into the FD arrays
    // Method: Apply all existing FD cuts only after raw-group sampling, then
    //         append the accepted shear, weight, magnitude, diagnostic,
    //         coordinate, and exposure fields without changing their formulas.
    // ==========================================
    auto flushDuplicates = [&]() -> bool {
        const int count = ndup > fc::MAX_DUP ? fc::MAX_DUP : ndup;
        for (int i = 0; i < count; ++i) {
            const CatalogRow& row = dup_buf[i];

            int ix = static_cast<int>(row[layout.source.pixx]);
            int iy = static_cast<int>(row[layout.source.pixy]);
            int ccd_val = static_cast<int>(std::lround(row[layout.ccd]));
            bool bad_ccd = false;
            for (int j = 0; j < fc::n_bad_ccds; ++j) {
                if (ccd_val == fc::bad_ccds[j]) {
                    bad_ccd = true;
                    break;
                }
            }
            if (bad_ccd) continue;
            if (ix < fc::chip_xmin || ix > fc::chip_xmax ||
                iy < fc::chip_ymin || iy > fc::chip_ymax) continue;
            if (row[layout.source.snr_f] < fc::snrfcut) continue;

            float snr = row[layout.source.h_flux]
                        / std::sqrt(row[layout.source.h_area]);
            if (fc::snrlow > 0.0 && snr < fc::snrlow) continue;
            if (fc::snrhigh > 0.0 && snr > fc::snrhigh) continue;
            if (fc::r_half_thresh > 0.0) {
                float r_half = static_cast<float>(
                    std::sqrt(row[layout.source.h_area] / LensingConfig::pi)
                    * pixel_size * 2.0);
                if (r_half <= fc::r_half_thresh * row[layout.source.psf]) continue;
            }
            if (row[magnitude_column] < fc::mag_min_val ||
                row[magnitude_column] > fc::mag_max_val) continue;
            if (row[layout.source.polychi2] > fc::psf_chi2_mltp) continue;
            if (row[layout.source.star] < fc::starcut) continue;
            if (!fc::FD_PER_EXPOSURE_STAR_BAR &&
                row[layout.source.chi2] > fc::chi2_thresh) continue;
            if (row[layout.source.flag] <= fc::flagcut) continue;
            if (row[layout.source.imax] >= fc::imaxcut) continue;
            if (row[layout.source.jmax] >= fc::jmaxcut) continue;
            if (row[layout.external.zp] <= fc::zplow) continue;
            if (row[layout.external.zp] >= fc::zphigh) continue;
            if (std::fabs(row[layout.source.gf1]) > fc::gf_lim) continue;
            if (std::fabs(row[layout.source.gf2]) > fc::gf_lim) continue;

            int source_exposure = 0;
            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                const double exposure_value = row[layout.expo];
                if (exposure_value < 1.0f
                    || exposure_value
                           > static_cast<double>(
                                 std::numeric_limits<int>::max())) {
                    continue;
                }
                source_exposure = static_cast<int>(
                    std::lround(exposure_value));
                if (source_exposure < 1) continue;
            }

            const int idx = data.ng;
            if (idx >= fc::nmax_per_core) {
                std::cerr << "nmax_per_core is too small!" << std::endl;
                return false;
            }

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
                data.iexpo[idx] = source_exposure;
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

        // Duplicate grouping intentionally precedes every scientific cut.
        const float current_dec = item[layout.external.dec];
        const float current_ra = item[layout.external.ra];
        const bool new_galaxy =
            have_group
            && (std::fabs(current_dec - last_dec) > multisam_thrsh
                || std::fabs(current_ra - last_ra) > multisam_thrsh);

        // Flush previous duplicates when a new galaxy is detected
        if (new_galaxy && ndup > 0) {
            if (!flushDuplicates()) return;
        }

        if (!have_group || new_galaxy) {
            last_dec = current_dec;
            last_ra = current_ra;
            have_group = true;
        }

        // ==========================================
        // Logic: Uniformly retain MAX_DUP measurements from the raw group
        // Method: Apply reservoir sampling before any scientific cut, using
        //         the per-rank ran1 stream initialized once by main.cpp.
        // ==========================================
        if (ndup < fc::MAX_DUP) {
            dup_buf[ndup] = item;
        } else {
            const int reservoir_index = static_cast<int>(
                NumericalRecipes::ran1() * static_cast<double>(ndup + 1));
            if (reservoir_index < fc::MAX_DUP) {
                dup_buf[reservoir_index] = item;
            }
        }
        ++ndup;
    }

    // Flush remaining duplicates
    if (ndup > 0) {
        if (!flushDuplicates()) return;
    }
    file.close();
}
