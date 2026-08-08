#include "SourceExtractor.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "Astrometry.hpp"
#include "ExternalCatalogReader.hpp"
#include "ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <system_error>

// Global/extern variables representing exposure filenames list (defined in main.cpp)
extern std::vector<std::string> EXPO_FILE;
extern int N_EXPO;

namespace SourceExtractor {

    // ==========================================
    // Function: Check external catalog path availability
    // Method: Use non-throwing filesystem existence check before opening optional catalog tiles.
    // ==========================================
    namespace {
        bool catalogFileExists(const std::string& path) {
            std::error_code ec;
            return std::filesystem::exists(path, ec) && !ec;
        }

    }

    void procSource(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = EXPO_FILE[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);
        
        int nchip = static_cast<int>(image_files.size());
        for (int ichip = 1; ichip <= nchip; ++ichip) {
            chipProcessSource(image_files, ichip, dir_output);
        }
    }


    // ==========================================
    // Function: Process one chip for source and star catalog generation
    // Method: Reject a failed Stage 1 normalized-map sentinel before reading coefficients or
    //         constructing a sigma map, then follow the lite external-catalog branch.
    // ==========================================
    void chipProcessSource(const std::vector<std::string>& imageFiles, int ichip, const std::string& dirOutput) {
        int proc_error = 0;
        int nstar = 0;
        int ngal = 0;

        int nx = 0, ny = 0;
        std::vector<float> array;

        if (!FitsIO::readImage(imageFiles[ichip - 1], nx, ny, array)) {
            std::cerr << "Error reading image: " << imageFiles[ichip - 1] << std::endl;
            return;
        }

        std::string raw_prefix = UniversalUtils::getPrefix(imageFiles[ichip - 1]);
        std::string PREFIX = raw_prefix;
        std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/Norm", PREFIX, "_norm.fits");

        std::vector<float> normap;
        int norm_nx = 0, norm_ny = 0;
        if (!FitsIO::readImage(filename, norm_nx, norm_ny, normap)) {
            std::cerr << "Error reading normalized map: " << filename << std::endl;
            return;
        }
        const size_t expected_size = static_cast<size_t>(nx) * static_cast<size_t>(ny);
        if (nx <= LensingConfig::CCD_split || ny < 3
            || norm_nx != nx || norm_ny != ny || normap.size() < expected_size || normap.empty()
            || !std::isfinite(normap[0]) || normap[0] >= 0.0f || normap[0] < -99990.0f) {
            std::cerr << "Error / proc_source rejected failed norm chip "
                      << imageFiles[ichip - 1] << std::endl;
            return;
        }

        float sigabc[2][3] = {0};
        for (int i = 0; i < LensingConfig::CCD_split; ++i) {
            for (int j = 0; j < 3; ++j) {
                sigabc[i][j] = normap[j * nx + (i + 1)];
            }
        }

        std::vector<int> weight(nx * ny, 1);
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int idx = y * nx + x;
                if (normap[idx] < -900.0f) {
                    weight[idx] = 0;
                }
            }
        }

        int nxc = nx / 2;
        std::vector<float> sigmap(nx * ny, 0.0f);
        if (LensingConfig::CCD_split == 2) {
            for (int x = 0; x < nxc; ++x) {
                int ii = x + nxc;
                for (int y = 0; y < ny; ++y) {
                    double val1 = sigabc[0][0] + sigabc[0][1] * (x + 1) + sigabc[0][2] * (y + 1);
                    sigmap[y * nx + x] = std::sqrt(0.5 * val1);

                    double val2 = sigabc[1][0] + sigabc[1][1] * (ii + 1) + sigabc[1][2] * (y + 1);
                    sigmap[y * nx + ii] = std::sqrt(0.5 * val2);
                }
            }
        } else {
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    double val = sigabc[0][0] + sigabc[0][1] * (x + 1) + sigabc[0][2] * (y + 1);
                    sigmap[y * nx + x] = std::sqrt(0.5 * val);
                }
            }
        }

        getExpoCatalog(dirOutput, PREFIX, nx, ny, sigmap, weight, normap, proc_error);

        std::string PREFIX_head = UniversalUtils::getPrefixExpo(imageFiles[0]);
        filename = dirOutput + "/astrometry/Head/" + PREFIX_head + ".head";

        double cRPIX[2] = {0.0, 0.0};
        double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double cRVAL[2] = {0.0, 0.0};
        double PU[2][LensingConfig::npd] = {0};

        Astrometry::readAstrometryPara(filename, ichip, cRPIX, cD, cRVAL, PU, LensingConfig::npd, proc_error);

        std::string catfile = LensingConfig::SOURCE_CAT;
        std::vector<std::string> sortfile(27);
        int sortnum = 0;

        if (proc_error == 0) {
            generateGalCatFileName(cRVAL, catfile, sortfile, sortnum);
        }

        deBlending(sortfile, sortnum, nx, ny, weight, cRPIX, cD, cRVAL, PU, proc_error);

        genSourceExtCatalog(dirOutput, sortfile, sortnum, PREFIX, nx, ny, array, weight, sigmap, cRPIX, cD, cRVAL, PU, ngal, proc_error);

        genStarCandidateDirect(dirOutput, PREFIX, nx, ny, array, weight, nstar, proc_error);

        if (proc_error != 0) {
            std::cout << "Error / proc_source " << imageFiles[ichip - 1] << " " << proc_error << " " << nstar << " " << ngal << std::endl;
        }
    }

    // ==========================================
    // Function: Deblend external-catalog sources by redshift consistency
    // Method: Read only configured RA, Dec, and ZP fields, mirror F77 missing-file behavior,
    //         and reject malformed or out-of-range photometric-redshift rows.
    // ==========================================
    void deBlending(const std::vector<std::string>& sortFile, int sortNum, int nx, int ny, std::vector<int>& weight,
                    const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                    const double PU[2][LensingConfig::npd], int& procError) {
        if (procError == 1) return;

        double ra_c = 0.0;
        double dec_bound[2] = {0.0, 0.0};
        double dra = 0.0;
        double astrometry_shift_ratio = 0.2;

        Astrometry::getRaDecRangeFine(nx, ny, ra_c, dec_bound, dra, cRPIX, cD, cRVAL, PU, LensingConfig::npd, astrometry_shift_ratio);

        for (int n = 0; n < sortNum; ++n) {
            if (!catalogFileExists(sortFile[n])) {
                continue;
            }

            std::ifstream fin(sortFile[n]);
            if (!fin.is_open()) {
                continue;
            }

            std::string header;
            std::getline(fin, header); // skip header line

            std::string line;
            while (std::getline(fin, line)) {
                if (line.empty()) continue;

                ExternalCatalogReader::Record record;
                if (!ExternalCatalogReader::parseRecord(line, record)) {
                    continue;
                }

                double ra = record.ra;
                double dec = record.dec;
                const double zp = record.zp;
                if (zp < 0.0 || zp > 5.0) continue;

                double z = zp;
                if (std::abs(Astrometry::diffra(ra, ra_c)) > dra * 0.5) continue;
                if (dec < dec_bound[0] || dec > dec_bound[1]) continue;

                double xx = 0.0, yy = 0.0;
                Astrometry::coordinateTransferPU(ra, dec, xx, yy, -1, cRPIX, cD, cRVAL, PU, LensingConfig::npd);

                int ix = static_cast<int>(xx + 0.5);
                int iy = static_cast<int>(yy + 0.5);

                int ix_0 = ix - 1;
                int iy_0 = iy - 1;

                if (ix_0 < 0 || ix_0 >= nx || iy_0 < 0 || iy_0 >= ny) continue;

                int idx = iy_0 * nx + ix_0;
                int val = weight[idx];

                if (val < 2) {
                    continue;
                } else if (val > 2) {
                    double z1 = 0.001 * (val - 10.0);
                    if (std::abs(z - z1) > LensingConfig::dz_thresh) {
                        int old_w = val;
                        int new_w = 0;
                        fillPatch(nx, ny, weight, ix_0, iy_0, old_w, new_w);
                    }
                } else {
                    int old_w = 2;
                    int new_w = static_cast<int>(z * 1000.0 + 10.0);
                    fillPatch(nx, ny, weight, ix_0, iy_0, old_w, new_w);
                }
            }
            fin.close();
        }

        for (int i = 0; i < nx * ny; ++i) {
            if (weight[i] > 2) {
                weight[i] = 2;
            }
        }
    }

    void fillPatch(int nx, int ny, std::vector<int>& map, int ix, int iy, int old_v, int new_v) {
        if (old_v == new_v) return;

        std::vector<std::pair<int, int>> q;
        q.push_back({ix, iy});
        map[iy * nx + ix] = new_v;

        size_t head = 0;
        while (head < q.size()) {
            auto [cx, cy] = q[head++];
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    int nx_x = cx + dx;
                    int ny_y = cy + dy;
                    if (nx_x >= 0 && nx_x < nx && ny_y >= 0 && ny_y < ny) {
                        int idx = ny_y * nx + nx_x;
                        if (map[idx] == old_v) {
                            map[idx] = new_v;
                            q.push_back({nx_x, ny_y});
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Build exposure source catalog
    // Method: Match F77 get_expo_catalog scan and neighbor traversal order.
    // ==========================================
    void getExpoCatalog(const std::string& dirOutput, const std::string& prefix, int nx, int ny, const std::vector<float>& sigmap,
                        std::vector<int>& weight, const std::vector<float>& normap, int& ierror) {
        if (ierror == 1) return;

        std::vector<int> mark(nx * ny, 0);
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (!normap.empty() && normap[idx] >= LensingConfig::source_thresh && weight[idx] >= 1) {
                    mark[idx] = 1;
                } else {
                    mark[idx] = 0;
                }
            }
        }

        std::string catname = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, ".cat");
        MainIO::OutputFile fout(catname);
        if (!fout.is_open()) {
            std::cerr << "Error: could not open catalog file for writing: " << catname << std::endl;
            ierror = 1;
            return;
        }

        fout << " xp yp total_area half_light_area sig sig_normed_total_flux sig_normed_half_light_flux sig_normed_peak rmax\n";

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int start_idx = y * nx + x;
                if (mark[start_idx] == 1) {
                    std::vector<std::pair<int, int>> buffer;
                    buffer.push_back({x, y});
                    mark[start_idx] = -1;
                    
                    bool toobig = false;
                    size_t head = 0;
                    
                    while (head < buffer.size() && !toobig) {
                        auto [cx, cy] = buffer[head++];
                        for (int dx = -1; dx <= 1 && !toobig; ++dx) {
                            for (int dy = -1; dy <= 1 && !toobig; ++dy) {
                                if (dx == 0 && dy == 0) continue;
                                int nx_x = cx + dx;
                                int ny_y = cy + dy;
                                if (nx_x >= 0 && nx_x < nx && ny_y >= 0 && ny_y < ny) {
                                    int n_idx = ny_y * nx + nx_x;
                                    if (mark[n_idx] == 1) {
                                        buffer.push_back({nx_x, ny_y});
                                        mark[n_idx] = -1;
                                        if (buffer.size() == static_cast<size_t>(LensingConfig::area_max)) {
                                            toobig = true;
                                        }
                                    } else if (mark[n_idx] > 1) {
                                        toobig = true;
                                    }
                                }
                            }
                        }
                    }
                    
                    int nb = static_cast<int>(buffer.size());
                    if (toobig) {
                        for (const auto& p : buffer) {
                            int idx = p.second * nx + p.first;
                            mark[idx] = LensingConfig::area_max;
                            weight[idx] = 2;
                        }
                    } else {
                        if (nb >= LensingConfig::area_thresh) {
                            double xc = 0.0;
                            double yc = 0.0;
                            double total_flux = 0.0;
                            double sig = 0.0;
                            double peak = -100000.0;
                            double xp = 0.0;
                            double yp = 0.0;
                            
                            for (const auto& p : buffer) {
                                int idx = p.second * nx + p.first;
                                mark[idx] = nb;
                                weight[idx] = 2;
                                
                                double ix_1 = p.first + 1;
                                double iy_1 = p.second + 1;
                                xc += ix_1 * normap[idx];
                                yc += iy_1 * normap[idx];
                                total_flux += normap[idx];
                                sig += sigmap[idx];
                                
                                if (normap[idx] > peak) {
                                    peak = normap[idx];
                                    xp = ix_1;
                                    yp = iy_1;
                                }
                            }
                            
                            xc /= total_flux;
                            yc /= total_flux;
                            int total_area = nb;
                            sig /= total_area;
                            
                            double thresh = peak * 0.5;
                            double rmax = 0.0;
                            int half_light_area = 0;
                            double half_light_flux = 0.0;
                            
                            for (const auto& p : buffer) {
                                int idx = p.second * nx + p.first;
                                double ix_1 = p.first + 1;
                                double iy_1 = p.second + 1;
                                double r2 = (ix_1 - xc) * (ix_1 - xc) + (iy_1 - yc) * (iy_1 - yc);
                                rmax = std::max(rmax, r2);
                                
                                if (normap[idx] >= thresh) {
                                    half_light_area++;
                                    half_light_flux += normap[idx];
                                }
                            }
                            
                            rmax = std::sqrt(rmax);
                            if (peak >= LensingConfig::core_thresh) {
                                fout << std::fixed << std::setprecision(10)
                                     << xp << " " << yp << " " << total_area << " " << half_light_area << " "
                                     << sig << " " << total_flux << " " << half_light_flux << " " << peak << " " << rmax << "\n";
                            }
                        } else {
                            for (const auto& p : buffer) {
                                int idx = p.second * nx + p.first;
                                mark[idx] = nb;
                            }
                        }
                    }
                }
            }
        }
        fout.close();
    }


    // ==========================================
    // Function: Extract source stamps from external catalog positions
    // Method: Read configured RA, Dec, and ZP fields while skipping missing tiles and malformed
    //         rows, preserving valid-row extraction and each accepted original catalog row.
    // ==========================================
    void genSourceExtCatalog(const std::string& dirOutput, const std::vector<std::string>& sortFile, int sortNum, const std::string& prefix,
                             int nx, int ny, const std::vector<float>& array, std::vector<int>& weight,
                             const std::vector<float>& sigmap, const double cRPIX[2], const double cD[2][2],
                             const double cRVAL[2], const double PU[2][LensingConfig::npd], int& ngal, int& procError) {
        ngal = 0;

        double ra_c = 0.0;
        double dec_bound[2] = {0.0, 0.0};
        double dra = 0.0;
        double astrometry_shift_ratio = 0.2;

        std::vector<float> source_collect(LensingConfig::ngal_max * LensingConfig::ns * LensingConfig::ns, 0.0f);
        std::vector<float> noise_collect(LensingConfig::ngal_max * LensingConfig::ns * LensingConfig::ns, 0.0f);
        std::vector<std::vector<float>> source_para(LensingConfig::ngal_max, std::vector<float>(LensingConfig::npara, 0.0f));
        std::vector<int> sid(LensingConfig::ngal_max, 0);
        std::vector<std::string> accepted_orig_lines;
        std::string orig_header = "";

        int ig = 0;
        bool header_captured = false;

        if (procError != 1) {
            Astrometry::getRaDecRangeFine(nx, ny, ra_c, dec_bound, dra, cRPIX, cD, cRVAL, PU, LensingConfig::npd, astrometry_shift_ratio);
        }

        for (int n = 0; procError != 1 && n < sortNum; ++n) {
            if (!catalogFileExists(sortFile[n])) {
                continue;
            }

            std::ifstream fin(sortFile[n]);
            if (!fin.is_open()) continue;

            std::string line;
            if (std::getline(fin, line)) {
                if (!header_captured) {
                    orig_header = line;
                    header_captured = true;
                }
            }

            while (std::getline(fin, line)) {
                if (line.empty()) continue;
                ExternalCatalogReader::Record record;
                if (!ExternalCatalogReader::parseRecord(line, record)) {
                    continue;
                }

                double ra = record.ra;
                double dec = record.dec;

                ig++;

                if (std::abs(Astrometry::diffra(ra, ra_c)) > dra * 0.5) continue;
                if (dec < dec_bound[0] || dec > dec_bound[1]) continue;

                double xx = 0.0, yy = 0.0;
                Astrometry::coordinateTransferPU(ra, dec, xx, yy, -1, cRPIX, cD, cRVAL, PU, LensingConfig::npd);

                double xp = xx;
                double yp = yy;

                if (xp - LensingConfig::nl_2 < LensingConfig::chip_margin || xp + LensingConfig::nl_2 > nx - LensingConfig::chip_margin ||
                    yp - LensingConfig::nl_2 < LensingConfig::chip_margin || yp + LensingConfig::nl_2 > ny - LensingConfig::chip_margin) {
                    continue;
                }

                int xp_idx = static_cast<int>(xp + 0.5) - 1;
                int yp_idx = static_cast<int>(yp + 0.5) - 1;
                if (xp_idx < 0 || xp_idx >= nx || yp_idx < 0 || yp_idx >= ny) continue;

                double sig = sigmap[yp_idx * nx + xp_idx];

                int flag = 0;
                std::vector<float> noise(LensingConfig::ns * LensingConfig::ns, 0.0f);
                int imax = 0, jmax = 0;

                findNoise(flag, noise, nx, ny, array, weight, xp, yp, sig, imax, jmax);
                if (flag < 0) continue;

                double peak = 0.0, half_light_flux = 0.0;
                int half_light_area = 0;
                std::vector<float> source(LensingConfig::ns * LensingConfig::ns, 0.0f);
                checkSource(flag, source, nx, ny, array, weight, xp, yp, sig, imax, jmax, peak, half_light_flux, half_light_area);
                if (flag < 0) continue;

                ngal++;

                int offset = (ngal - 1) * LensingConfig::ns * LensingConfig::ns;
                std::copy(source.begin(), source.end(), source_collect.begin() + offset);
                std::copy(noise.begin(), noise.end(), noise_collect.begin() + offset);

                source_para[ngal - 1][0] = ig;
                source_para[ngal - 1][1] = xp;
                source_para[ngal - 1][2] = yp;
                source_para[ngal - 1][3] = sig;
                source_para[ngal - 1][4] = peak;
                source_para[ngal - 1][5] = imax;
                source_para[ngal - 1][6] = jmax;
                source_para[ngal - 1][7] = half_light_flux;
                source_para[ngal - 1][8] = half_light_area;
                source_para[ngal - 1][9] = flag;

                sid[ngal - 1] = ig;
                accepted_orig_lines.push_back(line);

                if (ngal >= LensingConfig::ngal_max) {
                    fin.close();
                    break;
                }
            }
            fin.close();
            if (ngal >= LensingConfig::ngal_max) break;
        }

        if (ngal > 0) {
            int nn1 = LensingConfig::ns * LensingConfig::len_g;
            int nn2 = LensingConfig::ns * (ngal / LensingConfig::len_g + 1);
            
            std::string filename_src = OutputLayout::chipPath(
                dirOutput, "stamps/fits_Src", prefix, "_source.fits");
            FitsIO::writeStamps(LensingConfig::ngal_max, 1, ngal, LensingConfig::ns, LensingConfig::ns, source_collect, nn1, nn2, filename_src);

            std::string filename_noise = OutputLayout::chipPath(
                dirOutput, "stamps/fits_Noise", prefix, "_noise.fits");
            FitsIO::writeStamps(LensingConfig::ngal_max, 1, ngal, LensingConfig::ns, LensingConfig::ns, noise_collect, nn1, nn2, filename_noise);
        }

        std::string filename_info = OutputLayout::chipPath(
            dirOutput, "stamps/dat_SrcInfo", prefix, "_source_info.dat");
        MainIO::OutputFile fout(filename_info);
        if (fout.is_open()) {
            fout << "ig xp yp sigma peak imax jmax half_light_flux half_light_area flag\n";
            for (int i = 0; i < ngal; ++i) {
                for (int j = 0; j < LensingConfig::iflag + 1; ++j) {
                    fout << source_para[i][j] << (j == LensingConfig::iflag ? "" : " ");
                }
                fout << "\n";
            }
            fout.close();
        }

        std::string filename_orig = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
        MainIO::OutputFile fout_orig(filename_orig);
        if (fout_orig.is_open()) {
            if (procError == 1 || ngal == 0) {
                fout_orig << "No sources!!\n";
            } else {
                fout_orig << orig_header << "\n";
                for (const auto& line : accepted_orig_lines) {
                    fout_orig << line << "\n";
                }
            }
            fout_orig.close();
        }
    }

    void findNoise(int& flag, std::vector<float>& stamps, int nx, int ny, const std::vector<float>& array,
                   const std::vector<int>& weight, double xp, double yp, double sig, int& imax, int& jmax) {
        flag = 0;

        int cx_0 = static_cast<int>(xp + 0.5) - 1;
        int cy_0 = static_cast<int>(yp + 0.5) - 1;
        int x1_0 = cx_0 - LensingConfig::nl_2 - 1;
        int y1_0 = cy_0 - LensingConfig::nl_2 - 1;

        double temp = 100000.0;
        double pmax = temp;
        int jx = -1, jy = -1;

        int nsteph = 5;
        for (int i = 0; i < nsteph; ++i) {
            for (int j = 0; j < nsteph; ++j) {
                if (i == 0 || i == nsteph - 1 || j == 0 || j == nsteph - 1) {
                    int ix1_0 = x1_0 + LensingConfig::nl * (i - 2);
                    int iy1_0 = y1_0 + LensingConfig::nl * (j - 2);
                    int ix2_0 = ix1_0 + LensingConfig::nl - 1;
                    int iy2_0 = iy1_0 + LensingConfig::nl - 1;

                    if (ix1_0 < LensingConfig::chip_edge_margin - 1 || ix2_0 > nx - 1 - LensingConfig::chip_edge_margin ||
                        iy1_0 < LensingConfig::chip_edge_margin - 1 || iy2_0 > ny - 1 - LensingConfig::chip_edge_margin) {
                        continue;
                    }

                    double npeak = -temp;
                    int nd = 0;
                    for (int u = ix1_0; u <= ix2_0; ++u) {
                        for (int v = iy1_0; v <= iy2_0; ++v) {
                            int idx = v * nx + u;
                            if (array[idx] > npeak) {
                                npeak = array[idx];
                            }
                            if (weight[idx] == 0) {
                                nd++;
                            }
                        }
                    }

                    if (nd <= LensingConfig::nl && npeak < pmax) {
                        pmax = npeak;
                        jx = ix1_0;
                        jy = iy1_0;
                    }
                }
            }
        }

        if (pmax == temp) {
            flag = -1;
            return;
        }

        std::vector<float> stampl(LensingConfig::nl * LensingConfig::nl, 0.0f);
        std::vector<int> weightl(LensingConfig::nl * LensingConfig::nl, 0);

        for (int v_idx = 0; v_idx < LensingConfig::nl; ++v_idx) {
            int j = jy + v_idx;
            for (int u_idx = 0; u_idx < LensingConfig::nl; ++u_idx) {
                int i = jx + u_idx;
                int idx_src = j * nx + i;
                int idx_dest = v_idx * LensingConfig::nl + u_idx;
                stampl[idx_dest] = array[idx_src];
                weightl[idx_dest] = weight[idx_src];
            }
        }

        ImageProcessing::flattenStamp2D(LensingConfig::ns, LensingConfig::nl, stampl, weightl, flag);
        if (flag < 0) return;

        ImageProcessing::markNoise(LensingConfig::nl, stampl, weightl, sig, LensingConfig::source_thresh, LensingConfig::core_thresh);

        std::vector<int> weights(LensingConfig::ns * LensingConfig::ns, 0);
        int offset_xy = LensingConfig::nl_2 - LensingConfig::ns_2;
        for (int y = 0; y < LensingConfig::ns; ++y) {
            for (int x = 0; x < LensingConfig::ns; ++x) {
                stamps[y * LensingConfig::ns + x] = stampl[(y + offset_xy) * LensingConfig::nl + (x + offset_xy)];
                weights[y * LensingConfig::ns + x] = weightl[(y + offset_xy) * LensingConfig::nl + (x + offset_xy)];
            }
        }

        ImageProcessing::decorateStamp(LensingConfig::ns, sig, weights, stamps);
    }

    // ==========================================
    // Function: Extract and validate a source stamp
    // Method: Mirror F77 check_source defect counting and stamp decoration.
    // ==========================================
    void checkSource(int& flag, std::vector<float>& stamps, int nx, int ny, const std::vector<float>& array,
                     const std::vector<int>& weight, double xp, double yp, double sig, int& imax, int& jmax,
                     double& peak, double& half_light_flux, int& half_light_area) {
        flag = 0;

        int cx_0 = static_cast<int>(xp + 0.5) - 1;
        int cy_0 = static_cast<int>(yp + 0.5) - 1;

        if (cx_0 - LensingConfig::nl_2 < LensingConfig::chip_edge_margin - 1 ||
            cx_0 + LensingConfig::nl_2 > nx - LensingConfig::chip_edge_margin - 1 ||
            cy_0 - LensingConfig::nl_2 < LensingConfig::chip_edge_margin - 1 ||
            cy_0 + LensingConfig::nl_2 > ny - LensingConfig::chip_edge_margin - 1) {
            flag = -1;
            return;
        }

        int x1_0 = cx_0 - LensingConfig::nl_2;
        int y1_0 = cy_0 - LensingConfig::nl_2;

        std::vector<float> stampl(LensingConfig::nl * LensingConfig::nl, 0.0f);
        std::vector<int> weightl(LensingConfig::nl * LensingConfig::nl, 0);

        for (int v = 0; v < LensingConfig::nl; ++v) {
            for (int u = 0; u < LensingConfig::nl; ++u) {
                int idx_src = (y1_0 + v) * nx + (x1_0 + u);
                int idx_dest = v * LensingConfig::nl + u;
                stampl[idx_dest] = array[idx_src];
                weightl[idx_dest] = weight[idx_src];
            }
        }

        ImageProcessing::flattenStamp2D(LensingConfig::ns, LensingConfig::nl, stampl, weightl, flag);
        if (flag < 0) return;

        int boundx[2] = {0, 0};
        int boundy[2] = {0, 0};
        double total_flux = 0.0;
        int total_area = 0;
        double radius = 0.0;
        int xcenter = 0;
        int ycenter = 0;

        ImageProcessing::markSource(LensingConfig::nl, stampl, weightl, sig, LensingConfig::source_thresh, LensingConfig::core_thresh,
                                    boundx, boundy, total_flux, total_area, peak, half_light_flux, half_light_area, flag,
                                    radius, xcenter, ycenter);

        if (flag < 0) return;

        if (peak > LensingConfig::saturation_thresh / sig) {
            flag = -1;
            return;
        }

        double temp = LensingConfig::ns_2 - LensingConfig::flag_thresh;
        if (radius >= temp) {
            flag = -1;
            return;
        }

        int x1_cut_0 = xcenter - LensingConfig::ns_2;
        int y1_cut_0 = ycenter - LensingConfig::ns_2;

        if (x1_cut_0 < 0 || x1_cut_0 + LensingConfig::ns > LensingConfig::nl ||
            y1_cut_0 < 0 || y1_cut_0 + LensingConfig::ns > LensingConfig::nl) {
            flag = -1;
            return;
        }

        std::vector<int> weights(LensingConfig::ns * LensingConfig::ns, 0);

        for (int y_idx = 0; y_idx < LensingConfig::ns; ++y_idx) {
            for (int x_idx = 0; x_idx < LensingConfig::ns; ++x_idx) {
                int idx_dest = y_idx * LensingConfig::ns + x_idx;
                int idx_src = (y1_cut_0 + y_idx) * LensingConfig::nl + (x1_cut_0 + x_idx);
                stamps[idx_dest] = stampl[idx_src];
                weights[idx_dest] = weightl[idx_src];
            }
        }

        imax = 0;
        for (int x_idx = 0; x_idx < LensingConfig::ns; ++x_idx) {
            int u = 0;
            for (int y_idx = 0; y_idx < LensingConfig::ns; ++y_idx) {
                int idx_dest = y_idx * LensingConfig::ns + x_idx;
                if (weights[idx_dest] == 0) {
                    u++;
                }
            }
            if (u > imax) {
                imax = u;
            }
        }

        jmax = 0;
        for (int y_idx = 0; y_idx < LensingConfig::ns; ++y_idx) {
            int u = 0;
            for (int x_idx = 0; x_idx < LensingConfig::ns; ++x_idx) {
                int idx_dest = y_idx * LensingConfig::ns + x_idx;
                if (weights[idx_dest] == 0) {
                    u++;
                }
            }
            if (u > jmax) {
                jmax = u;
            }
        }

        ImageProcessing::decorateStamp(LensingConfig::ns, sig, weights, stamps);
    }


    // ==========================================
    // Function: Publish star candidates selected directly from detections
    // Method: Re-extract qualifying detection stamps and route candidate text
    //         and FITS products through checked main-process writers.
    // ==========================================
    void genStarCandidateDirect(const std::string& dirOutput, const std::string& prefix, int nx, int ny, const std::vector<float>& array,
                                const std::vector<int>& weight, int& nstar, int& procError) {
        nstar = 0;

        std::vector<float> star_source_collect(LensingConfig::ngal_max * LensingConfig::ns * LensingConfig::ns, 0.0f);
        std::vector<float> star_noise_collect(LensingConfig::ngal_max * LensingConfig::ns * LensingConfig::ns, 0.0f);
        std::vector<std::vector<float>> star_para;

        if (procError != 1) {
            std::string catname = OutputLayout::chipPath(
                dirOutput, "stamps/cat_Orig", prefix, ".cat");
            std::ifstream fin(catname);
            if (!fin.is_open()) {
                std::cerr << catname << std::endl;
                std::cerr << "Error / gen_star_candidate_direct catalog file error!!" << std::endl;
                std::exit(1);
            }

            std::string header;
            std::getline(fin, header); // skip header

            double xp, yp, sig, total_flux, half_light_flux, peak, rf;
            int total_area, half_light_area;

            while (fin >> xp >> yp >> total_area >> half_light_area >> sig >> total_flux >> half_light_flux >> peak >> rf) {
                double temp = half_light_area;
                double snr = half_light_flux / std::sqrt(temp);
                if (snr < LensingConfig::SNR_PSF * 0.5) continue;

                int flag = 0;
                std::vector<float> noise(LensingConfig::ns * LensingConfig::ns, 0.0f);
                int imax = 0, jmax = 0;

                findNoise(flag, noise, nx, ny, array, weight, xp, yp, sig, imax, jmax);
                if (flag < 0) continue;

                std::vector<float> source(LensingConfig::ns * LensingConfig::ns, 0.0f);
                checkSource(flag, source, nx, ny, array, weight, xp, yp, sig, imax, jmax, peak, half_light_flux, half_light_area);
                if (flag < 0) continue;

                temp = half_light_area;
                snr = half_light_flux / std::sqrt(temp);

                if (snr < LensingConfig::SNR_PSF) continue;

                nstar++;
                int offset_dest = (nstar - 1) * LensingConfig::ns * LensingConfig::ns;
                std::copy(source.begin(), source.end(), star_source_collect.begin() + offset_dest);
                std::copy(noise.begin(), noise.end(), star_noise_collect.begin() + offset_dest);

                std::vector<float> row(4, 0.0f);
                row[0] = nstar;
                row[1] = xp;
                row[2] = yp;
                row[3] = snr;
                star_para.push_back(row);

                if (nstar >= LensingConfig::nstar_max) break;
            }
            fin.close();
        }

        std::string filename_star_info = OutputLayout::chipPath(
            dirOutput, "stamps/dat_StarCanInfo", prefix, "_star_can_info.dat");
        MainIO::OutputFile fout(filename_star_info);
        if (fout.is_open()) {
            fout << "ig xp yp SNR\n";
            for (int i = 0; i < nstar; ++i) {
                fout << star_para[i][0] << " " << star_para[i][1] << " " << star_para[i][2] << " " << star_para[i][3] << "\n";
            }
            fout.close();
        }

        if (nstar > 0) {
            int nn1_s = LensingConfig::ns * LensingConfig::len_s;
            int nn2_s = LensingConfig::ns * (nstar / LensingConfig::len_s + 1);

            std::string filename_star_src = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCan", prefix, "_star_can.fits");
            FitsIO::writeStamps(LensingConfig::ngal_max, 1, nstar, LensingConfig::ns, LensingConfig::ns, star_source_collect, nn1_s, nn2_s, filename_star_src);

            std::string filename_star_noise = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanN", prefix, "_star_can_noise.fits");
            FitsIO::writeStamps(LensingConfig::ngal_max, 1, nstar, LensingConfig::ns, LensingConfig::ns, star_noise_collect, nn1_s, nn2_s, filename_star_noise);
        }
    }

    void generateGalCatFileName(const double cRVAL[2], std::string& filename, std::vector<std::string>& sortfile, int& sortnum) {
        std::vector<std::string> filenames = UniversalUtils::generateGalCatFileNames(filename, cRVAL);
        sortnum = static_cast<int>(filenames.size());
        sortfile.resize(filenames.size());
        for (size_t i = 0; i < filenames.size(); ++i) {
            sortfile[i] = filenames[i];
        }
    }
}
