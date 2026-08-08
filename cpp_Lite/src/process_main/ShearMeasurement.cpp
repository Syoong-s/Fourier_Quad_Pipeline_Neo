#include "ShearMeasurement.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "FitsIO.hpp"
#include "Astrometry.hpp"
#include "PSFModel.hpp"
#include "ExStar.hpp"
#include "ImageProcessing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iomanip>

extern std::vector<std::string> EXPO_FILE;

namespace ShearMeasurement {

void getWindowMinK(int ns, const std::vector<float>& psf_model, float thresh, float& k_win) {
    int c_pix = ns / 2;
    float k_min = 1e10f;
    for (int i = 0; i < ns; ++i) {
        float kx = static_cast<float>(i - c_pix);
        for (int j = 0; j < ns; ++j) {
            float ky = static_cast<float>(j - c_pix);
            if (psf_model[i * ns + j] > thresh) continue;
            float temp = kx * kx + ky * ky;
            if (temp < k_min) {
                k_min = temp;
            }
        }
    }
    if (k_min > 0.0f) {
        k_win = std::sqrt(k_min);
    } else {
        k_win = 0.0f;
    }
}

// ==========================================
// Function: Find alternate PSF low-k cutoff radius.
// Method: Mirrors F77 get_window_min_k_ver2 radial-bin stddev rule.
// ==========================================
void getWindowMinKVer2(int ns, const std::vector<float>& psf_model, float thresh, float& k_win) {
    int c_pix = ns / 2;
    int nbins = static_cast<int>(ns * std::sqrt(2.0) / 0.5) + 2;

    std::vector<int> cnt(nbins, 0);
    std::vector<double> sum_log(nbins, 0.0);
    std::vector<double> sum_sq(nbins, 0.0);
    std::vector<double> stddev(nbins, 0.0);
    double k_dent = 1e5;

    for (int i = 0; i < ns; ++i) {
        double kx = static_cast<double>(i - c_pix);
        for (int j = 0; j < ns; ++j) {
            double ky = static_cast<double>(j - c_pix);
            double kval = std::sqrt(kx * kx + ky * ky);
            float val = psf_model[i * ns + j];
            if (val <= 0.0f) {
                if (kval < k_dent) {
                    k_dent = kval;
                }
                continue;
            }
            int bin_idx_cpp = static_cast<int>(kval / 0.5);
            if (bin_idx_cpp < 0 || bin_idx_cpp >= nbins) continue;

            double logval = std::log(static_cast<double>(val));
            cnt[bin_idx_cpp]++;
            sum_log[bin_idx_cpp] += logval;
            sum_sq[bin_idx_cpp] += logval * logval;
        }
    }

    for (int b = 0; b < nbins; ++b) {
        if (cnt[b] > 1) {
            double var = sum_sq[b] / cnt[b] - std::pow(sum_log[b] / cnt[b], 2.0);
            stddev[b] = std::sqrt(std::max(0.0, var));
        }
    }

    double a = 0.0;
    int n_avg = 0;
    for (int b = 0; b < nbins; ++b) {
        double bin_k = b * 0.5;
        if (bin_k >= 10.0) break;
        if (cnt[b] > 1 && stddev[b] > 0.0) {
            a += stddev[b];
            n_avg++;
        }
    }

    if (n_avg > 0) {
        a /= n_avg;
    } else {
        k_win = 0.0f;
        return;
    }

    k_win = 0.0f;
    for (int b = 0; b < nbins; ++b) {
        double bin_k = b * 0.5;
        if (bin_k < 10.0) continue;
        if (cnt[b] <= 1) continue;
        if (stddev[b] > static_cast<double>(thresh) * a) {
            k_win = static_cast<float>(bin_k);
            return;
        }
    }

    if (k_dent < k_win) {
        k_win = static_cast<float>(k_dent);
    }
}

void getShear(int n, const float* gal, const float* psf, float& g1, float& g2, float& de, float& h1, float& h2) {
    float peak = psf[0];
    for (int i = 0; i < n * n; ++i) {
        if (psf[i] > peak) {
            peak = psf[i];
        }
    }

    float thresh = std::exp(-1.0f) * peak;

    float area = 0.0f;
    for (int i = 0; i < n * n; ++i) {
        if (psf[i] >= thresh) {
            area += 1.0f;
        }
    }

    float ks = std::sqrt(area / 3.1415926f);
    float PSFr_ratio = 0.75f;
    float ks_2 = std::pow(ks * PSFr_ratio, -2.0f);

    thresh = peak * 1e-4f;
    float r_win = 0.0f;

    std::vector<float> psf_vec(psf, psf + n * n);
    getWindowMinK(n, psf_vec, thresh, r_win);

    int cc = n / 2;

    g1 = 0.0f;
    g2 = 0.0f;
    de = 0.0f;
    h1 = 0.0f;
    h2 = 0.0f;

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            float kx = static_cast<float>(c - cc);
            float kx2 = kx * kx;
            float ky = static_cast<float>(r - cc);
            float ky2 = ky * ky;
            float k2 = kx2 + ky2;
            float k = std::sqrt(k2);
            if (k < r_win) {
                int idx = r * n + c;
                float ff = k2 * ks_2;
                float temp = std::exp(-ff) / psf[idx];
                float temp1 = temp * gal[idx];

                g1 -= temp1 * (kx2 - ky2);
                g2 -= temp1 * 2.0f * kx * ky;
                de += temp1 * k2 * (2.0f - ff);
                h1 += temp1 * ks_2 * (k2 * k2 - 8.0f * kx2 * ky2);
                h2 += temp1 * ks_2 * 4.0f * kx * ky * (kx2 - ky2);
            }
        }
    }
}

void getPSFArea(const float* model, float& FWHM) {
    int ns = LensingConfig::ns;
    float model_center = model[(ns / 2) * ns + (ns / 2)];
    float thresh = std::exp(-1.0f) * model_center;

    float area = -1e-5f;
    for (int i = 0; i < ns * ns; ++i) {
        if (model[i] >= thresh) {
            area += 1.0f;
        }
    }
    if (area <= 0.0f) {
        FWHM = -1.0f;
        return;
    }
    float r = std::sqrt(area / LensingConfig::pi);
    float beta = ns / (2.0f * LensingConfig::pi) / r;
    FWHM = beta * 2.0f * std::sqrt(2.0f * std::log(2.0f)) * 0.2628f;
}

// ==========================================
// Function: Measure shear for every source in one exposure.
// Method: Mirrors F77 expo_shear chip loop, PSF model evaluation, FQ estimator, and catalog write.
// ==========================================
void expoShear(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int chipnx, int chipny) {
    int ns = LensingConfig::ns;
    int npl = LensingConfig::npl;

    int proc_error = 0;
    int nstar = 0;
    int status = 0;
    float poly_ave = 0.0f;
    float poly_std = 0.0f;

    std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    std::string headname = dirOutput + "/astrometry/Head/" + prefix_expo + ".head";

    std::vector<double> local_coe(ns * ns * (npl + 1), 0.0);

    for (int ichip = 0; ichip < nchip; ++ichip) {
        proc_error = 0;
        std::string PREFIX = UniversalUtils::getPrefix(imageFiles[ichip]);
        {
            std::string filename = OutputLayout::chipPath(
                dirOutput, "stamps/dat_PsfFit", PREFIX, "_PSF_coe_local.dat");
            std::ifstream fin(filename);
            if (!fin.is_open()) {
                proc_error = 1;
            } else {
                if (fin >> nstar >> status >> poly_ave >> poly_std) {
                    if (status == -1) {
                        proc_error = 1;
                    } else {
                        for (int i = 0; i < ns; ++i) {
                            for (int j = 0; j < ns; ++j) {
                                for (int k = 0; k < npl + 1; ++k) {
                                    if (!(fin >> local_coe[(j * ns + i) * (npl + 1) + k])) {
                                        proc_error = 1;
                                        break;
                                    }
                                }
                                if (proc_error) break;
                            }
                            if (proc_error) break;
                        }
                    }
                } else {
                    proc_error = 1;
                }
                fin.close();
            }
        }

        std::string out_filename = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", PREFIX, "_shear.dat");
        MainIO::OutputFile fout10;

        auto write_empty_output = [&]() {
            fout10.open(out_filename);
            fout10 << std::setprecision(10);
            fout10 << "poly_chi2 xc yc sigma nstar imax jmax "
                   << "half_light_flux half_light_area flag psf_FWHM SNR_F "
                   << "ra dec gf1 gf2 g1 g2 de h1 h2 cos2 sin2 parity\n";
            fout10.close();
        };

        if (proc_error == 1) {
            write_empty_output();
            continue;
        }

        double cRPIX[2] = {0.0, 0.0};
        double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double cRVAL[2] = {0.0, 0.0};
        double PU[2][LensingConfig::npd] = {{0.0}, {0.0}};

        int astrometry_error = 0;
        Astrometry::readAstrometryPara(headname, ichip + 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd, astrometry_error);
        if (astrometry_error == 1) {
            write_empty_output();
            continue;
        }

        int ngal = 0;
        std::string info_filename = OutputLayout::chipPath(
            dirOutput, "stamps/dat_SrcInfo", PREFIX, "_source_info.dat");
        std::ifstream fin(info_filename);
        if (!fin.is_open()) {
            std::cerr << "Error / proc_shear source_info catalog file error!! " << info_filename << std::endl;
            std::exit(1);
        }

        std::string header;
        std::getline(fin, header); // skip header line

        std::vector<std::vector<float>> gal_para;
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::vector<float> row(LensingConfig::npara, 0.0f);
            bool success = true;
            for (int i = 0; i <= LensingConfig::iSNR_F; ++i) {
                if (!(ss >> row[i])) {
                    success = false;
                    break;
                }
            }
            if (success) {
                gal_para.push_back(row);
            }
        }
        fin.close();

        ngal = gal_para.size();
        if (ngal == 0) {
            write_empty_output();
            continue;
        }

        int len_g = LensingConfig::len_g;
        int ngal_max = LensingConfig::ngal_max;
        int nn1 = ns * len_g;
        int nn2 = ns * (ngal / len_g + 1);

        std::vector<float> gal_p_coll;
        std::string power_fits = OutputLayout::chipPath(
            dirOutput, "stamps/fits_SrcP", PREFIX, "_source_p.fits");
        if (!FitsIO::readStamps(ngal_max, 1, ngal, ns, ns, gal_p_coll, nn1, nn2, power_fits)) {
            std::cerr << "Error / proc_shear source_p FITS file error!! " << power_fits << std::endl;
            std::exit(1);
        }

        fout10.open(out_filename);
        fout10 << std::setprecision(10);
        fout10 << "poly_chi2 xc yc sigma nstar imax jmax "
               << "half_light_flux half_light_area flag psf_FWHM SNR_F "
               << "ra dec gf1 gf2 g1 g2 de h1 h2 cos2 sin2 parity\n";

        for (int i = 0; i < ngal; ++i) {
            std::vector<float> gal_p(ns * ns);
            std::copy(gal_p_coll.begin() + i * ns * ns, gal_p_coll.begin() + (i + 1) * ns * ns, gal_p.begin());

            double x = gal_para[i][1];
            double y = gal_para[i][2];

            std::vector<float> psf_model(ns * ns, 0.0f);
            std::vector<float> psf_model0(ns * ns, 0.0f);
            double xx = 2.0 * (x / static_cast<double>(chipnx)) - 1.0;
            double yy = 2.0 * (y / static_cast<double>(chipny)) - 1.0;
            PSFModel::getPSFModel(ns, npl, local_coe, xx, yy, psf_model, psf_model0);

            if (std::isnan(psf_model[0])) {
                for (int j = 0; j <= LensingConfig::iparity; ++j) {
                    fout10 << "-999.0 ";
                }
                fout10 << "\n";
                std::cerr << "Error / proc_shear PSF model layer1 for chip " << imageFiles[ichip] << std::endl;
                continue;
            }

            float poly_chi2 = 0.0f;
            ExStar::anaChi2Simple(ns, psf_model.data(), psf_model0.data(), poly_chi2);
            // poly_chi2 = MINVAL(psf_model); // alternative

            gal_para[i][0] = (poly_chi2 - poly_ave) / poly_std;

            float psf_FWHM = 0.0f;
            getPSFArea(psf_model.data(), psf_FWHM);

            gal_para[i][LensingConfig::iPSF] = psf_FWHM;
            gal_para[i][LensingConfig::istar] = static_cast<float>(nstar);

            double ra = 0.0, dec = 0.0;
            Astrometry::coordinateTransferPU(ra, dec, x, y, 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd);
            gal_para[i][LensingConfig::ira] = static_cast<float>(ra);
            gal_para[i][LensingConfig::idec] = static_cast<float>(dec);

            double gf1 = 0.0, gf2 = 0.0, cos2 = 0.0, sin2 = 0.0;
            int parity = 0;
            Astrometry::fieldDistortionPU(x, y, LensingConfig::npd, PU, cD, cRPIX, gf1, gf2, cos2, sin2, parity);

            gal_para[i][LensingConfig::igf1] = static_cast<float>(gf1);
            gal_para[i][LensingConfig::igf2] = static_cast<float>(gf2);

            float g1 = 0.0f, g2 = 0.0f, de = 0.0f, h1 = 0.0f, h2 = 0.0f;
            getShear(ns, gal_p.data(), psf_model.data(), g1, g2, de, h1, h2);

            gal_para[i][LensingConfig::ig1] = static_cast<float>(g1 * cos2 + g2 * sin2);
            gal_para[i][LensingConfig::ig2] = static_cast<float>(g2 * cos2 - g1 * sin2);
            gal_para[i][LensingConfig::ide] = de;

            double cos4 = cos2 * cos2 - sin2 * sin2;
            double sin4 = 2.0 * sin2 * cos2;
            gal_para[i][LensingConfig::ih1] = static_cast<float>(h1 * cos4 + h2 * sin4);
            gal_para[i][LensingConfig::ih2] = static_cast<float>(h2 * cos4 - h1 * sin4);

            if (parity == -1) {
                gal_para[i][LensingConfig::ig2] = -gal_para[i][LensingConfig::ig2];
                gal_para[i][LensingConfig::ih2] = -gal_para[i][LensingConfig::ih2];
            }
            gal_para[i][LensingConfig::icos2] = static_cast<float>(cos2);
            gal_para[i][LensingConfig::isin2] = static_cast<float>(sin2);
            gal_para[i][LensingConfig::iparity] = static_cast<float>(parity);

            for (int j = 0; j <= LensingConfig::iparity; ++j) {
                fout10 << gal_para[i][j] << " ";
            }
            fout10 << "\n";
        }
        fout10.close();
    }
}

void procShear(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    expoShear(image_files.size(), image_files, dir_output, LensingConfig::chipnx, LensingConfig::chipny);
}

} // namespace ShearMeasurement
