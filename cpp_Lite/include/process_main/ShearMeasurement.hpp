#ifndef SHEAR_MEASUREMENT_HPP
#define SHEAR_MEASUREMENT_HPP

#include <vector>
#include <string>

namespace ShearMeasurement {
    void getWindowMinK(int ns, const std::vector<float>& psf_model, float thresh, float& k_win);
    void getWindowMinKVer2(int ns, const std::vector<float>& psf_model, float thresh, float& k_win);
    void getShear(int n, const float* gal, const float* psf, float& g1, float& g2, float& de, float& h1, float& h2);
    void getPSFArea(const float* model, float& FWHM);
    void expoShear(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int chipnx, int chipny);
    void procShear(int iexpo);
}

#endif // SHEAR_MEASUREMENT_HPP
