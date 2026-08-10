#ifndef EXPOSURE_INFO_HPP
#define EXPOSURE_INFO_HPP

#include <vector>
#include <string>

namespace ExposureInfo {
    // Global exposure parameters array (originally /expo_para_pass/; sized to 6 * N_EXPO at runtime)
    extern std::vector<float> expo_para;

    void getExpoInfo(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput, float para[6]);
    void procInfo(int iexpo);
}

#endif // EXPOSURE_INFO_HPP
