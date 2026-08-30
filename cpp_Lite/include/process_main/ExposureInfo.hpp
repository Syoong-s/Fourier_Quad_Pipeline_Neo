#ifndef EXPOSURE_INFO_HPP
#define EXPOSURE_INFO_HPP

#include <vector>
#include <string>

namespace ExposureInfo {
    struct State {
        std::vector<float> parameters;

        // ==========================================
        // Function: Reset Stage-8 aggregate storage
        // Method: Allocate exactly the requested number of zeroed values.
        // ==========================================
        void reset(std::size_t count) {
            parameters.assign(count, 0.0f);
        }
    };

    extern State state;

    void getExpoInfo(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput, float para[6]);
    void procInfo(int iexpo);
}

#endif // EXPOSURE_INFO_HPP
