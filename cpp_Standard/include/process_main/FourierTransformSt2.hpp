#ifndef FOURIER_TRANSFORM_ST2_HPP
#define FOURIER_TRANSFORM_ST2_HPP

#include <string>

namespace FourierTransformSt2 {
    // ==========================================
    // Function: Run the Stage-6 galaxy FFT for one exposure
    // Method: Resolve its chip list and process each chip serially.
    // ==========================================
    void procFourierTSt2(int iexpo);

    // ==========================================
    // Function: Run the Stage-6 galaxy FFT for one chip
    // Method: Enforce the norm and required Stage-3 input contracts before transforming stamps.
    // ==========================================
    void chipProcessFourierTSt2(const std::string& imageFile,
                                const std::string& dirOutput);
}

#endif // FOURIER_TRANSFORM_ST2_HPP
