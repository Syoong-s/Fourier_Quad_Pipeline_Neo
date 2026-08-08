#ifndef FD_DATA_HPP
#define FD_DATA_HPP

#include "FDConfig.hpp"

#include <vector>

// ==========================================
// FDData - In-memory shear catalog data (replaces Fortran common blocks)
// Method: All per-source arrays collected in one struct so they can be passed
//         as a single reference instead of through global common blocks.
// ==========================================
struct FDData {
    // Shear variables (shear_data_pass)
    std::vector<float> x1, x2;    // gf1, gf2 (field distortion)
    std::vector<float> y1, y2;    // g1, g2  (galaxy shear)
    std::vector<float> de1, de2;  // de∓h1, de±h1 (response, PDF mode)
    std::vector<float> ww;      // Jackknife weight = ±1/sqrt((g1/de)²+(g2/de)²)
    int ng = 0;

    // Additional variables (additional_data_pass)
    std::vector<float> magg, magr, magi;  // magnitudes
    std::vector<float> sizerel;            // relative size
    std::vector<float> src_snr;            // source SNR
    std::vector<float> rra, ddec;          // RA, Dec (degrees)
    std::vector<int>   iexpo;              // exposure index (per-exposure mode)
    std::vector<float> snrf;               // SNR_F (per-exposure mode)

    // Jackknife region assignment
    std::vector<int>   labels;

    void reserve(int capacity) {
        x1.resize(capacity); x2.resize(capacity);
        y1.resize(capacity); y2.resize(capacity);
        de1.resize(capacity); de2.resize(capacity);
        ww.resize(capacity);
        magg.resize(capacity); magr.resize(capacity); magi.resize(capacity);
        sizerel.resize(capacity); src_snr.resize(capacity);
        rra.resize(capacity); ddec.resize(capacity);
        iexpo.resize(capacity); snrf.resize(capacity);
        labels.resize(capacity);
    }
};

#endif  // FD_DATA_HPP
