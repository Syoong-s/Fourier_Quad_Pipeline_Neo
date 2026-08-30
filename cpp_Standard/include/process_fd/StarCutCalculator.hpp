#ifndef STAR_CUT_CALCULATOR_HPP
#define STAR_CUT_CALCULATOR_HPP

#include "process_fd/FDData.hpp"

#include <vector>

// ==========================================
// StarCutCalculator - point-source (star) removal via star-bar fitting
// Method: Build a size-magnitude histogram, locate the stellar locus,
//         and compute a size cut threshold.  Two modes:
//   - Single: one S_cut for all exposures (Nto1 / DES style)
//   - Per-exposure: one S_cut per exposure (NtoN / HSC style)
// ==========================================
class StarCutCalculator {
public:
    // Single global star cut (all exposures share one bar)
    static void calculateGlobalStarCut(const FDData& data,
                                       float& S_mean, float& S_std,
                                       float& S_cut);

    // Per-exposure star cut (one bar per exposure)
    static void calculateGlobalStarCutAuto(
        const FDData& data,
        std::vector<float>& S_mean_arr,
        std::vector<float>& S_std_arr,
        std::vector<float>& S_cut_arr);

    // Apply advanced cuts (per-exposure star cut + SNR cuts)
    // Used only in per-exposure mode.
    static void applyAdvancedCuts(FDData& data,
                                  const std::vector<float>& S_cut_arr);

    // Apply single global star cut in-place
    static void applySingleStarCut(FDData& data, float S_cut);
};

#endif  // STAR_CUT_CALCULATOR_HPP
