#ifndef PSF_MODEL_STATE_HPP
#define PSF_MODEL_STATE_HPP

#include "LensingConfig.hpp"
#include "process_main/PSFStarSelection.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace PSFModel {
namespace Internal {

// ==========================================
// Structure: Store one candidate's explicit star-selection metadata
// Method: Keep scientific flags, cached Fourier window, Type-3 state, and
//         PRESS diagnostics outside the legacy StarRow column layout.
// ==========================================
struct StarSelectionState {
    bool gaia_matched = false;
    bool in_size_locus = false;
    bool selected_group = false;
    bool selected_press = false;
    double full_power_sum = 0.0;
    double chi_window_sum = 0.0;
    float min_chi = 0.0f;
    double bad_pair_fraction = 0.0;
    double press_raw_score = 0.0;
    double press_standardized_score = 0.0;
    double leverage = 0.0;
    std::vector<float> chi_window;
};

// ==========================================
// Structure: Cache one chip's final Stage-5 polynomial fit
// Method: Preserve the selected original indices, coefficients, and leverage so
//         output generation never repeats an unchanged fit.
// ==========================================
struct ChipPSFFitState {
    bool valid = false;
    bool press_removed_any = false;
    int initial_star_count = 0;
    std::vector<int> star_indices;
    std::vector<double> coefficients;
    std::vector<double> leverage;

    // ==========================================
    // Function: Reset cached PSF fitting products
    // Method: Clear all flags and vectors before processing a new chip state.
    // ==========================================
    void clear() {
        valid = false;
        press_removed_any = false;
        initial_star_count = 0;
        star_indices.clear();
        coefficients.clear();
        leverage.clear();
    }

    // ==========================================
    // Function: Commit a successful optional PRESS refit transaction
    // Method: Leave the cached first fit byte-for-byte unchanged unless the
    //         refit is valid and its final index/leverage dimensions agree.
    // ==========================================
    bool tryCommitPressRefit(
        bool refit_valid,
        std::vector<int> refit_star_indices,
        std::vector<double> refit_coefficients,
        std::vector<double> refit_leverage) {
        if (!refit_valid || refit_star_indices.empty()
            || refit_coefficients.empty()
            || refit_star_indices.size() != refit_leverage.size()) {
            return false;
        }
        valid = true;
        press_removed_any = true;
        star_indices = std::move(refit_star_indices);
        coefficients = std::move(refit_coefficients);
        leverage = std::move(refit_leverage);
        return true;
    }
};

// ==========================================
// Structure: Store one chip's dynamically sized PSF candidates and selection
// Method: Align explicit selection metadata with legacy parameter rows and keep
//         only O(N*window + N*K) grouping storage rather than a square matrix.
// ==========================================
struct ChipPSFState {
    static constexpr int star_area_index = 12;
    using StarRow = std::array<double, LensingConfig::src_npara + 1>;

    static_assert(LensingConfig::src_npara == star_area_index,
                  "external catalog width must end before star_area");
    static_assert(star_area_index < static_cast<int>(std::tuple_size<StarRow>::value),
                  "star_area index must fit the private PSF row");

    std::vector<StarRow> stars;
    std::vector<StarSelectionState> selection;
    ChipPSFFitState fit;
};

// ==========================================
// Structure: Own all per-chip PSF state for one exposure
// Method: Size only the live outer chip vector and derive star counts from dynamic storage.
// ==========================================
struct ExposurePSFState {
    std::vector<ChipPSFState> chips;

    // ==========================================
    // Function: Create empty dynamic state for the exposure's live chip count
    // Method: Avoid allocating any fixed per-chip candidate or chi matrix storage.
    // ==========================================
    explicit ExposurePSFState(int nchip)
        : chips(static_cast<std::size_t>(nchip)) {}

    // ==========================================
    // Function: Access one mutable star parameter
    // Method: Index the requested row in the chip's live candidate vector.
    // ==========================================
    double& getStarPara(int chip, int star, int para) {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Access one immutable star parameter
    // Method: Index the requested row in the chip's live candidate vector.
    // ==========================================
    const double& getStarPara(int chip, int star, int para) const {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Return one chip's live candidate count
    // Method: Derive the integer count from dynamic row storage.
    // ==========================================
    int getNStar(int chip) const {
        return static_cast<int>(chips[chip].stars.size());
    }
};

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_MODEL_STATE_HPP
