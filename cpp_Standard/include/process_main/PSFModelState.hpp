#ifndef PSF_MODEL_STATE_HPP
#define PSF_MODEL_STATE_HPP

#include "LensingConfig.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace PSFModel {
namespace Internal {

// ==========================================
// Structure: Store one chip's dynamically sized PSF candidates and chi grid
// Method: Reserve only metadata initially, then index a full row-major
//         chi matrix allocated from the actual candidate count.
// ==========================================
struct ChipPSFState {
    using StarRow = std::array<double, LensingConfig::src_npara>;

    std::vector<StarRow> stars;
    std::vector<float> chi_d;

    // ==========================================
    // Function: Allocate the full pairwise chi matrix for loaded candidates
    // Method: Zero one row-major n-by-n matrix using the live star count.
    // ==========================================
    void allocateChiD() {
        const std::size_t nstar = stars.size();
        chi_d.assign(nstar * nstar, 0.0f);
    }

    // ==========================================
    // Function: Access one mutable pairwise chi value
    // Method: Use the live star count as the row-major matrix stride.
    // ==========================================
    float& getChiD(int star1, int star2) {
        const std::size_t nstar = stars.size();
        return chi_d[static_cast<std::size_t>(star1) * nstar
                     + static_cast<std::size_t>(star2)];
    }

    // ==========================================
    // Function: Access one immutable pairwise chi value
    // Method: Use the live star count as the row-major matrix stride.
    // ==========================================
    const float& getChiD(int star1, int star2) const {
        const std::size_t nstar = stars.size();
        return chi_d[static_cast<std::size_t>(star1) * nstar
                     + static_cast<std::size_t>(star2)];
    }
};

// ==========================================
// Structure: Own all per-chip PSF state for one exposure
// Method: Create empty chip containers and preserve the legacy accessor
//         shape while deriving every star count from dynamic storage.
// ==========================================
struct ExposurePSFState {
    std::vector<ChipPSFState> chips;

    // ==========================================
    // Function: Create empty dynamic state for the exposure's live chip count
    // Method: Size only the outer chip vector without allocating star matrices.
    // ==========================================
    explicit ExposurePSFState(int nchip)
        : chips(static_cast<std::size_t>(nchip)) {}

    // ==========================================
    // Function: Access one mutable star parameter
    // Method: Delegate to the selected chip's actual-size star row.
    // ==========================================
    double& getStarPara(int chip, int star, int para) {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Access one immutable star parameter
    // Method: Delegate to the selected chip's actual-size star row.
    // ==========================================
    const double& getStarPara(int chip, int star, int para) const {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Access one mutable chip chi value
    // Method: Delegate to the chip's live-stride row-major matrix.
    // ==========================================
    float& getChiD(int chip, int star1, int star2) {
        return chips[chip].getChiD(star1, star2);
    }

    // ==========================================
    // Function: Access one immutable chip chi value
    // Method: Delegate to the chip's live-stride row-major matrix.
    // ==========================================
    const float& getChiD(int chip, int star1, int star2) const {
        return chips[chip].getChiD(star1, star2);
    }

    // ==========================================
    // Function: Return one chip's actual candidate count
    // Method: Derive the count from star-vector size to prevent divergence.
    // ==========================================
    int getNStar(int chip) const {
        return static_cast<int>(chips[chip].stars.size());
    }
};

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_MODEL_STATE_HPP
