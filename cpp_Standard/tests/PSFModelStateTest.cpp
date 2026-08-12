#include "PSFModelState.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

// ==========================================
// Function: Stop the PSF-state test when one invariant fails
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFModelState test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Populate one chip and allocate its production chi matrix
// Method: Reserve the configured hint, append the requested actual star count,
//         and allocate the full actual-size row-major matrix.
// ==========================================
void populateChip(PSFModel::Internal::ExposurePSFState& state,
                  int chip_index, int star_count) {
    auto& chip = state.chips[chip_index];
    chip.stars.reserve(LensingConfig::nstar_max);
    for (int star = 0; star < star_count; ++star) {
        PSFModel::Internal::ChipPSFState::StarRow row{};
        row[0] = star;
        chip.stars.push_back(row);
    }
    chip.allocateChiD();
}

// ==========================================
// Function: Verify dynamic candidate counts and full chi-matrix dimensions
// Method: Exercise empty, small, ordinary, and above-reservation chips while
//         checking wrapper access at the last above-2000 element.
// ==========================================
void testDynamicChipSizes() {
    PSFModel::Internal::ExposurePSFState state(4);
    require(state.chips.size() == 4, "outer state must match live chip count");
    for (const auto& chip : state.chips) {
        require(chip.stars.empty() && chip.chi_d.empty(),
                "construction must not preallocate per-chip matrices");
    }

    const int counts[4] = {0, 10, 300, 2301};
    for (int chip = 0; chip < 4; ++chip) {
        populateChip(state, chip, counts[chip]);
        const std::size_t expected =
            static_cast<std::size_t>(counts[chip]) * counts[chip];
        require(state.getNStar(chip) == counts[chip],
                "candidate count must follow dynamic star storage");
        require(state.chips[chip].chi_d.size() == expected,
                "chi matrix must use actual nstar squared");
    }

    state.getStarPara(3, 2300, 4) = 1.0;
    state.getChiD(3, 2300, 2299) = 7.5f;
    state.getChiD(3, 2299, 2300) = 7.5f;
    const auto& const_state = state;
    require(const_state.getStarPara(3, 2300, 4) == 1.0,
            "star access must extend beyond the reservation hint");
    require(const_state.getChiD(3, 2300, 2299) == 7.5f
                && const_state.getChiD(3, 2299, 2300) == 7.5f,
            "full chi matrix must preserve symmetric live-stride access");
}

}  // namespace

// ==========================================
// Function: Run focused Stage-5 dynamic-state regression cases
// Method: Execute all chip-size cases and report one success line.
// ==========================================
int main() {
    testDynamicChipSizes();
    std::cout << "PSFModelState tests passed\n";
    return EXIT_SUCCESS;
}
