#include "process_main/PSFModelState.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

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
// Function: Verify one actual candidate count and linear selection-state dimension
// Method: Append all rows, align explicit metadata, and cache bounded window/KNN data.
// ==========================================
void testChipSize(int star_count) {
    PSFModel::Internal::ExposurePSFState state(1);
    auto& chip = state.chips[0];
    chip.stars.reserve(LensingConfig::nstar_max);
    for (int star = 0; star < star_count; ++star) {
        PSFModel::Internal::ChipPSFState::StarRow row{};
        row[0] = star;
        chip.stars.push_back(row);
    }
    chip.selection.resize(static_cast<std::size_t>(star_count));
    for (int star = 0; star < star_count; ++star) {
        chip.selection[star].chi_window.push_back(static_cast<float>(star));
        chip.selection[star].knn.push_back({star, 0.0f});
    }

    require(state.getNStar(0) == star_count,
            "candidate count must follow dynamic storage");
    require(chip.selection.size() == static_cast<std::size_t>(star_count),
            "selection metadata must align one-to-one with candidate rows");
    if (star_count > 0) {
        require(chip.selection.back().chi_window.size() == 1
                    && chip.selection.back().knn.size() == 1,
                "last candidate must own bounded non-square cache vectors");
    }
}

// ==========================================
// Function: Run focused Stage-5 dynamic-state regression cases
// Method: Exercise zero, ordinary, boundary, and above-reservation star counts
//         without allocating any candidate-count-squared matrix.
// ==========================================
void testDynamicChipSizes() {
    const int counts[] = {0, 10, 1999, 2000, 2001, 2301};
    for (int count : counts) {
        testChipSize(count);
    }
}

// ==========================================
// Function: Verify failed PRESS refits preserve the cached first fit
// Method: Attempt an invalid transactional commit and compare every original
//         fit field before also checking a valid replacement can be committed.
// ==========================================
void testPressRefitTransaction() {
    PSFModel::Internal::ChipPSFFitState fit;
    fit.valid = true;
    fit.press_removed_any = false;
    fit.initial_star_count = 3;
    fit.star_indices = {1, 2, 3};
    fit.coefficients = {10.0, 20.0};
    fit.leverage = {0.1, 0.2, 0.3};
    require(!fit.tryCommitPressRefit(false, {1, 3}, {30.0, 40.0}, {0.15, 0.25}),
            "failed PRESS refit must not commit");
    require(fit.valid && !fit.press_removed_any && fit.initial_star_count == 3
                && fit.star_indices == std::vector<int>({1, 2, 3})
                && fit.coefficients == std::vector<double>({10.0, 20.0})
                && fit.leverage == std::vector<double>({0.1, 0.2, 0.3}),
            "failed PRESS refit must preserve all first-fit cache fields");
    require(fit.tryCommitPressRefit(true, {1, 3}, {30.0, 40.0}, {0.15, 0.25}),
            "valid PRESS refit must commit");
    require(fit.valid && fit.press_removed_any && fit.initial_star_count == 3
                && fit.star_indices == std::vector<int>({1, 3})
                && fit.coefficients == std::vector<double>({30.0, 40.0})
                && fit.leverage == std::vector<double>({0.15, 0.25}),
            "successful PRESS refit must retain the original first-fit count");
}

}  // namespace

// ==========================================
// Function: Run the PSF-state test suite
// Method: Execute all live-size cases and report one success line.
// ==========================================
int main() {
    testDynamicChipSizes();
    testPressRefitTransaction();
    std::cout << "PSFModelState tests passed\n";
    return EXIT_SUCCESS;
}
