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

}  // namespace

// ==========================================
// Function: Run the PSF-state test suite
// Method: Execute all live-size cases and report one success line.
// ==========================================
int main() {
    testDynamicChipSizes();
    std::cout << "PSFModelState tests passed\n";
    return EXIT_SUCCESS;
}
