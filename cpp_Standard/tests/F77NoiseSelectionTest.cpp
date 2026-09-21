#include "process_main/F77NoiseSelection.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using SourceExtractor::Internal::F77NoiseCandidate;
using SourceExtractor::Internal::selectF77NoiseCandidate;

// ==========================================
// Function: Stop the F77-noise test on one failed invariant
// Method: Emit a focused diagnostic and terminate with failure.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "F77 noise-selection test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Paint one square candidate with a controlled maximum
// Method: Fill every pixel below the requested peak and set one corner to it.
// ==========================================
void paintCandidate(
    std::vector<float>& image,
    int nx,
    int x0,
    int y0,
    int side,
    float peak) {
    for (int y = y0; y < y0 + side; ++y) {
        for (int x = x0; x < x0 + side; ++x) {
            image[static_cast<std::size_t>(y) * nx + x] = peak - 1.0f;
        }
    }
    image[static_cast<std::size_t>(y0) * nx + x0] = peak;
}

// ==========================================
// Function: Verify deterministic F77 candidate ranking and gates
// Method: Exercise global minimum peak, strict tie order, nd cutoff, and sentinel.
// ==========================================
void testF77CandidateSelection() {
    constexpr int nx = 80;
    constexpr int ny = 80;
    constexpr int side = 8;
    constexpr int half = 4;
    constexpr int margin = 2;
    constexpr double xp = 41.0;
    constexpr double yp = 41.0;
    std::vector<float> image(static_cast<std::size_t>(nx) * ny, 50.0f);
    std::vector<int> weight(static_cast<std::size_t>(nx) * ny, 1);

    const F77NoiseCandidate baseline = selectF77NoiseCandidate(
        nx, ny, image, weight, xp, yp, side, half, margin);
    require(baseline.found, "uniform outer-ring candidates must be eligible");

    paintCandidate(image, nx, baseline.x0, baseline.y0, side, 5.0f);
    const int second_x0 = baseline.x0;
    const int second_y0 = baseline.y0 + side;
    paintCandidate(image, nx, second_x0, second_y0, side, 5.0f);
    F77NoiseCandidate selected = selectF77NoiseCandidate(
        nx, ny, image, weight, xp, yp, side, half, margin);
    require(selected.x0 == baseline.x0 && selected.y0 == baseline.y0,
            "equal peaks must preserve the first Fortran scan position");

    for (int offset = 0; offset <= side; ++offset) {
        const int x = baseline.x0 + (offset % side);
        const int y = baseline.y0 + (offset / side);
        weight[static_cast<std::size_t>(y) * nx + x] = 0;
    }
    selected = selectF77NoiseCandidate(
        nx, ny, image, weight, xp, yp, side, half, margin);
    require(selected.x0 == second_x0 && selected.y0 == second_y0,
            "a minimum-peak candidate with nd>nl must be excluded");

    std::fill(image.begin(), image.end(), 100000.0f);
    selected = selectF77NoiseCandidate(
        nx, ny, image, weight, xp, yp, side, half, margin);
    require(!selected.found,
            "the historical 100000 sentinel must reject equal-or-larger peaks");
}

}  // namespace

// ==========================================
// Function: Run the focused F77 noise-selection regression
// Method: Execute deterministic candidate geometry and boundary tests.
// ==========================================
int main() {
    testF77CandidateSelection();
    std::cout << "F77 noise-selection tests passed\n";
    return EXIT_SUCCESS;
}
