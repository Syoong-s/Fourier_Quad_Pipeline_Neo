#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace PSFModel::Internal;

// ==========================================
// Function: Stop the star-selection test on one failed invariant
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF star-selection test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Verify shared chi-window bounds and normalized distance
// Method: Check the locked 64-pixel bounds and exact zero/nonzero distances.
// ==========================================
void testChiWindowAndDistance() {
    const PSFChiWindow window = getPSFChiWindow(64);
    require(window.first == 15 && window.last == 47
                && window.pixelCount() == 1089,
            "64x64 chi window must remain inclusive 15..47");
    std::vector<float> first(9, 0.1f);
    std::vector<float> second = first;
    require(normalizedChiDistance(first, second) == 0.0f,
            "identical normalized windows must have zero chi distance");
    second[4] += 0.05f;
    require(normalizedChiDistance(first, second) > 0.0f,
            "different normalized windows must have positive chi distance");
}

// ==========================================
// Function: Verify histogram FWHM peak selection and discretization guards
// Method: Exercise density selection, Gaia-supported alternate peak, repeated
//         values, and the configured minimum-sample failure.
// ==========================================
void testFWHMLocus() {
    std::vector<FWHMSample> density_samples;
    for (int index = 0; index < 80; ++index) {
        density_samples.push_back({1.00 + 0.01 * (index % 5), false});
    }
    for (int index = 0; index < 120; ++index) {
        density_samples.push_back({1.70 + 0.02 * (index % 25), false});
    }
    FWHMLocus locus;
    require(estimateFWHMLocus(density_samples, 128, 4.0, 30, 10, locus),
            "two-population FWHM fixture must produce a locus");
    require(std::abs(locus.center - 1.02) < 0.08,
            "highest-density narrow stellar peak must be selected without Gaia");

    std::vector<FWHMSample> gaia_samples;
    for (int index = 0; index < 45; ++index) {
        gaia_samples.push_back({0.90 + 0.01 * (index % 5), index < 12});
    }
    for (int index = 0; index < 100; ++index) {
        gaia_samples.push_back({1.80 + 0.005 * (index % 5), false});
    }
    require(estimateFWHMLocus(gaia_samples, 128, 4.0, 30, 10, locus),
            "Gaia-supported two-peak fixture must produce a locus");
    require(std::abs(locus.center - 0.92) < 0.08,
            "Gaia median must select the supported smaller peak");

    std::vector<FWHMSample> repeated(40, {1.25, false});
    require(estimateFWHMLocus(repeated, 128, 4.0, 30, 10, locus)
                && locus.width > 0.0 && locus.lower < 1.25 && locus.upper > 1.25,
            "repeated FWHM values must receive a finite positive width floor");
    repeated.resize(29);
    require(!estimateFWHMLocus(repeated, 128, 4.0, 30, 10, locus),
            "FWHM locus must reject fewer than the configured samples");
}

// ==========================================
// Function: Verify astro parsing and nearest-position Gaia labels
// Method: Parse finite matched rows, test radius boundaries and shared labels,
//         and reject a malformed matched-source row.
// ==========================================
void testGaiaParsingAndMatching() {
    std::istringstream valid(
        "1 2 3 4\n"
        "0.1 0.0 0.0 0.1\n"
        "2 20 30\n"
        "10 20 100.0 200.0\n"
        "11 21 150.0 250.0\n");
    std::vector<std::array<double, 2>> gaia_xy;
    std::string error;
    require(parseAstrometryGaiaPositions(valid, gaia_xy, error)
                == AstrometryGaiaReadStatus::Accepted,
            "well-formed astro rows must parse");
    require(gaia_xy.size() == 2
                && hasNearestGaiaMatch(101.0, 201.0, gaia_xy, 2.5)
                && !hasNearestGaiaMatch(104.0, 204.0, gaia_xy, 2.5),
            "nearest Gaia radius must be applied in image pixels");
    require(hasNearestGaiaMatch(100.5, 200.5, gaia_xy, 2.5)
                && hasNearestGaiaMatch(99.5, 199.5, gaia_xy, 2.5),
            "multiple candidates may share one supporting Gaia position");

    std::istringstream malformed(
        "1 2 3 4\n0.1 0.0 0.0 0.1\n1 2 3\n10 20 bad 40\n");
    require(parseAstrometryGaiaPositions(malformed, gaia_xy, error)
                == AstrometryGaiaReadStatus::Malformed,
            "malformed astro matched-source rows must be rejected");
}

// ==========================================
// Function: Verify streaming top-K, mutual edges, and shared group selection
// Method: Compare deterministic top-K order, build mutual components, and
//         exercise the secondary size-and-Gaia conjunction.
// ==========================================
void testGrouping() {
    std::vector<NeighborEdge> top_k;
    updateTopK(top_k, 3, 3.0f, 2);
    updateTopK(top_k, 1, 1.0f, 2);
    updateTopK(top_k, 2, 2.0f, 2);
    require(top_k.size() == 2 && top_k[0].star_index == 1
                && top_k[1].star_index == 2,
            "streaming top-K must match sorted full-distance reference");

    std::vector<std::vector<NeighborEdge>> neighbours(4);
    neighbours[0] = {{1, 1.0f}, {2, 2.0f}};
    neighbours[1] = {{0, 1.0f}, {2, 1.5f}};
    neighbours[2] = {{1, 1.5f}, {0, 2.0f}};
    neighbours[3] = {{2, 0.5f}};
    const std::vector<int> active = {0, 1, 2, 3};
    const std::vector<GraphEdge> mutual = buildMutualKNNEdges(active, neighbours);
    const std::vector<bool> gaia = {false, true, false, true};
    const std::vector<StarGroup> connected =
        buildConnectedGroups(active, mutual, gaia);
    require(connected.size() == 2,
            "mutual-KNN fixture must form one triple and one singleton");

    const std::vector<StarGroup> groups = {
        {{0, 1, 2, 3}, 0},
        {{4, 5}, 2},
        {{6, 7}, 1},
        {{8}, 3}
    };
    const std::vector<int> selected =
        selectMainAndSecondaryGroups(groups, 0.40, 2);
    require(selected == std::vector<int>({0, 1, 2, 3, 4, 5}),
            "secondary groups must pass both relative size and Gaia count");
}

// ==========================================
// Structure: Provide the minimal cached fields required by active-KNN rebuilds
// Method: Mirror production chi-window and neighbour storage without PSF I/O.
// ==========================================
struct SyntheticKNNSelectionState {
    std::vector<float> chi_window;
    std::vector<NeighborEdge> knn;
};

// ==========================================
// Function: Verify KNN slots are refilled after the minChi survivor cut
// Method: Build an initial top-2 containing rejected close stars, rebuild on
//         survivors only, and require the next two valid neighbours to replace them.
// ==========================================
void testKNNRebuiltAfterMinChiCut() {
    std::vector<SyntheticKNNSelectionState> candidates(5);
    for (int index = 0; index < static_cast<int>(candidates.size()); ++index) {
        candidates[index].chi_window = {
            1.0f + static_cast<float>(index) * 0.1f};
    }

    rebuildActiveKNN(std::vector<int>({0, 1, 2, 3, 4}), candidates, 2);
    require(candidates[0].knn.size() == 2
                && candidates[0].knn[0].star_index == 1
                && candidates[0].knn[1].star_index == 2,
            "pre-cut top-K fixture must be occupied by the closest rejected stars");

    rebuildActiveKNN(std::vector<int>({0, 3, 4}), candidates, 2);
    require(candidates[0].knn.size() == 2
                && candidates[0].knn[0].star_index == 3
                && candidates[0].knn[1].star_index == 4,
            "survivor-only rebuild must refill every vacated top-K slot");
    require(candidates[1].knn.empty() && candidates[2].knn.empty(),
            "minChi-rejected candidates must retain no stale neighbour state");
}

// ==========================================
// Function: Verify capped reference selection and one-pass minChi semantics
// Method: Exercise exposure top-fraction ranking, per-chip caps, deterministic
//         ties, non-reference nearest pairs, and unordered-pair deduplication.
// ==========================================
void testMinChiReferencesAndPairs() {
    std::vector<MinChiReferenceCandidate> candidates;
    for (int star = 0; star < 6; ++star) candidates.push_back({0, star, 100.0 - star});
    candidates.push_back({1, 0, 95.0});
    for (int star = 1; star < 15; ++star) candidates.push_back({1, star, 94.0 - star});
    const std::vector<MinChiReferenceCandidate> selected =
        selectMinChiReferenceStars(candidates, 1.0 / 3.0, 5);
    require(selected.size() == 6,
            "top-third reference pool must apply the cap without backfilling");
    for (int star = 0; star < 5; ++star) {
        require(selected[star].chip_index == 0 && selected[star].star_index == star,
                "largest five candidates on one chip must be retained in order");
    }
    require(selected[5].chip_index == 1 && selected[5].star_index == 0,
            "a capped top-pool entry must not suppress another chip reference");

    const std::vector<MinChiReferenceCandidate> tied = {
        {1, 2, 10.0}, {0, 5, 10.0}, {0, 3, 10.0}, {1, 1, 10.0},
        {2, 0, 9.0}, {2, 1, 8.0}, {2, 2, 7.0}, {2, 3, 6.0}, {2, 4, 5.0}
    };
    const std::vector<MinChiReferenceCandidate> tied_selected =
        selectMinChiReferenceStars(tied, 1.0 / 3.0, 5);
    require(tied_selected.size() == 3
                && tied_selected[0].chip_index == 0 && tied_selected[0].star_index == 3
                && tied_selected[1].chip_index == 0 && tied_selected[1].star_index == 5
                && tied_selected[2].chip_index == 1 && tied_selected[2].star_index == 1,
            "size ties must resolve by chip index and then star index");

    const std::vector<std::vector<float>> windows = {{1.0f}, {1.2f}, {5.0f}, {5.01f}};
    std::vector<MinChiCandidateView> views;
    for (int index = 0; index < 4; ++index) views.push_back({&windows[index], true, index < 2});
    const MinChiPairResult pair_result = computeMinChiAndThresholdPairs(views);
    const float nonreference_distance = normalizedChiDistance(windows[2], windows[3]);
    require(std::abs(pair_result.min_chi[2] - nonreference_distance) < 1.0e-7f
                && std::abs(pair_result.min_chi[3] - nonreference_distance) < 1.0e-7f,
            "a non-reference pair must still update both nearest distances");
    require(pair_result.threshold_pair_chi.size() == 5,
            "two references among four stars must contribute five unique pairs");
    const float reference_pair_distance = normalizedChiDistance(windows[0], windows[1]);
    int reference_pair_occurrences = 0;
    for (float chi : pair_result.threshold_pair_chi) {
        if (std::abs(chi - reference_pair_distance) < 1.0e-7f) reference_pair_occurrences++;
    }
    require(reference_pair_occurrences == 1,
            "the reference-reference pair must enter the threshold sample once");
}

// ==========================================
// Function: Verify leverage-standardized PRESS and removal safeguards
// Method: Check the exact score formula, denominator guard, removal cap, and
//         inclusive retained-star minimum boundaries.
// ==========================================
void testPressStandardizationAndDecision() {
    double standardized = 0.0;
    require(computeLeverageStandardizedPress(4.0, 0.75, 1.0e-6, standardized)
                && std::abs(standardized - 2.0) < 1.0e-12,
            "standardized PRESS must equal raw PRESS times sqrt(1-h)");
    require(!computeLeverageStandardizedPress(4.0, 1.0 - 5.0e-7, 1.0e-6, standardized),
            "standardized PRESS must reject a denominator below its floor");
    require(!computeLeverageStandardizedPress(
                std::numeric_limits<double>::infinity(), 0.1, 1.0e-6, standardized),
            "standardized PRESS must reject non-finite raw scores");
    require(decidePressRemoval(false, 30, 2, 16, 5) == PressRemovalDecision::Disabled,
            "disabled rejection must preserve the first fit");
    require(decidePressRemoval(true, 30, 0, 16, 5) == PressRemovalDecision::NoOutliers,
            "zero outliers must preserve the first fit");
    require(decidePressRemoval(true, 21, 5, 16, 5) == PressRemovalDecision::Apply,
            "five removals leaving sixteen stars must be allowed");
    require(decidePressRemoval(true, 30, 6, 16, 5) == PressRemovalDecision::TooManyOutliers,
            "six removals must trigger the configured cap");
    require(decidePressRemoval(true, 20, 5, 16, 5)
                == PressRemovalDecision::WouldUnderrunMinimum,
            "fifteen retained stars must trigger the local-minimum guard");
}

// ==========================================
// Function: Invert one finite nonsingular three-by-three matrix
// Method: Apply pivoted Gauss-Jordan elimination to a small augmented matrix.
// ==========================================
bool invertThreeByThree(
    const std::array<std::array<double, 3>, 3>& matrix,
    std::array<std::array<double, 3>, 3>& inverse) {
    std::array<std::array<double, 6>, 3> augmented{};
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) augmented[row][column] = matrix[row][column];
        augmented[row][row + 3] = 1.0;
    }
    for (int pivot = 0; pivot < 3; ++pivot) {
        int best_row = pivot;
        for (int row = pivot + 1; row < 3; ++row) {
            if (std::abs(augmented[row][pivot]) > std::abs(augmented[best_row][pivot])) best_row = row;
        }
        if (std::abs(augmented[best_row][pivot]) < 1.0e-14) return false;
        std::swap(augmented[pivot], augmented[best_row]);
        const double scale = augmented[pivot][pivot];
        for (double& value : augmented[pivot]) value /= scale;
        for (int row = 0; row < 3; ++row) {
            if (row == pivot) continue;
            const double factor = augmented[row][pivot];
            for (int column = 0; column < 6; ++column) {
                augmented[row][column] -= factor * augmented[pivot][column];
            }
        }
    }
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) inverse[row][column] = augmented[row][column + 3];
    }
    return true;
}

// ==========================================
// Function: Fit a quadratic while optionally omitting one sample
// Method: Form the three-term normal equations, invert them, and multiply by
//         the observation cross-products for the dependency-free LOO fixture.
// ==========================================
bool fitQuadratic(
    const std::vector<double>& x,
    const std::vector<double>& observed,
    int omitted,
    std::array<double, 3>& coefficients,
    std::array<std::array<double, 3>, 3>* covariance = nullptr) {
    std::array<std::array<double, 3>, 3> normal{};
    std::array<double, 3> rhs{};
    for (int row = 0; row < static_cast<int>(x.size()); ++row) {
        if (row == omitted) continue;
        const std::array<double, 3> basis = {1.0, x[row], x[row] * x[row]};
        for (int first = 0; first < 3; ++first) {
            rhs[first] += basis[first] * observed[row];
            for (int second = 0; second < 3; ++second) {
                normal[first][second] += basis[first] * basis[second];
            }
        }
    }
    std::array<std::array<double, 3>, 3> inverse{};
    if (!invertThreeByThree(normal, inverse)) return false;
    for (int first = 0; first < 3; ++first) {
        for (int second = 0; second < 3; ++second) coefficients[first] += inverse[first][second] * rhs[second];
    }
    if (covariance != nullptr) *covariance = inverse;
    return true;
}

// ==========================================
// Function: Verify analytic LOO against explicit leave-one-out fits
// Method: Fit a small quadratic design, use each hat diagonal in the production
//         formula, and compare with a brute-force refit that removes that row.
// ==========================================
void testAnalyticLOO() {
    constexpr int sample_count = 20;
    std::vector<double> x(static_cast<std::size_t>(sample_count));
    std::vector<double> observed(static_cast<std::size_t>(sample_count));
    for (int row = 0; row < sample_count; ++row) {
        x[row] = -1.0 + 2.0 * row / static_cast<double>(sample_count - 1);
        observed[row] = 2.0 - 0.4 * x[row] + 0.7 * x[row] * x[row]
            + 0.02 * std::sin(3.0 * row);
    }
    std::array<double, 3> coefficients{};
    std::array<std::array<double, 3>, 3> covariance{};
    require(fitQuadratic(x, observed, -1, coefficients, &covariance),
            "full quadratic fixture must be nonsingular");

    for (int omitted = 0; omitted < sample_count; ++omitted) {
        const std::array<double, 3> basis = {1.0, x[omitted], x[omitted] * x[omitted]};
        double fitted = 0.0;
        double leverage = 0.0;
        for (int first = 0; first < 3; ++first) {
            fitted += basis[first] * coefficients[first];
            for (int second = 0; second < 3; ++second) {
                leverage += basis[first] * covariance[first][second] * basis[second];
            }
        }
        double loo_residual = 0.0;
        double loo_model = 0.0;
        require(computeAnalyticLOO(
                    observed[omitted], fitted, leverage,
                    1.0e-6, loo_residual, loo_model),
                "well-conditioned leverage must produce an analytic LOO model");

        std::array<double, 3> reduced_coefficients{};
        require(fitQuadratic(x, observed, omitted, reduced_coefficients, nullptr),
                "leave-one-out quadratic fixture must remain nonsingular");
        double explicit_model = 0.0;
        for (int term = 0; term < 3; ++term) explicit_model += basis[term] * reduced_coefficients[term];
        require(std::abs(loo_model - explicit_model) < 1.0e-10,
                "analytic and explicit leave-one-out predictions must agree");
        require(std::abs(loo_residual - (observed[omitted] - explicit_model))
                    < 1.0e-10,
                "analytic LOO residual must equal observed minus explicit model");
    }

    double residual = 0.0;
    double model = 0.0;
    require(!computeAnalyticLOO(1.0, 0.0, 1.0 - 5.0e-7,
                                1.0e-6, residual, model),
            "near-unit leverage must fail the configured denominator guard");
}

}  // namespace

// ==========================================
// Function: Run the focused PSF star-selection regression suite
// Method: Execute quality-independent locus, Gaia, grouping, and LOO cases.
// ==========================================
int main() {
    testChiWindowAndDistance();
    testFWHMLocus();
    testGaiaParsingAndMatching();
    testGrouping();
    testKNNRebuiltAfterMinChiCut();
    testMinChiReferencesAndPairs();
    testAnalyticLOO();
    testPressStandardizationAndDecision();
    std::cout << "PSF star-selection tests passed\n";
    return EXIT_SUCCESS;
}
