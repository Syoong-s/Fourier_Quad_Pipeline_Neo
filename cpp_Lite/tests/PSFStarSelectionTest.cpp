#include "PSFStarSelection.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
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
// Function: Verify analytic LOO against explicit leave-one-out fits
// Method: Fit a small quadratic design, use each hat diagonal in the production
//         formula, and compare with a brute-force refit that removes that row.
// ==========================================
void testAnalyticLOO() {
    constexpr int sample_count = 20;
    Eigen::MatrixXd design(sample_count, 3);
    Eigen::VectorXd observed(sample_count);
    for (int row = 0; row < sample_count; ++row) {
        const double x = -1.0 + 2.0 * row / static_cast<double>(sample_count - 1);
        design(row, 0) = 1.0;
        design(row, 1) = x;
        design(row, 2) = x * x;
        observed(row) = 2.0 - 0.4 * x + 0.7 * x * x
            + 0.02 * std::sin(3.0 * row);
    }
    const Eigen::VectorXd coefficients = design.colPivHouseholderQr().solve(observed);
    const Eigen::VectorXd fitted = design * coefficients;
    const Eigen::MatrixXd covariance =
        (design.transpose() * design).inverse();

    for (int omitted = 0; omitted < sample_count; ++omitted) {
        const double leverage = (
            design.row(omitted) * covariance
            * design.row(omitted).transpose())(0, 0);
        double loo_residual = 0.0;
        double loo_model = 0.0;
        require(computeAnalyticLOO(
                    observed(omitted), fitted(omitted), leverage,
                    1.0e-6, loo_residual, loo_model),
                "well-conditioned leverage must produce an analytic LOO model");

        Eigen::MatrixXd reduced_design(sample_count - 1, 3);
        Eigen::VectorXd reduced_observed(sample_count - 1);
        int target = 0;
        for (int row = 0; row < sample_count; ++row) {
            if (row == omitted) continue;
            reduced_design.row(target) = design.row(row);
            reduced_observed(target) = observed(row);
            target++;
        }
        const Eigen::VectorXd reduced_coefficients =
            reduced_design.colPivHouseholderQr().solve(reduced_observed);
        const double explicit_model =
            (design.row(omitted) * reduced_coefficients)(0, 0);
        require(std::abs(loo_model - explicit_model) < 1.0e-10,
                "analytic and explicit leave-one-out predictions must agree");
        require(std::abs(loo_residual - (observed(omitted) - explicit_model))
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
    testAnalyticLOO();
    std::cout << "PSF star-selection tests passed\n";
    return EXIT_SUCCESS;
}
