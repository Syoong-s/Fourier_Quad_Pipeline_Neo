#ifndef PSF_STAR_SELECTION_HPP
#define PSF_STAR_SELECTION_HPP

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace PSFModel {
namespace Internal {

// ==========================================
// Structure: Describe the inclusive central Fourier window used by PSF chi
// Method: Derive one shared bound pair from the square stamp side length.
// ==========================================
struct PSFChiWindow {
    int first = 0;
    int last = -1;

    // ==========================================
    // Function: Return the inclusive chi-window side length
    // Method: Convert valid bounds to a positive size and invalid bounds to zero.
    // ==========================================
    int side() const { return last >= first ? last - first + 1 : 0; }

    // ==========================================
    // Function: Return the square chi-window pixel count
    // Method: Square the derived inclusive side length.
    // ==========================================
    int pixelCount() const { return side() * side(); }
};

// ==========================================
// Function: Derive the shared PSF chi-window bounds
// Method: Preserve the legacy inclusive n/4-1 through 3n/4-1 definition.
// ==========================================
PSFChiWindow getPSFChiWindow(int n);

// ==========================================
// Structure: Supply one quality-valid FWHM measurement to locus estimation
// Method: Carry only the scalar FWHM and its optional Gaia supporting label.
// ==========================================
struct FWHMSample {
    double fwhm = 0.0;
    bool gaia_matched = false;
};

// ==========================================
// Structure: Publish a robust exposure-wide stellar FWHM locus
// Method: Store the median center, robust width, cut bounds, and histogram scale.
// ==========================================
struct FWHMLocus {
    bool valid = false;
    double center = 0.0;
    double width = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    double histogram_bin_width = 0.0;
};

// ==========================================
// Function: Estimate an exposure-wide stellar FWHM locus
// Method: Smooth a fixed-bin robust histogram, optionally choose the peak
//         nearest the Gaia median, then apply median/MAD clipping with a
//         discretization floor.
// ==========================================
bool estimateFWHMLocus(
    const std::vector<FWHMSample>& samples,
    int histogram_bins,
    double sigma_cut,
    int minimum_samples,
    int minimum_gaia_matches,
    FWHMLocus& locus);

enum class AstrometryGaiaReadStatus {
    Accepted,
    Empty,
    Malformed
};

// ==========================================
// Function: Parse matched Gaia image positions from one astro stream
// Method: Read the two WCS lines, the matched/user/reference counts, and exactly
//         n_matched finite RA/Dec/x/y rows while retaining only x/y.
// ==========================================
AstrometryGaiaReadStatus parseAstrometryGaiaPositions(
    std::istream& input,
    std::vector<std::array<double, 2>>& gaia_xy,
    std::string& error);

// ==========================================
// Function: Match one candidate to the nearest same-chip Gaia position
// Method: Accept when the minimum squared image-plane distance is within the
//         configured radius; do not impose one-to-one assignment.
// ==========================================
bool hasNearestGaiaMatch(
    double x,
    double y,
    const std::vector<std::array<double, 2>>& gaia_xy,
    double radius_pixels);

// ==========================================
// Structure: Store one exact same-chip nearest-neighbour relation
// Method: Retain the original candidate index and its normalized Fourier chi.
// ==========================================
struct NeighborEdge {
    int star_index = -1;
    float chi = 0.0f;
};

// ==========================================
// Function: Compute the exact normalized PSF chi distance
// Method: Apply the legacy sqrt(sum squared difference / mean signed flux)
//         directly to two cached central windows.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second);

// ==========================================
// Function: Maintain one exact sorted top-K neighbour list
// Method: Insert or improve the candidate edge, sort by chi/index, and truncate.
// ==========================================
void updateTopK(
    std::vector<NeighborEdge>& neighbours,
    int star_index,
    float chi,
    int k);

// ==========================================
// Function: Rebuild exact top-K neighbours on one active survivor set
// Method: Clear every stale list, compare only active same-chip cached windows,
//         and refill bounded neighbours in exact distance/index order.
// ==========================================
template <typename CandidateSelectionState>
void rebuildActiveKNN(
    const std::vector<int>& active_indices,
    std::vector<CandidateSelectionState>& selection,
    int k) {
    for (CandidateSelectionState& candidate : selection) {
        candidate.knn.clear();
    }
    if (k <= 0) return;

    for (std::size_t first = 0; first + 1 < active_indices.size(); ++first) {
        const int first_index = active_indices[first];
        if (first_index < 0
            || first_index >= static_cast<int>(selection.size())) {
            continue;
        }
        for (std::size_t second = first + 1;
             second < active_indices.size(); ++second) {
            const int second_index = active_indices[second];
            if (second_index < 0 || second_index == first_index
                || second_index >= static_cast<int>(selection.size())) {
                continue;
            }
            const float chi = normalizedChiDistance(
                selection[first_index].chi_window,
                selection[second_index].chi_window);
            updateTopK(
                selection[first_index].knn, second_index, chi, k);
            updateTopK(
                selection[second_index].knn, first_index, chi, k);
        }
    }
}

// ==========================================
// Structure: Describe one undirected grouping-graph edge
// Method: Address both endpoints by their original per-chip candidate indices.
// ==========================================
struct GraphEdge {
    int first = -1;
    int second = -1;
};

// ==========================================
// Structure: Return one connected stellar component
// Method: Preserve original candidate indices and count Gaia supporting labels.
// ==========================================
struct StarGroup {
    std::vector<int> members;
    int gaia_count = 0;
};

// ==========================================
// Function: Extract mutual-KNN graph edges among active candidates
// Method: Keep an undirected edge only when both retained top-K lists contain
//         the opposite endpoint and both endpoints survived the shared cut.
// ==========================================
std::vector<GraphEdge> buildMutualKNNEdges(
    const std::vector<int>& active_indices,
    const std::vector<std::vector<NeighborEdge>>& neighbours_by_star);

// ==========================================
// Function: Convert one same-chip graph into connected components
// Method: Use disjoint sets over active original indices and attach Gaia counts.
// ==========================================
std::vector<StarGroup> buildConnectedGroups(
    const std::vector<int>& active_indices,
    const std::vector<GraphEdge>& edges,
    const std::vector<bool>& gaia_matched);

// ==========================================
// Function: Select the shared main and eligible secondary stellar groups
// Method: Always keep the largest component and keep another component only
//         when both its relative size and Gaia-count requirements pass.
// ==========================================
std::vector<int> selectMainAndSecondaryGroups(
    const std::vector<StarGroup>& groups,
    double minimum_size_ratio,
    int minimum_gaia_count);

// ==========================================
// Function: Compute one analytic leave-one-out residual and prediction
// Method: Divide the ordinary residual by 1-h after finite leverage guards.
// ==========================================
bool computeAnalyticLOO(
    double observed,
    double fitted,
    double leverage,
    double minimum_denominator,
    double& loo_residual,
    double& loo_model);

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_STAR_SELECTION_HPP
