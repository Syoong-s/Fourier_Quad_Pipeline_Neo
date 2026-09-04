#ifndef PSF_STAR_SELECTION_HPP
#define PSF_STAR_SELECTION_HPP

#include <array>
#include <cstddef>
#include <istream>
#include <string>
#include <vector>

namespace PSFModel {
namespace Internal {

inline constexpr int PSFCountHistogramBinWidth = 2;

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
// Structure: Supply one quality-valid integer star area to locus estimation
// Method: Carry the exact exp(-1) Fourier-pixel count and Gaia support label.
// ==========================================
struct PSFCountSample {
    int star_area = 0;
    bool gaia_matched = false;
};

// ==========================================
// Structure: Configure exposure-wide integer star-area locus estimation
// Method: Group pilot clipping, zero-MAD quantiles, and final-cut controls.
// ==========================================
struct PSFCountLocusConfig {
    double pilot_clip_sigma = 0.0;
    int pilot_clip_iterations = 0;
    double zero_mad_quantile = -1.0;
    double histogram_range_sigma = 0.0;
    double locus_sigma = 0.0;
    int minimum_samples = 0;
    int minimum_gaia_matches = 0;
};

// ==========================================
// Structure: Publish a robust exposure-wide integer star-area locus
// Method: Store the median center, side-specific widths, and strict-cut bounds.
// ==========================================
struct PSFCountLocus {
    bool valid = false;
    double center = 0.0;
    double lower_width = 0.0;
    double upper_width = 0.0;
    double lower = 0.0;
    double upper = 0.0;
};

// ==========================================
// Structure: Publish exact integer star-area locus diagnostics
// Method: Keep raw counts immutable, expose the hole-filled working histogram,
//         and retain pilot, range, Gaia, minChi, and group accounting.
// ==========================================
struct PSFCountLocusDiagnostics {
    int sample_count = 0;
    int gaia_match_count = 0;
    bool pilot_uses_gaia = false;
    int pilot_input_count = 0;
    int pilot_retained_count = 0;
    double pilot_center = 0.0;
    double pilot_width = 0.0;
    double pilot_lower = 0.0;
    double pilot_upper = 0.0;
    bool pilot_uses_quantile_range = false;
    bool pilot_rejected_zero_mad_clip = false;
    int histogram_sample_count = 0;
    int histogram_below_count = 0;
    int histogram_above_count = 0;
    int gaia_histogram_sample_count = 0;
    int gaia_histogram_below_count = 0;
    int gaia_histogram_above_count = 0;
    int minchi_survivor_count = 0;
    int selected_group_count = 0;
    bool has_gaia_median = false;
    double gaia_median = 0.0;
    int histogram_first_count = 0;
    int histogram_last_count = -1;
    int peak_bin = -1;
    double mad_lower = 0.0;
    double mad_upper = 0.0;
    int left_elbow_bin = -1;
    int right_elbow_bin = -1;
    bool left_elbow_guard_applied = false;
    bool right_elbow_guard_applied = false;
    std::vector<double> histogram;
    std::vector<double> working_histogram;
    std::vector<double> smoothed_histogram;
    std::vector<double> gaia_histogram;
    std::vector<double> minchi_survivor_histogram;
    std::vector<double> selected_group_histogram;
};

// ==========================================
// Structure: Publish an inclusive PSF count-histogram bin range
// Method: Use invalid negative bounds when no deterministic range is available.
// ==========================================
struct PSFCountBinRange {
    int first = -1;
    int last = -1;
};

// ==========================================
// Structure: Publish independently detected outer histogram elbows
// Method: Keep each side unavailable until a positive-curvature candidate wins.
// ==========================================
struct PSFCountElbows {
    int left = -1;
    int right = -1;
};

// ==========================================
// Structure: Publish one re-absorbing asymmetric-MAD refinement
// Method: Return final center, side widths, and retained real-sample count.
// ==========================================
struct PSFCountRefinement {
    bool valid = false;
    double center = 0.0;
    double lower_width = 0.0;
    double upper_width = 0.0;
    int sample_count = 0;
};

// ==========================================
// Function: Return the nominal center of one fixed two-count histogram bin
// Method: Offset the first allowed count by two per bin and one half count.
// ==========================================
double psfCountHistogramBinCenter(
    int histogram_first_count,
    int bin);

// ==========================================
// Function: Find the selected peak's complete significant peak complex
// Method: Include every local peak strictly above H_selected/e, then descend
//         outward from the outermost qualifying peaks until the next rise.
// ==========================================
PSFCountBinRange findPSFCountPeakComplexBasin(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    int selected_peak);

// ==========================================
// Function: Find independent outer elbows around the selected count peak
// Method: Cross below ten percent of peak height on each side, then maximize
//         positive signed curvature outward with nearest-bin tie breaking.
// ==========================================
PSFCountElbows findPSFCountOuterElbows(
    const std::vector<double>& smoothed_histogram,
    int selected_peak);

// ==========================================
// Function: Refine a peak-basin seed with re-absorbing asymmetric MAD
// Method: Rebuild each pass from every real sample in the pilot histogram domain.
// ==========================================
PSFCountRefinement refinePSFCountPopulation(
    const std::vector<double>& seed_values,
    const std::vector<double>& domain_values,
    double locus_sigma,
    int iterations);

// ==========================================
// Function: Select one integer star-area histogram peak deterministically
// Method: Anchor Gaia eligibility to nominal bin-center distance plus one count,
//         then rank by Gaia count, density, exact distance, and lower bin.
// ==========================================
int selectPSFCountPeak(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    const std::vector<double>& gaia_histogram,
    int histogram_first_count,
    double pilot_center,
    bool pilot_uses_gaia);

// ==========================================
// Function: Fill bounded one- or two-bin holes in a working histogram
// Method: Detect zero runs only in the immutable raw histogram and linearly
//         interpolate between positive endpoints without chaining.
// ==========================================
std::vector<double> interpolateShortInternalHoles(
    const std::vector<double>& histogram);

// ==========================================
// Function: Estimate an exposure-wide integer star-area locus
// Method: Build fixed two-count bins, refine the significant peak complex with
//         re-absorbing asymmetric MAD, then apply independent outer elbows.
// ==========================================
bool estimatePSFCountLocus(
    const std::vector<PSFCountSample>& samples,
    const PSFCountLocusConfig& config,
    PSFCountLocus& locus,
    PSFCountLocusDiagnostics* diagnostics = nullptr);

// ==========================================
// Function: Count exp(-1)-threshold pixels in one square Fourier-power stamp
// Method: Compare every finite stamp value with the central value times exp(-1).
// ==========================================
int countPSFStarArea(
    const std::vector<float>& power,
    int stamp_side);

// ==========================================
// Function: Convert an exp(-1) star area to the historical PSF FWHM
// Method: Preserve the legacy area-minus-1e-5 formula and explicit pixel scale.
// ==========================================
double fwhmFromStarArea(
    double star_area,
    int stamp_side,
    double pixel_size);

// ==========================================
// Function: Populate the minChi-survivor integer-area histogram
// Method: Count all actual grouping inputs and bin only values on the science
//         grid without changing upstream count-locus diagnostics.
// ==========================================
void populateMinChiSurvivorCountHistogram(
    const std::vector<int>& minchi_star_areas,
    PSFCountLocusDiagnostics& diagnostics);

// ==========================================
// Function: Populate the historical pre-PRESS integer-area histogram
// Method: Count all selected stars for any grouping type on the science grid
//         without changing upstream count-locus diagnostics.
// ==========================================
void populateSelectedGroupCountHistogram(
    const std::vector<int>& selected_star_areas,
    PSFCountLocusDiagnostics& diagnostics);

// ==========================================
// Enumeration: Identify one adaptive upper-elbow histogram outcome
// Method: Preserve deterministic failure reasons for fail-open production logs.
// ==========================================
enum class PSFUpperElbowStatus {
    NoFiniteValues,
    InvalidConfig,
    InvalidInput,
    NoFDSamples,
    NonPositiveWidth,
    UnsafeBinCount,
    AllocationFailure,
    NoPeaks,
    NoElbow,
    Valid
};

// ==========================================
// Structure: Configure one Freedman-Diaconis upper-elbow histogram
// Method: Separate peak validity, FD sampling, zero-IQR fallback, and origin rules.
// ==========================================
struct PSFUpperElbowHistogramConfig {
    double valid_peak_fraction = 0.0;
    bool exclude_zero_from_fd = false;
    bool zero_iqr_uses_min_positive = false;
    bool force_zero_origin = false;
    // Zero uses the actual FD-distribution sample count.
    std::size_t fd_scale_sample_count = 0;
    // Candidates must be strictly below this smoothed main-peak fraction.
    double elbow_search_height_fraction = 0.0;
};

// ==========================================
// Structure: Publish one adaptive histogram and upper-elbow decision
// Method: Retain FD statistics, topology indices, and raw/smoothed diagnostics.
// ==========================================
struct PSFUpperElbowHistogramResult {
    bool valid = false;
    PSFUpperElbowStatus status = PSFUpperElbowStatus::NoFiniteValues;
    std::size_t finite_value_count = 0;
    std::size_t fd_sample_count = 0;
    std::size_t fd_scale_sample_count = 0;
    double fd_iqr = 0.0;
    double bin_origin = 0.0;
    double bin_width = 0.0;
    double cut = 0.0;
    double elbow_search_height = 0.0;
    int elbow_search_first_bin = -1;
    int elbow_search_last_bin = -1;
    std::size_t elbow_candidate_count = 0;
    int main_peak_bin = -1;
    int rightmost_valid_peak_bin = -1;
    int first_invalid_peak_bin = -1;
    int elbow_bin = -1;
    std::vector<int> peaks;
    std::vector<int> valid_peaks;
    std::vector<int> invalid_peaks;
    std::vector<double> histogram;
    std::vector<double> smoothed_histogram;
};

// ==========================================
// Function: Analyze one already-binned upper-elbow histogram
// Method: Expose smoothing, strict peak classes, the outer-tail height gate,
//         and the right positive-curvature elbow as a deterministic test seam.
// ==========================================
bool analyzePSFUpperElbowHistogram(
    const std::vector<double>& histogram,
    double bin_origin,
    double bin_width,
    double valid_peak_fraction,
    double elbow_search_height_fraction,
    PSFUpperElbowHistogramResult& result);

// ==========================================
// Function: Estimate an FD-histogram upper elbow from double samples
// Method: Apply strict peak validity and the candidate-height gate, then
//         maximize positive curvature before the first invalid peak.
// ==========================================
bool estimatePSFUpperElbowCut(
    const std::vector<double>& values,
    const PSFUpperElbowHistogramConfig& config,
    PSFUpperElbowHistogramResult& result);

// ==========================================
// Function: Estimate an FD-histogram upper elbow from float samples
// Method: Preserve compact pair-chi storage while sharing the double algorithm.
// ==========================================
bool estimatePSFUpperElbowCut(
    const std::vector<float>& values,
    const PSFUpperElbowHistogramConfig& config,
    PSFUpperElbowHistogramResult& result);

// ==========================================
// Function: Return a stable adaptive-histogram status label
// Method: Map every public status enumerator to one diagnostic token.
// ==========================================
const char* psfUpperElbowStatusName(PSFUpperElbowStatus status);

// ==========================================
// Structure: Publish one chip's Type-3 fraction-gate decision
// Method: Align Boolean survivors with active inputs and report minimum rejection.
// ==========================================
struct PSFType3ChipSelection {
    std::vector<bool> selected;
    std::size_t finite_pair_count = 0;
    std::size_t retained_count = 0;
    bool rejected_by_minimum = false;
};

// ==========================================
// Function: Classify one finite pair against the Type-3 upper cut
// Method: Mark bad only when both values are finite and chi is strictly greater.
// ==========================================
bool isPSFType3BadPair(double chi, double pair_chi_cut);

// ==========================================
// Function: Apply the Type-3 bad-pair-fraction gate to one chip
// Method: Require finite denominators, reject only fractions strictly above an
//         optional cut, then atomically enforce the retained-star minimum.
// ==========================================
PSFType3ChipSelection selectPSFType3FractionSurvivors(
    const std::vector<double>& bad_pair_fractions,
    const std::vector<bool>& has_finite_pair_denominator,
    bool apply_fraction_cut,
    double fraction_cut,
    int minimum_retained);

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
// Structure: Identify one exposure-level minChi reference candidate
// Method: Carry deterministic chip/star keys and the legacy index-7 size rank.
// ==========================================
struct MinChiReferenceCandidate {
    int chip_index = -1;
    int star_index = -1;
    double size = 0.0;
};

// ==========================================
// Function: Select capped exposure-wide large-size minChi references
// Method: Sort the top fraction by size/chip/star, then apply a per-chip cap
//         without filling from candidates below the exposure-wide pool.
// ==========================================
std::vector<MinChiReferenceCandidate> selectMinChiReferenceStars(
    const std::vector<MinChiReferenceCandidate>& locus_candidates,
    double reference_fraction,
    int maximum_per_chip);

// ==========================================
// Structure: View one candidate during the same-chip minChi pair pass
// Method: Borrow its normalized chi window and carry locus/reference labels.
// ==========================================
struct MinChiCandidateView {
    const std::vector<float>* chi_window = nullptr;
    bool in_size_locus = false;
    bool is_reference = false;
};

// ==========================================
// Structure: Return nearest-neighbour distances and reference-pair samples
// Method: Align minChi with input candidates and retain every qualifying
//         unordered pair distance exactly once for exposure thresholding.
// ==========================================
struct MinChiPairResult {
    std::vector<float> min_chi;
    std::vector<float> threshold_pair_chi;
};

// ==========================================
// Function: Compute one chip's minChi values and threshold-pair sample
// Method: Visit every unordered locus-locus pair once, update both endpoints,
//         and sample the distance when either endpoint is a reference.
// ==========================================
MinChiPairResult computeMinChiAndThresholdPairs(
    const std::vector<MinChiCandidateView>& candidates);

// ==========================================
// Function: Compute the exact normalized PSF chi distance
// Method: Apply the legacy sqrt(sum squared difference / mean signed flux)
//         directly to two cached central windows.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second);

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

// ==========================================
// Function: Convert raw analytic PRESS to its leverage-standardized score
// Method: Multiply by sqrt(1-h) after the same finite denominator guard used
//         by the analytic leave-one-out calculation.
// ==========================================
bool computeLeverageStandardizedPress(
    double raw_press,
    double leverage,
    double minimum_denominator,
    double& standardized_press);

enum class PressRemovalDecision {
    Disabled,
    NoOutliers,
    TooManyOutliers,
    WouldUnderrunMinimum,
    Apply
};

// ==========================================
// Function: Decide whether an optional PRESS removal may be attempted
// Method: Apply the rejection switch, outlier count, removal cap, and minimum
//         retained-star guard without mutating any fitting or selection state.
// ==========================================
PressRemovalDecision decidePressRemoval(
    bool rejection_enabled,
    int initial_count,
    int flagged_count,
    int minimum_count,
    int maximum_removals);

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_STAR_SELECTION_HPP
