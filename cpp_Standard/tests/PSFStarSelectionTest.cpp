#include "process_main/PSFStarSelection.hpp"
#include "process_main/PSFModelState.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
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
// Function: Build the production-shaped integer star-area locus configuration
// Method: Keep the current pilot, quantile, final-cut, and Gaia controls fixed.
// ==========================================
PSFCountLocusConfig countLocusConfig() {
    return {3.0, 3, 0.05, 5.0, 4.0, 30, 5};
}

// ==========================================
// Function: Verify exact star-area measurement and historical FWHM conversion
// Method: Use a controlled stamp, check index 12, and preserve prior row fields.
// ==========================================
void testStarAreaMeasurementAndStorage() {
    std::vector<float> power(25, 0.0f);
    power[12] = 10.0f;
    power[6] = 4.0f;
    power[7] = 5.0f;
    power[11] = 6.0f;
    require(countPSFStarArea(power, 5) == 4,
            "star_area must count exact central exp(-1) threshold pixels");
    const double fwhm = fwhmFromStarArea(4.0, 5, 0.2628);
    constexpr double pi = 3.14159265358979323846;
    const double legacy_area = 4.0 - 1.0e-5;
    const double expected_fwhm =
        (5.0 / (2.0 * pi) / std::sqrt(legacy_area / pi))
        * 2.0 * std::sqrt(2.0 * std::log(2.0)) * 0.2628;
    require(std::abs(fwhm - expected_fwhm) < 1.0e-14
                && fwhmFromStarArea(9.0, 5, 0.2628) < fwhm,
            "star_area conversion must preserve the historical FWHM exactly");

    using ChipState = PSFModel::Internal::ChipPSFState;
    static_assert(ChipState::star_area_index == 12);
    static_assert(ChipState::star_area_index
        < static_cast<int>(std::tuple_size<ChipState::StarRow>::value));
    ChipState::StarRow row{};
    row[7] = 17.0;
    row[10] = fwhm;
    row[11] = 0.25;
    row[ChipState::star_area_index] = 4.0;
    require(row[7] == 17.0 && row[10] == fwhm && row[11] == 0.25
                && row[ChipState::star_area_index] == 4.0,
            "index-12 star_area storage must not alter existing PSF fields");
}

// ==========================================
// Function: Verify bounded non-chaining interpolation of short count-bin holes
// Method: Exercise one/two-bin fills, long/edge gaps, and raw-run detection.
// ==========================================
void testCountHoleInterpolation() {
    require(interpolateShortInternalHoles({100.0, 0.0, 80.0})
                == std::vector<double>({100.0, 90.0, 80.0}),
            "one two-count-bin internal hole must be linearly interpolated");
    require(interpolateShortInternalHoles({100.0, 0.0, 0.0, 70.0})
                == std::vector<double>({100.0, 90.0, 80.0, 70.0}),
            "two two-count-bin internal holes must be linearly interpolated");
    require(interpolateShortInternalHoles({100.0, 0.0, 0.0, 0.0, 70.0})
                == std::vector<double>({100.0, 0.0, 0.0, 0.0, 70.0}),
            "three-bin gaps must remain raw zeros");
    require(interpolateShortInternalHoles({0.0, 0.0, 50.0, 0.0})
                == std::vector<double>({0.0, 0.0, 50.0, 0.0}),
            "edge zeros must never be extrapolated");
    require(interpolateShortInternalHoles(
                {10.0, 0.0, 8.0, 0.0, 0.0, 0.0, 4.0})
                == std::vector<double>(
                    {10.0, 9.0, 8.0, 0.0, 0.0, 0.0, 4.0}),
            "a filled short hole must not chain into a raw long gap");
}

// ==========================================
// Function: Verify fixed one-count Gaia near-tie peak selection on two-count bins
// Method: Cover every ranking tier, global anchoring, order invariance, and the
//         unchanged no-Gaia maximum-density route.
// ==========================================
void testCountGaiaPeakTieBreaks() {
    std::vector<double> smoothed(8, 0.0);
    std::vector<double> gaia(8, 0.0);
    smoothed[2] = 6.0;
    smoothed[4] = 5.0;
    gaia[2] = 2.0;
    gaia[4] = 4.0;
    require(selectPSFCountPeak(
                {2, 4}, smoothed, gaia, 30, 36.5, true) == 4,
            "equal-distance count peaks must prefer higher raw Gaia count");
    smoothed[6] = 20.0;
    gaia[6] = 100.0;
    require(selectPSFCountPeak(
                {2, 4, 6}, smoothed, gaia, 30, 37.9, true) == 4
                && selectPSFCountPeak(
                    {6, 4, 2}, smoothed, gaia, 30, 37.9, true) == 4,
            "one-count eligibility must use a global anchor without chaining");
    gaia[4] = gaia[2];
    smoothed[4] = 7.0;
    require(selectPSFCountPeak(
                {2, 4}, smoothed, gaia, 30, 36.5, true) == 4,
            "Gaia ties must prefer higher smoothed density");
    smoothed[4] = smoothed[2];
    require(selectPSFCountPeak(
                {2, 4}, smoothed, gaia, 30, 36.25, true) == 2,
            "density ties must prefer exact pilot distance");
    require(selectPSFCountPeak(
                {4, 2}, smoothed, gaia, 30, 36.5, true) == 2,
            "complete ties must prefer the lower count level");
    require(selectPSFCountPeak(
                {2, 6}, smoothed, gaia, 30, 33.0, false) == 6,
            "no-Gaia selection must retain maximum smoothed density");
}

// ==========================================
// Function: Verify integer pilot, histogram, locus, and width-floor behavior
// Method: Exercise zero MAD, exact two-count bins, immutable raw diagnostics,
//         single-level support, strict cuts, and insufficient-sample failure.
// ==========================================
void testPSFCountLocus() {
    const PSFCountLocusConfig config = countLocusConfig();
    PSFCountLocus locus;
    PSFCountLocusDiagnostics diagnostics;
    std::vector<PSFCountSample> zero_mad;
    zero_mad.insert(zero_mad.end(), 60, {37, true});
    zero_mad.insert(zero_mad.end(), 30, {40, true});
    zero_mad.insert(zero_mad.end(), 10, {44, true});
    require(estimatePSFCountLocus(zero_mad, config, locus, &diagnostics),
            "integer zero-MAD Gaia pilot must produce a count locus");
    require(diagnostics.pilot_uses_gaia
                && diagnostics.pilot_retained_count == 100
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range
                && diagnostics.pilot_lower == 37.0
                && diagnostics.pilot_upper == 44.0,
            "zero-MAD pilot must retain all samples and use unpadded Q05-Q95");
    require(diagnostics.histogram_first_count == 37
                && diagnostics.histogram_last_count == 44
                && diagnostics.histogram.size() == 4
                && diagnostics.histogram[0] == 60.0
                && diagnostics.histogram[1] == 30.0
                && diagnostics.histogram[3] == 10.0
                && diagnostics.working_histogram != diagnostics.histogram,
            "count histogram must map adjacent integer levels into width-two bins");

    PSFCountLocusConfig custom_quantile_config = config;
    custom_quantile_config.zero_mad_quantile = 0.20;
    std::vector<PSFCountSample> custom_quantile;
    custom_quantile.insert(custom_quantile.end(), 60, {37, false});
    custom_quantile.insert(custom_quantile.end(), 20, {40, false});
    custom_quantile.insert(custom_quantile.end(), 20, {44, false});
    require(estimatePSFCountLocus(
                custom_quantile,
                custom_quantile_config,
                locus,
                &diagnostics)
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_lower == 37.0
                && std::abs(diagnostics.pilot_upper - 40.8) < 1.0e-12,
            "one zero-MAD quantile must control Q(q) and Q(1-q) bounds");

    std::vector<PSFCountSample> repeated(40, {37, false});
    require(estimatePSFCountLocus(repeated, config, locus, &diagnostics)
                && diagnostics.histogram_first_count == 37
                && diagnostics.histogram_last_count == 37
                && diagnostics.histogram.size() == 1
                && diagnostics.histogram[0] == 40.0
                && locus.center == 37.0
                && locus.lower_width == 1.0
                && locus.upper_width == 1.0
                && locus.lower == 33.0
                && locus.upper == 41.0,
            "single count-level support must keep deterministic one-count MAD floors");
    require(!(33.0 > locus.lower && 33.0 < locus.upper)
                && (37.0 > locus.lower && 37.0 < locus.upper)
                && !(41.0 > locus.lower && 41.0 < locus.upper),
            "final production star-area selection must remain strict");

    std::vector<PSFCountSample> guarded;
    guarded.insert(guarded.end(), 6, {30, false});
    guarded.insert(guarded.end(), 88, {40, false});
    guarded.insert(guarded.end(), 6, {50, false});
    require(estimatePSFCountLocus(guarded, config, locus, &diagnostics)
                && diagnostics.histogram_first_count == 30
                && diagnostics.histogram_last_count == 50
                && diagnostics.histogram.size() == 11
                && diagnostics.left_elbow_bin == 2
                && diagnostics.right_elbow_bin == 8
                && diagnostics.mad_lower == 36.0
                && diagnostics.mad_upper == 44.0
                && diagnostics.left_elbow_guard_applied
                && diagnostics.right_elbow_guard_applied
                && locus.lower == 34.5
                && locus.upper == 46.5
                && locus.center == 40.0,
            "outer elbows must widen final cuts without changing MAD statistics");

    std::vector<PSFCountSample> right_skew;
    std::vector<PSFCountSample> left_skew;
    for (int step = 1; step <= 8; ++step) {
        const int copies = 9 - step;
        for (int copy = 0; copy < copies; ++copy) {
            right_skew.push_back({40 - step, false});
            right_skew.push_back({40 + 2 * step, false});
            left_skew.push_back({40 - 2 * step, false});
            left_skew.push_back({40 + step, false});
        }
    }
    right_skew.insert(right_skew.end(), 30, {40, false});
    left_skew.insert(left_skew.end(), 30, {40, false});
    PSFCountLocus right_locus;
    PSFCountLocus left_locus;
    require(estimatePSFCountLocus(right_skew, config, right_locus)
                && estimatePSFCountLocus(left_skew, config, left_locus)
                && right_locus.upper_width > right_locus.lower_width
                && left_locus.lower_width > left_locus.upper_width,
            "integer asymmetric MAD must broaden only the populated tail side");

    repeated.resize(29);
    require(!estimatePSFCountLocus(repeated, config, locus, &diagnostics)
                && diagnostics.sample_count == 29
                && diagnostics.histogram.empty(),
            "count locus must reject fewer than its configured samples");
}

// ==========================================
// Function: Verify peak-complex, elbow, and re-absorbing refinement helpers
// Method: Lock strict height/crossing rules, signed curvature, nearest ties,
//         unavailable sides, nominal centers, and pilot-domain re-entry.
// ==========================================
void testPSFCountTopologyAndRefinement() {
    require(psfCountHistogramBinCenter(30, 0) == 30.5
                && psfCountHistogramBinCenter(30, 4) == 38.5,
            "two-count bins must use their nominal half-count centers");

    const PSFCountBinRange complex = findPSFCountPeakComplexBasin(
        {2, 4, 6},
        {5.0, 4.0, 40.0, 10.0, 100.0, 20.0, 50.0, 4.0, 6.0},
        4);
    require(complex.first == 1 && complex.last == 7,
            "all peaks above H_selected/e must form one valley-agnostic complex");

    const double exact_floor = 100.0 * std::exp(-1.0);
    const PSFCountBinRange strict = findPSFCountPeakComplexBasin(
        {1, 3, 5},
        {5.0, exact_floor, 1.0, 100.0, 1.0, 20.0, 5.0},
        3);
    require(strict.first == 2 && strict.last == 4,
            "a peak exactly at H_selected/e must be excluded from the complex");

    const PSFCountElbows elbows = findPSFCountOuterElbows(
        {0.0, 1.0, 8.0, 20.0, 100.0, 20.0, 8.0, 1.0, 0.0},
        4);
    require(elbows.left == 1 && elbows.right == 7,
            "elbow search must retain candidates from the crossing to each edge");

    const PSFCountElbows tied = findPSFCountOuterElbows(
        {2.0, 2.0, 8.0, 20.0, 100.0, 20.0, 8.0, 2.0, 2.0},
        4);
    require(tied.left == 2 && tied.right == 6,
            "equal positive curvature must prefer the candidate nearest the peak");

    const PSFCountElbows unavailable = findPSFCountOuterElbows(
        {9.0, 9.9, 10.0, 100.0, 10.0, 9.9, 9.0},
        3);
    require(unavailable.left == -1 && unavailable.right == -1,
            "nonpositive curvature after a strict crossing must leave elbows unavailable");
    const PSFCountElbows no_crossing = findPSFCountOuterElbows(
        {20.0, 30.0, 100.0, 30.0, 20.0},
        2);
    require(no_crossing.left == -1 && no_crossing.right == -1,
            "a side without a below-ten-percent crossing must stay unavailable");
    const PSFCountElbows edge_crossing = findPSFCountOuterElbows(
        {0.0, 100.0, 20.0},
        1);
    require(edge_crossing.left == -1,
            "an edge crossing without an interior curvature bin must stay unavailable");

    const PSFCountRefinement refinement = refinePSFCountPopulation(
        {10.0, 10.0, 11.0, 12.0, 12.0},
        {10.0, 10.0, 11.0, 12.0, 12.0, 13.0},
        2.0,
        2);
    require(refinement.valid && refinement.sample_count == 6
                && refinement.center == 11.5,
            "MAD passes must re-absorb eligible real samples from the domain");
}

// ==========================================
// Function: Verify post-minChi and selected diagnostics use the science grid
// Method: Check exact bins, nested subset bounds, out-of-range accounting, and
//         that every upstream count-locus diagnostic remains unchanged.
// ==========================================
void testCountOverlayHistograms() {
    PSFCountLocusDiagnostics diagnostics;
    diagnostics.sample_count = 14;
    diagnostics.gaia_match_count = 3;
    diagnostics.pilot_center = 31.0;
    diagnostics.pilot_lower = 30.0;
    diagnostics.pilot_upper = 33.0;
    diagnostics.histogram_sample_count = 14;
    diagnostics.histogram_first_count = 30;
    diagnostics.histogram_last_count = 37;
    diagnostics.peak_bin = 1;
    diagnostics.histogram = {4.0, 5.0, 3.0, 2.0};
    diagnostics.working_histogram = {4.0, 5.0, 3.0, 2.0};
    diagnostics.smoothed_histogram = {4.5, 4.0, 3.0, 2.5};
    diagnostics.gaia_histogram = {0.0, 2.0, 1.0, 0.0};
    const PSFCountLocusDiagnostics baseline = diagnostics;

    populateMinChiSurvivorCountHistogram({30, 31, 31, 37}, diagnostics);
    require(diagnostics.minchi_survivor_count == 4
                && diagnostics.minchi_survivor_histogram
                    == std::vector<double>({3.0, 0.0, 0.0, 1.0}),
            "minChi survivors must use the fixed two-count science grid");
    for (std::size_t bin = 0; bin < diagnostics.histogram.size(); ++bin) {
        require(diagnostics.minchi_survivor_histogram[bin]
                    <= diagnostics.histogram[bin],
                "minChi bins must remain subsets of candidate bins");
    }
    require(diagnostics.sample_count == baseline.sample_count
                && diagnostics.gaia_match_count == baseline.gaia_match_count
                && diagnostics.pilot_center == baseline.pilot_center
                && diagnostics.pilot_lower == baseline.pilot_lower
                && diagnostics.pilot_upper == baseline.pilot_upper
                && diagnostics.histogram_sample_count
                    == baseline.histogram_sample_count
                && diagnostics.histogram_first_count
                    == baseline.histogram_first_count
                && diagnostics.histogram_last_count
                    == baseline.histogram_last_count
                && diagnostics.peak_bin == baseline.peak_bin
                && diagnostics.histogram == baseline.histogram
                && diagnostics.working_histogram
                    == baseline.working_histogram
                && diagnostics.smoothed_histogram
                    == baseline.smoothed_histogram
                && diagnostics.gaia_histogram == baseline.gaia_histogram,
            "minChi histogram completion must not alter count science");

    populateSelectedGroupCountHistogram({31, 37}, diagnostics);
    require(diagnostics.selected_group_count == 2
                && diagnostics.selected_group_histogram
                    == std::vector<double>({1.0, 0.0, 0.0, 1.0}),
            "selected stars must use the fixed two-count science grid");
    for (std::size_t bin = 0; bin < diagnostics.histogram.size(); ++bin) {
        require(diagnostics.selected_group_histogram[bin]
                    <= diagnostics.minchi_survivor_histogram[bin],
                "selected count bins must remain subsets of minChi bins");
    }

    populateMinChiSurvivorCountHistogram({29, 30, 31, 38}, diagnostics);
    require(diagnostics.minchi_survivor_count == 4
                && diagnostics.minchi_survivor_histogram
                    == std::vector<double>({2.0, 0.0, 0.0, 0.0}),
            "out-of-grid stars must count as minChi survivors without SVG bins");

    populateSelectedGroupCountHistogram({29, 30, 31, 38}, diagnostics);
    require(diagnostics.selected_group_count == 4
                && diagnostics.selected_group_histogram
                    == std::vector<double>({2.0, 0.0, 0.0, 0.0}),
            "out-of-grid stars must count as selected without entering SVG bins");
}

// ==========================================
// Function: Verify adaptive FD grids and gated upper-elbow topology
// Method: Cover FD geometry, strict peak classes, the ten-percent boundary,
//         positive curvature, first-invalid bounds, ties, and fail-open states.
// ==========================================
void testAdaptiveUpperElbowHistogram() {
    const PSFUpperElbowHistogramConfig pair_config = {
        std::exp(-1.0), false, false, false, 0U, 0.10};
    PSFUpperElbowHistogramResult result;
    const std::vector<double> ordinary = {
        0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0,
        std::numeric_limits<double>::quiet_NaN()};
    estimatePSFUpperElbowCut(ordinary, pair_config, result);
    require(result.finite_value_count == 8
                && result.fd_sample_count == 8
                && result.fd_scale_sample_count == 8
                && std::abs(result.fd_iqr - 1.5) < 1.0e-12
                && std::abs(result.bin_width - 1.5) < 1.0e-12
                && result.bin_origin == 0.0
                && result.histogram.size() == 3
                && std::accumulate(
                    result.histogram.begin(), result.histogram.end(), 0.0)
                    == 8.0
                && result.histogram.back() == 2.0,
            "ordinary FD bins must filter nonfinite input and include the maximum");

    PSFUpperElbowHistogramConfig invalid_gate_config = pair_config;
    invalid_gate_config.elbow_search_height_fraction = 1.0;
    require(!estimatePSFUpperElbowCut(
                ordinary, invalid_gate_config, result)
                && result.status == PSFUpperElbowStatus::InvalidConfig,
            "the elbow search fraction must lie strictly between zero and one");

    PSFUpperElbowHistogramConfig explicit_scale_config = pair_config;
    explicit_scale_config.fd_scale_sample_count = 4U;
    estimatePSFUpperElbowCut(ordinary, explicit_scale_config, result);
    const double explicit_scale_width = 3.0 / std::cbrt(4.0);
    require(result.fd_sample_count == 8
                && result.fd_scale_sample_count == 4
                && std::abs(result.fd_iqr - 1.5) < 1.0e-12
                && std::abs(result.bin_width - explicit_scale_width)
                    < 1.0e-12,
            "explicit FD scale count must change only the cube-root factor");

    PSFUpperElbowHistogramConfig scale_1000_config = pair_config;
    scale_1000_config.fd_scale_sample_count = 1000U;
    estimatePSFUpperElbowCut(ordinary, scale_1000_config, result);
    const double width_1000 = result.bin_width;
    PSFUpperElbowHistogramConfig scale_10000_config = pair_config;
    scale_10000_config.fd_scale_sample_count = 10000U;
    estimatePSFUpperElbowCut(ordinary, scale_10000_config, result);
    require(result.fd_sample_count == 8
                && result.fd_scale_sample_count == 10000
                && std::abs(
                    result.bin_width / width_1000 - std::cbrt(0.1))
                    < 1.0e-12,
            "FD width must scale with the explicit star count to power -1/3");

    std::vector<double> outlier = {
        0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 100.0};
    estimatePSFUpperElbowCut(outlier, pair_config, result);
    require(result.finite_value_count == 9
                && result.fd_iqr == 2.0
                && result.histogram.size() > 3
                && std::accumulate(
                    result.histogram.begin(), result.histogram.end(), 0.0)
                    == 9.0
                && result.histogram.back() == 1.0,
            "an outlier must extend the true FD grid without clipping its bin");

    require(!estimatePSFUpperElbowCut(
                std::vector<double>({1.0, 1.0, 1.0}),
                pair_config,
                result)
                && result.status == PSFUpperElbowStatus::NonPositiveWidth,
            "pair-chi zero IQR must fail open without an invented width");

    const PSFUpperElbowHistogramConfig fraction_config = {
        0.10, true, true, true, 0U, 0.10};
    require(!estimatePSFUpperElbowCut(
                std::vector<double>({0.0, 0.0, 0.0}),
                fraction_config,
                result)
                && result.status == PSFUpperElbowStatus::NoFDSamples,
            "all-zero fractions must remain distinguishable from estimator failure");
    estimatePSFUpperElbowCut(
        std::vector<double>({
            0.0, 0.0, 0.25, 0.25, 0.25, 0.25,
            std::numeric_limits<double>::infinity()}),
        fraction_config,
        result);
    require(result.finite_value_count == 6
                && result.fd_sample_count == 4
                && result.fd_scale_sample_count == 4
                && result.fd_iqr == 0.0
                && result.bin_width == 0.25
                && result.bin_origin == 0.0
                && result.histogram == std::vector<double>({2.0, 4.0}),
            "fraction FD must exclude zeros only from width estimation");

    std::vector<double> unsafe(80, 1.0e-300);
    unsafe.push_back(1.0);
    require(!estimatePSFUpperElbowCut(
                unsafe, fraction_config, result)
                && result.status == PSFUpperElbowStatus::UnsafeBinCount,
            "an unrepresentable FD grid must take the deterministic fail-open path");

    require(!analyzePSFUpperElbowHistogram(
                {5.0, 5.0, 5.0, 5.0}, 0.0, 1.0, 0.5, 0.10, result)
                && result.peaks == std::vector<int>({1})
                && result.status == PSFUpperElbowStatus::NoElbow,
            "an even plateau must collapse to its lower middle bin");
    require(analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 10.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 5.0, 0.0, 0.0, 0.0},
                0.0, 1.0, 0.5, 0.10, result)
                && result.main_peak_bin == 2
                && result.valid_peaks == std::vector<int>({2})
                && result.invalid_peaks == std::vector<int>({8})
                && result.first_invalid_peak_bin == 8
                && result.elbow_search_first_bin == 5
                && result.elbow_search_last_bin == 5
                && result.elbow_candidate_count == 1
                && result.elbow_bin == 5
                && result.cut == 5.5,
            "strict peak equality must retain the first-invalid search bound");
    require(!analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0},
                0.0, 1.0, 0.5, 0.10, result)
                && result.main_peak_bin == 2
                && result.rightmost_valid_peak_bin == 6
                && result.elbow_candidate_count == 0
                && result.status == PSFUpperElbowStatus::NoElbow,
            "a search domain without bins below ten percent must fail open");
    require(analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 100.0, 0.0, 0.0, 15.0,
                 0.0, 0.0, 5.0, 0.0, 0.0, 0.0},
                0.0, 1.0, 0.5, 0.10, result)
                && result.elbow_bin == 8
                && result.elbow_search_first_bin == 7
                && result.elbow_search_last_bin == 10
                && result.elbow_candidate_count == 4
                && result.smoothed_histogram[5]
                    > result.elbow_search_height
                && result.smoothed_histogram[result.elbow_bin]
                    < result.elbow_search_height,
            "the strongest high-count curvature must yield to a gated tail elbow");
    require(analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 30.0, 7.0, 6.0, 9.0, 3.0,
                 4.0, 6.0, 2.0, 0.0, 0.0, 1.0, 4.0},
                0.0, 1.0, 0.5, 0.10, result)
                && result.first_invalid_peak_bin == 13
                && result.elbow_search_first_bin == 11
                && result.elbow_search_last_bin == 11
                && result.elbow_candidate_count == 1
                && result.elbow_bin == 11
                && result.smoothed_histogram[10]
                    == result.elbow_search_height,
            "a bin exactly at ten percent must be excluded from the search");
    require(!analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 100.0, 19.0, 15.0, 22.0, 12.0, 2.0,
                 11.0, 21.0, 2.0, 1.0, 9.0, 6.0, 0.0, 0.0},
                0.0, 1.0, 0.5, 0.10, result)
                && result.elbow_candidate_count == 1
                && result.elbow_search_first_bin == 14
                && result.elbow_search_last_bin == 14
                && result.status == PSFUpperElbowStatus::NoElbow,
            "gated bins without positive curvature must not invent an elbow");
    require(analyzePSFUpperElbowHistogram(
                {2.0, 0.0, 310.0, 0.0, 5.0, 7.0,
                 6.0, 7.0, 1.0, 0.0, 7.0, 5.0},
                0.0, 1.0, 0.05, 0.10, result)
                && result.rightmost_valid_peak_bin == 6
                && result.first_invalid_peak_bin == 11
                && result.elbow_search_first_bin == 7
                && result.elbow_search_last_bin == 10
                && result.elbow_candidate_count == 4
                && result.elbow_bin == 8,
            "equal gated curvatures must keep the nearer lower bin");
    require(analyzePSFUpperElbowHistogram(
                {0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 5.0, 0.0, 0.0, 0.0},
                0.0, 1.0, 0.10, 0.10, result)
                && result.valid_peaks == std::vector<int>({2})
                && result.invalid_peaks == std::vector<int>({8})
                && result.elbow_bin == 5
                && result.smoothed_histogram[result.elbow_bin]
                    < result.elbow_search_height,
            "fraction peak validity and elbow height gates must remain independent");
    require(std::string(psfUpperElbowStatusName(
                PSFUpperElbowStatus::UnsafeBinCount)) == "UNSAFE_BIN_COUNT",
            "adaptive histogram status labels must remain stable for logs");
}

// ==========================================
// Function: Verify the pure per-chip Type-3 fraction gate
// Method: Lock finite-denominator eligibility, strict-above rejection, equality
//         retention, fail-open cuts, and atomic minimum-star enforcement.
// ==========================================
void testType3FractionSelection() {
    const std::vector<double> fractions = {
        0.0, 0.20, 0.30, std::numeric_limits<double>::quiet_NaN()};
    const std::vector<bool> denominators = {true, true, true, false};
    require(!isPSFType3BadPair(0.20, 0.20)
                && isPSFType3BadPair(
                    std::nextafter(0.20, 1.0), 0.20)
                && !isPSFType3BadPair(
                    std::numeric_limits<double>::infinity(), 0.20),
            "Type-3 pair rejection must be finite and strictly above its cut");
    PSFType3ChipSelection selected = selectPSFType3FractionSurvivors(
        fractions, denominators, true, 0.20, 2);
    require(selected.selected
                == std::vector<bool>({true, true, false, false})
                && selected.finite_pair_count == 3
                && selected.retained_count == 2
                && !selected.rejected_by_minimum,
            "Type-3 must retain a fraction equal to the cut and reject only above");

    selected = selectPSFType3FractionSurvivors(
        fractions, denominators, false, 0.0, 3);
    require(selected.selected
                == std::vector<bool>({true, true, true, false})
                && selected.retained_count == 3,
            "failed fraction estimation must retain every finite-denominator star");

    selected = selectPSFType3FractionSurvivors(
        fractions, denominators, true, 0.20, 3);
    require(selected.selected
                == std::vector<bool>({false, false, false, false})
                && selected.retained_count == 0
                && selected.rejected_by_minimum,
            "a sub-minimum Type-3 result must reject the whole chip atomically");
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
    testStarAreaMeasurementAndStorage();
    testCountHoleInterpolation();
    testCountGaiaPeakTieBreaks();
    testPSFCountLocus();
    testPSFCountTopologyAndRefinement();
    testCountOverlayHistograms();
    testAdaptiveUpperElbowHistogram();
    testType3FractionSelection();
    testGaiaParsingAndMatching();
    testGrouping();
    testKNNRebuiltAfterMinChiCut();
    testMinChiReferencesAndPairs();
    testAnalyticLOO();
    testPressStandardizationAndDecision();
    std::cout << "PSF star-selection tests passed\n";
    return EXIT_SUCCESS;
}
