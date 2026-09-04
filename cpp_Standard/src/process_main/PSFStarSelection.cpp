#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>

namespace PSFModel {
namespace Internal {
namespace {

// ==========================================
// Function: Return the median of a sorted finite vector
// Method: Select the middle element or average the two central elements.
// ==========================================
double sortedMedian(const std::vector<double>& sorted) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t middle = sorted.size() / 2;
    if (sorted.size() % 2 == 1) return sorted[middle];
    return 0.5 * (sorted[middle - 1] + sorted[middle]);
}

// ==========================================
// Function: Interpolate one quantile from a sorted finite vector
// Method: Map the bounded fraction onto [0, N-1] and linearly blend neighbours.
// ==========================================
double sortedQuantile(
    const std::vector<double>& sorted,
    double fraction) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    fraction = std::clamp(fraction, 0.0, 1.0);
    const double position = fraction
        * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    if (lower == upper) return sorted[lower];
    const double weight = position - static_cast<double>(lower);
    return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
}

// ==========================================
// Function: Estimate median and scaled median absolute deviation
// Method: Sort finite samples and their deviations and multiply MAD by 1.4826.
// ==========================================
bool medianAndMad(
    const std::vector<double>& values,
    double& median,
    double& scaled_mad) {
    if (values.empty()) return false;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    median = sortedMedian(sorted);
    std::vector<double> deviations;
    deviations.reserve(sorted.size());
    for (double value : sorted) deviations.push_back(std::abs(value - median));
    std::sort(deviations.begin(), deviations.end());
    scaled_mad = 1.4826 * sortedMedian(deviations);
    return std::isfinite(median) && std::isfinite(scaled_mad);
}

// ==========================================
// Function: Estimate median and side-specific scaled median deviations
// Method: Exclude center duplicates from both side MADs and apply the fixed
//         one-count measurement-resolution floor independently to both sides.
// ==========================================
bool medianAndAsymmetricMad(
    const std::vector<double>& values,
    double& median,
    double& lower_width,
    double& upper_width) {
    if (values.empty()) return false;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    median = sortedMedian(sorted);

    std::vector<double> lower_deviations;
    std::vector<double> upper_deviations;
    lower_deviations.reserve(sorted.size());
    upper_deviations.reserve(sorted.size());
    for (double value : sorted) {
        if (value < median) {
            lower_deviations.push_back(median - value);
        } else if (value > median) {
            upper_deviations.push_back(value - median);
        }
    }
    std::sort(lower_deviations.begin(), lower_deviations.end());
    std::sort(upper_deviations.begin(), upper_deviations.end());
    lower_width = lower_deviations.empty()
        ? 0.0
        : 1.4826 * sortedMedian(lower_deviations);
    upper_width = upper_deviations.empty()
        ? 0.0
        : 1.4826 * sortedMedian(upper_deviations);
    lower_width = std::max(lower_width, 1.0);
    upper_width = std::max(upper_width, 1.0);
    return std::isfinite(median)
        && std::isfinite(lower_width) && lower_width > 0.0
        && std::isfinite(upper_width) && upper_width > 0.0;
}

// ==========================================
// Structure: Publish one first-layer integer star-area pilot and its bounds
// Method: Keep center/width, retained count, and zero-MAD branch decisions atomic.
// ==========================================
struct PSFCountPilotEstimate {
    double center = 0.0;
    double width = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    int retained_count = 0;
    bool uses_quantile_range = false;
    bool rejected_zero_mad_clip = false;
};

// ==========================================
// Function: Estimate a quantization-safe first-layer integer star-area pilot
// Method: Use configurable symmetric Q(q)--Q(1-q) bounds for an initially zero
//         MAD; otherwise accept clips only while their next MAD remains positive.
// ==========================================
bool estimatePSFCountPilot(
    const std::vector<double>& sorted_input,
    double clip_sigma,
    int max_iterations,
    double zero_mad_quantile,
    double range_sigma,
    PSFCountPilotEstimate& pilot) {
    pilot = {};
    if (sorted_input.empty() || !std::isfinite(clip_sigma)
        || clip_sigma <= 0.0 || max_iterations <= 0
        || !std::isfinite(zero_mad_quantile)
        || zero_mad_quantile < 0.0 || zero_mad_quantile >= 0.5
        || !std::isfinite(range_sigma) || range_sigma <= 0.0) {
        return false;
    }

    std::vector<double> population = sorted_input;
    double center = 0.0;
    double width = 0.0;
    if (!medianAndMad(population, center, width)) return false;

    if (width <= 0.0) {
        pilot.center = center;
        pilot.width = 0.0;
        pilot.lower = sortedQuantile(population, zero_mad_quantile);
        pilot.upper = sortedQuantile(population, 1.0 - zero_mad_quantile);
        pilot.retained_count = static_cast<int>(population.size());
        pilot.uses_quantile_range = true;
        return std::isfinite(pilot.center)
            && std::isfinite(pilot.width)
            && pilot.width >= 0.0
            && std::isfinite(pilot.lower)
            && std::isfinite(pilot.upper)
            && pilot.upper >= pilot.lower
            && pilot.retained_count > 0;
    }

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (!medianAndMad(population, center, width)) return false;
        std::vector<double> clipped;
        clipped.reserve(population.size());
        for (double value : population) {
            if (std::abs(value - center) <= clip_sigma * width) {
                clipped.push_back(value);
            }
        }
        if (clipped.empty()) return false;
        if (clipped.size() == population.size()) break;

        double next_center = 0.0;
        double next_width = 0.0;
        if (!medianAndMad(clipped, next_center, next_width)) return false;
        if (next_width <= 0.0) {
            pilot.rejected_zero_mad_clip = true;
            break;
        }
        population.swap(clipped);
    }

    if (!medianAndMad(population, pilot.center, pilot.width)) return false;
    pilot.retained_count = static_cast<int>(population.size());
    if (pilot.width <= 0.0) {
        pilot.uses_quantile_range = true;
        pilot.lower = sortedQuantile(population, zero_mad_quantile);
        pilot.upper = sortedQuantile(population, 1.0 - zero_mad_quantile);
    } else {
        pilot.lower = pilot.center - range_sigma * pilot.width;
        pilot.upper = pilot.center + range_sigma * pilot.width;
    }
    return std::isfinite(pilot.center)
        && std::isfinite(pilot.width)
        && pilot.width >= 0.0
        && std::isfinite(pilot.lower)
        && std::isfinite(pilot.upper)
        && pilot.upper >= pilot.lower
        && pilot.retained_count > 0;
}

// ==========================================
// Function: Test whether a top-K list contains one candidate index
// Method: Perform a bounded linear scan over the short neighbour vector.
// ==========================================
bool containsNeighbour(
    const std::vector<NeighborEdge>& neighbours,
    int star_index) {
    return std::any_of(
        neighbours.begin(), neighbours.end(),
        [star_index](const NeighborEdge& edge) {
            return edge.star_index == star_index;
        });
}

// ==========================================
// Class: Maintain connected components over a compact active-index domain
// Method: Apply path compression and union by rank.
// ==========================================
class DisjointSet {
public:
    // ==========================================
    // Function: Initialize singleton disjoint-set components
    // Method: Assign every compact active index as its own parent with zero rank.
    // ==========================================
    explicit DisjointSet(int size)
        : parent_(static_cast<std::size_t>(size)),
          rank_(static_cast<std::size_t>(size), 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    // ==========================================
    // Function: Find one disjoint-set root
    // Method: Compress every traversed parent link recursively.
    // ==========================================
    int find(int value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    // ==========================================
    // Function: Unite two disjoint-set members
    // Method: Attach the lower-rank root and increment rank on a tie.
    // ==========================================
    void unite(int first, int second) {
        int root_first = find(first);
        int root_second = find(second);
        if (root_first == root_second) return;
        if (rank_[root_first] < rank_[root_second]) {
            std::swap(root_first, root_second);
        }
        parent_[root_second] = root_first;
        if (rank_[root_first] == rank_[root_second]) rank_[root_first]++;
    }

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

// ==========================================
// Function: Smooth one continuous FD histogram without filling empty bins
// Method: Apply the shared 1-2-3-2-1 kernel with edge renormalization.
// ==========================================
std::vector<double> smoothUpperElbowHistogram(
    const std::vector<double>& histogram) {
    const int offsets[] = {-2, -1, 0, 1, 2};
    const double weights[] = {1.0, 2.0, 3.0, 2.0, 1.0};
    std::vector<double> smoothed(histogram.size(), 0.0);
    for (int bin = 0; bin < static_cast<int>(histogram.size()); ++bin) {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (int index = 0; index < 5; ++index) {
            const int neighbour = bin + offsets[index];
            if (neighbour < 0
                || neighbour >= static_cast<int>(histogram.size())) {
                continue;
            }
            weighted_sum += weights[index] * histogram[neighbour];
            weight_sum += weights[index];
        }
        smoothed[bin] = weight_sum > 0.0
            ? weighted_sum / weight_sum
            : 0.0;
    }
    return smoothed;
}

// ==========================================
// Function: Collapse positive local-maximum plateaus to deterministic peaks
// Method: Treat each exact-height run atomically and use its lower middle bin.
// ==========================================
std::vector<int> findUpperElbowPeaks(
    const std::vector<double>& smoothed_histogram) {
    std::vector<int> peaks;
    int first = 0;
    while (first < static_cast<int>(smoothed_histogram.size())) {
        int last = first;
        while (last + 1 < static_cast<int>(smoothed_histogram.size())
               && smoothed_histogram[last + 1]
                    == smoothed_histogram[first]) {
            ++last;
        }
        const double height = smoothed_histogram[first];
        const double left = first == 0
            ? -std::numeric_limits<double>::infinity()
            : smoothed_histogram[first - 1];
        const double right = last + 1
                == static_cast<int>(smoothed_histogram.size())
            ? -std::numeric_limits<double>::infinity()
            : smoothed_histogram[last + 1];
        if (std::isfinite(height) && height > 0.0
            && height >= left && height >= right) {
            peaks.push_back(first + (last - first) / 2);
        }
        first = last + 1;
    }
    return peaks;
}

// ==========================================
// Function: Classify peaks and select one gated right-side curvature elbow
// Method: Search below the configured main-peak fraction between the rightmost
//         valid peak and first invalid peak; ties keep the nearer lower bin.
// ==========================================
bool analyzeUpperElbowTopology(
    double valid_peak_fraction,
    double elbow_search_height_fraction,
    PSFUpperElbowHistogramResult& result) {
    result.valid = false;
    result.main_peak_bin = -1;
    result.rightmost_valid_peak_bin = -1;
    result.first_invalid_peak_bin = -1;
    result.elbow_bin = -1;
    result.cut = 0.0;
    result.elbow_search_height = 0.0;
    result.elbow_search_first_bin = -1;
    result.elbow_search_last_bin = -1;
    result.elbow_candidate_count = 0;
    result.peaks.clear();
    result.valid_peaks.clear();
    result.invalid_peaks.clear();
    if (!std::isfinite(elbow_search_height_fraction)
        || elbow_search_height_fraction <= 0.0
        || elbow_search_height_fraction >= 1.0) {
        result.status = PSFUpperElbowStatus::InvalidConfig;
        return false;
    }
    result.smoothed_histogram =
        smoothUpperElbowHistogram(result.histogram);
    result.peaks = findUpperElbowPeaks(result.smoothed_histogram);
    if (result.peaks.empty()) {
        result.status = PSFUpperElbowStatus::NoPeaks;
        return false;
    }

    result.main_peak_bin = result.peaks.front();
    for (int peak : result.peaks) {
        if (result.smoothed_histogram[peak]
            > result.smoothed_histogram[result.main_peak_bin]) {
            result.main_peak_bin = peak;
        }
    }
    const double valid_height =
        result.smoothed_histogram[result.main_peak_bin]
        * valid_peak_fraction;
    for (int peak : result.peaks) {
        if (result.smoothed_histogram[peak] > valid_height) {
            result.valid_peaks.push_back(peak);
        } else {
            result.invalid_peaks.push_back(peak);
        }
    }
    if (result.valid_peaks.empty()) {
        result.status = PSFUpperElbowStatus::NoPeaks;
        return false;
    }
    result.rightmost_valid_peak_bin = result.valid_peaks.back();
    for (int peak : result.invalid_peaks) {
        if (peak > result.rightmost_valid_peak_bin) {
            result.first_invalid_peak_bin = peak;
            break;
        }
    }

    const int first_candidate = std::max(
        1, result.rightmost_valid_peak_bin + 1);
    const int last_candidate = result.first_invalid_peak_bin >= 0
        ? result.first_invalid_peak_bin - 1
        : static_cast<int>(result.smoothed_histogram.size()) - 2;
    const double main_peak_height =
        result.smoothed_histogram[result.main_peak_bin];
    result.elbow_search_height =
        main_peak_height * elbow_search_height_fraction;
    double best_curvature = 0.0;
    for (int bin = first_candidate; bin <= last_candidate; ++bin) {
        const double height = result.smoothed_histogram[bin];
        if (!std::isfinite(height)
            || !(height < result.elbow_search_height)) {
            continue;
        }
        if (result.elbow_search_first_bin < 0) {
            result.elbow_search_first_bin = bin;
        }
        result.elbow_search_last_bin = bin;
        ++result.elbow_candidate_count;
        const double curvature = result.smoothed_histogram[bin - 1]
            - 2.0 * height
            + result.smoothed_histogram[bin + 1];
        if (std::isfinite(curvature) && curvature > best_curvature) {
            best_curvature = curvature;
            result.elbow_bin = bin;
        }
    }
    if (result.elbow_bin < 0) {
        result.status = PSFUpperElbowStatus::NoElbow;
        return false;
    }
    result.cut = result.bin_origin
        + (static_cast<double>(result.elbow_bin) + 0.5)
            * result.bin_width;
    if (!std::isfinite(result.cut)) {
        result.status = PSFUpperElbowStatus::InvalidInput;
        return false;
    }
    result.valid = true;
    result.status = PSFUpperElbowStatus::Valid;
    return true;
}

// ==========================================
// Function: Estimate one generic adaptive upper-elbow histogram
// Method: Filter samples, derive an FD grid, collapse peaks, then maximize
//         gated positive curvature right of the rightmost valid peak.
// ==========================================
template <typename Sample>
bool estimatePSFUpperElbowCutImpl(
    const std::vector<Sample>& input,
    const PSFUpperElbowHistogramConfig& config,
    PSFUpperElbowHistogramResult& result) {
    result = {};
    if (!std::isfinite(config.valid_peak_fraction)
        || config.valid_peak_fraction <= 0.0
        || config.valid_peak_fraction >= 1.0
        || !std::isfinite(config.elbow_search_height_fraction)
        || config.elbow_search_height_fraction <= 0.0
        || config.elbow_search_height_fraction >= 1.0) {
        result.status = PSFUpperElbowStatus::InvalidConfig;
        return false;
    }

    try {
        std::vector<double> values;
        values.reserve(input.size());
        for (Sample sample : input) {
            const double value = static_cast<double>(sample);
            if (std::isfinite(value)) values.push_back(value);
        }
        result.finite_value_count = values.size();
        if (values.empty()) return false;
        if (config.force_zero_origin
            && *std::min_element(values.begin(), values.end()) < 0.0) {
            result.status = PSFUpperElbowStatus::InvalidInput;
            return false;
        }

        std::vector<double> fd_values;
        fd_values.reserve(values.size());
        for (double value : values) {
            if (!config.exclude_zero_from_fd || value > 0.0) {
                fd_values.push_back(value);
            }
        }
        result.fd_sample_count = fd_values.size();
        result.fd_scale_sample_count = config.fd_scale_sample_count > 0U
            ? config.fd_scale_sample_count
            : result.fd_sample_count;
        if (fd_values.empty() || result.fd_scale_sample_count == 0U) {
            result.status = PSFUpperElbowStatus::NoFDSamples;
            return false;
        }
        std::sort(fd_values.begin(), fd_values.end());
        const double first_quartile = sortedQuantile(fd_values, 0.25);
        const double third_quartile = sortedQuantile(fd_values, 0.75);
        result.fd_iqr = third_quartile - first_quartile;
        if (std::isfinite(result.fd_iqr) && result.fd_iqr > 0.0) {
            result.bin_width = 2.0 * result.fd_iqr
                / std::cbrt(static_cast<double>(
                    result.fd_scale_sample_count));
        } else if (config.zero_iqr_uses_min_positive
                   && fd_values.front() > 0.0) {
            result.bin_width = fd_values.front();
        }
        if (!std::isfinite(result.bin_width) || result.bin_width <= 0.0) {
            result.status = PSFUpperElbowStatus::NonPositiveWidth;
            return false;
        }

        const auto minimum_and_maximum =
            std::minmax_element(values.begin(), values.end());
        const double minimum = *minimum_and_maximum.first;
        const double maximum = *minimum_and_maximum.second;
        if (config.force_zero_origin) {
            result.bin_origin = 0.0;
        } else {
            const long double origin = std::floor(
                static_cast<long double>(minimum)
                    / static_cast<long double>(result.bin_width))
                * static_cast<long double>(result.bin_width);
            result.bin_origin = static_cast<double>(origin);
        }
        if (!std::isfinite(result.bin_origin)
            || result.bin_origin > minimum) {
            result.status = PSFUpperElbowStatus::InvalidInput;
            return false;
        }

        const long double scaled_range =
            (static_cast<long double>(maximum)
             - static_cast<long double>(result.bin_origin))
            / static_cast<long double>(result.bin_width);
        if (!std::isfinite(scaled_range) || scaled_range < 0.0L) {
            result.status = PSFUpperElbowStatus::UnsafeBinCount;
            return false;
        }
        const long double last_bin_value = std::floor(scaled_range);
        const long double maximum_last_bin = static_cast<long double>(
            std::numeric_limits<int>::max() - 1);
        if (last_bin_value > maximum_last_bin) {
            result.status = PSFUpperElbowStatus::UnsafeBinCount;
            return false;
        }
        const std::size_t bin_count =
            static_cast<std::size_t>(last_bin_value) + 1U;
        if (bin_count == 0
            || bin_count > result.histogram.max_size()) {
            result.status = PSFUpperElbowStatus::UnsafeBinCount;
            return false;
        }
        result.histogram.assign(bin_count, 0.0);
        for (double value : values) {
            long double scaled =
                (static_cast<long double>(value)
                 - static_cast<long double>(result.bin_origin))
                / static_cast<long double>(result.bin_width);
            const long double tolerance = 64.0L
                * std::numeric_limits<long double>::epsilon()
                * std::max(1.0L, std::abs(scaled));
            if (!std::isfinite(scaled) || scaled < -tolerance) {
                result.status = PSFUpperElbowStatus::InvalidInput;
                return false;
            }
            if (scaled < 0.0L) scaled = 0.0L;
            std::size_t bin = static_cast<std::size_t>(std::floor(scaled));
            if (bin >= bin_count) bin = bin_count - 1U;
            result.histogram[bin] += 1.0;
        }
        return analyzeUpperElbowTopology(
            config.valid_peak_fraction,
            config.elbow_search_height_fraction,
            result);
    } catch (const std::bad_alloc&) {
        result = {};
        result.status = PSFUpperElbowStatus::AllocationFailure;
        return false;
    } catch (const std::length_error&) {
        result = {};
        result.status = PSFUpperElbowStatus::AllocationFailure;
        return false;
    }
}

}  // namespace

// ==========================================
// Function: Analyze one already-binned upper-elbow histogram
// Method: Validate nonnegative finite bins, then expose the production topology.
// ==========================================
bool analyzePSFUpperElbowHistogram(
    const std::vector<double>& histogram,
    double bin_origin,
    double bin_width,
    double valid_peak_fraction,
    double elbow_search_height_fraction,
    PSFUpperElbowHistogramResult& result) {
    result = {};
    if (histogram.empty() || !std::isfinite(bin_origin)
        || !std::isfinite(bin_width) || bin_width <= 0.0
        || !std::isfinite(valid_peak_fraction)
        || valid_peak_fraction <= 0.0 || valid_peak_fraction >= 1.0
        || std::any_of(
            histogram.begin(), histogram.end(),
            [](double count) {
                return !std::isfinite(count) || count < 0.0;
            })) {
        result.status = PSFUpperElbowStatus::InvalidInput;
        return false;
    }
    try {
        result.bin_origin = bin_origin;
        result.bin_width = bin_width;
        result.histogram = histogram;
        return analyzeUpperElbowTopology(
            valid_peak_fraction,
            elbow_search_height_fraction,
            result);
    } catch (const std::bad_alloc&) {
        result = {};
        result.status = PSFUpperElbowStatus::AllocationFailure;
        return false;
    } catch (const std::length_error&) {
        result = {};
        result.status = PSFUpperElbowStatus::AllocationFailure;
        return false;
    }
}

// ==========================================
// Function: Estimate an FD-histogram upper elbow from double samples
// Method: Delegate to the type-generic finite-filtering implementation.
// ==========================================
bool estimatePSFUpperElbowCut(
    const std::vector<double>& values,
    const PSFUpperElbowHistogramConfig& config,
    PSFUpperElbowHistogramResult& result) {
    return estimatePSFUpperElbowCutImpl(values, config, result);
}

// ==========================================
// Function: Estimate an FD-histogram upper elbow from float samples
// Method: Delegate without first expanding the compact production pair vector.
// ==========================================
bool estimatePSFUpperElbowCut(
    const std::vector<float>& values,
    const PSFUpperElbowHistogramConfig& config,
    PSFUpperElbowHistogramResult& result) {
    return estimatePSFUpperElbowCutImpl(values, config, result);
}

// ==========================================
// Function: Return a stable adaptive-histogram status label
// Method: Map every outcome to one uppercase diagnostic token.
// ==========================================
const char* psfUpperElbowStatusName(PSFUpperElbowStatus status) {
    switch (status) {
        case PSFUpperElbowStatus::NoFiniteValues: return "NO_FINITE_VALUES";
        case PSFUpperElbowStatus::InvalidConfig: return "INVALID_CONFIG";
        case PSFUpperElbowStatus::InvalidInput: return "INVALID_INPUT";
        case PSFUpperElbowStatus::NoFDSamples: return "NO_FD_SAMPLES";
        case PSFUpperElbowStatus::NonPositiveWidth: return "NONPOSITIVE_WIDTH";
        case PSFUpperElbowStatus::UnsafeBinCount: return "UNSAFE_BIN_COUNT";
        case PSFUpperElbowStatus::AllocationFailure: return "ALLOCATION_FAILURE";
        case PSFUpperElbowStatus::NoPeaks: return "NO_PEAKS";
        case PSFUpperElbowStatus::NoElbow: return "NO_ELBOW";
        case PSFUpperElbowStatus::Valid: return "VALID";
    }
    return "UNKNOWN";
}

// ==========================================
// Function: Classify one finite pair against the Type-3 upper cut
// Method: Preserve the scientific strict-greater-than boundary exactly.
// ==========================================
bool isPSFType3BadPair(double chi, double pair_chi_cut) {
    return std::isfinite(chi) && std::isfinite(pair_chi_cut)
        && chi > pair_chi_cut;
}

// ==========================================
// Function: Apply the Type-3 bad-pair-fraction gate to one chip
// Method: Require finite denominators, apply a strict upper rejection only when
//         requested, and clear the complete result below the chip minimum.
// ==========================================
PSFType3ChipSelection selectPSFType3FractionSurvivors(
    const std::vector<double>& bad_pair_fractions,
    const std::vector<bool>& has_finite_pair_denominator,
    bool apply_fraction_cut,
    double fraction_cut,
    int minimum_retained) {
    PSFType3ChipSelection result;
    result.selected.assign(bad_pair_fractions.size(), false);
    if (bad_pair_fractions.size() != has_finite_pair_denominator.size()
        || minimum_retained < 0
        || (apply_fraction_cut
            && (!std::isfinite(fraction_cut) || fraction_cut < 0.0))) {
        return result;
    }
    for (std::size_t index = 0; index < bad_pair_fractions.size(); ++index) {
        const double fraction = bad_pair_fractions[index];
        if (!has_finite_pair_denominator[index]
            || !std::isfinite(fraction) || fraction < 0.0
            || fraction > 1.0) {
            continue;
        }
        result.finite_pair_count++;
        if (!apply_fraction_cut || fraction <= fraction_cut) {
            result.selected[index] = true;
            result.retained_count++;
        }
    }
    if (result.retained_count < static_cast<std::size_t>(minimum_retained)) {
        std::fill(result.selected.begin(), result.selected.end(), false);
        result.retained_count = 0;
        result.rejected_by_minimum = true;
    }
    return result;
}

// ==========================================
// Function: Derive the shared PSF chi-window bounds
// Method: Preserve the legacy inclusive n/4-1 through 3n/4-1 definition.
// ==========================================
PSFChiWindow getPSFChiWindow(int n) {
    if (n <= 0) return {};
    PSFChiWindow window;
    window.first = std::max(0, n / 4 - 1);
    window.last = std::min(n - 1, (n / 4) * 3 - 1);
    return window;
}

// ==========================================
// Function: Return the nominal center of one fixed two-count histogram bin
// Method: Offset the first allowed count by two per bin and one half count.
// ==========================================
double psfCountHistogramBinCenter(
    int histogram_first_count,
    int bin) {
    return static_cast<double>(histogram_first_count)
        + static_cast<double>(PSFCountHistogramBinWidth * bin) + 0.5;
}

// ==========================================
// Function: Select one integer star-area histogram peak deterministically
// Method: Use density without Gaia; with Gaia, anchor eligibility to the global
//         nominal-center distance plus one count and rank deterministically.
// ==========================================
int selectPSFCountPeak(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    const std::vector<double>& gaia_histogram,
    int histogram_first_count,
    double pilot_center,
    bool pilot_uses_gaia) {
    if (peaks.empty() || smoothed_histogram.empty()
        || gaia_histogram.size() != smoothed_histogram.size()
        || histogram_first_count <= 0
        || !std::isfinite(pilot_center)) {
        return -1;
    }
    for (int peak : peaks) {
        if (peak < 0
            || peak >= static_cast<int>(smoothed_histogram.size())) {
            return -1;
        }
    }

    int selected = peaks.front();
    if (!pilot_uses_gaia) {
        for (int candidate : peaks) {
            if (smoothed_histogram[candidate]
                > smoothed_histogram[selected]) {
                selected = candidate;
            }
        }
        return selected;
    }

    double minimum_distance = std::numeric_limits<double>::infinity();
    for (int peak : peaks) {
        const double distance = std::abs(
            psfCountHistogramBinCenter(histogram_first_count, peak)
                - pilot_center);
        minimum_distance = std::min(minimum_distance, distance);
    }

    selected = -1;
    double selected_distance = std::numeric_limits<double>::infinity();
    for (int candidate : peaks) {
        const double candidate_distance = std::abs(
            psfCountHistogramBinCenter(histogram_first_count, candidate)
                - pilot_center);
        if (candidate_distance > minimum_distance + 1.0) {
            continue;
        }
        if (selected < 0) {
            selected = candidate;
            selected_distance = candidate_distance;
            continue;
        }

        bool replace = false;
        if (gaia_histogram[candidate] != gaia_histogram[selected]) {
            replace = gaia_histogram[candidate]
                > gaia_histogram[selected];
        } else if (smoothed_histogram[candidate]
                   != smoothed_histogram[selected]) {
            replace = smoothed_histogram[candidate]
                > smoothed_histogram[selected];
        } else if (candidate_distance != selected_distance) {
            replace = candidate_distance < selected_distance;
        } else {
            replace = candidate < selected;
        }
        if (replace) {
            selected = candidate;
            selected_distance = candidate_distance;
        }
    }
    return selected;
}

// ==========================================
// Function: Find the selected peak's complete significant peak complex
// Method: Include local peaks strictly above H_selected/e and descend outward
//         from their extrema until the smoothed histogram rises again.
// ==========================================
PSFCountBinRange findPSFCountPeakComplexBasin(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    int selected_peak) {
    PSFCountBinRange basin;
    if (selected_peak < 0
        || selected_peak >= static_cast<int>(smoothed_histogram.size())
        || !std::isfinite(smoothed_histogram[selected_peak])
        || smoothed_histogram[selected_peak] <= 0.0) {
        return basin;
    }

    const double minimum_height =
        smoothed_histogram[selected_peak] * std::exp(-1.0);
    for (int peak : peaks) {
        if (peak < 0 || peak >= static_cast<int>(smoothed_histogram.size())
            || !std::isfinite(smoothed_histogram[peak])) {
            return {};
        }
        if (smoothed_histogram[peak] > minimum_height) {
            if (basin.first < 0) basin.first = peak;
            basin.first = std::min(basin.first, peak);
            basin.last = std::max(basin.last, peak);
        }
    }
    if (basin.first < 0 || basin.last < 0
        || selected_peak < basin.first || selected_peak > basin.last) {
        return {};
    }
    while (basin.first > 0
           && smoothed_histogram[basin.first - 1]
                <= smoothed_histogram[basin.first]) {
        --basin.first;
    }
    while (basin.last + 1 < static_cast<int>(smoothed_histogram.size())
           && smoothed_histogram[basin.last + 1]
                <= smoothed_histogram[basin.last]) {
        ++basin.last;
    }
    return basin;
}

// ==========================================
// Function: Find independent outer elbows around the selected count peak
// Method: Cross below ten percent of peak height, then maximize positive signed
//         curvature toward each domain edge with nearest-bin tie breaking.
// ==========================================
PSFCountElbows findPSFCountOuterElbows(
    const std::vector<double>& smoothed_histogram,
    int selected_peak) {
    PSFCountElbows elbows;
    if (selected_peak < 0
        || selected_peak >= static_cast<int>(smoothed_histogram.size())
        || !std::isfinite(smoothed_histogram[selected_peak])
        || smoothed_histogram[selected_peak] <= 0.0) {
        return elbows;
    }
    const double threshold = 0.10 * smoothed_histogram[selected_peak];

    int left_crossing = -1;
    for (int bin = selected_peak - 1; bin >= 0; --bin) {
        if (smoothed_histogram[bin] < threshold) {
            left_crossing = bin;
            break;
        }
    }
    double best_curvature = 0.0;
    for (int bin = left_crossing; bin >= 1; --bin) {
        const double curvature = smoothed_histogram[bin - 1]
            - 2.0 * smoothed_histogram[bin]
            + smoothed_histogram[bin + 1];
        if (std::isfinite(curvature) && curvature > best_curvature) {
            best_curvature = curvature;
            elbows.left = bin;
        }
    }

    int right_crossing = -1;
    for (int bin = selected_peak + 1;
         bin < static_cast<int>(smoothed_histogram.size()); ++bin) {
        if (smoothed_histogram[bin] < threshold) {
            right_crossing = bin;
            break;
        }
    }
    best_curvature = 0.0;
    for (int bin = right_crossing;
         bin >= 0 && bin + 1 < static_cast<int>(smoothed_histogram.size());
         ++bin) {
        const double curvature = smoothed_histogram[bin - 1]
            - 2.0 * smoothed_histogram[bin]
            + smoothed_histogram[bin + 1];
        if (std::isfinite(curvature) && curvature > best_curvature) {
            best_curvature = curvature;
            elbows.right = bin;
        }
    }
    return elbows;
}

// ==========================================
// Function: Fill bounded one- or two-bin holes in a working histogram
// Method: Detect every zero run from the immutable raw input, then interpolate
//         only short internal runs between positive raw endpoints.
// ==========================================
std::vector<double> interpolateShortInternalHoles(
    const std::vector<double>& histogram) {
    std::vector<double> working = histogram;
    std::size_t index = 0;
    while (index < histogram.size()) {
        if (histogram[index] != 0.0) {
            ++index;
            continue;
        }
        const std::size_t first = index;
        while (index < histogram.size() && histogram[index] == 0.0) ++index;
        const std::size_t length = index - first;
        if (first == 0 || index >= histogram.size() || length > 2
            || histogram[first - 1] <= 0.0 || histogram[index] <= 0.0) {
            continue;
        }
        const double left = histogram[first - 1];
        const double right = histogram[index];
        for (std::size_t offset = 0; offset < length; ++offset) {
            const double fraction = static_cast<double>(offset + 1)
                / static_cast<double>(length + 1);
            working[first + offset] = left + fraction * (right - left);
        }
    }
    return working;
}

// ==========================================
// Function: Refine a peak-basin seed with re-absorbing asymmetric MAD
// Method: Rebuild each pass from all real pilot-domain samples inside the current
//         asymmetric bounds, allowing previously excluded samples to return.
// ==========================================
PSFCountRefinement refinePSFCountPopulation(
    const std::vector<double>& seed_values,
    const std::vector<double>& domain_values,
    double locus_sigma,
    int iterations) {
    PSFCountRefinement result;
    if (seed_values.empty() || domain_values.empty()
        || !std::isfinite(locus_sigma) || locus_sigma <= 0.0
        || iterations <= 0) {
        return result;
    }
    std::vector<double> population = seed_values;
    std::vector<double> domain = domain_values;
    if (std::any_of(population.begin(), population.end(),
                    [](double value) { return !std::isfinite(value); })
        || std::any_of(domain.begin(), domain.end(),
                       [](double value) { return !std::isfinite(value); })) {
        return result;
    }
    std::sort(population.begin(), population.end());
    std::sort(domain.begin(), domain.end());

    double center = 0.0;
    double lower_width = 0.0;
    double upper_width = 0.0;
    for (int iteration = 0; iteration < iterations; ++iteration) {
        if (!medianAndAsymmetricMad(
                population, center, lower_width, upper_width)) {
            return {};
        }
        const double clip_lower = center - locus_sigma * lower_width;
        const double clip_upper = center + locus_sigma * upper_width;
        std::vector<double> refined;
        refined.reserve(domain.size());
        for (double value : domain) {
            if (value >= clip_lower && value <= clip_upper) {
                refined.push_back(value);
            }
        }
        if (refined.empty()) return {};
        population.swap(refined);
    }
    if (!medianAndAsymmetricMad(
            population, center, lower_width, upper_width)) {
        return {};
    }
    result.valid = true;
    result.center = center;
    result.lower_width = lower_width;
    result.upper_width = upper_width;
    result.sample_count = static_cast<int>(population.size());
    return result;
}

// ==========================================
// Function: Count exp(-1)-threshold pixels in one square Fourier-power stamp
// Method: Compare every stamp value with the finite central value times exp(-1).
// ==========================================
int countPSFStarArea(
    const std::vector<float>& power,
    int stamp_side) {
    if (stamp_side <= 0
        || power.size() < static_cast<std::size_t>(stamp_side * stamp_side)) {
        return 0;
    }
    const float center = power[static_cast<std::size_t>(stamp_side / 2)
        * static_cast<std::size_t>(stamp_side)
        + static_cast<std::size_t>(stamp_side / 2)];
    if (!std::isfinite(center)) return 0;
    const float threshold = center * std::exp(-1.0f);
    int star_area = 0;
    for (int index = 0; index < stamp_side * stamp_side; ++index) {
        if (std::isfinite(power[static_cast<std::size_t>(index)])
            && power[static_cast<std::size_t>(index)] >= threshold) {
            ++star_area;
        }
    }
    return star_area;
}

// ==========================================
// Function: Convert an exp(-1) star area to the historical PSF FWHM
// Method: Preserve the legacy area-minus-1e-5 formula and explicit pixel scale.
// ==========================================
double fwhmFromStarArea(
    double star_area,
    int stamp_side,
    double pixel_size) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double area = star_area - 1.0e-5;
    if (!std::isfinite(star_area) || area <= 0.0 || stamp_side <= 0
        || !std::isfinite(pixel_size) || pixel_size <= 0.0) {
        return 0.0;
    }
    const double beta = static_cast<double>(stamp_side) / (2.0 * pi)
        / std::sqrt(area / pi);
    return beta * 2.0 * std::sqrt(2.0 * std::log(2.0)) * pixel_size;
}

// ==========================================
// Function: Estimate an exposure-wide integer star-area locus
// Method: Build fixed two-count bins, refine the significant peak complex from
//         all pilot-domain samples, and widen only with independent outer elbows.
// ==========================================
bool estimatePSFCountLocus(
    const std::vector<PSFCountSample>& samples,
    const PSFCountLocusConfig& config,
    PSFCountLocus& locus,
    PSFCountLocusDiagnostics* diagnostics) {
    locus = {};
    if (diagnostics != nullptr) *diagnostics = {};
    if (!std::isfinite(config.pilot_clip_sigma)
        || config.pilot_clip_sigma <= 0.0
        || config.pilot_clip_iterations <= 0
        || !std::isfinite(config.zero_mad_quantile)
        || config.zero_mad_quantile < 0.0
        || config.zero_mad_quantile >= 0.5
        || !std::isfinite(config.histogram_range_sigma)
        || config.histogram_range_sigma <= 0.0
        || !std::isfinite(config.locus_sigma)
        || config.locus_sigma <= 0.0
        || config.minimum_samples <= 0
        || config.minimum_gaia_matches <= 0) {
        return false;
    }

    std::vector<double> values;
    std::vector<double> gaia_values;
    values.reserve(samples.size());
    gaia_values.reserve(samples.size());
    for (const PSFCountSample& sample : samples) {
        if (sample.star_area <= 0) continue;
        values.push_back(static_cast<double>(sample.star_area));
        if (sample.gaia_matched) {
            gaia_values.push_back(static_cast<double>(sample.star_area));
        }
    }
    if (diagnostics != nullptr) {
        diagnostics->sample_count = static_cast<int>(values.size());
        diagnostics->gaia_match_count = static_cast<int>(gaia_values.size());
    }
    if (static_cast<int>(values.size()) < config.minimum_samples) return false;

    std::sort(values.begin(), values.end());
    std::sort(gaia_values.begin(), gaia_values.end());
    const bool has_gaia_support =
        static_cast<int>(gaia_values.size()) >= config.minimum_gaia_matches;
    if (diagnostics != nullptr && has_gaia_support) {
        diagnostics->has_gaia_median = true;
        diagnostics->gaia_median = sortedMedian(gaia_values);
    }

    bool pilot_uses_gaia = has_gaia_support;
    int pilot_input_count = pilot_uses_gaia
        ? static_cast<int>(gaia_values.size())
        : static_cast<int>(values.size());
    PSFCountPilotEstimate pilot;
    bool pilot_valid = estimatePSFCountPilot(
        pilot_uses_gaia ? gaia_values : values,
        config.pilot_clip_sigma,
        config.pilot_clip_iterations,
        config.zero_mad_quantile,
        config.histogram_range_sigma,
        pilot);
    if (pilot_uses_gaia
        && (!pilot_valid
            || pilot.retained_count < config.minimum_gaia_matches)) {
        pilot_uses_gaia = false;
        pilot_input_count = static_cast<int>(values.size());
        pilot = {};
        pilot_valid = estimatePSFCountPilot(
            values,
            config.pilot_clip_sigma,
            config.pilot_clip_iterations,
            config.zero_mad_quantile,
            config.histogram_range_sigma,
            pilot);
    }
    if (!pilot_valid) return false;

    if (diagnostics != nullptr) {
        diagnostics->pilot_uses_gaia = pilot_uses_gaia;
        diagnostics->pilot_input_count = pilot_input_count;
        diagnostics->pilot_retained_count = pilot.retained_count;
        diagnostics->pilot_center = pilot.center;
        diagnostics->pilot_width = pilot.width;
        diagnostics->pilot_lower = pilot.lower;
        diagnostics->pilot_upper = pilot.upper;
        diagnostics->pilot_uses_quantile_range = pilot.uses_quantile_range;
        diagnostics->pilot_rejected_zero_mad_clip =
            pilot.rejected_zero_mad_clip;
    }

    int histogram_first_count = std::max(
        1, static_cast<int>(std::ceil(pilot.lower)));
    int histogram_last_count = static_cast<int>(std::floor(pilot.upper));
    if (histogram_first_count > histogram_last_count) {
        const std::vector<double>& source = pilot_uses_gaia
            ? gaia_values : values;
        if (source.empty()) return false;
        int nearest = static_cast<int>(std::llround(source.front()));
        double nearest_distance = std::abs(
            static_cast<double>(nearest) - pilot.center);
        for (double value : source) {
            const int observed = static_cast<int>(std::llround(value));
            const double distance = std::abs(value - pilot.center);
            if (distance < nearest_distance
                || (distance == nearest_distance && observed < nearest)) {
                nearest = observed;
                nearest_distance = distance;
            }
        }
        histogram_first_count = nearest;
        histogram_last_count = nearest;
    }
    const int histogram_bin_count =
        (histogram_last_count - histogram_first_count)
            / PSFCountHistogramBinWidth + 1;
    if (histogram_bin_count <= 0) return false;

    std::vector<double> histogram(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    std::vector<double> domain_values;
    domain_values.reserve(values.size());
    int histogram_sample_count = 0;
    int histogram_below_count = 0;
    int histogram_above_count = 0;
    for (double value : values) {
        const int star_area = static_cast<int>(std::llround(value));
        if (star_area < histogram_first_count) {
            histogram_below_count++;
            continue;
        }
        if (star_area > histogram_last_count) {
            histogram_above_count++;
            continue;
        }
        const int bin = (star_area - histogram_first_count)
            / PSFCountHistogramBinWidth;
        histogram[bin] += 1.0;
        domain_values.push_back(value);
        histogram_sample_count++;
    }

    std::vector<double> gaia_histogram(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    int gaia_histogram_sample_count = 0;
    int gaia_histogram_below_count = 0;
    int gaia_histogram_above_count = 0;
    for (double value : gaia_values) {
        const int star_area = static_cast<int>(std::llround(value));
        if (star_area < histogram_first_count) {
            gaia_histogram_below_count++;
            continue;
        }
        if (star_area > histogram_last_count) {
            gaia_histogram_above_count++;
            continue;
        }
        const int bin = (star_area - histogram_first_count)
            / PSFCountHistogramBinWidth;
        gaia_histogram[bin] += 1.0;
        gaia_histogram_sample_count++;
    }

    const std::vector<double> working =
        interpolateShortInternalHoles(histogram);
    const int offsets[] = {-2, -1, 0, 1, 2};
    const double weights[] = {1.0, 2.0, 3.0, 2.0, 1.0};
    std::vector<double> smoothed(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    for (int bin = 0; bin < histogram_bin_count; ++bin) {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (int offset_index = 0; offset_index < 5; ++offset_index) {
            const int neighbour = bin + offsets[offset_index];
            if (neighbour < 0 || neighbour >= histogram_bin_count) continue;
            weighted_sum += weights[offset_index] * working[neighbour];
            weight_sum += weights[offset_index];
        }
        smoothed[bin] = weighted_sum / weight_sum;
    }

    std::vector<int> peaks;
    for (int bin = 0; bin < histogram_bin_count; ++bin) {
        const double left = bin == 0 ? -1.0 : smoothed[bin - 1];
        const double right = bin + 1 == histogram_bin_count
            ? -1.0
            : smoothed[bin + 1];
        if (smoothed[bin] > 0.0 && smoothed[bin] >= left && smoothed[bin] >= right) {
            peaks.push_back(bin);
        }
    }
    if (peaks.empty()) return false;

    const int peak_bin = selectPSFCountPeak(
        peaks,
        smoothed,
        gaia_histogram,
        histogram_first_count,
        pilot.center,
        pilot_uses_gaia);
    if (peak_bin < 0) return false;
    if (diagnostics != nullptr) {
        diagnostics->histogram_sample_count = histogram_sample_count;
        diagnostics->histogram_below_count = histogram_below_count;
        diagnostics->histogram_above_count = histogram_above_count;
        diagnostics->gaia_histogram_sample_count =
            gaia_histogram_sample_count;
        diagnostics->gaia_histogram_below_count =
            gaia_histogram_below_count;
        diagnostics->gaia_histogram_above_count =
            gaia_histogram_above_count;
        diagnostics->histogram_first_count = histogram_first_count;
        diagnostics->histogram_last_count = histogram_last_count;
        diagnostics->peak_bin = peak_bin;
        diagnostics->histogram = histogram;
        diagnostics->working_histogram = working;
        diagnostics->smoothed_histogram = smoothed;
        diagnostics->gaia_histogram = gaia_histogram;
    }

    const PSFCountBinRange basin = findPSFCountPeakComplexBasin(
        peaks, smoothed, peak_bin);
    if (basin.first < 0 || basin.last < basin.first) return false;
    const int basin_low = histogram_first_count
        + PSFCountHistogramBinWidth * basin.first;
    const int basin_high = std::min(
        histogram_last_count,
        histogram_first_count
            + PSFCountHistogramBinWidth * (basin.last + 1) - 1);
    std::vector<double> population;
    for (double value : values) {
        if (value >= static_cast<double>(basin_low)
            && value <= static_cast<double>(basin_high)) {
            population.push_back(value);
        }
    }
    if (population.size() < 3) {
        const int fallback_first = std::max(0, peak_bin - 2);
        const int fallback_last = std::min(
            histogram_bin_count - 1, peak_bin + 2);
        const int fallback_low = histogram_first_count
            + PSFCountHistogramBinWidth * fallback_first;
        const int fallback_high = std::min(
            histogram_last_count,
            histogram_first_count
                + PSFCountHistogramBinWidth * (fallback_last + 1) - 1);
        population.clear();
        for (double value : values) {
            if (value >= static_cast<double>(fallback_low)
                && value <= static_cast<double>(fallback_high)) {
                population.push_back(value);
            }
        }
    }
    if (population.empty()) return false;

    const PSFCountRefinement refinement = refinePSFCountPopulation(
        population, domain_values, config.locus_sigma, 2);
    if (!refinement.valid) return false;
    const double mad_lower = refinement.center
        - config.locus_sigma * refinement.lower_width;
    const double mad_upper = refinement.center
        + config.locus_sigma * refinement.upper_width;
    const PSFCountElbows elbows = findPSFCountOuterElbows(
        smoothed, peak_bin);
    const double left_elbow = elbows.left >= 0
        ? psfCountHistogramBinCenter(histogram_first_count, elbows.left)
        : mad_lower;
    const double right_elbow = elbows.right >= 0
        ? psfCountHistogramBinCenter(histogram_first_count, elbows.right)
        : mad_upper;

    locus.valid = true;
    locus.center = refinement.center;
    locus.lower_width = refinement.lower_width;
    locus.upper_width = refinement.upper_width;
    locus.lower = std::min(mad_lower, left_elbow);
    locus.upper = std::max(mad_upper, right_elbow);
    if (diagnostics != nullptr) {
        diagnostics->mad_lower = mad_lower;
        diagnostics->mad_upper = mad_upper;
        diagnostics->left_elbow_bin = elbows.left;
        diagnostics->right_elbow_bin = elbows.right;
        diagnostics->left_elbow_guard_applied =
            elbows.left >= 0 && left_elbow < mad_lower;
        diagnostics->right_elbow_guard_applied =
            elbows.right >= 0 && right_elbow > mad_upper;
    }
    return std::isfinite(locus.lower) && std::isfinite(locus.upper);
}

// ==========================================
// Function: Populate one integer-area histogram on the science count grid
// Method: Reset the output bins and count only in-range integer star areas.
// ==========================================
static void populateCountHistogramOnScienceGrid(
    const std::vector<int>& star_areas,
    int histogram_first_count,
    int histogram_last_count,
    std::size_t histogram_bin_count,
    std::vector<double>& output_histogram) {
    output_histogram.assign(histogram_bin_count, 0.0);
    if (histogram_bin_count == 0) return;

    const int first_count = histogram_first_count;
    const int last_count = histogram_last_count;
    if (last_count < first_count) return;
    for (int star_area : star_areas) {
        if (star_area < first_count || star_area > last_count) continue;
        const int bin = (star_area - first_count)
            / PSFCountHistogramBinWidth;
        if (bin < 0 || bin >= static_cast<int>(histogram_bin_count)) continue;
        output_histogram[static_cast<std::size_t>(bin)] += 1.0;
    }
}

// ==========================================
// Function: Populate the minChi-survivor integer-area histogram
// Method: Count all actual grouping inputs and bin only values on the science
//         grid without changing upstream count-locus diagnostics.
// ==========================================
void populateMinChiSurvivorCountHistogram(
    const std::vector<int>& minchi_star_areas,
    PSFCountLocusDiagnostics& diagnostics) {
    diagnostics.minchi_survivor_count =
        static_cast<int>(minchi_star_areas.size());
    populateCountHistogramOnScienceGrid(
        minchi_star_areas,
        diagnostics.histogram_first_count,
        diagnostics.histogram_last_count,
        diagnostics.histogram.size(),
        diagnostics.minchi_survivor_histogram);
}

// ==========================================
// Function: Populate the historical pre-PRESS integer-area histogram
// Method: Count all selected stars for any grouping type on the science grid
//         without changing upstream count-locus diagnostics.
// ==========================================
void populateSelectedGroupCountHistogram(
    const std::vector<int>& selected_star_areas,
    PSFCountLocusDiagnostics& diagnostics) {
    diagnostics.selected_group_count =
        static_cast<int>(selected_star_areas.size());
    populateCountHistogramOnScienceGrid(
        selected_star_areas,
        diagnostics.histogram_first_count,
        diagnostics.histogram_last_count,
        diagnostics.histogram.size(),
        diagnostics.selected_group_histogram);
}

// ==========================================
// Function: Parse matched Gaia image positions from one astro stream
// Method: Read the two WCS lines, the matched/user/reference counts, and exactly
//         n_matched finite RA/Dec/x/y rows while retaining only x/y.
// ==========================================
AstrometryGaiaReadStatus parseAstrometryGaiaPositions(
    std::istream& input,
    std::vector<std::array<double, 2>>& gaia_xy,
    std::string& error) {
    gaia_xy.clear();
    error.clear();
    std::string first_line;
    std::string second_line;
    if (!std::getline(input, first_line) || !std::getline(input, second_line)) {
        error = "missing astrometry WCS lines";
        return AstrometryGaiaReadStatus::Malformed;
    }

    int matched = 0;
    int user = 0;
    int reference = 0;
    if (!(input >> matched >> user >> reference) || matched < 0
        || user < 0 || reference < 0) {
        error = "invalid astrometry match-count row";
        return AstrometryGaiaReadStatus::Malformed;
    }

    gaia_xy.reserve(static_cast<std::size_t>(matched));
    for (int row = 0; row < matched; ++row) {
        double ra = 0.0;
        double dec = 0.0;
        double x = 0.0;
        double y = 0.0;
        if (!(input >> ra >> dec >> x >> y)
            || !std::isfinite(ra) || !std::isfinite(dec)
            || !std::isfinite(x) || !std::isfinite(y)) {
            error = "invalid astrometry matched-source row " + std::to_string(row);
            gaia_xy.clear();
            return AstrometryGaiaReadStatus::Malformed;
        }
        gaia_xy.push_back({x, y});
    }
    return gaia_xy.empty()
        ? AstrometryGaiaReadStatus::Empty
        : AstrometryGaiaReadStatus::Accepted;
}

// ==========================================
// Function: Match one candidate to the nearest same-chip Gaia position
// Method: Accept when the minimum squared image-plane distance is within the
//         configured radius; do not impose one-to-one assignment.
// ==========================================
bool hasNearestGaiaMatch(
    double x,
    double y,
    const std::vector<std::array<double, 2>>& gaia_xy,
    double radius_pixels) {
    if (!std::isfinite(x) || !std::isfinite(y)
        || !std::isfinite(radius_pixels) || radius_pixels < 0.0) {
        return false;
    }
    const double radius_squared = radius_pixels * radius_pixels;
    double minimum_squared = std::numeric_limits<double>::infinity();
    for (const std::array<double, 2>& gaia : gaia_xy) {
        const double dx = x - gaia[0];
        const double dy = y - gaia[1];
        minimum_squared = std::min(minimum_squared, dx * dx + dy * dy);
    }
    return minimum_squared <= radius_squared;
}

// ==========================================
// Function: Compute the exact normalized PSF chi distance
// Method: Apply the legacy sqrt(sum squared difference / mean signed flux)
//         directly to two cached central windows.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second) {
    if (first.empty() || first.size() != second.size()) {
        return std::numeric_limits<float>::infinity();
    }
    double flux = 0.0;
    double squared_difference = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        const double first_value = first[index];
        const double second_value = second[index];
        if (!std::isfinite(first_value) || !std::isfinite(second_value)) {
            return std::numeric_limits<float>::infinity();
        }
        flux += 0.5 * (first_value + second_value);
        const double difference = first_value - second_value;
        squared_difference += difference * difference;
    }
    if (!std::isfinite(flux) || flux <= 0.0
        || !std::isfinite(squared_difference)) {
        return std::numeric_limits<float>::infinity();
    }
    const double distance = std::sqrt(squared_difference / flux);
    return std::isfinite(distance)
        ? static_cast<float>(distance)
        : std::numeric_limits<float>::infinity();
}

// ==========================================
// Function: Maintain one exact sorted top-K neighbour list
// Method: Insert or improve the candidate edge, sort by chi/index, and truncate.
// ==========================================
void updateTopK(
    std::vector<NeighborEdge>& neighbours,
    int star_index,
    float chi,
    int k) {
    if (star_index < 0 || k <= 0 || !std::isfinite(chi)) return;
    auto existing = std::find_if(
        neighbours.begin(), neighbours.end(),
        [star_index](const NeighborEdge& edge) {
            return edge.star_index == star_index;
        });
    if (existing == neighbours.end()) {
        neighbours.push_back({star_index, chi});
    } else if (chi < existing->chi) {
        existing->chi = chi;
    }
    std::sort(
        neighbours.begin(), neighbours.end(),
        [](const NeighborEdge& first, const NeighborEdge& second) {
            if (first.chi != second.chi) return first.chi < second.chi;
            return first.star_index < second.star_index;
        });
    if (static_cast<int>(neighbours.size()) > k) {
        neighbours.resize(static_cast<std::size_t>(k));
    }
}

// ==========================================
// Function: Select capped exposure-wide large-size minChi references
// Method: Sort the top fraction by size/chip/star, then apply a per-chip cap
//         without filling from candidates below the exposure-wide pool.
// ==========================================
std::vector<MinChiReferenceCandidate> selectMinChiReferenceStars(
    const std::vector<MinChiReferenceCandidate>& locus_candidates,
    double reference_fraction,
    int maximum_per_chip) {
    std::vector<MinChiReferenceCandidate> ranked;
    if (!std::isfinite(reference_fraction) || reference_fraction <= 0.0
        || reference_fraction > 1.0 || maximum_per_chip <= 0) {
        return ranked;
    }
    ranked.reserve(locus_candidates.size());
    for (const MinChiReferenceCandidate& candidate : locus_candidates) {
        if (candidate.chip_index >= 0 && candidate.star_index >= 0
            && std::isfinite(candidate.size)) {
            ranked.push_back(candidate);
        }
    }
    std::sort(
        ranked.begin(), ranked.end(),
        [](const MinChiReferenceCandidate& first,
           const MinChiReferenceCandidate& second) {
            if (first.size != second.size) return first.size > second.size;
            if (first.chip_index != second.chip_index) {
                return first.chip_index < second.chip_index;
            }
            return first.star_index < second.star_index;
        });

    const std::size_t pool_count = std::min(
        ranked.size(),
        static_cast<std::size_t>(std::ceil(
            reference_fraction * static_cast<double>(ranked.size()))));
    std::vector<MinChiReferenceCandidate> selected;
    selected.reserve(pool_count);
    std::vector<int> chip_counts;
    for (std::size_t index = 0; index < pool_count; ++index) {
        const MinChiReferenceCandidate& candidate = ranked[index];
        if (candidate.chip_index >= static_cast<int>(chip_counts.size())) {
            chip_counts.resize(
                static_cast<std::size_t>(candidate.chip_index + 1), 0);
        }
        if (chip_counts[candidate.chip_index] >= maximum_per_chip) continue;
        chip_counts[candidate.chip_index]++;
        selected.push_back(candidate);
    }
    return selected;
}

// ==========================================
// Function: Compute one chip's minChi values and threshold-pair sample
// Method: Visit every unordered locus-locus pair once, update both endpoints,
//         and sample the distance when either endpoint is a reference.
// ==========================================
MinChiPairResult computeMinChiAndThresholdPairs(
    const std::vector<MinChiCandidateView>& candidates) {
    MinChiPairResult result;
    result.min_chi.assign(
        candidates.size(), std::numeric_limits<float>::infinity());
    for (std::size_t first = 0; first + 1 < candidates.size(); ++first) {
        if (!candidates[first].in_size_locus
            || candidates[first].chi_window == nullptr) {
            continue;
        }
        for (std::size_t second = first + 1;
             second < candidates.size(); ++second) {
            if (!candidates[second].in_size_locus
                || candidates[second].chi_window == nullptr) {
                continue;
            }
            const float chi = normalizedChiDistance(
                *candidates[first].chi_window,
                *candidates[second].chi_window);
            if (!std::isfinite(chi)) continue;
            result.min_chi[first] = std::min(result.min_chi[first], chi);
            result.min_chi[second] = std::min(result.min_chi[second], chi);
            if (candidates[first].is_reference
                || candidates[second].is_reference) {
                result.threshold_pair_chi.push_back(chi);
            }
        }
    }
    return result;
}

// ==========================================
// Function: Extract mutual-KNN graph edges among active candidates
// Method: Keep an undirected edge only when both retained top-K lists contain
//         the opposite endpoint and both endpoints survived the shared cut.
// ==========================================
std::vector<GraphEdge> buildMutualKNNEdges(
    const std::vector<int>& active_indices,
    const std::vector<std::vector<NeighborEdge>>& neighbours_by_star) {
    std::vector<GraphEdge> edges;
    if (active_indices.empty()) return edges;
    const int maximum_index = *std::max_element(
        active_indices.begin(), active_indices.end());
    std::vector<bool> active(static_cast<std::size_t>(maximum_index + 1), false);
    for (int index : active_indices) {
        if (index >= 0) active[index] = true;
    }
    for (int first : active_indices) {
        if (first < 0 || first >= static_cast<int>(neighbours_by_star.size())) continue;
        for (const NeighborEdge& neighbour : neighbours_by_star[first]) {
            const int second = neighbour.star_index;
            if (second <= first || second > maximum_index || !active[second]
                || second >= static_cast<int>(neighbours_by_star.size())) {
                continue;
            }
            if (containsNeighbour(neighbours_by_star[second], first)) {
                edges.push_back({first, second});
            }
        }
    }
    return edges;
}

// ==========================================
// Function: Convert one same-chip graph into connected components
// Method: Use disjoint sets over active original indices and attach Gaia counts.
// ==========================================
std::vector<StarGroup> buildConnectedGroups(
    const std::vector<int>& active_indices,
    const std::vector<GraphEdge>& edges,
    const std::vector<bool>& gaia_matched) {
    std::vector<StarGroup> groups;
    if (active_indices.empty()) return groups;

    const int maximum_index = *std::max_element(
        active_indices.begin(), active_indices.end());
    std::vector<int> local_index(static_cast<std::size_t>(maximum_index + 1), -1);
    for (int local = 0; local < static_cast<int>(active_indices.size()); ++local) {
        const int original = active_indices[local];
        if (original >= 0) local_index[original] = local;
    }
    DisjointSet disjoint_set(static_cast<int>(active_indices.size()));
    for (const GraphEdge& edge : edges) {
        if (edge.first < 0 || edge.second < 0
            || edge.first > maximum_index || edge.second > maximum_index) {
            continue;
        }
        const int first_local = local_index[edge.first];
        const int second_local = local_index[edge.second];
        if (first_local >= 0 && second_local >= 0) {
            disjoint_set.unite(first_local, second_local);
        }
    }

    std::vector<int> root_to_group(active_indices.size(), -1);
    for (int local = 0; local < static_cast<int>(active_indices.size()); ++local) {
        const int root = disjoint_set.find(local);
        if (root_to_group[root] < 0) {
            root_to_group[root] = static_cast<int>(groups.size());
            groups.push_back({});
        }
        StarGroup& group = groups[root_to_group[root]];
        const int original = active_indices[local];
        group.members.push_back(original);
        if (original >= 0 && original < static_cast<int>(gaia_matched.size())
            && gaia_matched[original]) {
            group.gaia_count++;
        }
    }
    return groups;
}

// ==========================================
// Function: Select the shared main and eligible secondary stellar groups
// Method: Always keep the largest component and keep another component only
//         when both its relative size and Gaia-count requirements pass.
// ==========================================
std::vector<int> selectMainAndSecondaryGroups(
    const std::vector<StarGroup>& groups,
    double minimum_size_ratio,
    int minimum_gaia_count) {
    std::vector<int> selected;
    if (groups.empty() || minimum_size_ratio < 0.0 || minimum_gaia_count < 0) {
        return selected;
    }
    std::size_t main_group = 0;
    for (std::size_t group = 1; group < groups.size(); ++group) {
        if (groups[group].members.size() > groups[main_group].members.size()) {
            main_group = group;
        }
    }
    const double main_size = static_cast<double>(groups[main_group].members.size());
    for (std::size_t group = 0; group < groups.size(); ++group) {
        const bool keep_main = group == main_group;
        const bool size_ok = main_size > 0.0
            && static_cast<double>(groups[group].members.size())
                >= minimum_size_ratio * main_size;
        const bool gaia_ok = groups[group].gaia_count >= minimum_gaia_count;
        if (keep_main || (size_ok && gaia_ok)) {
            selected.insert(
                selected.end(), groups[group].members.begin(), groups[group].members.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    return selected;
}

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
    double& loo_model) {
    constexpr double leverage_tolerance = 1.0e-10;
    if (!std::isfinite(observed) || !std::isfinite(fitted)
        || !std::isfinite(leverage) || !std::isfinite(minimum_denominator)
        || minimum_denominator <= 0.0 || leverage < -leverage_tolerance
        || leverage >= 1.0 - minimum_denominator) {
        return false;
    }
    if (leverage < 0.0) leverage = 0.0;
    const double denominator = 1.0 - leverage;
    loo_residual = (observed - fitted) / denominator;
    loo_model = observed - loo_residual;
    return std::isfinite(loo_residual) && std::isfinite(loo_model);
}

// ==========================================
// Function: Convert raw analytic PRESS to its leverage-standardized score
// Method: Multiply by sqrt(1-h) after the same finite denominator guard used
//         by the analytic leave-one-out calculation.
// ==========================================
bool computeLeverageStandardizedPress(
    double raw_press,
    double leverage,
    double minimum_denominator,
    double& standardized_press) {
    const double denominator = 1.0 - leverage;
    if (!std::isfinite(raw_press) || !std::isfinite(denominator)
        || !std::isfinite(minimum_denominator)
        || minimum_denominator <= 0.0
        || denominator <= minimum_denominator) {
        return false;
    }
    standardized_press = raw_press * std::sqrt(denominator);
    return std::isfinite(standardized_press);
}

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
    int maximum_removals) {
    if (!rejection_enabled) return PressRemovalDecision::Disabled;
    if (flagged_count <= 0) return PressRemovalDecision::NoOutliers;
    if (maximum_removals < 0 || flagged_count > maximum_removals) {
        return PressRemovalDecision::TooManyOutliers;
    }
    if (initial_count - flagged_count < minimum_count) {
        return PressRemovalDecision::WouldUnderrunMinimum;
    }
    return PressRemovalDecision::Apply;
}

}  // namespace Internal
}  // namespace PSFModel
