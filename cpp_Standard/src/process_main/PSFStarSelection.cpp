#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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
// Function: Return a linearly interpolated quantile from sorted samples
// Method: Interpolate between the two enclosing zero-based order statistics.
// ==========================================
double sortedQuantile(const std::vector<double>& sorted, double quantile) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = std::clamp(quantile, 0.0, 1.0)
        * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction;
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
// Function: Find the smallest positive sample spacing
// Method: Scan adjacent sorted values and ignore exact duplicates.
// ==========================================
double minimumPositiveSpacing(const std::vector<double>& sorted) {
    double spacing = std::numeric_limits<double>::infinity();
    for (std::size_t index = 1; index < sorted.size(); ++index) {
        const double difference = sorted[index] - sorted[index - 1];
        if (difference > 0.0 && std::isfinite(difference)) {
            spacing = std::min(spacing, difference);
        }
    }
    return spacing;
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

}  // namespace

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
    FWHMLocus& locus) {
    locus = {};
    if (histogram_bins < 3 || sigma_cut <= 0.0 || minimum_samples <= 0) {
        return false;
    }

    std::vector<double> values;
    std::vector<double> gaia_values;
    values.reserve(samples.size());
    gaia_values.reserve(samples.size());
    for (const FWHMSample& sample : samples) {
        if (!std::isfinite(sample.fwhm) || sample.fwhm <= 0.0) continue;
        values.push_back(sample.fwhm);
        if (sample.gaia_matched) gaia_values.push_back(sample.fwhm);
    }
    if (static_cast<int>(values.size()) < minimum_samples) return false;

    std::sort(values.begin(), values.end());
    std::sort(gaia_values.begin(), gaia_values.end());
    const double minimum_spacing = minimumPositiveSpacing(values);
    double range_low = sortedQuantile(values, 0.01);
    double range_high = sortedQuantile(values, 0.99);

    if (!(range_high > range_low)) {
        const double center = sortedMedian(values);
        const double floor = std::isfinite(minimum_spacing)
            ? minimum_spacing
            : std::max(std::abs(center) * 1.0e-6, 1.0e-6);
        locus.valid = true;
        locus.center = center;
        locus.width = floor;
        locus.lower = center - sigma_cut * floor;
        locus.upper = center + sigma_cut * floor;
        locus.histogram_bin_width = floor;
        return true;
    }

    const double bin_width = (range_high - range_low)
        / static_cast<double>(histogram_bins);
    if (!std::isfinite(bin_width) || bin_width <= 0.0) return false;

    std::vector<double> histogram(static_cast<std::size_t>(histogram_bins), 0.0);
    for (double value : values) {
        if (value < range_low || value > range_high) continue;
        int bin = static_cast<int>((value - range_low) / bin_width);
        bin = std::clamp(bin, 0, histogram_bins - 1);
        histogram[bin] += 1.0;
    }

    const int offsets[] = {-2, -1, 0, 1, 2};
    const double weights[] = {1.0, 2.0, 3.0, 2.0, 1.0};
    std::vector<double> smoothed(static_cast<std::size_t>(histogram_bins), 0.0);
    for (int bin = 0; bin < histogram_bins; ++bin) {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (int offset_index = 0; offset_index < 5; ++offset_index) {
            const int neighbour = bin + offsets[offset_index];
            if (neighbour < 0 || neighbour >= histogram_bins) continue;
            weighted_sum += weights[offset_index] * histogram[neighbour];
            weight_sum += weights[offset_index];
        }
        smoothed[bin] = weighted_sum / weight_sum;
    }

    std::vector<int> peaks;
    for (int bin = 0; bin < histogram_bins; ++bin) {
        const double left = bin == 0 ? -1.0 : smoothed[bin - 1];
        const double right = bin + 1 == histogram_bins ? -1.0 : smoothed[bin + 1];
        if (smoothed[bin] > 0.0 && smoothed[bin] >= left && smoothed[bin] >= right) {
            peaks.push_back(bin);
        }
    }
    if (peaks.empty()) return false;

    int peak_bin = peaks.front();
    if (static_cast<int>(gaia_values.size()) >= minimum_gaia_matches) {
        const double gaia_median = sortedMedian(gaia_values);
        double best_distance = std::numeric_limits<double>::infinity();
        for (int candidate_peak : peaks) {
            const double peak_center = range_low
                + (static_cast<double>(candidate_peak) + 0.5) * bin_width;
            const double distance = std::abs(peak_center - gaia_median);
            if (distance < best_distance) {
                best_distance = distance;
                peak_bin = candidate_peak;
            }
        }
    } else {
        for (int candidate_peak : peaks) {
            if (smoothed[candidate_peak] > smoothed[peak_bin]) {
                peak_bin = candidate_peak;
            }
        }
    }

    int basin_first = peak_bin;
    while (basin_first > 0
           && smoothed[basin_first - 1] <= smoothed[basin_first]) {
        basin_first--;
    }
    int basin_last = peak_bin;
    while (basin_last + 1 < histogram_bins
           && smoothed[basin_last + 1] <= smoothed[basin_last]) {
        basin_last++;
    }

    const double basin_low = range_low + static_cast<double>(basin_first) * bin_width;
    const double basin_high = range_low
        + static_cast<double>(basin_last + 1) * bin_width;
    std::vector<double> population;
    for (double value : values) {
        if (value >= basin_low && value <= basin_high) population.push_back(value);
    }
    if (population.size() < 3) {
        const int fallback_first = std::max(0, peak_bin - 2);
        const int fallback_last = std::min(histogram_bins - 1, peak_bin + 2);
        const double fallback_low = range_low
            + static_cast<double>(fallback_first) * bin_width;
        const double fallback_high = range_low
            + static_cast<double>(fallback_last + 1) * bin_width;
        population.clear();
        for (double value : values) {
            if (value >= fallback_low && value <= fallback_high) {
                population.push_back(value);
            }
        }
    }
    if (population.empty()) return false;

    const double spacing_floor = std::isfinite(minimum_spacing)
        ? std::max(bin_width, minimum_spacing)
        : bin_width;
    double center = 0.0;
    double width = 0.0;
    for (int iteration = 0; iteration < 2; ++iteration) {
        if (!medianAndMad(population, center, width)) return false;
        width = std::max(width, spacing_floor);
        std::vector<double> clipped;
        clipped.reserve(population.size());
        for (double value : population) {
            if (std::abs(value - center) <= sigma_cut * width) {
                clipped.push_back(value);
            }
        }
        if (clipped.empty() || clipped.size() == population.size()) break;
        population.swap(clipped);
    }
    if (!medianAndMad(population, center, width)) return false;
    width = std::max(width, spacing_floor);

    locus.valid = true;
    locus.center = center;
    locus.width = width;
    locus.lower = center - sigma_cut * width;
    locus.upper = center + sigma_cut * width;
    locus.histogram_bin_width = bin_width;
    return std::isfinite(locus.lower) && std::isfinite(locus.upper);
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
        if (!candidates[first].in_fwhm_locus
            || candidates[first].chi_window == nullptr) {
            continue;
        }
        for (std::size_t second = first + 1;
             second < candidates.size(); ++second) {
            if (!candidates[second].in_fwhm_locus
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
