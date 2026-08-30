#include "process_fd/StarCutCalculator.hpp"
#include "FDConfig.hpp"
#include "LensingConfig.hpp"
#include "general/MPIScheduler.hpp"
#include "process_main/MPIFailure.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fc = FDConfig;

namespace {

// ==========================================
// Function: Compute a checked flattened FD array size
// Method: Multiply runtime exposure and fixed histogram dimensions in size_t
//         and abort the MPI world before overflow could under-allocate storage.
// ==========================================
std::size_t checkedElementCount(
    std::initializer_list<std::size_t> dimensions,
    const std::string& operation) {
    std::size_t count = 1;
    for (const std::size_t dimension : dimensions) {
        if (dimension != 0
            && count > std::numeric_limits<std::size_t>::max() / dimension) {
            MPIFailure::abortWorld(
                operation, "runtime dimensions overflow size_t");
        }
        count *= dimension;
    }
    return count;
}

// ==========================================
// Function: Sum one runtime-sized vector across MPI ranks
// Method: Split the reduction into INT_MAX-bounded collectives so the dynamic
//         exposure count is not constrained by one MPI count argument.
// ==========================================
template <typename T>
void allreduceSum(const std::vector<T>& local, std::vector<T>& global,
                  MPI_Datatype datatype) {
    if (local.size() != global.size()) {
        MPIFailure::abortWorld(
            "reduce FD runtime vectors", "local/global sizes differ");
    }
    std::size_t offset = 0;
    const std::size_t max_chunk =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    while (offset < local.size()) {
        const int chunk = static_cast<int>(
            std::min(max_chunk, local.size() - offset));
        MPI_Allreduce(local.data() + offset, global.data() + offset, chunk,
                      datatype, MPI_SUM, MPIScheduler::state.communicator);
        offset += static_cast<std::size_t>(chunk);
    }
}

}  // namespace

// ==========================================
// Function: Calculate one global star cut
// Method: Build the size-selected-magnitude histogram across all ranks and
//         derive a shared stellar-locus size threshold.
// ==========================================
void StarCutCalculator::calculateGlobalStarCut(const FDData& data,
                                               float& S_mean, float& S_std,
                                               float& S_cut) {
    const int rank = MPIScheduler::state.rank;
    const MPI_Comm communicator = MPIScheduler::state.communicator;
    S_mean = 0.0; S_std = -1.0; S_cut = 0.0;

    if (fc::star_bar_mltp <= 0.0) { S_cut = 0.0; return; }
    float k_sigma = fc::star_bar_mltp;

    const int ns = fc::n_size_bins, nm = fc::n_mag_bins;
    float size_bin_w = (fc::size_max - fc::size_min) / float(ns);
    float mag_bin_w = (fc::mag_max_val - fc::mag_min_val) / float(nm);
    int size_limit_idx = int((-0.5 - fc::size_min) / size_bin_w);  // 0-based
    if (size_limit_idx < 0) size_limit_idx = 0;
    if (size_limit_idx >= ns) size_limit_idx = ns - 1;

    // Local 2D histogram
    std::vector<int> hist2d(ns * nm, 0), global_hist(ns * nm, 0);
    std::vector<int> mag_count(nm, 0), global_mag_count(nm, 0);

    for (int idx = 0; idx < data.ng; ++idx) {
        int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
        int i = int((data.sizerel[idx] - fc::size_min) / size_bin_w);
        if (j >= 0 && j < nm) {
            mag_count[j]++;
            if (i >= 0 && i < ns)
                hist2d[j * ns + i]++;
        }
    }

    MPI_Allreduce(hist2d.data(), global_hist.data(), ns * nm, MPI_INT,
                  MPI_SUM, communicator);
    MPI_Allreduce(mag_count.data(), global_mag_count.data(), nm, MPI_INT,
                  MPI_SUM, communicator);

    // Find peak concentration
    float max_concentration = 0.0, S_init = 0.5;
    int best_j = -1;
    for (int j = 0; j < nm; ++j) {
        if (global_mag_count[j] < fc::min_bin_count) continue;
        int peak = size_limit_idx;
        for (int i = size_limit_idx + 1; i < ns; ++i)
            if (global_hist[j * ns + i] > global_hist[j * ns + peak]) peak = i;
        int sum_peak = global_hist[j * ns + peak];
        if (peak > 0) sum_peak += global_hist[j * ns + peak - 1];
        if (peak < ns - 1) sum_peak += global_hist[j * ns + peak + 1];
        float concentration = float(sum_peak) / float(global_mag_count[j]);
        if (concentration > max_concentration) {
            max_concentration = concentration;
            best_j = j;
            S_init = fc::size_min + (peak + 0.5) * size_bin_w;
        }
    }

    // Determine active magnitude bins
    std::vector<bool> active_mag_bins(nm, false);
    if (best_j >= 0) {
        for (int j = 0; j < nm; ++j) {
            if (global_mag_count[j] < fc::min_bin_count) continue;
            int peak = size_limit_idx;
            for (int i = size_limit_idx + 1; i < ns; ++i)
                if (global_hist[j * ns + i] > global_hist[j * ns + peak]) peak = i;
            float peak_size = fc::size_min + (peak + 0.5) * size_bin_w;
            int sum_peak = global_hist[j * ns + peak];
            if (peak > 0) sum_peak += global_hist[j * ns + peak - 1];
            if (peak < ns - 1) sum_peak += global_hist[j * ns + peak + 1];
            float concentration = float(sum_peak) / float(global_mag_count[j]);
            if (std::fabs(peak_size - S_init) <= fc::peak_match_tol &&
                concentration >= fc::min_concentration)
                active_mag_bins[j] = true;
        }
    }

    // Compute global mean and std of star sequence size
    float local_sum = 0.0;
    int local_count = 0;
    for (int idx = 0; idx < data.ng; ++idx) {
        int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm && active_mag_bins[j] &&
            std::fabs(data.sizerel[idx] - S_init) < 0.1) {
            local_sum += data.sizerel[idx];
            local_count++;
        }
    }
    float global_sum = 0.0;
    int global_count = 0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, communicator);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, communicator);

    if (global_count > 1) {
        S_mean = global_sum / float(global_count);
        float local_sq_diff = 0.0;
        for (int idx = 0; idx < data.ng; ++idx) {
            int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm && active_mag_bins[j] &&
                std::fabs(data.sizerel[idx] - S_init) < 0.1)
                local_sq_diff += (data.sizerel[idx] - S_mean) * (data.sizerel[idx] - S_mean);
        }
        float global_sq_diff = 0.0;
        MPI_Allreduce(&local_sq_diff, &global_sq_diff, 1, MPI_FLOAT, MPI_SUM,
                      communicator);
        S_std = std::sqrt(global_sq_diff / float(global_count - 1));
        S_cut = S_mean + k_sigma * S_std;
    } else {
        if (best_j >= 0) S_cut = S_init + 0.05f * k_sigma;
        else S_cut = 0.6;
        S_mean = S_init;
        S_std = 0.05;
    }

    if (rank == 0)
        std::cout << "Global star cut: S_mean=" << S_mean
                  << " S_std=" << S_std << " S_cut=" << S_cut << std::endl;
}

// ==========================================
// Function: Calculate per-exposure global star cuts
// Method: Build size-selected-magnitude histograms per exposure, combine them
//         across ranks, and iteratively estimate each stellar-locus threshold.
// ==========================================
void StarCutCalculator::calculateGlobalStarCutAuto(
    const FDData& data,
    std::vector<float>& S_mean_arr,
    std::vector<float>& S_std_arr,
    std::vector<float>& S_cut_arr) {
    const int rank = MPIScheduler::state.rank;
    const MPI_Comm communicator = MPIScheduler::state.communicator;

    if (fc::star_bar_mltp <= 0.0) return;
    float k_sigma = fc::star_bar_mltp;

    const int ns = fc::n_size_bins, nm = fc::n_mag_bins;
    float size_bin_w = (fc::size_max - fc::size_min) / float(ns);
    float mag_bin_w = (fc::mag_max_val - fc::mag_min_val) / float(nm);
    int star_min_idx = int((fc::star_phy_min - fc::size_min) / size_bin_w);
    int star_max_idx = int((fc::star_phy_max - fc::size_min) / size_bin_w);
    if (star_min_idx < 0) star_min_idx = 0;
    if (star_max_idx >= ns) star_max_idx = ns - 1;

    // Find max exposure index across all nodes
    int local_max_iex = 1;
    for (int idx = 0; idx < data.ng; ++idx)
        if (data.iexpo[idx] > local_max_iex) local_max_iex = data.iexpo[idx];
    int global_max_iex = 1;
    MPI_Allreduce(&local_max_iex, &global_max_iex, 1, MPI_INT, MPI_MAX,
                  communicator);
    if (global_max_iex < 1) global_max_iex = 1;
    const std::size_t exposure_count =
        static_cast<std::size_t>(global_max_iex);
    S_mean_arr.assign(exposure_count, 0.0f);
    S_std_arr.assign(exposure_count, 0.0f);
    S_cut_arr.assign(exposure_count, 0.0f);

    // 3D histogram (ns × nm × global_max_iex) — flattened
    const std::size_t histogram_elements = checkedElementCount(
        {exposure_count, static_cast<std::size_t>(nm),
         static_cast<std::size_t>(ns)},
        "allocate FD exposure histograms");
    const std::size_t magnitude_elements = checkedElementCount(
        {exposure_count, static_cast<std::size_t>(nm)},
        "allocate FD exposure magnitude counts");
    std::vector<int> hist3d(histogram_elements, 0);
    std::vector<int> global_hist3d(histogram_elements, 0);
    std::vector<int> mag_count3d(magnitude_elements, 0);
    std::vector<int> global_mag_count3d(magnitude_elements, 0);
    const auto magnitudeIndex = [nm](int iex, int magnitude_bin) {
        return static_cast<std::size_t>(iex - 1)
             * static_cast<std::size_t>(nm)
             + static_cast<std::size_t>(magnitude_bin);
    };
    const auto histogramIndex = [ns, &magnitudeIndex](
                                    int iex, int magnitude_bin,
                                    int size_bin) {
        return magnitudeIndex(iex, magnitude_bin)
             * static_cast<std::size_t>(ns)
             + static_cast<std::size_t>(size_bin);
    };

    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex) continue;
        int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
        int i = int((data.sizerel[idx] - fc::size_min) / size_bin_w);
        if (j >= 0 && j < nm) {
            mag_count3d[magnitudeIndex(iex, j)]++;
            if (i >= 0 && i < ns)
                hist3d[histogramIndex(iex, j, i)]++;
        }
    }

    allreduceSum(hist3d, global_hist3d, MPI_INT);
    allreduceSum(mag_count3d, global_mag_count3d, MPI_INT);

    // Per-exposure peak analysis
    std::vector<float> S_init_arr(exposure_count, fc::default_s_init);
    std::vector<bool> active_mag_bins3d(magnitude_elements, false);
    std::vector<bool> use_fallback(exposure_count, false);
    std::vector<bool> skip_iter(exposure_count, false);

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        float max_concentration = 0.0;
        int best_j = -1;
        for (int j = 0; j < nm; ++j) {
            int mc = global_mag_count3d[magnitudeIndex(iex, j)];
            if (mc < fc::min_bin_count) continue;
            int peak = star_min_idx;
            for (int i = star_min_idx + 1; i <= star_max_idx; ++i)
                if (global_hist3d[histogramIndex(iex, j, i)] >
                    global_hist3d[histogramIndex(iex, j, peak)]) peak = i;
            int sum_peak = global_hist3d[histogramIndex(iex, j, peak)];
            if (peak > 0) {
                sum_peak += global_hist3d[histogramIndex(iex, j, peak - 1)];
            }
            if (peak < ns - 1) {
                sum_peak += global_hist3d[histogramIndex(iex, j, peak + 1)];
            }
            float concentration = float(sum_peak) / float(mc);
            if (concentration > max_concentration) {
                max_concentration = concentration;
                best_j = j;
                S_init_arr[eidx] = fc::size_min + (peak + 0.5) * size_bin_w;
            }
        }
        // Determine active mag bins for this exposure
        if (best_j >= 0) {
            for (int j = 0; j < nm; ++j) {
                int mc = global_mag_count3d[magnitudeIndex(iex, j)];
                if (mc < fc::min_bin_count) continue;
                int peak = star_min_idx;
                for (int i = star_min_idx + 1; i <= star_max_idx; ++i)
                    if (global_hist3d[histogramIndex(iex, j, i)] >
                        global_hist3d[histogramIndex(iex, j, peak)]) peak = i;
                float peak_size = fc::size_min + (peak + 0.5) * size_bin_w;
                int sum_peak = global_hist3d[histogramIndex(iex, j, peak)];
                if (peak > 0) {
                    sum_peak += global_hist3d[
                        histogramIndex(iex, j, peak - 1)];
                }
                if (peak < ns - 1) {
                    sum_peak += global_hist3d[
                        histogramIndex(iex, j, peak + 1)];
                }
                float concentration = float(sum_peak) / float(mc);
                if (std::fabs(peak_size - S_init_arr[eidx]) <= fc::peak_match_tol &&
                    concentration >= fc::min_concentration)
                    active_mag_bins3d[magnitudeIndex(iex, j)] = true;
            }
        }
    }

    // Initialize mean using stage1_snr
    std::vector<float> local_sum_arr(exposure_count, 0.0f);
    std::vector<float> global_sum_arr(exposure_count, 0.0f);
    std::vector<int> local_count_arr(exposure_count, 0);
    std::vector<int> global_count_arr(exposure_count, 0);

    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex) continue;
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm
            && active_mag_bins3d[magnitudeIndex(iex, j)] &&
            std::fabs(data.sizerel[idx] - S_init_arr[eidx]) < fc::init_win_active) {
            local_sum_arr[eidx] += data.sizerel[idx];
            local_count_arr[eidx]++;
        }
    }
    allreduceSum(local_sum_arr, global_sum_arr, MPI_FLOAT);
    allreduceSum(local_count_arr, global_count_arr, MPI_INT);

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        use_fallback[eidx] = (global_count_arr[eidx] <= 1);
    }

    // Fallback with wider window
    std::fill(local_sum_arr.begin(), local_sum_arr.end(), 0.0);
    std::fill(local_count_arr.begin(), local_count_arr.end(), 0);
    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex) continue;
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        if (!use_fallback[eidx]) continue;
        int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm &&
            std::fabs(data.sizerel[idx] - S_init_arr[eidx]) < fc::init_win_fallback) {
            local_sum_arr[eidx] += data.sizerel[idx];
            local_count_arr[eidx]++;
        }
    }
    std::vector<float> fb_sum(exposure_count, 0.0f);
    std::vector<int> fb_count(exposure_count, 0);
    allreduceSum(local_sum_arr, fb_sum, MPI_FLOAT);
    allreduceSum(local_count_arr, fb_count, MPI_INT);
    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        if (use_fallback[eidx]) {
            global_sum_arr[eidx] = fb_sum[eidx];
            global_count_arr[eidx] = fb_count[eidx];
        }
    }

    // Initialize S_mean_temp and S_std_temp
    std::vector<float> S_mean_t(exposure_count, fc::default_s_init);
    std::vector<float> S_std_t(exposure_count, fc::default_s_std);
    std::vector<float> clip_limit(exposure_count, fc::init_win_active);

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        if (global_count_arr[eidx] > 1) {
            skip_iter[eidx] = false;
            S_mean_t[eidx] = global_sum_arr[eidx] / float(global_count_arr[eidx]);
            S_std_t[eidx] = use_fallback[eidx] ? fc::init_win_fallback : fc::init_win_active;
        } else {
            skip_iter[eidx] = true;
            S_mean_t[eidx] = S_init_arr[eidx];
            S_std_t[eidx] = fc::default_s_std;
        }
    }

    // Iterative sigma-clipping (3 iterations)
    for (int i_iter = 0; i_iter < 3; ++i_iter) {
        for (int iex = 1; iex <= global_max_iex; ++iex) {
            const std::size_t eidx = static_cast<std::size_t>(iex - 1);
            if (skip_iter[eidx]) continue;
            float iw = use_fallback[eidx] ? fc::init_win_fallback : fc::init_win_active;
            if (i_iter == 0) clip_limit[eidx] = iw;
            else if (i_iter == 1) clip_limit[eidx] = 0.5f * iw;
            else {
                clip_limit[eidx] = fc::clip_nsigma * S_std_t[eidx];
                if (clip_limit[eidx] > 0.5f * iw) clip_limit[eidx] = 0.5f * iw;
                if (clip_limit[eidx] < fc::min_clip_limit) clip_limit[eidx] = fc::min_clip_limit;
            }
        }

        // A. Compute mean
        std::fill(local_sum_arr.begin(), local_sum_arr.end(), 0.0);
        std::fill(local_count_arr.begin(), local_count_arr.end(), 0);
        for (int idx = 0; idx < data.ng; ++idx) {
            if (data.src_snr[idx] <= fc::stage1_snr) continue;
            int iex = data.iexpo[idx];
            if (iex < 1 || iex > global_max_iex) continue;
            const std::size_t eidx = static_cast<std::size_t>(iex - 1);
            if (skip_iter[eidx]) continue;
            int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm &&
                (use_fallback[eidx]
                 || active_mag_bins3d[magnitudeIndex(iex, j)]) &&
                std::fabs(data.sizerel[idx] - S_mean_t[eidx]) < clip_limit[eidx]) {
                local_sum_arr[eidx] += data.sizerel[idx];
                local_count_arr[eidx]++;
            }
        }
        allreduceSum(local_sum_arr, global_sum_arr, MPI_FLOAT);
        allreduceSum(local_count_arr, global_count_arr, MPI_INT);
        for (int iex = 1; iex <= global_max_iex; ++iex) {
            const std::size_t eidx = static_cast<std::size_t>(iex - 1);
            if (!skip_iter[eidx]) {
                if (global_count_arr[eidx] > 1)
                    S_mean_t[eidx] = global_sum_arr[eidx] / float(global_count_arr[eidx]);
                else skip_iter[eidx] = true;
            }
        }

        // B. Compute std
        std::vector<float> local_sq(exposure_count, 0.0f);
        std::vector<float> global_sq(exposure_count, 0.0f);
        std::vector<int> local_cnt_s(exposure_count, 0);
        std::vector<int> global_cnt_s(exposure_count, 0);
        for (int idx = 0; idx < data.ng; ++idx) {
            if (data.src_snr[idx] <= fc::stage2_snr) continue;
            int iex = data.iexpo[idx];
            if (iex < 1 || iex > global_max_iex) continue;
            const std::size_t eidx = static_cast<std::size_t>(iex - 1);
            if (skip_iter[eidx]) continue;
            int j = int((data.star_mag[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm &&
                (use_fallback[eidx]
                 || active_mag_bins3d[magnitudeIndex(iex, j)]) &&
                std::fabs(data.sizerel[idx] - S_mean_t[eidx]) < clip_limit[eidx]) {
                local_sq[eidx] += (data.sizerel[idx] - S_mean_t[eidx]) * (data.sizerel[idx] - S_mean_t[eidx]);
                local_cnt_s[eidx]++;
            }
        }
        allreduceSum(local_sq, global_sq, MPI_FLOAT);
        allreduceSum(local_cnt_s, global_cnt_s, MPI_INT);
        for (int iex = 1; iex <= global_max_iex; ++iex) {
            const std::size_t eidx = static_cast<std::size_t>(iex - 1);
            if (!skip_iter[eidx]) {
                if (global_cnt_s[eidx] > 1)
                    S_std_t[eidx] = std::sqrt(global_sq[eidx] / float(global_cnt_s[eidx] - 1));
                else S_std_t[eidx] = fc::default_s_std;
                if (S_std_t[eidx] < fc::min_clip_limit / fc::clip_nsigma)
                    S_std_t[eidx] = fc::min_clip_limit / fc::clip_nsigma;
            }
        }
    }

    // Finalize: average successful exposures for fallback
    int n_success = 0;
    float sum_mean = 0.0, sum_std = 0.0;
    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        if (!skip_iter[eidx]) {
            n_success++;
            sum_mean += S_mean_t[eidx];
            sum_std += S_std_t[eidx];
        }
    }
    float avg_mean = n_success > 0 ? sum_mean / n_success : fc::default_s_init;
    float avg_std = n_success > 0 ? sum_std / n_success : fc::default_s_std;

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        const std::size_t eidx = static_cast<std::size_t>(iex - 1);
        if (!skip_iter[eidx]) {
            S_mean_arr[eidx] = S_mean_t[eidx];
            S_std_arr[eidx] = S_std_t[eidx];
        } else {
            S_mean_arr[eidx] = avg_mean;
            S_std_arr[eidx] = avg_std;
        }
        S_cut_arr[eidx] = S_mean_arr[eidx] + k_sigma * S_std_arr[eidx];
    }
    if (rank == 0)
        std::cout << "Per-exposure star cut done for " << global_max_iex
                  << " exposures." << std::endl;
}

// ------------------------------------------------------------------
// Apply single global star cut in-place
// ------------------------------------------------------------------
void StarCutCalculator::applySingleStarCut(FDData& data, float S_cut) {
    int write_idx = 0;
    for (int idx = 0; idx < data.ng; ++idx) {
        // Discard stars: small sizerel AND high SNR (star-like)
        if (data.sizerel[idx] <= S_cut && data.src_snr[idx] > 20.0) continue;
        data.x1[write_idx]  = data.x1[idx];
        data.x2[write_idx]  = data.x2[idx];
        data.rra[write_idx]  = data.rra[idx];
        data.ddec[write_idx] = data.ddec[idx];

        if constexpr (!fc::FD_USE_SWSE_DATA) {
            // PDF modes (1 & 2): keep raw shear and response variables
            data.y1[write_idx]  = data.y1[idx];
            data.de1[write_idx] = data.de1[idx];
            data.y2[write_idx]  = data.y2[idx];
            data.de2[write_idx] = data.de2[idx];
        } else {
            // SWSE mode (3): y = (g1/de)*ww, de = ww
            float de_val = (data.de1[idx] + data.de2[idx]) * 0.5f;
            data.y1[write_idx]  = (data.y1[idx] / de_val) * data.ww[idx];
            data.de1[write_idx] = data.ww[idx];
            data.y2[write_idx]  = (data.y2[idx] / de_val) * data.ww[idx];
            data.de2[write_idx] = data.ww[idx];
        }
        write_idx++;
    }
    data.ng = write_idx;
}

// ==========================================
// Function: Apply advanced per-exposure star and SNR cuts
// Method: Map serialized 1-based exposure IDs to zero-based cut-vector entries,
//         falling back to the largest cut only for an unavailable identity.
// ==========================================
void StarCutCalculator::applyAdvancedCuts(FDData& data,
                                          const std::vector<float>& S_cut_arr) {
    int write_idx = 0;
    float max_scut = 0.0;
    for (const float value : S_cut_arr)
        if (value > max_scut) max_scut = value;

    for (int idx = 0; idx < data.ng; ++idx) {
        // SNRF cut
        if (data.snrf[idx] < fc::snrfcut) continue;
        // SNR cuts
        if (fc::snrlow > 0.0 && data.src_snr[idx] < fc::snrlow) continue;
        if (fc::snrhigh > 0.0 && data.src_snr[idx] > fc::snrhigh) continue;

        // Per-exposure size cut
        int iex = data.iexpo[idx];
        float scut = max_scut;
        if (iex >= 1 && static_cast<std::size_t>(iex) <= S_cut_arr.size()) {
            scut = S_cut_arr[static_cast<std::size_t>(iex - 1)];
        }
        if (data.sizerel[idx] <= scut && data.src_snr[idx] > 20.0) continue;

        // Keep this galaxy
        data.x1[write_idx]  = data.x1[idx];
        data.x2[write_idx]  = data.x2[idx];
        data.rra[write_idx]  = data.rra[idx];
        data.ddec[write_idx] = data.ddec[idx];

        if constexpr (!fc::FD_USE_SWSE_DATA) {
            data.y1[write_idx]  = data.y1[idx];
            data.de1[write_idx] = data.de1[idx];
            data.y2[write_idx]  = data.y2[idx];
            data.de2[write_idx] = data.de2[idx];
        } else {
            float de_val = (data.de1[idx] + data.de2[idx]) * 0.5f;
            data.y1[write_idx]  = (data.y1[idx] / de_val) * data.ww[idx];
            data.de1[write_idx] = data.ww[idx];
            data.y2[write_idx]  = (data.y2[idx] / de_val) * data.ww[idx];
            data.de2[write_idx] = data.ww[idx];
        }
        write_idx++;
    }
    data.ng = write_idx;
    std::cout << "Combined cuts applied: kept " << data.ng << " sources." << std::endl;
}
