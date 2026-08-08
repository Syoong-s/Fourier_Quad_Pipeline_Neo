#include "process_fd/StarCutCalculator.hpp"
#include "FDConfig.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace fc = FDConfig;

// ------------------------------------------------------------------
// Single global star cut (all exposures share one bar)
// ------------------------------------------------------------------
void StarCutCalculator::calculateGlobalStarCut(const FDData& data, int rank,
                                               int num_procs,
                                               float& S_mean, float& S_std,
                                               float& S_cut) {
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
        int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
        int i = int((data.sizerel[idx] - fc::size_min) / size_bin_w);
        if (j >= 0 && j < nm) {
            mag_count[j]++;
            if (i >= 0 && i < ns)
                hist2d[j * ns + i]++;
        }
    }

    MPI_Allreduce(hist2d.data(), global_hist.data(), ns * nm, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(mag_count.data(), global_mag_count.data(), nm, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);

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
        int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm && active_mag_bins[j] &&
            std::fabs(data.sizerel[idx] - S_init) < 0.1) {
            local_sum += data.sizerel[idx];
            local_count++;
        }
    }
    float global_sum = 0.0;
    int global_count = 0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (global_count > 1) {
        S_mean = global_sum / float(global_count);
        float local_sq_diff = 0.0;
        for (int idx = 0; idx < data.ng; ++idx) {
            int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm && active_mag_bins[j] &&
                std::fabs(data.sizerel[idx] - S_init) < 0.1)
                local_sq_diff += (data.sizerel[idx] - S_mean) * (data.sizerel[idx] - S_mean);
        }
        float global_sq_diff = 0.0;
        MPI_Allreduce(&local_sq_diff, &global_sq_diff, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
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

// ------------------------------------------------------------------
// Per-exposure star cut (one bar per exposure)
// ------------------------------------------------------------------
void StarCutCalculator::calculateGlobalStarCutAuto(
    const FDData& data, int rank, int num_procs,
    std::vector<float>& S_mean_arr,
    std::vector<float>& S_std_arr,
    std::vector<float>& S_cut_arr) {

    const int NMAX_E = LensingConfig::NMAX_EXPO;
    S_mean_arr.assign(NMAX_E, 0.0);
    S_std_arr.assign(NMAX_E, 0.0);
    S_cut_arr.assign(NMAX_E, 0.0);

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
    MPI_Allreduce(&local_max_iex, &global_max_iex, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_max_iex > NMAX_E) global_max_iex = NMAX_E;
    if (global_max_iex < 1) global_max_iex = 1;

    // 3D histogram (ns × nm × global_max_iex) — flattened
    std::vector<int> hist3d(ns * nm * global_max_iex, 0);
    std::vector<int> global_hist3d(ns * nm * global_max_iex, 0);
    std::vector<int> mag_count3d(nm * global_max_iex, 0);
    std::vector<int> global_mag_count3d(nm * global_max_iex, 0);

    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex) continue;
        int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
        int i = int((data.sizerel[idx] - fc::size_min) / size_bin_w);
        if (j >= 0 && j < nm) {
            mag_count3d[(iex - 1) * nm + j]++;
            if (i >= 0 && i < ns)
                hist3d[((iex - 1) * nm + j) * ns + i]++;
        }
    }

    MPI_Allreduce(hist3d.data(), global_hist3d.data(),
                  ns * nm * global_max_iex, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(mag_count3d.data(), global_mag_count3d.data(),
                  nm * global_max_iex, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Per-exposure peak analysis
    std::vector<float> S_init_arr(NMAX_E, fc::default_s_init);
    std::vector<bool> active_mag_bins3d(nm * global_max_iex, false);
    std::vector<bool> use_fallback(NMAX_E, false);
    std::vector<bool> skip_iter(NMAX_E, false);

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        float max_concentration = 0.0;
        int best_j = -1;
        for (int j = 0; j < nm; ++j) {
            int mc = global_mag_count3d[(iex - 1) * nm + j];
            if (mc < fc::min_bin_count) continue;
            int peak = star_min_idx;
            for (int i = star_min_idx + 1; i <= star_max_idx; ++i)
                if (global_hist3d[((iex - 1) * nm + j) * ns + i] >
                    global_hist3d[((iex - 1) * nm + j) * ns + peak]) peak = i;
            int sum_peak = global_hist3d[((iex - 1) * nm + j) * ns + peak];
            if (peak > 0) sum_peak += global_hist3d[((iex - 1) * nm + j) * ns + peak - 1];
            if (peak < ns - 1) sum_peak += global_hist3d[((iex - 1) * nm + j) * ns + peak + 1];
            float concentration = float(sum_peak) / float(mc);
            if (concentration > max_concentration) {
                max_concentration = concentration;
                best_j = j;
                S_init_arr[iex] = fc::size_min + (peak + 0.5) * size_bin_w;
            }
        }
        // Determine active mag bins for this exposure
        if (best_j >= 0) {
            for (int j = 0; j < nm; ++j) {
                int mc = global_mag_count3d[(iex - 1) * nm + j];
                if (mc < fc::min_bin_count) continue;
                int peak = star_min_idx;
                for (int i = star_min_idx + 1; i <= star_max_idx; ++i)
                    if (global_hist3d[((iex - 1) * nm + j) * ns + i] >
                        global_hist3d[((iex - 1) * nm + j) * ns + peak]) peak = i;
                float peak_size = fc::size_min + (peak + 0.5) * size_bin_w;
                int sum_peak = global_hist3d[((iex - 1) * nm + j) * ns + peak];
                if (peak > 0) sum_peak += global_hist3d[((iex - 1) * nm + j) * ns + peak - 1];
                if (peak < ns - 1) sum_peak += global_hist3d[((iex - 1) * nm + j) * ns + peak + 1];
                float concentration = float(sum_peak) / float(mc);
                if (std::fabs(peak_size - S_init_arr[iex]) <= fc::peak_match_tol &&
                    concentration >= fc::min_concentration)
                    active_mag_bins3d[(iex - 1) * nm + j] = true;
            }
        }
    }

    // Initialize mean using stage1_snr
    std::vector<float> local_sum_arr(NMAX_E, 0.0), global_sum_arr(NMAX_E, 0.0);
    std::vector<int> local_count_arr(NMAX_E, 0), global_count_arr(NMAX_E, 0);

    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex) continue;
        int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm && active_mag_bins3d[(iex - 1) * nm + j] &&
            std::fabs(data.sizerel[idx] - S_init_arr[iex]) < fc::init_win_active) {
            local_sum_arr[iex] += data.sizerel[idx];
            local_count_arr[iex]++;
        }
    }
    MPI_Allreduce(local_sum_arr.data(), global_sum_arr.data(), global_max_iex,
                  MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_count_arr.data(), global_count_arr.data(), global_max_iex,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int iex = 1; iex <= global_max_iex; ++iex)
        use_fallback[iex] = (global_count_arr[iex] <= 1);

    // Fallback with wider window
    std::fill(local_sum_arr.begin(), local_sum_arr.end(), 0.0);
    std::fill(local_count_arr.begin(), local_count_arr.end(), 0);
    for (int idx = 0; idx < data.ng; ++idx) {
        if (data.src_snr[idx] <= fc::stage1_snr) continue;
        int iex = data.iexpo[idx];
        if (iex < 1 || iex > global_max_iex || !use_fallback[iex]) continue;
        int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
        if (j >= 0 && j < nm &&
            std::fabs(data.sizerel[idx] - S_init_arr[iex]) < fc::init_win_fallback) {
            local_sum_arr[iex] += data.sizerel[idx];
            local_count_arr[iex]++;
        }
    }
    std::vector<float> fb_sum(NMAX_E, 0.0);
    std::vector<int> fb_count(NMAX_E, 0);
    MPI_Allreduce(local_sum_arr.data(), fb_sum.data(), global_max_iex,
                  MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_count_arr.data(), fb_count.data(), global_max_iex,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    for (int iex = 1; iex <= global_max_iex; ++iex)
        if (use_fallback[iex]) {
            global_sum_arr[iex] = fb_sum[iex];
            global_count_arr[iex] = fb_count[iex];
        }

    // Initialize S_mean_temp and S_std_temp
    std::vector<float> S_mean_t(NMAX_E, fc::default_s_init);
    std::vector<float> S_std_t(NMAX_E, fc::default_s_std);
    std::vector<float> clip_limit(NMAX_E, fc::init_win_active);

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        if (global_count_arr[iex] > 1) {
            skip_iter[iex] = false;
            S_mean_t[iex] = global_sum_arr[iex] / float(global_count_arr[iex]);
            S_std_t[iex] = use_fallback[iex] ? fc::init_win_fallback : fc::init_win_active;
        } else {
            skip_iter[iex] = true;
            S_mean_t[iex] = S_init_arr[iex];
            S_std_t[iex] = fc::default_s_std;
        }
    }

    // Iterative sigma-clipping (3 iterations)
    for (int i_iter = 0; i_iter < 3; ++i_iter) {
        for (int iex = 1; iex <= global_max_iex; ++iex) {
            if (skip_iter[iex]) continue;
            float iw = use_fallback[iex] ? fc::init_win_fallback : fc::init_win_active;
            if (i_iter == 0) clip_limit[iex] = iw;
            else if (i_iter == 1) clip_limit[iex] = 0.5f * iw;
            else {
                clip_limit[iex] = fc::clip_nsigma * S_std_t[iex];
                if (clip_limit[iex] > 0.5f * iw) clip_limit[iex] = 0.5f * iw;
                if (clip_limit[iex] < fc::min_clip_limit) clip_limit[iex] = fc::min_clip_limit;
            }
        }

        // A. Compute mean
        std::fill(local_sum_arr.begin(), local_sum_arr.end(), 0.0);
        std::fill(local_count_arr.begin(), local_count_arr.end(), 0);
        for (int idx = 0; idx < data.ng; ++idx) {
            if (data.src_snr[idx] <= fc::stage1_snr) continue;
            int iex = data.iexpo[idx];
            if (iex < 1 || iex > global_max_iex || skip_iter[iex]) continue;
            int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm &&
                (use_fallback[iex] || active_mag_bins3d[(iex - 1) * nm + j]) &&
                std::fabs(data.sizerel[idx] - S_mean_t[iex]) < clip_limit[iex]) {
                local_sum_arr[iex] += data.sizerel[idx];
                local_count_arr[iex]++;
            }
        }
        MPI_Allreduce(local_sum_arr.data(), global_sum_arr.data(), global_max_iex,
                      MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local_count_arr.data(), global_count_arr.data(), global_max_iex,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for (int iex = 1; iex <= global_max_iex; ++iex)
            if (!skip_iter[iex]) {
                if (global_count_arr[iex] > 1)
                    S_mean_t[iex] = global_sum_arr[iex] / float(global_count_arr[iex]);
                else skip_iter[iex] = true;
            }

        // B. Compute std
        std::vector<float> local_sq(NMAX_E, 0.0), global_sq(NMAX_E, 0.0);
        std::vector<int> local_cnt_s(NMAX_E, 0), global_cnt_s(NMAX_E, 0);
        for (int idx = 0; idx < data.ng; ++idx) {
            if (data.src_snr[idx] <= fc::stage2_snr) continue;
            int iex = data.iexpo[idx];
            if (iex < 1 || iex > global_max_iex || skip_iter[iex]) continue;
            int j = int((data.magi[idx] - fc::mag_min_val) / mag_bin_w);
            if (j >= 0 && j < nm &&
                (use_fallback[iex] || active_mag_bins3d[(iex - 1) * nm + j]) &&
                std::fabs(data.sizerel[idx] - S_mean_t[iex]) < clip_limit[iex]) {
                local_sq[iex] += (data.sizerel[idx] - S_mean_t[iex]) * (data.sizerel[idx] - S_mean_t[iex]);
                local_cnt_s[iex]++;
            }
        }
        MPI_Allreduce(local_sq.data(), global_sq.data(), global_max_iex,
                      MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local_cnt_s.data(), global_cnt_s.data(), global_max_iex,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for (int iex = 1; iex <= global_max_iex; ++iex)
            if (!skip_iter[iex]) {
                if (global_cnt_s[iex] > 1)
                    S_std_t[iex] = std::sqrt(global_sq[iex] / float(global_cnt_s[iex] - 1));
                else S_std_t[iex] = fc::default_s_std;
                if (S_std_t[iex] < fc::min_clip_limit / fc::clip_nsigma)
                    S_std_t[iex] = fc::min_clip_limit / fc::clip_nsigma;
            }
    }

    // Finalize: average successful exposures for fallback
    int n_success = 0;
    float sum_mean = 0.0, sum_std = 0.0;
    for (int iex = 1; iex <= global_max_iex; ++iex)
        if (!skip_iter[iex]) {
            n_success++;
            sum_mean += S_mean_t[iex];
            sum_std += S_std_t[iex];
        }
    float avg_mean = n_success > 0 ? sum_mean / n_success : fc::default_s_init;
    float avg_std = n_success > 0 ? sum_std / n_success : fc::default_s_std;

    for (int iex = 1; iex <= global_max_iex; ++iex) {
        if (!skip_iter[iex]) {
            S_mean_arr[iex] = S_mean_t[iex];
            S_std_arr[iex] = S_std_t[iex];
        } else {
            S_mean_arr[iex] = avg_mean;
            S_std_arr[iex] = avg_std;
        }
        S_cut_arr[iex] = S_mean_arr[iex] + k_sigma * S_std_arr[iex];
    }
    for (int iex = global_max_iex + 1; iex < NMAX_E; ++iex) {
        S_mean_arr[iex] = avg_mean;
        S_std_arr[iex] = avg_std;
        S_cut_arr[iex] = avg_mean + k_sigma * avg_std;
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

// ------------------------------------------------------------------
// Apply advanced cuts (per-exposure star cut + SNR cuts)
// ------------------------------------------------------------------
void StarCutCalculator::applyAdvancedCuts(FDData& data,
                                          const std::vector<float>& S_cut_arr) {
    const int NMAX_E = LensingConfig::NMAX_EXPO;
    int write_idx = 0;
    float max_scut = 0.0;
    for (int i = 0; i < NMAX_E; ++i)
        if (S_cut_arr[i] > max_scut) max_scut = S_cut_arr[i];

    for (int idx = 0; idx < data.ng; ++idx) {
        // SNRF cut
        if (data.snrf[idx] < fc::snrfcut) continue;
        // SNR cuts
        if (fc::snrlow > 0.0 && data.src_snr[idx] < fc::snrlow) continue;
        if (fc::snrhigh > 0.0 && data.src_snr[idx] > fc::snrhigh) continue;

        // Per-exposure size cut
        int iex = data.iexpo[idx];
        float scut = max_scut;
        if (iex >= 1 && iex < NMAX_E) scut = S_cut_arr[iex];
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
