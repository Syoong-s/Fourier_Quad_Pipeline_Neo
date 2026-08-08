#include "process_fd/KMeansClusterer.hpp"
#include "FDConfig.hpp"
#include "process_main/NumericalRecipes.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace fc = FDConfig;

namespace {

// Reservoir sampling: select n_select elements from ra/dec
void randomSelect(int ng, int n_select,
                  const std::vector<float>& ra,
                  const std::vector<float>& dec,
                  std::vector<float>& ras,
                  std::vector<float>& decs) {
    ras.assign(n_select, 0.0);
    decs.assign(n_select, 0.0);
    int i_current = 0;
    for (int i = 0; i < ng; ++i) {
        i_current = i + 1;
        if (i < n_select) {
            ras[i] = ra[i];
            decs[i] = dec[i];
        } else {
            float rand_val = static_cast<float>(NumericalRecipes::ran1());
            if (rand_val < float(n_select) / float(i_current)) {
                int replace = static_cast<int>(NumericalRecipes::ran1() * n_select);
                if (replace >= n_select) replace = n_select - 1;
                ras[replace] = ra[i];
                decs[replace] = dec[i];
            }
        }
    }
}

// Spherical k-means (maximize dot product)
void simpleKmeans(int n_gal, const std::vector<float>& pos,
                  int n_k, int max_iter,
                  std::vector<float>& centers /* 3×n_k */) {
    for (int k = 0; k < n_k; ++k) {
        int idx = static_cast<int>(NumericalRecipes::ran1() * n_gal);
        if (idx >= n_gal) idx = n_gal - 1;
        if (idx < 0) idx = 0;
        for (int d = 0; d < 3; ++d)
            centers[d * n_k + k] = pos[d * n_gal + idx];
    }

    std::vector<float> new_centers(3 * n_k, 0.0);
    std::vector<int> counts(n_k, 0);

    for (int iter = 0; iter < max_iter; ++iter) {
        std::fill(new_centers.begin(), new_centers.end(), 0.0);
        std::fill(counts.begin(), counts.end(), 0);

        for (int i = 0; i < n_gal; ++i) {
            float max_dot = -2.0;
            int best_k = 0;
            for (int k = 0; k < n_k; ++k) {
                float dot = pos[0 * n_gal + i] * centers[0 * n_k + k]
                          + pos[1 * n_gal + i] * centers[1 * n_k + k]
                          + pos[2 * n_gal + i] * centers[2 * n_k + k];
                if (dot > max_dot) { max_dot = dot; best_k = k; }
            }
            for (int d = 0; d < 3; ++d)
                new_centers[d * n_k + best_k] += pos[d * n_gal + i];
            counts[best_k]++;
        }

        for (int k = 0; k < n_k; ++k) {
            if (counts[k] > 0) {
                float norm = 0.0;
                for (int d = 0; d < 3; ++d)
                    norm += new_centers[d * n_k + k] * new_centers[d * n_k + k];
                norm = std::sqrt(norm);
                if (norm > 0.0) {
                    for (int d = 0; d < 3; ++d)
                        centers[d * n_k + k] = new_centers[d * n_k + k] / norm;
                }
            } else {
                int idx = static_cast<int>(NumericalRecipes::ran1() * n_gal);
                if (idx >= n_gal) idx = n_gal - 1;
                if (idx < 0) idx = 0;
                for (int d = 0; d < 3; ++d)
                    centers[d * n_k + k] = pos[d * n_gal + idx];
            }
        }
    }
}

}  // namespace

// ==========================================
// Function: runMPI
// Method: Gather random position samples from all ranks, run spherical
//         k-means on rank 0, broadcast the resulting cluster centers.
// ==========================================
void KMeansClusterer::runMPI(int ng, const std::vector<float>& ra,
                            const std::vector<float>& dec, int rank,
                            int num_procs,
                            std::vector<float>& centers) {
    centers.assign(3 * fc::N_jack, 0.0);

    // Gather per-rank source counts
    std::vector<int> ng_all(num_procs, 0);
    MPI_Allgather(&ng, 1, MPI_INT, ng_all.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Determine sample size (min across non-zero ranks)
    int min_ng = ng_all[0];
    for (int p = 1; p < num_procs; ++p)
        if (ng_all[p] > 0 && (min_ng == 0 || ng_all[p] < min_ng))
            min_ng = ng_all[p];
    int n_select = std::min(min_ng, fc::nmax_total / std::max(1, num_procs - 1));
    if (n_select < 1) n_select = 1;
    int n_tot = n_select * num_procs;

    std::vector<float> ra_select(n_select, 0.0), dec_select(n_select, 0.0);
    std::vector<float> ra_tot, dec_tot;

    if (rank != 0 && ng > 0)
        randomSelect(ng, n_select, ra, dec, ra_select, dec_select);

    if (rank == 0) {
        ra_tot.resize(n_tot, 0.0);
        dec_tot.resize(n_tot, 0.0);
    }

    MPI_Gather(ra_select.data(), n_select, MPI_FLOAT,
               ra_tot.data(), n_select, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(dec_select.data(), n_select, MPI_FLOAT,
               dec_tot.data(), n_select, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Convert to 3D unit vectors (skip rank-0 samples which are zeros)
        int n_end = n_tot - n_select;
        std::vector<float> pos(3 * n_end);
        for (int i = 0; i < n_end; ++i) {
            int src = i + n_select;
            float ra_rad = ra_tot[src] * LensingConfig::pi / 180.0;
            float dec_rad = dec_tot[src] * LensingConfig::pi / 180.0;
            pos[0 * n_end + i] = std::cos(dec_rad) * std::cos(ra_rad);
            pos[1 * n_end + i] = std::cos(dec_rad) * std::sin(ra_rad);
            pos[2 * n_end + i] = std::sin(dec_rad);
        }
        simpleKmeans(n_end, pos, fc::N_jack, fc::Km_iter, centers);
        std::cout << "simple_kmeans_MPI done!" << std::endl;
    }

    MPI_Bcast(centers.data(), 3 * fc::N_jack, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}
