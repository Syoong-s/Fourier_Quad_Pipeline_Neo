#ifndef KMEANS_CLUSTERER_HPP
#define KMEANS_CLUSTERER_HPP

#include <vector>

// ==========================================
// KMeansClusterer - spherical k-means for jackknife region assignment
// Method: Randomly select sources from all MPI nodes, run spherical k-means
//         (maximizing dot product) on rank 0, then broadcast the cluster
//         centers.  Equivalent to Fortran simple_kmeans_MPI + simple_kmeans.
// ==========================================
class KMeansClusterer {
public:
    // Run distributed k-means and return N_jack cluster centers (3D unit vectors)
    static void runMPI(int ng, const std::vector<float>& ra,
                       const std::vector<float>& dec, int rank,
                       int num_procs,
                       std::vector<float>& centers /* 3×N_jack */);
};

#endif  // KMEANS_CLUSTERER_HPP
