#include <Eigen/Dense>
#include <fitsio.h>
#include <fftw3.h>
#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

extern "C" {
void dgesv_(const int* n, const int* nrhs, double* matrix, const int* leading_dim,
            int* pivots, double* rhs, const int* rhs_leading_dim, int* info);
}

// ==========================================
// Function: Validate the C++ MPI and scientific-library stack
// Method: Exercise Eigen, LAPACK, FFTW, and CFITSIO from every initialized MPI rank
// ==========================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Eigen::Matrix2d eigen_matrix;
    eigen_matrix << 2.0, 0.0, 0.0, 4.0;
    const Eigen::Vector2d eigen_rhs(4.0, 8.0);
    const Eigen::Vector2d eigen_solution =
        eigen_matrix.colPivHouseholderQr().solve(eigen_rhs);

    int n = 1;
    int nrhs = 1;
    int leading_dim = 1;
    int rhs_leading_dim = 1;
    int pivot = 0;
    int info = 0;
    double lapack_matrix = 2.0;
    double lapack_rhs = 4.0;
    dgesv_(&n, &nrhs, &lapack_matrix, &leading_dim, &pivot, &lapack_rhs,
           &rhs_leading_dim, &info);

    double fft_input[4] = {1.0, 0.0, 0.0, 0.0};
    fftw_complex fft_output[3] = {};
    fftw_plan plan = fftw_plan_dft_r2c_1d(4, fft_input, fft_output, FFTW_ESTIMATE);
    if (plan == nullptr) {
        MPI_Abort(MPI_COMM_WORLD, 4);
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    float cfitsio_version = 0.0F;
    fits_get_version(&cfitsio_version);

    const bool eigen_ok = (eigen_solution.array() - 2.0).abs().maxCoeff() < 1.0e-12;
    const bool lapack_ok = info == 0 && std::abs(lapack_rhs - 2.0) < 1.0e-12;
    const bool fftw_ok = std::abs(fft_output[0][0] - 1.0) < 1.0e-12;
    const bool cfitsio_ok = cfitsio_version >= 4.0F;

    if (!(eigen_ok && lapack_ok && fftw_ok && cfitsio_ok)) {
        MPI_Abort(MPI_COMM_WORLD, 5);
    }

    if (rank == 0) {
        std::cout << "science stack passed: CFITSIO " << cfitsio_version
                  << ", Eigen/LAPACK/FFTW OK" << std::endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
