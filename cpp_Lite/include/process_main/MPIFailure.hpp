#ifndef MPI_FAILURE_HPP
#define MPI_FAILURE_HPP

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <string>

namespace MPIFailure {

// ==========================================
// Function: Abort the complete MPI world after an existing fatal pipeline error
// Method: Report rank, operation, and detail when MPI is active, then call
//         MPI_Abort so peer ranks cannot remain blocked in schedulers or collectives.
// ==========================================
[[noreturn]] inline void abortWorld(const std::string& operation,
                                    const std::string& detail,
                                    int error_code = EXIT_FAILURE) {
    int mpi_initialized = 0;
    int mpi_finalized = 0;
    int rank = -1;

    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized != 0) {
        MPI_Finalized(&mpi_finalized);
    }
    if (mpi_initialized != 0 && mpi_finalized == 0) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    std::cerr << "Fatal pipeline error"
              << " [rank=" << rank << "]"
              << " operation=" << operation
              << " detail=" << detail
              << std::endl;

    if (mpi_initialized != 0 && mpi_finalized == 0) {
        MPI_Abort(MPI_COMM_WORLD, error_code);
    }
    std::_Exit(error_code);
}

}  // namespace MPIFailure

#endif  // MPI_FAILURE_HPP
