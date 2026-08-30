#ifndef MPI_SCHEDULER_HPP
#define MPI_SCHEDULER_HPP

#include <mpi.h>

#include <string>
#include <functional>

namespace MPIScheduler {
    struct State {
        MPI_Comm communicator = MPI_COMM_NULL;
        int rank = 0;
        int size = 1;
    };

    extern State state;

    void init(int& argc, char**& argv);
    void finalize();
    void barrier();

    // Standard MPI master-worker dynamic distribution (1-based job indexes)
    void distribute(int num_jobs, const std::function<void(int)>& job_func, const std::string& message);

    // Dynamic distribution with rank throttling per node (for PCA/Covariance fitting)
    void forcecov(int ppn, int work_pn, int num_jobs, const std::function<void(int, int)>& job_func, const std::string& message, int nexpo);
}

#endif // MPI_SCHEDULER_HPP
