#ifndef MPI_SCHEDULER_HPP
#define MPI_SCHEDULER_HPP

#include <string>
#include <functional>

namespace MPIScheduler {
    extern int my_id;
    extern int num_procs;

    void init(int& argc, char**& argv);
    void finalize();
    void barrier();

    // Standard MPI master-worker dynamic distribution (1-based job indexes)
    void distribute(int num_jobs, const std::function<void(int)>& job_func, const std::string& message);

    // cpp_lite: forcecov() (rank-throttled distribution for the PCA/covariance fit) was only
    // used by the PSF_Ms=1 reconstruction stage and has been removed with it.
}

#endif // MPI_SCHEDULER_HPP
