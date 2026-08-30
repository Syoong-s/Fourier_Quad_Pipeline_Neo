#include "general/MPIScheduler.hpp"
#include <iostream>
#include <algorithm>

namespace MPIScheduler {
    State state;

    void init(int& argc, char**& argv) {
        MPI_Init(&argc, &argv);
        state.communicator = MPI_COMM_WORLD;
        MPI_Comm_rank(state.communicator, &state.rank);
        MPI_Comm_size(state.communicator, &state.size);
    }

    void finalize() {
        MPI_Finalize();
    }

    void barrier() {
        MPI_Barrier(state.communicator);
    }

    void distribute(int num_jobs, const std::function<void(int)>& job_func, const std::string& message) {
        if (state.size <= 1) {
            for (int i = 1; i <= num_jobs; ++i) {
                job_func(i);
            }
            return;
        }

        MPI_Barrier(state.communicator);
        int complete = 0;
        int i = 0;
        int j = 0;
        if (state.rank == 0) {
            i = 1;
            j = num_jobs;
        }

        while (complete == 0) {
            if (state.rank != 0) {
                // Workers request a job by sending their status (initially 0, later the job index they completed)
                MPI_Send(&i, 1, MPI_INT, 0, 0, state.communicator);
                MPI_Recv(&i, 1, MPI_INT, 0, 0, state.communicator, MPI_STATUS_IGNORE);
                if (i == 0) {
                    complete = 1;
                } else {
                    job_func(i);
                }
            } else {
                // Master coordinates job distribution
                int k = 0;
                MPI_Status status;
                MPI_Recv(&k, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, state.communicator, &status);
                int source = status.MPI_SOURCE;
                int tag = status.MPI_TAG;

                // Send current job index i to worker
                MPI_Send(&i, 1, MPI_INT, source, tag, state.communicator);
                if (k > 0) {
                    j--;
                }

                if (i != 0) {
                    i++;
                    if (i % 100 == 1) {
                        std::cout << "Master: Process " << source << " assigned job " << i
                                  << ", remaining jobs to complete: " << j << " (" << message << ")" << std::endl;
                        std::cout.flush();
                    }
                } else {
                    // Send termination signal to remaining workers requesting jobs after queue is empty
                    // Do not increment i, keep it as 0
                }

                if (i > num_jobs) {
                    i = 0;
                }

                if (j == 0) {
                    complete = 1;
                    std::cout << "Master: Process " << source << " assigned job " << i
                              << ", remaining jobs to complete: " << j << " (" << message << ")" << std::endl;
                    std::cout.flush();
                }
            }
        }
        MPI_Barrier(state.communicator);
    }
}
