#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <string>

// ==========================================
// Function: Verify MPI rank identity and placement
// Method: Initialize MPI, reject singleton execution, and print each rank and processor
// ==========================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int world_size = 0;
    int name_length = 0;
    char processor_name[MPI_MAX_PROCESSOR_NAME] = {};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Get_processor_name(processor_name, &name_length);

    if (world_size < 2) {
        MPI_Abort(MPI_COMM_WORLD, 3);
    }

    std::cout << "rank " << rank << " of " << world_size << " on "
              << std::string(processor_name, name_length) << std::endl;

    MPI_Finalize();
    return EXIT_SUCCESS;
}
