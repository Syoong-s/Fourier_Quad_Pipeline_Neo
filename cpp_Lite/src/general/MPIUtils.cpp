#include "general/MPIUtils.hpp"

#include "general/MPIScheduler.hpp"

#include <limits>

namespace MPIUtils {

bool broadcastString(std::string& value,
                     int root,
                     std::string& error) {
    const int rank = MPIScheduler::state.rank;
    const MPI_Comm communicator = MPIScheduler::state.communicator;

    int length = 0;
    if (rank == root) {
        length = value.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
                     ? -1
                     : static_cast<int>(value.size());
    }
    if (MPI_Bcast(&length, 1, MPI_INT, root, communicator) != MPI_SUCCESS) {
        error = "failed to broadcast string length";
        return false;
    }
    if (length < 0) {
        error = "string length exceeds the MPI int range";
        return false;
    }
    if (rank != root) {
        value.resize(static_cast<std::size_t>(length));
    }
    if (length > 0
        && MPI_Bcast(value.data(), length, MPI_CHAR, root, communicator)
               != MPI_SUCCESS) {
        error = "failed to broadcast string payload";
        return false;
    }
    error.clear();
    return true;
}

bool broadcastStrings(std::vector<std::string>& values,
                      int root,
                      std::string& error) {
    const int rank = MPIScheduler::state.rank;
    const MPI_Comm communicator = MPIScheduler::state.communicator;

    int count = 0;
    if (rank == root) {
        count = values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
                    ? -1
                    : static_cast<int>(values.size());
    }
    if (MPI_Bcast(&count, 1, MPI_INT, root, communicator) != MPI_SUCCESS) {
        error = "failed to broadcast string-vector count";
        return false;
    }
    if (count < 0) {
        error = "string-vector count exceeds the MPI int range";
        return false;
    }
    if (rank != root) {
        values.resize(static_cast<std::size_t>(count));
    }
    for (std::string& value : values) {
        if (!broadcastString(value, root, error)) {
            return false;
        }
    }
    error.clear();
    return true;
}

bool allRanksSucceeded(bool local_success,
                       bool& global_success,
                       std::string& error) {
    const int local_value = local_success ? 1 : 0;
    int global_value = 0;
    if (MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_MIN,
                      MPIScheduler::state.communicator) != MPI_SUCCESS) {
        error = "failed to reduce MPI success flags";
        global_success = false;
        return false;
    }
    global_success = global_value != 0;
    error.clear();
    return true;
}

}  // namespace MPIUtils
