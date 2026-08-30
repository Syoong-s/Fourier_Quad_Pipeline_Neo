#ifndef GENERAL_MPI_UTILS_HPP
#define GENERAL_MPI_UTILS_HPP

#include <string>
#include <vector>

namespace MPIUtils {

// ==========================================
// Function: Broadcast one dynamically sized string
// Method: Send an int-safe byte length followed by the exact payload.
// ==========================================
bool broadcastString(std::string& value,
                     int root,
                     std::string& error);

// ==========================================
// Function: Broadcast an ordered vector of strings
// Method: Send an int-safe item count and reuse the length-prefixed string
//         transport for every element.
// ==========================================
bool broadcastStrings(std::vector<std::string>& values,
                      int root,
                      std::string& error);

// ==========================================
// Function: Combine rank-local success flags
// Method: Reduce integer flags with MPI_MIN so every rank receives the same
//         collective outcome.
// ==========================================
bool allRanksSucceeded(bool local_success,
                       bool& global_success,
                       std::string& error);

}  // namespace MPIUtils

#endif  // GENERAL_MPI_UTILS_HPP
