#ifndef OUTPUT_FILE_HPP
#define OUTPUT_FILE_HPP

#include <mpi.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace MainIO {

// ==========================================
// Function: Terminate the complete MPI program after an output-file failure
// Method: Print rank, operation, path, and root cause, then abort MPI so peer
//         ranks cannot remain blocked in schedulers or collectives.
// ==========================================
[[noreturn]] inline void failOutput(const std::string& operation,
                                    const std::string& path,
                                    const std::string& reason) {
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

    std::cerr << "Output creation failed"
              << " [rank=" << rank << "]"
              << " operation=" << operation
              << " path=\"" << path << "\""
              << " reason=" << reason << std::endl;

    if (mpi_initialized != 0 && mpi_finalized == 0) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    std::_Exit(EXIT_FAILURE);
}

// ==========================================
// Function: Describe one failed standard-library file operation
// Method: Preserve errno captured immediately after open when available and
//         otherwise return a deterministic stream-state explanation.
// ==========================================
inline std::string systemFailureReason(int error_number) {
    if (error_number != 0) {
        return std::strerror(error_number);
    }
    return "standard output stream entered a failure state";
}

// ==========================================
// Class: Fail-fast text output stream
// Method: Validate file creation on open and validate buffered writes on every
//         explicit or scope-driven close before allowing pipeline execution to continue.
// ==========================================
class OutputFile final : public std::ofstream {
public:
    OutputFile() = default;

    explicit OutputFile(const std::string& filename,
                        std::ios_base::openmode mode = std::ios_base::out) {
        open(filename, mode);
    }

    ~OutputFile() {
        if (is_open()) {
            close();
        }
    }

    // ==========================================
    // Function: Open one checked text output file
    // Method: Capture errno adjacent to std::ofstream::open and abort all MPI
    //         ranks immediately when the path cannot be created.
    // ==========================================
    void open(const std::string& filename,
              std::ios_base::openmode mode = std::ios_base::out) {
        path_ = filename;
        errno = 0;
        std::ofstream::open(filename, mode);
        const int open_error = errno;
        if (!is_open()) {
            failOutput("open text output", path_, systemFailureReason(open_error));
        }
    }

    // ==========================================
    // Function: Close one checked text output file
    // Method: Flush and test the stream before and after close so delayed write
    //         and filesystem-close failures terminate the complete MPI program.
    // ==========================================
    void close() {
        if (!is_open()) {
            return;
        }
        errno = 0;
        flush();
        const int flush_error = errno;
        if (!good()) {
            failOutput("write text output", path_, systemFailureReason(flush_error));
        }

        errno = 0;
        std::ofstream::close();
        const int close_error = errno;
        if (fail()) {
            failOutput("close text output", path_, systemFailureReason(close_error));
        }
    }

private:
    std::string path_;
};

}  // namespace MainIO

#endif  // OUTPUT_FILE_HPP
