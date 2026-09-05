#include "process_main/FitsIO.hpp"
#include "process_main/PSFRecons.hpp"
#include "process_main/Universalblock.hpp"

#include <fitsio.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

namespace {

// ==========================================
// Function: Stop the lookup regression when one requirement is not met
// Method: Print a focused failure message and return a non-zero process code.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFRecons CCD lookup test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Class: Own synthetic sequentially named Science FITS inputs
// Method: Create an isolated temporary directory and remove it on scope exit.
// ==========================================
class TemporaryScienceImages {
public:
    TemporaryScienceImages()
        : root_(std::filesystem::temp_directory_path()
                / ("fq_psfrecons_ccd_" + std::to_string(::getpid()))) {
        std::filesystem::create_directories(root_);
    }

    ~TemporaryScienceImages() {
        std::error_code error;
        std::filesystem::remove_all(root_, error);
    }

    // ==========================================
    // Function: Create one sequentially named Science FITS
    // Method: Write a minimal image and optionally attach its physical CCDNUM.
    // ==========================================
    std::string add(int sequence_number, int ccd_num,
                    bool write_ccdnum = true) const {
        const std::string filename =
            (root_ / ("exposure_" + std::to_string(sequence_number)
                      + ".fits")).string();
        require(FitsIO::writeImage(filename, 2, 2,
                                   {0.0f, 0.0f, 0.0f, 0.0f}),
                "cannot write synthetic Science FITS");
        if (!write_ccdnum) {
            return filename;
        }

        fitsfile* file = nullptr;
        int status = 0;
        fits_open_file(&file, filename.c_str(), READWRITE, &status);
        if (status == 0) {
            fits_update_key(file, TINT, "CCDNUM", &ccd_num, nullptr, &status);
        }
        int close_status = 0;
        if (file != nullptr) {
            fits_close_file(file, &close_status);
        }
        require(status == 0 && close_status == 0,
                "cannot write synthetic CCDNUM");
        return filename;
    }

    // ==========================================
    // Function: Write the Norm associated with one continuous Science basename
    // Method: Derive the production path and store a selectable first-pixel sentinel.
    // ==========================================
    void writeNorm(const std::string& science_file, float sentinel) const {
        const std::string filename = Universalblock::normFilename(
            science_file, root_.string());
        std::filesystem::create_directories(
            std::filesystem::path(filename).parent_path());
        require(FitsIO::writeImage(filename, 1, 1, {sentinel}),
                "cannot write synthetic Norm FITS");
    }

    // ==========================================
    // Function: Return the synthetic pipeline output root
    // Method: Expose the root used by production Norm path derivation.
    // ==========================================
    std::string outputRoot() const {
        return root_.string();
    }

private:
    std::filesystem::path root_;
};

// ==========================================
// Function: Observe a fatal lookup-index build in an isolated process
// Method: Fork before invoking the production helper and inspect child status.
// ==========================================
bool buildFails(const std::vector<std::string>& image_files, int max_chip_id) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        (void)PSFRecons::Internal::buildChipImageIndex(
            image_files, max_chip_id);
        std::_Exit(EXIT_SUCCESS);
    }
    int status = 0;
    require(::waitpid(child, &status, 0) == child, "waitpid failed");
    return !WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS;
}

// ==========================================
// Function: Verify physical CCDNUM maps to the sequential Science basename
// Method: Build a reordered list with CCD 2 absent and inspect every index slot.
// ==========================================
void testGapAndReorderedList(const TemporaryScienceImages& images) {
    const std::string sequence1 = images.add(1, 1);
    const std::string sequence2 = images.add(2, 3);
    const std::string sequence3 = images.add(3, 4);
    const std::vector<std::string> image_files = {
        sequence3, sequence1, sequence2
    };

    const std::vector<int> index =
        PSFRecons::Internal::buildChipImageIndex(image_files, 4);
    require(index.size() == 4, "mapping has the wrong CCD range");
    require(index[0] == 1, "CCD 1 did not map to sequence 1");
    require(index[1] == -1, "missing CCD 2 must remain absent");
    require(index[2] == 2, "CCD 3 did not map to sequence 2");
    require(index[3] == 0, "CCD 4 did not map to sequence 3");
    const std::string* ccd3 = PSFRecons::Internal::indexedChipImage(
        image_files, index, 0, 3, 4);
    require(ccd3 != nullptr && *ccd3 == sequence2,
            "CCD 3 resolved the wrong continuous-basename Science FITS");
    images.writeNorm(sequence2, -1.0f);
    images.writeNorm(sequence3, 1.0f);
    require(Universalblock::checkNorm(*ccd3, images.outputRoot())
                == Universalblock::NormStatus::Valid,
            "CCD 3 did not validate the Norm for sequence 2");
    require(PSFRecons::Internal::indexedChipImage(
                image_files, index, 0, 2, 4) == nullptr,
            "missing CCD 2 must resolve as absent");

    const std::vector<std::string> second_image_files = {
        sequence2, sequence3
    };
    const std::vector<int> second_index =
        PSFRecons::Internal::buildChipImageIndex(second_image_files, 4);
    std::vector<int> flattened_indices = index;
    flattened_indices.insert(flattened_indices.end(),
                             second_index.begin(), second_index.end());
    const std::string* exposure2_ccd4 =
        PSFRecons::Internal::indexedChipImage(
            second_image_files, flattened_indices, 1, 4, 4);
    require(exposure2_ccd4 != nullptr && *exposure2_ccd4 == sequence3,
            "second exposure used the wrong flattened mapping offset");
}

// ==========================================
// Function: Verify malformed CCDNUM mappings fail before PCA scheduling
// Method: Exercise missing, duplicate, and out-of-range physical identities.
// ==========================================
void testMalformedMappings(const TemporaryScienceImages& images) {
    const std::string missing_header = images.add(10, 0, false);
    require(buildFails({missing_header}, 62),
            "missing CCDNUM must fail");

    const std::string duplicate1 = images.add(11, 3);
    const std::string duplicate2 = images.add(12, 3);
    require(buildFails({duplicate1, duplicate2}, 62),
            "duplicate CCDNUM must fail");

    const std::string out_of_range = images.add(13, 63);
    require(buildFails({out_of_range}, 62),
            "out-of-range CCDNUM must fail");
}

}  // namespace

// ==========================================
// Function: Run physical-CCD to sequential-Science lookup regressions
// Method: Exercise valid gaps/reordering and malformed mapping contracts.
// ==========================================
int main() {
    TemporaryScienceImages images;
    testGapAndReorderedList(images);
    testMalformedMappings(images);
    std::cout << "PSFRecons CCD lookup tests passed\n";
    return EXIT_SUCCESS;
}
