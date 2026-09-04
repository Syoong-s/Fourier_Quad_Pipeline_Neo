#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace Initialize {

// Instrument / CCD
inline constexpr int N_CCD = 62;
inline constexpr double pixel_size = 0.2628;
inline constexpr int chipnx = 2046;
inline constexpr int chipny = 4094;
inline constexpr double saturation_thresh = 25000.0;
inline constexpr int CCD_split = 2;

// Pipeline stages
inline constexpr int PROCESS_stage =
    2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23;

// Workflow
inline constexpr bool RUN_PROCESS_ASTROCAT = false;
inline constexpr bool RUN_PROCESS_EXTCAT = false;
inline constexpr bool RUN_PROCESS_INIT = true;
inline constexpr bool RUN_PROCESS_MAIN = true;
inline constexpr bool RUN_PROCESS_REARR = false;
inline constexpr bool RUN_PROCESS_FD = false;

// ==========================================
// Structure: Describe one input dataset
// Method: Pair the output target directory with its archive basename prefix.
// ==========================================
struct DatasetSpec {
    std::string target;
    std::string prefix;
};

inline const std::vector<DatasetSpec> DATASETS = {
    {"gband", "c4d_"}
};

inline const std::vector<std::string> CONTAINS = {
    "v1"
};

inline constexpr const char* EXPO_LIST = "";

// Main paths
inline constexpr const char* ASTROMETRY_CAT =
    "/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted";
inline constexpr const char* SOURCE_CAT_DEFAULT =
    "/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat";
inline constexpr const char* SCIENCE_ROOT =
    "/lustre/home/acct-phyzj/share/DES/g";
inline constexpr const char* DQ_ROOT =
    "/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask";
inline constexpr const char* OUTPUT_ROOT =
    "/lustre/home/acct-phyzj/share/DES/g_band_v1";

// Astrometry catalog
inline constexpr const char* ASTROMETRY_TILE_PREFIX = "astra_";
inline constexpr const char* ASTROCAT_INPUT_DIRECTORY = "";
inline constexpr bool ASTROCAT_ADD_HEADER = true;

// External source catalog
inline constexpr const char* SOURCE_CAT_TILE_PREFIX = "extern_";
inline constexpr const char* EXTCAT_INPUT_DIRECTORY = "";
inline constexpr std::size_t EXTCAT_TOTAL_COLUMNS = 18;
inline constexpr std::size_t EXTCAT_RA_COLUMN_ONE_BASED = 5;
inline constexpr std::size_t EXTCAT_DEC_COLUMN_ONE_BASED = 6;
inline constexpr std::size_t EXTCAT_ZP_COLUMN_ONE_BASED = 17;

// FITS / archive conventions
inline constexpr const char* ARCHIVE_SUFFIX = ".fits.fz";
inline constexpr const char* CCDNUM_KEYWORD = "CCDNUM";
inline constexpr const char* DQ_STEM_REPLACE_FROM = "ood";
inline constexpr const char* DQ_STEM_REPLACE_TO = "ooi";

}  // namespace Initialize

#endif  // INITIALIZE_HPP
