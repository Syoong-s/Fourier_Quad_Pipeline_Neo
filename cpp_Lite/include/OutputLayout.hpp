#ifndef OUTPUT_LAYOUT_HPP
#define OUTPUT_LAYOUT_HPP

#include <array>
#include <stdexcept>
#include <string>

namespace OutputLayout {

// Base directories that never receive process_main chip products. Together
// with CHIP_PRODUCT_DIRECTORIES, this is the complete process_init directory
// contract. std::filesystem::create_directories makes creation idempotent.
inline constexpr std::array<const char*, 14> NON_CHIP_BASE_DIRECTORIES = {
    "science",
    "dqmask",
    "stamps",
    "result",
    "stamps/dat_StarInfo",
    "stamps/fits_StarP",
    "stamps/fits_PsfSrc",
    "stamps/dat_ExpoInfo",
    "stamps/dat_StarComp",
    "stamps/dat_Rescale",
    "stamps/dat_Pcs",
    "stamps/dat_StarCompV2",
    "astrometry/Head",
    "astrometry/dat_Chk",
};

// Union of every Standard/Lite process_main directory whose products have one
// file per chip. Lite intentionally reserves the Standard-only PSF directories
// so both variants expose one stable initializer layout. After publishing
// expo_<target>.list, process_init uses the same list to materialize one
// exposure subdirectory beneath each product directory.
inline constexpr std::array<const char*, 16> CHIP_PRODUCT_DIRECTORIES = {
    "stamps/Norm",
    "stamps/cat_Orig",
    "stamps/dat_StarCanInfo",
    "stamps/fits_StarCan",
    "stamps/fits_StarCanN",
    "stamps/fits_StarCanP",
    "stamps/dat_SrcInfo",
    "stamps/fits_Src",
    "stamps/fits_Noise",
    "stamps/fits_SrcP",
    "stamps/dat_PsfFit",
    "stamps/fits_PsfLocal",
    "stamps/dat_Shear",
    "stamps/dat_StarXY",
    "stamps/fits_PsfResi",
    "astrometry/dat_Astro",
};

// Extract the exposure component from the established <exposure>_<chip>
// prefix contract used by process_init and UniversalUtils::getPrefix.
inline std::string exposureFromChipPrefix(const std::string& chip_prefix) {
    const std::size_t delimiter = chip_prefix.find_last_of('_');
    if (delimiter == std::string::npos || delimiter == 0
        || delimiter + 1 >= chip_prefix.size()) {
        throw std::invalid_argument("invalid chip prefix: " + chip_prefix);
    }
    return chip_prefix.substr(0, delimiter);
}

// Build <output-root>/<product-directory>/<exposure>/<chip><suffix>.
inline std::string chipPath(const std::string& output_root,
                            const std::string& product_directory,
                            const std::string& chip_prefix,
                            const std::string& suffix) {
    return output_root + "/" + product_directory + "/"
           + exposureFromChipPrefix(chip_prefix) + "/"
           + chip_prefix + suffix;
}

}  // namespace OutputLayout

#endif  // OUTPUT_LAYOUT_HPP
