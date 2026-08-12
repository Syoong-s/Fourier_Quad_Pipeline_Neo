#ifndef SHEAR_CATALOG_READER_HPP
#define SHEAR_CATALOG_READER_HPP

#include "CatalogLayout.hpp"
#include "FDData.hpp"

#include <cstddef>
#include <string>
#include <vector>

// ==========================================
// ShearCatalogReader - reads one exposure's _all.cat shear catalog
// Method: Parse whitespace-separated rows, apply quality cuts, deduplicate
//         overlapping detections, and append valid sources to the shared
//         FDData arrays.  Equivalent to Fortran read_shear_cat_v2.
// ==========================================
class ShearCatalogReader {
public:
    // ==========================================
    // Function: Parse one exact-width numeric FD catalog row
    // Method: Require the runtime layout width, finite values, and no trailing columns.
    // ==========================================
    static bool parseCatalogRow(const std::string& line,
                                std::size_t column_count,
                                std::vector<float>& values);

    // ==========================================
    // Function: Read one exposure catalog and append accepted sources
    // Method: Use the shared layout for row width/fields and one FD-selected
    //         magnitude column for both range filtering and star-bar storage.
    // ==========================================
    static void readExposure(int iexpo, FDData& data,
                             const std::vector<std::string>& expo_files,
                             const PipelineCatalog::CatalogLayout& layout,
                             std::size_t magnitude_column,
                             int rank);
};

#endif  // SHEAR_CATALOG_READER_HPP
