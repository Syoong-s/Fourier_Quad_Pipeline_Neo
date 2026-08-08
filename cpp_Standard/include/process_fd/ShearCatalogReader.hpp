#ifndef SHEAR_CATALOG_READER_HPP
#define SHEAR_CATALOG_READER_HPP

#include "FDData.hpp"

#include <string>

// ==========================================
// ShearCatalogReader - reads one exposure's _all.cat shear catalog
// Method: Parse whitespace-separated rows, apply quality cuts, deduplicate
//         overlapping detections, and append valid sources to the shared
//         FDData arrays.  Equivalent to Fortran read_shear_cat_v2.
// ==========================================
class ShearCatalogReader {
public:
    // Read one exposure (1-based index into EXPO_FILE) and append to data.
    // The caller is responsible for MPI distribution of exposures.
    static void readExposure(int iexpo, FDData& data,
                             const std::vector<std::string>& expo_files,
                             int rank);
};

#endif  // SHEAR_CATALOG_READER_HPP
