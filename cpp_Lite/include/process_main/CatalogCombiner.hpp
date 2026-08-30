#ifndef CATALOG_COMBINER_HPP
#define CATALOG_COMBINER_HPP

#include <vector>
#include <string>

namespace CatalogCombiner {
    void combineExpoCatalog(int nchip,
                            const std::vector<std::string>& imageFiles,
                            const std::string& dirOutput,
                            int expo_index,
                            float chi2);
    void procComb(int iexpo);
}

#endif // CATALOG_COMBINER_HPP
