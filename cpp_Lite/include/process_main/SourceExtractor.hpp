#ifndef SOURCE_EXTRACTOR_HPP
#define SOURCE_EXTRACTOR_HPP

#include <string>
#include <vector>
#include "LensingConfig.hpp"
#include "FitsIO.hpp"

namespace SourceExtractor {
    void procSource(int iexpo);
    void chipProcessSource(const std::vector<std::string>& imageFiles, int ichip, const std::string& dirOutput);

    void deBlending(const std::vector<std::string>& sortFile, int sortNum, int nx, int ny, std::vector<int>& weight,
                    const double cRPIX[2], const double cD[2][2], const double cRVAL[2],
                    const double PU[2][LensingConfig::npd], int& procError);
                    
    void fillPatch(int nx, int ny, std::vector<int>& map, int ix, int iy, int old_v, int new_v);
    
    void getExpoCatalog(const std::string& dirOutput, const std::string& prefix, int nx, int ny, const std::vector<float>& sigmap,
                        std::vector<int>& weight, const std::vector<float>& normap, int& ierror);
                        
    void genSourceExtCatalog(const std::string& dirOutput, const std::vector<std::string>& sortFile, int sortNum, const std::string& prefix,
                             int nx, int ny, const std::vector<float>& array, std::vector<int>& weight,
                             const std::vector<float>& sigmap, const double cRPIX[2], const double cD[2][2],
                             const double cRVAL[2], const double PU[2][LensingConfig::npd], int& ngal, int& procError);
                             
    void checkSourceAndEstimateNoisePower(
        int& flag, std::vector<float>& sourceStamp, std::vector<float>& noisePower,
        int nx, int ny, const std::vector<float>& array, const std::vector<int>& weight,
        double xp, double yp, double sig, int& imax, int& jmax,
        double& peak, double& half_light_flux, int& half_light_area);
                     
    void genStarCandidateDirect(const std::string& dirOutput, const std::string& prefix, int nx, int ny, const std::vector<float>& array,
                                const std::vector<int>& weight, int& nstar, int& procError);

    void generateGalCatFileName(const double cRVAL[2], std::string& filename,
                                std::vector<std::string>& sortfile, int& sortnum);
}

#endif // SOURCE_EXTRACTOR_HPP
