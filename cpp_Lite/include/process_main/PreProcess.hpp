#ifndef PRE_PROCESS_HPP
#define PRE_PROCESS_HPP

#include <string>
#include <vector>
#include "FitsIO.hpp"

namespace PreProcess {

    // Stage 1 driver
    void preProcess(int iexpo);

    // Individual chip preprocessing
    void chipPreProcess(const std::string& imageFile, const std::string& dirOutput, int cid);

    // Helper functions for preprocessing
    void setBackground(int x_start, int x_end, int y_start, int y_end, int nx, int ny,
                       std::vector<float>& image, const std::vector<int>& weight,
                       int blocksize, int nct, int ncx, int& ierror);

    void flattenChip(int x_start, int x_end, int y_start, int y_end, int nx, int ny, std::vector<float>& array,
                     int nct, int ncx, int& ierror);

    // ==========================================
    // Function: Estimate, validate, and apply one amplifier's noise-sigma plane
    // Method: Use the caller's immutable base-validity map and explicit named sig_scale.
    // ==========================================
    void setSig(int x_start, int x_end, int y_start, int y_end, int nx, int ny, std::vector<float>& image,
                const std::vector<int>& weight, double& aa, double& bb, double& cc, int& ierror,
                double sig_scale);

    void locateDefects(int nx, int ny, const std::vector<float>& array, std::vector<float>& normap,
                       std::vector<int>& weight, int area_max, int area_thresh, int& ierror);

    void mergeDefects(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                      int area_max, double source_thresh, int area_thresh, int& ierror);

    void detectArtificialStripes(int nx, int ny, std::vector<int>& weight,
                                 const std::vector<float>& diffx, const std::vector<float>& diffy,
                                 float sigx, float sigy, float medx, float medy);

    void maskSourceRegions(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                           int area_max, double source_thresh, int area_thresh);

    void detectStripes(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                       int x_smooth, int y_smooth);

    void detectStellarHalo(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                           int npmax, double defect_halo_thresh);

    void detectDent(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                    int npmax, double defect_halo_thresh);

}

#endif // PRE_PROCESS_HPP
