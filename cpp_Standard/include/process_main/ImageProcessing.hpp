#ifndef IMAGE_PROCESSING_HPP
#define IMAGE_PROCESSING_HPP

#include <vector>
#include <array>
#include <complex>
#include <functional>

namespace ImageProcessing {
    // Source detection & masking
    void markSource(int n, const std::vector<float>& stamp, std::vector<int>& weight, double sig, 
                    double source_thresh, double core_thresh, int boundx[2], int boundy[2], 
                    double& total_flux, int& total_area, double& peak, double& half_light_flux, 
                    int& half_light_area, int& flag, double& radius, int& xp, int& yp);
                    
    void markNoise(int n, const std::vector<float>& stamp, std::vector<int>& weight, double sig, 
                   double source_thresh, double core_thresh);

    // Stamp flattening
    void flattenStamp2D(int ns, int nl, std::vector<float>& stamp, const std::vector<int>& weight, int& ierror);
    void flattenStampNew(int ns, int nl, std::vector<float>& stamp, const std::vector<int>& weight, int& ierror);

    // Grid decorating (replacing masked pixels with noise)
    void decorateStamp(int ns, double sig, const std::vector<int>& weights, std::vector<float>& stamp);

    // Image smoothing (3x3 grid)
    void smoothGrid33(std::vector<float>& f);
    void smoothImage33(int nx, int ny, std::vector<float>& map);

    // Image smoothing (5x5 grid)
    void smoothGrid55(std::vector<float>& f);
    void smoothGrid55_3rd_order(std::vector<float>& f);
    void smoothImage55(int nx, int ny, std::vector<float>& map, int ord);

    // Image smoothing with holes (excludes central pixel)
    void smoothGrid55WithHole(const std::vector<float>& f, int xh, int yh, float& fc);
    void smoothImage55Hole(int nx, int ny, std::vector<float>& map);
    void smoothImage55HoleLn(int nx, int ny, std::vector<float>& map);

    // Power spectrum utilities
    void processPowers(int n, std::vector<float>& sourcep, const std::vector<float>& noisep);
    void regularizePower(int nx, int ny, std::vector<float>& power, int star_smooth);
    void getPowerShape(int nx, int ny, const std::vector<float>& power, double& e, double thresh_ratio);

    // Background statistics
    void getSigMed(int nx, int ny, const std::vector<float>& image, float& sig, float& med);
    void getEntropy(int nx, int ny, const std::vector<float>& image, double sig, double med, int r, std::vector<float>& entropy);

    // Background continuous removal
    void removeContinuous(int nx, int ny, int npx, int npy, std::vector<float>& map, 
                          const std::function<double(double, int)>& func, int ord);

    // 2D FFT wrappers using FFTW3
    void FFT2D(int n1, int n2, std::vector<std::complex<float>>& arr, int direction);
    void getPower(int n1, int n2, const std::vector<float>& map, std::vector<float>& power, int smooth, double& pc);

    // Drawing utilities for diagnostic shape catalogs
    void drawDot(int nx, int ny, std::vector<float>& map, double x, double y, double intensity, double thickness);
    void drawLine(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity, double thickness);
    void drawRectangle(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity, double thickness);
    void drawBoxFill(int nx, int ny, std::vector<float>& map, double x1, double y1, double x2, double y2, double intensity);
    void reverseColor(int nx, int ny, std::vector<float>& map);
    void drawShearExpo(int n, std::vector<float>& map, 
                       const std::vector<std::array<double, 4>>& pc, 
                       const std::vector<std::array<double, 5>>& sk, 
                       double intensity, double thickness);
}

#endif // IMAGE_PROCESSING_HPP
