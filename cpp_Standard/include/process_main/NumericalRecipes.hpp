#ifndef NUMERICAL_RECIPES_HPP
#define NUMERICAL_RECIPES_HPP

#include <vector>

namespace NumericalRecipes {
    // Sorting and indexing
    void sort(std::vector<float>& arr);
    void sort(int n, std::vector<float>& arr);
    void sortDoub(std::vector<double>& arr);
    void sortDoub(int n, std::vector<double>& arr);
    void indexx(const std::vector<float>& arr, std::vector<int>& indx);
    void indexx(int n, const std::vector<float>& arr, std::vector<int>& indx);
    void sort2i(std::vector<float>& arr, std::vector<int>& brr);
    void sort2i(int n, std::vector<float>& arr, std::vector<int>& brr);

    // Random numbers (thread-local with one F77-equivalent initialization per MPI rank)
    unsigned int initializeRan1Seed(int rank, int numProcs);
    void seedRandom(unsigned int seed);
    double ran1();
    double gasdev();
    void gasdev2(double& x, double& y);

    // Distribution statistics (replacements for get_peak_width and get_peak_width_low_side)
    void getPeakWidth(const std::vector<float>& arr, float& p, float& sig, int& status, int direc);
    void getPeakWidthLowSide(const std::vector<float>& arr, float& p, float& sig);

    // Special functions
    double gammln(double xx);
    double gammq(double a, double x);
}

#endif // NUMERICAL_RECIPES_HPP
