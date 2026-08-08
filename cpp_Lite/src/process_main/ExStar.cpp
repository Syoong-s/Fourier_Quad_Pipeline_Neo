#include "ExStar.hpp"
#include <cmath>
#include <iostream>

namespace ExStar {

void anaChi2Simple(int n, const float* map1, const float* map2, float& p) {
    int n1 = n / 4;
    int n2 = (n / 4) * 3;
    p = 0.0f;

    int start = n1 - 1;
    int end = n2 - 1;

    for (int i = start; i <= end; ++i) {
        for (int j = start; j <= end; ++j) {
            float diff = map1[i * n + j] - map2[i * n + j];
            p += diff * diff;
        }
    }
}

void getArrayAveStd(const std::vector<float>& arr, float& mean, float& std_dev) {
    int n = arr.size();
    if (n <= 1) {
        mean = 0.0f;
        std_dev = 0.0f;
        return;
    }

    double sum = 0.0;
    double sumsq = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += arr[i];
        sumsq += static_cast<double>(arr[i]) * arr[i];
    }

    mean = static_cast<float>(sum / n);
    double variance = (sumsq - (sum * sum) / n) / (n - 1);
    if (variance < 0.0) {
        variance = 0.0;
    }
    std_dev = static_cast<float>(std::sqrt(variance));
}

} // namespace ExStar
