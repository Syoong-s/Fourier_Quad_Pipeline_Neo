#ifndef EX_STAR_HPP
#define EX_STAR_HPP

#include <vector>

namespace ExStar {
    void anaChi2Simple(int n, const float* map1, const float* map2, float& p);
    void getArrayAveStd(const std::vector<float>& arr, float& mean, float& std);
}

#endif // EX_STAR_HPP
