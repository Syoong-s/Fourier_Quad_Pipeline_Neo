#ifndef QUADRATIC_FITTING_HPP
#define QUADRATIC_FITTING_HPP

#include <vector>

// ==========================================
// QuadraticFitting - chi2 minimization utilities
// Method: Least-squares quadratic fit via Cramer's rule (3×3 determinant)
//         to extract best-fit c and uncertainty sigma from chi2(c) samples.
// ==========================================
class QuadraticFitting {
public:
    // Fit x(2,i) = a1*x(1,i)^2 + a2*x(1,i) + a3  (returns a[3])
    static void fit(int n, const std::vector<float>& x_data,
                    std::vector<float>& a);

    // 3×3 determinant
    static float determinant3(const float c[3][3]);
};

#endif  // QUADRATIC_FITTING_HPP
