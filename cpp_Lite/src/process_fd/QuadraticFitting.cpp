#include "process_fd/QuadraticFitting.hpp"
#include "general/NumericalRecipes.hpp"

#include <algorithm>

// ==========================================
// Function: determinant3
// Method: Direct cofactor expansion of a 3×3 matrix.
// ==========================================
float QuadraticFitting::determinant3(const float c[3][3]) {
    return c[0][0] * c[1][1] * c[2][2]
         + c[1][0] * c[2][1] * c[0][2]
         + c[2][0] * c[0][1] * c[1][2]
         - c[0][2] * c[1][1] * c[2][0]
         - c[0][1] * c[1][0] * c[2][2]
         - c[0][0] * c[1][2] * c[2][1];
}

// ==========================================
// Function: fit
// Method: Least-squares quadratic fit of x(2,i) = a1*u^2 + a2*u + a3
//         where u = (x1 - cc) / b is a shifted-scaled variable for stability.
//         Solve the 3×3 normal equations via Cramer's rule, then transform
//         coefficients back to the original coordinate.  Faithful translation
//         of Fortran simple_quadratic_fitting.
// ==========================================
void QuadraticFitting::fit(int n, const std::vector<float>& x_data,
                           std::vector<float>& a) {
    // x_data layout: x_data[2*i] = c_i, x_data[2*i+1] = chi2_i
    // Extract and sort the c values to find center and half-range
    std::vector<float> u_vals(n);
    for (int i = 0; i < n; ++i)
        u_vals[i] = x_data[2 * i];

    std::sort(u_vals.begin(), u_vals.end());
    float cc = u_vals[n / 2];
    float b = u_vals[n - 1] - u_vals[0];
    if (b == 0.0) b = 1.0;  // guard against division by zero

    std::vector<float> u(n);
    for (int i = 0; i < n; ++i)
        u[i] = (x_data[2 * i] - cc) / b;

    // Build normal-equation matrix c[3][4] (augmented)
    float c[3][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    for (int i = 0; i < n; ++i) {
        float ui = u[i];
        float yi = x_data[2 * i + 1];
        c[0][0] += ui * ui * ui * ui;
        c[1][0] += ui * ui * ui;
        c[2][0] += ui * ui;

        c[0][1] += ui * ui * ui;
        c[1][1] += ui * ui;
        c[2][1] += ui;

        c[0][2] += ui * ui;
        c[1][2] += ui;
        c[2][2] += 1.0;

        c[0][3] += ui * ui * yi;
        c[1][3] += ui * yi;
        c[2][3] += yi;
    }

    float temp[3][3];
    // Denominator
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            temp[i][j] = c[i][j];
    float denomi = determinant3(temp);
    if (denomi == 0.0) denomi = 1e-30;

    // a1 (replace column 0 with RHS)
    for (int i = 0; i < 3; ++i) {
        for (int j = 1; j < 3; ++j)
            temp[i][j] = c[i][j];
        temp[i][0] = c[i][3];
    }
    a[0] = determinant3(temp) / denomi;

    // a2 (replace column 1 with RHS)
    for (int i = 0; i < 3; ++i) {
        temp[i][0] = c[i][0];
        temp[i][2] = c[i][2];
        temp[i][1] = c[i][3];
    }
    a[1] = determinant3(temp) / denomi;

    // a3 (replace column 2 with RHS)
    for (int i = 0; i < 3; ++i) {
        temp[i][0] = c[i][0];
        temp[i][1] = c[i][1];
        temp[i][2] = c[i][3];
    }
    a[2] = determinant3(temp) / denomi;

    // Transform back from u-space to original c coordinate
    a[0] = a[0] / (b * b);
    a[1] = a[1] / b - 2.0f * a[0] * cc;
    a[2] = a[2] - a[0] * cc * cc - a[1] * cc;
}
