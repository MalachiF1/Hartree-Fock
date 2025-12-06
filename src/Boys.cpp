#include "Boys.hpp"

#include "Utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

/**
 * This implementation is based on:
 * Weiss, A. K., & Ochsenfeld, C. (2015). A rigorous and optimized strategy for the evaluation of the Boys function
 * kernel in molecular electronic structure theory. Journal of Computational Chemistry, 36(18), 1390-1398.
 */
void Boys::calculateBoys(unsigned m_max, double T, std::span<double> F)
{
    // For large T, use the asymptotic expansion.
    if (T > 30.0)
    {
        F[0]              = SQRT_PI_BY_2 / std::sqrt(T);
        const double invT = 1.0 / T;
        for (unsigned m = 0; m < m_max; ++m) { F[m + 1] = (m + 0.5) * F[m] * invT; }
        return;
    }

    // Initialize grid
    const auto& gridManager          = GridManager::getInstance();
    const auto gridIndex             = static_cast<size_t>(T / GRID_STEP_DELTA);
    const double T_grid              = gridIndex * GRID_STEP_DELTA;
    const double T_diff              = T - T_grid;
    const double half_T_diff_squared = 0.5 * T_diff * T_diff;


    // m_max is 0 or 1, handled directly.
    if (m_max < 2)
    {
        // Pointer for m = 0, index = gridIndex, the following are m = 1, 2, 3... for the same gridIndex, until m = MAX_L
        // + 2. i.e gridManager.getValue(1, gridIndex) is equivalent to grid_ptr[1]
        const double* const grid_ptr = &gridManager.getValue(0, gridIndex);

        F[0] = grid_ptr[0] - T_diff * grid_ptr[1] + half_T_diff_squared * grid_ptr[2];

        if (m_max == 0)
            return;

        F[1] = grid_ptr[1] - T_diff * grid_ptr[2] + half_T_diff_squared * grid_ptr[3];

        return;
    }

    const double expMinusT = std::exp(-T);

    // Decide whether to use upward or downward recursion based on the value of T.
    if (T <= TRANSITION_VALUES[m_max])
    {
        // Downward recursion
        const double* const grid_ptr = &gridManager.getValue(m_max, gridIndex);
        F[m_max]                     = grid_ptr[0] - T_diff * grid_ptr[1] + half_T_diff_squared * grid_ptr[2];

        const double twoT = 2.0 * T;
        for (unsigned m = m_max; m >= 1; --m) { F[m - 1] = (twoT * F[m] + expMinusT) / (2 * (m - 1) + 1); }
    }
    else
    {
        // Upward recursion
        const double* const grid_ptr = &gridManager.getValue(0, gridIndex);
        F[0]                         = grid_ptr[0] - T_diff * grid_ptr[1] + half_T_diff_squared * grid_ptr[2];

        const double inv_twoT = 1 / (2 * T);
        for (unsigned m = 0; m < m_max; ++m) { F[m + 1] = ((2.0 * m + 1.0) * F[m] - expMinusT) * inv_twoT; }
    }
}


double Boys::analyticalBoys(unsigned m, double T)
{
    double sum  = 0.0;
    double term = 1.0 / (m + 0.5); // Initial term for the iterative series.
    sum += term;

    for (int i = 1; i < 200; ++i)
    {
        term *= T / (m + i + 0.5);
        sum += term;
        if (std::abs(term) < 1e-15) // Convergence check.
            break;
    }

    return sum * std::exp(-T) / 2.0;
}


GridManager::GridManager()
{
    std::ifstream in("boys_grid.data", std::ios_base::binary);
    if (!in.read((char*)grid.data(), totalSize * sizeof(double)))
    {
        for (size_t i = 0; i <= numGridPointsPerOrder; ++i)
        {
            double T_i = i * GRID_STEP_DELTA;
            for (size_t j = 0; j < gridMaxOrder; ++j) { grid[i * gridMaxOrder + j] = Boys::analyticalBoys(j, T_i); }
        }
        std::ofstream out("boys_grid.data", std::ios_base::binary);
        out.write((char*)grid.data(), totalSize * sizeof(double));
    }
}
