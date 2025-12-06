#pragma once
#include <array>
#include <span>

constexpr static double SQRT_PI_BY_2    = 0.8862269254527579;
constexpr static double GRID_STEP_DELTA = 5.0e-5;
constexpr static double GRID_MAX        = 30.0;
constexpr static size_t MAX_L           = 32;
constexpr static double EXP_CUTOFF      = 14.39;

class Boys
{
  public:
    static void calculateBoys(unsigned m_max, double T, std::span<double> F);

    static double analyticalBoys(unsigned m, double T);

  private:
    constexpr static std::array<double, 83> TRANSITION_VALUES = {
        0.00,  2.35,  2.35,  2.35,  2.35,  2.38,  2.82,  3.05,  3.72,  3.76,  4.01,  4.61,  4.92,  5.15,
        5.63,  6.02,  6.33,  6.83,  7.15,  7.43,  7.94,  8.27,  8.59,  8.88,  9.21,  9.84,  9.99,  10.40,
        10.71, 11.15, 11.54, 11.88, 12.37, 12.69, 13.21, 13.42, 13.81, 14.25, 14.65, 14.97, 15.31, 15.84,
        15.96, 16.22, 16.69, 17.15, 17.49, 17.92, 18.28, 18.61, 18.85, 19.18, 19.68, 19.94, 20.53, 20.73,
        21.25, 21.54, 21.89, 22.29, 22.72, 23.07, 23.44, 23.90, 24.17, 24.56, 24.85, 25.23, 25.63, 25.96,
        26.43, 26.78, 27.09, 27.56, 27.89, 28.27, 28.71, 29.00, 29.38, 29.84, 30.22, 30.56, 30.85
    };
};

class GridManager
{
  private:
    static constexpr size_t gridMaxOrder          = MAX_L + 2;
    static constexpr size_t numGridPointsPerOrder = static_cast<size_t>(GRID_MAX / GRID_STEP_DELTA) + 1;
    static constexpr size_t totalSize             = (numGridPointsPerOrder + 1) * gridMaxOrder;
    std::array<double, GridManager::totalSize> grid;

    GridManager();

  public:
    static GridManager& getInstance()
    {
        static GridManager instance;
        return instance;
    }

    const std::array<double, GridManager::totalSize>& getGrid() const { return grid; }
    const double& getValue(size_t i, size_t j) const
    {
        constexpr size_t gridMaxOrder = MAX_L + 2;
        return grid[j * gridMaxOrder + i];
    }

    GridManager(const GridManager&)    = delete;
    void operator=(const GridManager&) = delete;
};
