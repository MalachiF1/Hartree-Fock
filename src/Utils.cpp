#include "Utils.hpp"

#include <cmath>
#include <ranges>
#include <string>
#include <unordered_map>

namespace Utils
{

extern const std::unordered_map<int, double> atomicMasses = {
    {1, 1.008},   {2, 4.0026},  {3, 6.94},    {4, 9.0122},  {5, 10.81},   {6, 12.011},  {7, 14.007},  {8, 15.999},
    {9, 18.998},  {10, 20.180}, {11, 22.990}, {12, 24.305}, {13, 26.982}, {14, 28.085}, {15, 30.974}, {16, 32.06},
    {17, 35.45},  {18, 39.948}, {19, 39.098}, {20, 40.078}, {21, 44.956}, {22, 47.867}, {23, 50.942}, {24, 51.996},
    {25, 54.938}, {26, 55.845}, {27, 58.933}, {28, 58.693}, {29, 63.546}, {30, 65.38},  {31, 69.723}, {32, 72.630},
    {33, 74.922}, {34, 78.971}, {35, 79.904}, {36, 83.798}, {37, 85.468}, {38, 87.62},  {39, 88.906}, {40, 91.224},
    {41, 92.906}, {42, 95.95},  {43, 98},     {44, 101.07}, {45, 102.91}, {46, 106.42}, {47, 107.87}, {48, 112.41},
    {49, 114.82}, {50, 118.71}, {51, 121.76}, {52, 127.60}, {53, 126.90}, {54, 131.29}, {55, 132.91}, {56, 137.33},
    {57, 138.91}, {58, 140.12}, {59, 140.91}, {60, 144.24}, {61, 145},    {62, 150.36}, {63, 151.96}, {64, 157.25},
    {65, 158.93}, {66, 162.50}, {67, 164.93}, {68, 167.26}, {69, 168.93}, {70, 173.05}, {71, 174.97}, {72, 178.49},
    {73, 180.95}, {74, 183.84}, {75, 186.21}, {76, 190.23}, {77, 192.22}, {78, 195.08}, {79, 196.97}, {80, 200.59},
    {81, 204.38}, {82, 207.2},  {83, 208.98}, {84, 209},    {85, 210},    {86, 222},    {87, 223},    {88, 226},
    {89, 227},    {90, 232.04}, {91, 231.04}, {92, 238.03},
    // I seriously doubt we need anything above Uranium.
};

extern const std::unordered_map<unsigned, std::string> atomicNumberToName = {
    {1, "H"},   {2, "He"},  {3, "Li"},  {4, "Be"},  {5, "B"},   {6, "C"},   {7, "N"},   {8, "O"},   {9, "F"},
    {10, "Ne"}, {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"}, {15, "P"},  {16, "S"},  {17, "Cl"}, {18, "Ar"},
    {19, "K"},  {20, "Ca"}, {21, "Sc"}, {22, "Ti"}, {23, "V"},  {24, "Cr"}, {25, "Mn"}, {26, "Fe"}, {27, "Co"},
    {28, "Ni"}, {29, "Cu"}, {30, "Zn"}, {31, "Ga"}, {32, "Ge"}, {33, "As"}, {34, "Se"}, {35, "Br"}, {36, "Kr"},
    {37, "Rb"}, {38, "Sr"}, {39, "Y"},  {40, "Zr"}, {41, "Nb"}, {42, "Mo"}, {43, "Tc"}, {44, "Ru"}, {45, "Rh"},
    {46, "Pd"}, {47, "Ag"}, {48, "Cd"}, {49, "In"}, {50, "Sn"}, {51, "Sb"}, {52, "Te"}, {53, "I"},  {54, "Xe"},
    {55, "Cs"}, {56, "Ba"}, {57, "La"}, {58, "Ce"}, {59, "Pr"}, {60, "Nd"}, {61, "Pm"}, {62, "Sm"}, {63, "Eu"},
    {64, "Gd"}, {65, "Tb"}, {66, "Dy"}, {67, "Ho"}, {68, "Er"}, {69, "Tm"}, {70, "Yb"}, {71, "Lu"}, {72, "Hf"},
    {73, "Ta"}, {74, "W"},  {75, "Re"}, {76, "Os"}, {77, "Ir"}, {78, "Pt"}, {79, "Au"}, {80, "Hg"}, {81, "Tl"},
    {82, "Pb"}, {83, "Bi"}, {84, "Po"}, {85, "At"}, {86, "Rn"}, {87, "Fr"}, {88, "Ra"}, {89, "Ac"}, {90, "Th"},
    {91, "Pa"}, {92, "U"},

};

size_t CaseInsensitiveHash::operator()(const std::string& str) const
{
    auto toLowerString = [](std::string_view sv) -> std::string
    {
        auto view = sv | std::ranges::views::transform(::tolower);
        return std::string(view.begin(), view.end());
    };
    std::string lowerStr = toLowerString(str);

    return std::hash<std::string>()(lowerStr);
}

bool CaseInsensitiveEqual::operator()(const std::string& a, const std::string& b) const
{
    if (a.length() != b.length())
    {
        return false;
    }
    return std::equal(
        a.begin(),
        a.end(),
        b.begin(),
        [](unsigned char c1, unsigned char c2) { return std::tolower(c1) == std::tolower(c2); }
    );
}

extern const std::unordered_map<std::string, unsigned, CaseInsensitiveHash, CaseInsensitiveEqual> nameToAtomicNumber = {
    {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},   {"O", 8},   {"F", 9},
    {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17}, {"Ar", 18},
    {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27},
    {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"I", 53},
    {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39},  {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44},
    {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"Xe", 54},
    {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63},
    {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
    {"Ta", 73}, {"W", 74},  {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Tl", 81},
    {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92},
};

} // namespace Utils

namespace
{
double SQRT_PI_BY_2 = std::sqrt(M_PI) / 2.0;
}

double boys(int m, double T)
{
    // boys_call_count++;

    if (T < 1e-9)
    {
        // For small T, use the series expansion for stability.
        // F_m(T) = 1/(2m+1) - T/(2m+3) + ...
        return 1.0 / (2.0 * m + 1.0);
    }

    if (T > 30.0)
    {
        // For large T, use the asymptotic expansion.
        double F_i = SQRT_PI_BY_2 / std::sqrt(T);

        for (int i = 0; i < m; ++i) { F_i = (i + 0.5) * F_i / T; }
        return F_i;
    }

    // Use the upward recursion relation.
    double sqrt_T = std::sqrt(T);
    double F_i    = SQRT_PI_BY_2 * std::erf(sqrt_T) / sqrt_T;
    for (int i = 0; i < m; ++i) { F_i = ((2.0 * i + 1.0) * F_i - std::exp(-T)) / (2.0 * T); }
    return F_i;
}


Eigen::MatrixXd inverseSqrtMatrix(const Eigen::MatrixXd& S)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    const Eigen::VectorXd& eigenvalues  = solver.eigenvalues();
    const Eigen::MatrixXd& eigenvectors = solver.eigenvectors();

    Eigen::VectorXd sqrtInvEigenvalues(eigenvalues.size());
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        if (eigenvalues(i) > 0)
        {
            sqrtInvEigenvalues(i) = 1.0 / std::sqrt(eigenvalues(i));
        }
        else
        {
            // Handle very small or zero eigenvalues if they occur
            sqrtInvEigenvalues(i) = 0;
        }
    }

    return eigenvectors * sqrtInvEigenvalues.asDiagonal() * eigenvectors.transpose();
}
