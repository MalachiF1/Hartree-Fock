#include "Symmetry/SymmetryElement.hpp"

#include "Symmetry/Helpers.hpp"

#include <Eigen/Core>
#include <numeric>

namespace Symmetry
{


std::vector<SymmetryElement> generateCiSymels()
{
    std::vector<SymmetryElement> elements;

    // Identity element
    elements.push_back(SymmetryElement {.symbol = "E", .matrix = Eigen::Matrix3d::Identity(), .vector = Vec3::Zero()});

    // Inversion element
    elements.push_back(SymmetryElement {.symbol = "i", .matrix = -Eigen::Matrix3d::Identity(), .vector = Vec3::Zero()});

    return elements;
}

std::vector<SymmetryElement> generateCsSymels()
{
    std::vector<SymmetryElement> elements;

    // Identity element
    elements.emplace_back("E", Eigen::Matrix3d::Identity(), Vec3::Zero());

    // Reflection element (σ)
    elements.emplace_back("s", createReflectionMatrix(Vec3(0, 0, 1)), Vec3(0, 0, 1));

    return elements;
}

std::vector<SymmetryElement> generateCnSymels(size_t n)
{
    std::vector<SymmetryElement> elements;

    // Identity element
    elements.emplace_back("E", Eigen::Matrix3d::Identity(), Vec3::Zero());

    // Rotation elements (C_n)
    for (size_t k = 2; k <= n; ++k)
    {
        if (n % k == 0)
            elements.emplace_back(fmt::format("C{}", n / k), Cn(Vec3(0, 0, 1), n / k), Vec3(0, 0, 1));
        else
        {
            double angle      = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
            Eigen::Matrix3d R = createRotationMatrix(Vec3(0, 0, 1), angle);
            size_t gcd        = std::gcd(n, k);
            size_t order      = n / gcd;
            size_t power      = k / gcd;
            elements.emplace_back(fmt::format("C{}^{}", order, power), R, Vec3(0, 0, 1));
        }
    }

    return elements;
}

std::vector<SymmetryElement> generateS2nSymmetryElements(size_t n)
{
    std::vector<SymmetryElement> elements = generateCnSymels(n);
    elements.emplace_back(fmt::format("S{}", 2 * n), Sn(Vec3(0, 0, 1), 2 * n), Vec3(0, 0, 1));
    size_t order = 2 * n;
    for (size_t k = 3; k < order; ++k)
    {
        if (order % k == 0 || k % 2 == 0)
            continue; // Just a Cn operation
        //
        size_t gcd      = std::gcd(order, k);
        size_t subOrder = order / gcd;
        size_t power    = k / gcd;

        double angle      = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
        Eigen::Matrix3d S = createImproperRotationMatrix(Vec3(0, 0, 1), angle);

        if (power == 1)
            elements.emplace_back(fmt::format("S{}", subOrder), S, Vec3(0, 0, 1));
        else
            elements.emplace_back(fmt::format("S{}^{}", subOrder, power), S, Vec3(0, 0, 1));
    }

    return elements;
}

std::vector<SymmetryElement> generateCnhSymels(size_t n)
{
    std::vector<SymmetryElement> elements = generateCnSymels(n);
    elements.emplace_back("sh", createReflectionMatrix(Vec3(0, 0, 1)), Vec3(0, 0, 1));

    if (n % 2 == 0)
        elements.emplace_back("i", -Eigen::Matrix3d::Identity(), Vec3::Zero());

    if (n <= 2)
        return elements;

    elements.emplace_back(fmt::format("S{}", n), Sn(Vec3(0, 0, 1), n), Vec3(0, 0, 1));
    for (size_t k = 3; k < (n % 2 == 0 ? n : 2 * n); ++k)
    {
        if (n % k == 0 || k % 2 == 0)
            continue; // Just a Cn operation
        //
        size_t gcd   = std::gcd(n, k);
        size_t order = n / gcd;
        size_t power = k / gcd;

        double angle      = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
        Eigen::Matrix3d S = createImproperRotationMatrix(Vec3(0, 0, 1), angle);

        if (power == 1)
            elements.emplace_back(fmt::format("S{}", order), S, Vec3(0, 0, 1));
        else
            elements.emplace_back(fmt::format("S{}^{}", order, power), S, Vec3(0, 0, 1));
    }

    return elements;
}

std::vector<SymmetryElement> generateCnvSymels(size_t n)
{
    std::vector<SymmetryElement> elements = generateCnSymels(n);
    for (size_t k = 0; k < n; ++k)
    {
        double angle = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
        Vec3 normal  = createRotationMatrix(Vec3(0, 0, 1), angle) * Vec3(0, 1, 0);
        if (n % 2 == 0 && k % 2 == 1)
            elements.emplace_back("sd", createReflectionMatrix(normal), normal);
        else
            elements.emplace_back("sv", createReflectionMatrix(normal), normal);
    }

    return elements;
}


std::vector<SymmetryElement> generateDnSymels(size_t n)
{
    std::vector<SymmetryElement> elements = generateCnSymels(n);

    for (size_t k = 0; k < n; ++k)
    {
        double angle = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
        Vec3 axis    = createRotationMatrix(Vec3(0, 0, 1), angle) * Vec3(1, 0, 0);
        if (n % 2 == 0 && k % 2 == 1)
            elements.emplace_back("C2''", Cn(axis, 2), axis);
        else
            elements.emplace_back("C2'", Cn(axis, 2), axis);
    }

    return elements;
}

// std::vector<SymmetryElement> generateDnhSymels(size_t n);
// std::vector<SymmetryElement> generateDndSymels(size_t n);
// std::vector<SymmetryElement> generateTSymels();
// std::vector<SymmetryElement> generateThSymels();
// std::vector<SymmetryElement> generateTdSymels();
// std::vector<SymmetryElement> generateOhSymels();
// std::vector<SymmetryElement> generateOSymels();
// std::vector<SymmetryElement> generateIhSymels();
// std::vector<SymmetryElement> generateISymels();


} // namespace Symmetry
