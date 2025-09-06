#pragma once
#include "Utils.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

namespace Symmetry
{

struct SymmetryElement
{
    std::string symbol;     // Symbol for the symmetry operation, e.g., "C2", "σv", "i", etc.
    Eigen::Matrix3d matrix; // 3x3 Cartesian matrix representation of the symmetry operation.
    Vec3 vector;            // Axis of rotation or reflection plane normal (zero vector for inversion and identity).

    std::string toString() const
    {
        return fmt::format(
            "{}:\nVector: ({:.3f}, {:.3f}, {:.3f})\nMatrix:\n{}", symbol, vector.x(), vector.y(), vector.z(), fmt::streamed(matrix)
        );
    }
};


std::vector<SymmetryElement> generateCiSymels();
std::vector<SymmetryElement> generateCsSymels();
std::vector<SymmetryElement> generateCnSymels(size_t n);
std::vector<SymmetryElement> generateS2nSymmetryElements(size_t n);
std::vector<SymmetryElement> generateCnhSymels(size_t n);
std::vector<SymmetryElement> generateCnvSymels(size_t n);
std::vector<SymmetryElement> generateDnSymels(size_t n);
std::vector<SymmetryElement> generateDnhSymels(size_t n);
std::vector<SymmetryElement> generateDndSymels(size_t n);
std::vector<SymmetryElement> generateTSymels();
std::vector<SymmetryElement> generateThSymels();
std::vector<SymmetryElement> generateTdSymels();
std::vector<SymmetryElement> generateOhSymels();
std::vector<SymmetryElement> generateOSymels();
std::vector<SymmetryElement> generateIhSymels();
std::vector<SymmetryElement> generateISymels();


} // namespace Symmetry
