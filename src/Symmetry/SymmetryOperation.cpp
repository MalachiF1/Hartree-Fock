#include "Symmetry/SymmetryOperation.hpp"

#include <Eigen/Dense>


SymmetryOperation::SymmetryOperation(
    OperationName name, const Eigen::Matrix3d& matrix, const Vec3& element, size_t n, size_t power
) :
    name(name), matrix(matrix), element(element), n(n), power(power)
{
    if (name == OperationName::E)
    {
        this->element = Vec3(0, 0, 0);
    }
    else if (name == OperationName::i)
    {
        this->element = Vec3(0, 0, 0);
    }
}

std::string SymmetryOperation::toString() const
{
    std::string nStr       = std::to_string(n);
    std::string powerStr   = (power > 1 ? "^" + std::to_string(power) : "");
    std::string elementStr = "[" + std::to_string(element.x()) + ", " + std::to_string(element.y()) + ", "
                           + std::to_string(element.z()) + "]";
    switch (name)
    {
        case OperationName::E: return "Identity symmetry operation";
        case OperationName::i: return "Inversion symmetry operation";
        case OperationName::Cn: return "C" + nStr + powerStr + ", Axis: " + elementStr;
        case OperationName::Sn: return "S" + nStr + powerStr + ", Axis: " + elementStr;
        case OperationName::sigma: return "σ, Normal: " + elementStr;
        case OperationName::sigma_v: return "σv, Normal: " + elementStr;
        case OperationName::sigma_d: return "σd, Normal: " + elementStr;
        case OperationName::sigma_h: return "σh, Normal: " + elementStr;
        default: return "Unknown symmetry operation";
    }
}

std::string SymmetryOperation::getName() const
{
    std::string nStr     = std::to_string(n);
    std::string powerStr = (power > 1 ? "^" + std::to_string(power) : "");
    switch (name)
    {
        case OperationName::E: return "E";
        case OperationName::i: return "i";
        case OperationName::Cn: return "C" + nStr + powerStr;
        case OperationName::Sn: return "S" + nStr + powerStr;
        case OperationName::sigma: return "sigma";
        case OperationName::sigma_v: return "sigmav";
        case OperationName::sigma_d: return "sigmad";
        case OperationName::sigma_h: return "sigmah";
        default: return "Unknown";
    }
};

SymmetryOperation SymmetryOperation::Identity()
{
    return SymmetryOperation(OperationName::E, Eigen::Matrix3d::Identity(), Vec3(0, 0, 0));
}

SymmetryOperation SymmetryOperation::Inversion()
{
    return SymmetryOperation(OperationName::i, -Eigen::Matrix3d::Identity(), Vec3(0, 0, 0));
}
