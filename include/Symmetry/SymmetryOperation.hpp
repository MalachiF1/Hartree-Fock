#pragma once

#include "Utils.hpp"

/**
 * @brief Represents a symmetry operation as a 3x3 matrix acting on position vectors.
 */
struct SymmetryOperation
{
    enum OperationName
    {
        E,
        i,
        Cn,
        Sn,
        sigma,
        sigma_v,
        sigma_d,
        sigma_h,
    };

    std::string toString() const;
    std::string getName() const;

    SymmetryOperation(OperationName name, const Eigen::Matrix3d& matrix, const Vec3& element, size_t n = 1, size_t power = 1);

    /**
     * @return the identity (E) symmetry operation as a SymmetryOperation.
     */
    static SymmetryOperation Identity();

    /**
     * @return the inversion (i) symmetry operation as a SymmetryOperation.
     */
    static SymmetryOperation Inversion();

    OperationName name;
    Eigen::Matrix3d matrix;
    Vec3 element;
    size_t n;     // only relevent to Sn and Cn, 1 otherwise.
    size_t power; // only relevent to Sn and Cn, 1 otherwise.
};
