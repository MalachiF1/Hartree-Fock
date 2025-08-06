#pragma once

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

#include <Eigen/Dense>
#include <unordered_map>

using Vec3 = Eigen::Vector3d;

/**
 * Represents an atom with its atomic number and coordinates.
 */
struct Atom
{
    int atomicNumber;
    Vec3 coords;
};

namespace Utils
{
extern const std::unordered_map<int, double> atomicMasses;
extern const std::unordered_map<unsigned, std::string> atomicNumberToName;
extern const std::unordered_map<std::string, unsigned> nameToAtomicNumber;
} // namespace Utils

/**
 * Struct to hold shell information from a basis set
 */
struct Shell
{
    Shell() = default;
    Shell(int angularMomentum, const std::vector<double>& exponents, const std::vector<double>& coefficients) :
        angularMomentum(angularMomentum), exponents(exponents), coefficients(coefficients)
    {
    }
    int angularMomentum;
    std::vector<double> exponents;
    std::vector<double> coefficients;
};

/**
 * Computes the double factorial of an integer n.
 * If n is -1, returns 1.
 *
 * @param n The integer for which to compute the double factorial.
 * @return The double factorial of n, or 1 if n is -1.
 */
int dfact(int n);

/**
 * Computes the Boys function F_m(T) using the hypergeometric function.
 *
 * @param m The order of the Boys function.
 * @param T The argument of the Boys function.
 * @return The value of the Boys function F_m(T).
 */
double boys(int m, double T);

/*
 * Compute the inverse square root of a symmetric positive definite matrix S.
 * This is used for symmetric orthogonalization (Löwdin orthogonalization).
 *
 * @param S The symmetric positive definite matrix.
 * @return The inverse square root matrix S^(-1/2).
 */
Eigen::MatrixXd inverseSqrtMatrix(const Eigen::MatrixXd& S);
