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
    Atom() = default;
    Atom(unsigned atomicNumber, const Vec3& coords) : atomicNumber(atomicNumber), coords(coords) {}
    unsigned atomicNumber;
    Vec3 coords;
};

namespace Utils
{
struct CaseInsensitiveHash
{
    size_t operator()(const std::string& str) const;
};

struct CaseInsensitiveEqual
{
    bool operator()(const std::string& a, const std::string& b) const;
};

extern const std::unordered_map<int, double> atomicMasses;
extern const std::unordered_map<unsigned, std::string> atomicNumberToName;
extern const std::unordered_map<std::string, unsigned, CaseInsensitiveHash, CaseInsensitiveEqual> nameToAtomicNumber;
} // namespace Utils

/**
 * Computes the double factorial of an integer n.
 * If n is -1, returns 1.
 *
 * @param n The integer for which to compute the double factorial.
 * @return The double factorial of n, or 1 if n is -1.
 */
int dfact(int n);

/**
 * Computes the factorial of an integer n. This is simple itterative implementation, mosty for small n.
 *
 * @param n The integer for which to compute the factorial.
 * @return The factorial of n
 */
int fact(int n);

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
