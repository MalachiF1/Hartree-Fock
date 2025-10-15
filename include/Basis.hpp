#pragma once

#include "Utils.hpp"

#include <map>
#include <string>
#include <vector>

/**
 * Struct to hold shell information from a basis set
 */
struct RawShell
{
    RawShell() = default;
    RawShell(int angularMomentum, const std::vector<double>& exponents, const std::vector<double>& coefficients) :
        angularMomentum(angularMomentum), exponents(exponents), coefficients(coefficients)
    {
    }
    unsigned angularMomentum;
    std::vector<double> exponents;
    std::vector<double> coefficients;
};

struct Primitive
{
    double exponent;
    double coefficient;
};

struct Shell
{
    unsigned l;
    unsigned nPrimitives;
    unsigned nAOs;
    unsigned nNormFactors;
    size_t atomIndex;
    std::vector<Primitive> primitives;
    Eigen::Matrix<int, 3, Eigen::Dynamic> angularMomentum;
    Eigen::VectorXd normalizationFactors;
    Eigen::Vector3d center;
};

struct BasisSet
{
    BasisSet() = default;
    BasisSet(const std::string& name, const std::vector<Atom>& geometry);
    std::string name;
    size_t nShells;
    size_t nPrimitives;
    size_t nNormFactors;
    size_t nAOs;
    std::vector<Shell> shells;

    static std::map<int, std::vector<RawShell>> readBasis(const std::string& name, const std::vector<int>& elements);
};
