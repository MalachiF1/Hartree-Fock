#pragma once
#include "Symmetry/AlignGeometry.hpp"
#include "Utils.hpp"

namespace Symmetry
{

struct SEA
{
    std::vector<size_t> indices;
    Vec3 center;
    InertiaTensor inertiaTensor;
};

std::vector<SEA> findSEAs(const std::vector<Atom>& geometry, double tol);

struct MolecularPG
{
    std::string symbol;
    Vec3 primaryAxis;
    Vec3 secondaryAxis;
};

MolecularPG findPointGroup(const std::vector<Atom>& geometry, double tol);

bool areGeometriesSuperimposable(const std::vector<Atom>& geometry1, const std::vector<Atom>& geometry2, double tol);

// Finds all rotations on the principal axis of rotation (i.e the perpindicular C2s for D groups are nnot found).
// Returns empty vector for linear and spherical top molecules.
// perpindicular C2s and spherical tops are dealt with seperately.
bool isSymmetryOperation(const std::vector<Atom>& geometry, const Eigen::Matrix3d& operation, double tol);

std::vector<std::pair<Vec3, size_t>> findRotations(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol);

std::vector<Vec3> findC2s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol);

std::vector<Vec3> findIcosahedralC3s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol);

std::vector<Vec3> findOctahedralC4s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol);

std::vector<Vec3> findIcosahedralC5s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol);

std::vector<Vec3> findSigmaVs(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, const Vec3& paxis, double tol);

bool hasSigmaH(const std::vector<Atom>& geometry, const Vec3& paxis, double tol);

Vec3 findC2Prime(const std::vector<Atom>& geometry, const std::vector<Vec3>& C2s, const Vec3& pAxis, double tol);

Vec3 findSigmaVPrime(const std::vector<Atom>& geometry, const std::vector<Vec3>& sigmaVs, double tol);

} // namespace Symmetry
