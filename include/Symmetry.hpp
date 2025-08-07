#pragma once

#include "Utils.hpp"

#include <string>

/**
 * @brief Represents a symmetry operation in a molecule.
 *
 * This class encapsulates the type of symmetry operation (proper rotation, reflection, inversion, or improper
 * rotation), the order of the operation, and the axis of rotation or reflection plane normal.
 */
struct SymmetryOperation
{
    enum SymmetryType
    {
        Cn,       // Proper rotation
        Cn_power, // Proper rotation power (e.g, C_3^2)
        sigma,    // Reflection
        i,        // Inversion
        Sn,       // Improper rotation
        Sn_power, // Improper rotation power (e.g, S_6^5)
        E         // Identity operation
    };

    SymmetryOperation() = default;
    SymmetryOperation(SymmetryType type, size_t order, const Vec3& axis, const Eigen::Matrix3d& matrix, size_t power = 1);
    std::string toString() const;

    SymmetryType type;
    size_t order;
    Vec3 axis;              // Axis of rotation or reflection plane normal
    Eigen::Matrix3d matrix; // Matrix representation of the symmetry operation
    size_t power;
};

/**
 * @brief Represents a point group classification of a molecule.
 *
 * This class encapsulates the type of point group and its order.
 */
struct PointGroup
{
    enum PointGroupClass
    {
        Dinfh,
        Cinfh,
        I,
        Ih,
        O,
        Oh,
        T,
        Td,
        Th,
        Dn,
        Dnh,
        Dnd,
        Cn,
        Cnh,
        Cnv,
        Cs,
        Ci,
        S2n,
        C1,
    };

    PointGroup();
    PointGroup(PointGroupClass pointGroupClass, size_t n, const std::vector<SymmetryOperation>& symmetryOperations);

    std::string toString() const;

    PointGroupClass pointGroupClass;
    size_t n;
    std::vector<SymmetryOperation> symmetryOperations;
    size_t order;
};

namespace Symmetry
{

/**
 * Diagonalizes the inertia tensor of a molecule to find principal moments and axes.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @return A pair containing the principal moments (Eigen::Vec3d) and the principal axes (Eigen::Matrix3d) as columns of a 3x3 matrix.
 */
std::pair<Vec3, Eigen::Matrix3d> diagonalizeInertiaTensor(const std::vector<Atom>& geometry);

/**
 * Translates the molecule so the center of mass is at the origin and rotates it to align with the principal axes.
 * The z axis is aligned with the principal axis of rotation, if it exists.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param principalAxes The principal axes of the molecule as a 3x3 matrix.
 * @param tol The tolerance for considering two distances or coordinates identical, used for determining the principle axis of rotation.
 * @return A new vector of atoms representing the translated and rotated geometry.
 */
std::vector<Atom> translateAndRotate(const std::vector<Atom>& geometry, const Eigen::Matrix3d& principalAxes, double tol);

/**
 * Finds the symmetry equivalent atom (SEA) groups in the molecule's geometry.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return A vector of vectors, where each inner vector contains the indices of atoms that are symmetry equivalent.
 */
std::vector<std::vector<size_t>> findSEAGroups(const std::vector<Atom>& geometry, double tol);

/**
 * Checks if two molecular geometries are superimposable within a given tolerance.
 *
 * @param geometry1
 * @param geometry2
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return True if the molecules are superimposable, false otherwise.
 */
bool areMoleculesSuperimposable(const std::vector<Atom>& geometry1, const std::vector<Atom>& geometry2, double tol);

/**
 * Checks if a given symmetry operation is valid for the molecule's geometry.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param operation The symmetry operation represented as a 3x3 matrix.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return True if the operation is a valid symmetry operation, false otherwise.
 */
bool isSymmetryOperation(const std::vector<Atom>& geometry, const Eigen::Matrix3d& operation, double tol);

/**
 * Checks if the molecule has an inversion center. This function assumes that the Center of mass is centered at the
 * origin and that the geometry is aligned with the principal axes.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return A vector containing a single symmetry operation representing the inversion center if it exists, otherwise an empty vector.
 */
std::vector<SymmetryOperation> findInvertion(const std::vector<Atom>& geometry, double tol);

/**
 * Finds proper rotations in the molecule's geometry. This function assumes that the Center of mass is centered at the
 * origin and that the geometry is aligned with the principal axes.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param SEAIndecies The symmetry equivalent atom groups.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return A vector of symmetry operations representing the proper rotations.
 */
std::vector<SymmetryOperation> findProperRotations(
    const std::vector<Atom>& geometry, const std::vector<std::vector<size_t>>& SEAIndecies, double tol
);

/**
 * Finds the reflection planes in the molecule's geometry. This function assumes that the Center of mass is centered
 * at the origin and that the geometry is aligned with the principal axes.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param SEAIndecies The symmetry equivalent atom groups.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return A vector of symmetry operations representing the reflection planes.
 */
std::vector<SymmetryOperation> findReflectionPlanes(
    const std::vector<Atom>& geometry, const std::vector<std::vector<size_t>>& SEAIndecies, double tol
);

/**
 * Finds improper rotations in the molecule's geometry. This function assumes that the Center of mass is centered at
 * the origin and that the geometry is aligned with the principal axes. Inversion operation is not returned by this
 * function, see findInvertion() instead.
 *
 * @param geometry The geometry of the molecule as a vector of atoms.
 * @param SEAIndecies The symmetry equivalent atom groups.
 * @param tol The tolerance for considering two distances or coordiantes identical.
 * @return A vector of symmetry operations representing the improper rotations.
 */
std::vector<SymmetryOperation> findImproperRotations(
    const std::vector<Atom>& geometry, const std::vector<std::vector<size_t>>& SEAIndecies, double tol
);

/**
 * Classifies the point group of a molecule based on its symmetry operations and principal moments.
 * Assumes that the molecule is aligned with the principal axes such that the z-axis is the principal axis of rotation, if it exists.
 *
 * @param symmetryOperations The symmetry operations found in the molecule's geometry.
 * @param principalMoments The principal moments of inertia of the molecule.
 * @param tol The tolerance for considering a principal moment as zero, or two principal moments as identical.
 * @return A PointGroup object representing the classified point group.
 */
PointGroup classifyPointGroup(const std::vector<SymmetryOperation>& symmetryOperations, const Vec3& principalMoments, double tol);

/**
 * Creates a rotation matrix for a given axis and angle.
 *
 * @param axis The axis of rotation.
 * @param angle The angle of rotation in radians.
 * @return A 3x3 rotation matrix as an Eigen::Matrix3d.
 */
Eigen::Matrix3d createRotationMatrix(const Vec3& axis, double angle);

/**
 * Creates a reflection matrix for a given normal vector.
 *
 * @param normal The normal vector of the reflection plane.
 * @return A 3x3 reflection matrix as an Eigen::Matrix3d.
 */
Eigen::Matrix3d createReflectionMatrix(const Vec3& normal);

}; // namespace Symmetry
