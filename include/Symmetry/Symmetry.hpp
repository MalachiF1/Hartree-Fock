#pragma once

#include "AtomicOrbital.hpp"
#include "Symmetry/PointGroup.hpp"
#include "Symmetry/SymmetryOperation.hpp"
#include "Utils.hpp"

#include <functional>

namespace Symmetry
{

/**
 * Diagonalizes the inertia tensor of a molecular geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @return A pair containing the principal moments (eigenvalues) and the principal axes (eigenvectors).
 */
std::pair<Vec3, Eigen::Matrix3d> diagonalizeInertiaTensor(const std::vector<Atom>& geometry);

/**
 * Translates the center of mass of a molecular geometry to the origin.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @return A new vector of Atoms with the center of mass at the origin.
 */
std::vector<Atom> translateCOMToOrigin(const std::vector<Atom>& geometry);

/**
 * Aligns the molecular geometry canonically based on the point group and principal moments of inertia, according to conventions.
 * This function assumes that the geometry is already translated to the origin.
 * This function modifies the geometry in place and updates the point group operations (in place) accordingly.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param pointGroup The point group object of the molecule.
 * @param principalMoments The principal moments of inertia.
 * @param principalAxes The principal axes of inertia.
 * @param tol Tolerance for numerical comparisons.
 */
void alignCanonically(
    std::vector<Atom>& geometry, PointGroup& pointGroup, const Vec3& principalMoments, const Eigen::Matrix3d& principalAxes, double tol
);

/**
 * Checks if two geometries are identical within a given tolerance.
 *
 * @param geometry1 The first geometry as a vector of Atoms.
 * @param geometry2 The second geometry as a vector of Atoms.
 * @param tol The tolerance for numerical comparisons.
 */
bool areGeometriesSuperimposable(const std::vector<Atom>& geometry1, const std::vector<Atom>& geometry2, double tol);

/**
 * Checks if a given symmetry operation is valid for the provided geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param operation The symmetry operation represented as a 3x3 matrix.
 * @param tol The tolerance for numerical comparisons.
 * @return True if the operation is a valid symmetry operation, false otherwise.
 */
bool isSymmetryOperation(const std::vector<Atom>& geometry, const Eigen::Matrix3d& operation, double tol);

/**
 * Finds groups of symmetrically equivalent atoms (SEA) in the molecular geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of vectors, where each inner vector contains references to symmetrically equivalent Atoms.
 */
std::vector<std::vector<std::reference_wrapper<const Atom>>> findSEAGroups(const std::vector<Atom>& geometry, double tol);

/**
 * Creates a rotation matrix for a given axis and angle.
 *
 * @param axis The axis of rotation as a Vec3.
 * @param angle The angle of rotation in radians.
 * @return A 3x3 rotation matrix.
 */
Eigen::Matrix3d createRotationMatrix(const Vec3& axis, double angle);

/**
 * Creates a reflection matrix for a given normal vector.
 *
 * @param normal The normal vector of the reflection plane as a Vec3.
 * @return A 3x3 reflection matrix.
 */
Eigen::Matrix3d createReflectionMatrix(const Vec3& normal);

/**
 * Creates an improper rotation matrix for a given axis and angle.
 * An improper rotation is a combination of a proper rotation and a reflection.
 *
 * @param axis The axis of rotation as a Vec3.
 * @param angle The angle of rotation in radians.
 * @return A 3x3 improper rotation matrix.
 */
Eigen::Matrix3d createImproperRotationMatrix(const Vec3& axis, double angle);

/**
 * Checks if a given geometry has an inversion center.
 * Assumes that the geometry is already translated to the origin.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param tol The tolerance for numerical comparisons.
 * @return True if the geometry has an inversion center, false otherwise.
 */
bool hasInversion(const std::vector<Atom>& geometry, double tol);

/**
 * Finds candidate rotation axes based on the symmetrically equivalent atom groups.
 * Candidate axes include axes going between pairs of symmetrically equivalent atoms, or through
 * regular polygons formed by three or more symmetrically equivalent atoms.
 *
 * @param SEAGroups The groups of symmetrically equivalent atoms.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of candidate rotation axes as Vec3.
 */
std::vector<Vec3> getRotationAxesCandidates(
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups, double tol
);

/**
 * Finds candidate reflection plane normals based on the symmetrically equivalent atom groups.
 * Candidate normals include normals to planes bisecting pairs of symmetrically equivalent atoms,
 * and normals to planes bisecting triplets of atoms.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param SEAGroups The groups of symmetrically equivalent atoms.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of candidate reflection plane normals as Vec3.
 */
std::vector<Vec3> getReflectionPlaneNormalsCandidates(
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups, double tol
);

/**
 * Checks the given candidate axes for proper rotations.
 * Validates the candidate axes by checking if the rotation operation is a symmetry operation for the geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param SEAGroups The groups of symmetrically equivalent atoms.
 * @param candidateAxes The candidate axes for proper rotations.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of valid proper rotation symmetry operations.
 */
std::vector<SymmetryOperation> findProperRotations(
    const std::vector<Atom>& geometry,
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups,
    const std::vector<Vec3>& candidateAxes,
    double tol
);

/**
 * Checks the given candidate axes for improper rotations.
 * Validates the candidate axes by checking if the improper rotation operation is a symmetry operation for the geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param SEAGroups The groups of symmetrically equivalent atoms.
 * @param candidateAxes The candidate axes for improper rotations.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of valid improper rotation symmetry operations.
 */
std::vector<SymmetryOperation> findImproperRotations(
    const std::vector<Atom>& geometry,
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups,
    const std::vector<Vec3>& candidateAxes,
    double tol
);

/**
 * Checks the given candidate normals for reflection planes.
 * Validates the candidate normals by checking if the reflection operation is a symmetry operation for the geometry.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param candidateNormals The candidate normals for reflection planes.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of valid reflection symmetry operations.
 */
std::vector<SymmetryOperation> findReflectionPlanes(
    const std::vector<Atom>& geometry, const std::vector<Vec3>& candidateNormals, double tol
);

/**
 * Classifies the reflection planes (sigma, sigma_v, sigma_h, sigma_d) based on the proper rotations and reflection
 * planes found in the geometry. This function modifies the reflectionPlanes vector in place.
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param reflectionPlanes The vector of reflection symmetry operations to be classified.
 * @param properRotations The vector of proper rotation symmetry operations.
 * @param tol The tolerance for numerical comparisons.
 */
void classifyReflectionPlanes(
    const std::vector<Atom>& geometry,
    std::vector<SymmetryOperation>& reflectionPlanes,
    const std::vector<SymmetryOperation>& properRotations,
    const std::vector<SymmetryOperation>& improperRotations,
    double tol
);

/**
 * Returns the NxN matrix (where N is the number of basis functions) representation of a 3x3 symmetry operation matrix
 * in the AO basis.
 *
 * @param op the SymmetryOperation object to get the AO representation of.
 * @aos The atomic orbital basis
 * @tol The tolerance for numerical comparisons
 */
Eigen::MatrixXd getAOBasisRep(const SymmetryOperation& op, const std::vector<AtomicOrbital>& aos, double tol);

/**
 * Classifies the point group of a molecular geometry based on its symmetry operations.
 * This function assumes that the geometry is already translated to the origin (and not necessarilly aligned canonically).
 *
 * @param geometry The molecular geometry as a vector of Atoms.
 * @param principalMoments The principal moments of inertia.
 * @param tol The tolerance for numerical comparisons.
 * @return A PointGroup object representing the classified point group.
 */
PointGroup classifyPointGroup(const std::vector<Atom>& geometry, const Vec3& principalMoments, double tol);

/**
 * Classifies a given set of symmetryOperations into classes, where each class contains symmetry operations that
 * are similarity transformations of one another by another symmetry operation.
 *
 * @param operations The symmetry oprations to be classified.
 * @tol The tolerance for numerical comparisons
 * @returns A vector of vectors where each inner vector contains SymmetryOperations of the same class.
 */
std::vector<std::vector<SymmetryOperation>> getClasses(const std::vector<SymmetryOperation>& operations, double tol);

/**
 * Given a set of classes of symmetry operations, returns a canonical name for each class.
 *
 * @param classes A vector of vectors where each inner vector contains SymmetryOperations of the same class.
 * @param pointGroup A PointGroup object representing the point group of the molecule.
 * @param tol The tolerance for numerical comparisons.
 * @return A vector of strings where each string is the canonical name of the class of the matching index.
 */
std::vector<std::string> nameClasses(
    const std::vector<std::vector<SymmetryOperation>>& classes, const PointGroup& pointGroup, double tol
);

} // namespace Symmetry
