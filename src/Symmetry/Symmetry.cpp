#include "Symmetry/Symmetry.hpp"

#include "Eigen/Core"
#include "Symmetry/PointGroup.hpp"
#include "Symmetry/SymmetryOperation.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <unordered_set>


// Anonymous namespace for helper functions.
namespace
{

bool areCollinear(const Vec3& v1, const Vec3& v2, double tol)
{
    return (v1.normalized().cross(v2.normalized())).squaredNorm() < (tol * tol);
}


bool areOrthogonal(const Vec3& v1, const Vec3& v2, double tol)
{
    return std::abs(v1.dot(v2)) < tol;
}


// Check if two vectors are equal within a tolerance
bool areEqual(const Vec3& v1, const Vec3& v2, double tol)
{
    return (v1 - v2).squaredNorm() < (tol * tol);
}


// Check if two matrices are equal within a tolerance
bool areEqual(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2, double tol)
{
    return (m1 - m2).squaredNorm() < tol * tol;
}


// Check if two doubles are equal within a tolerance
bool areEqual(double d1, double d2, double tol)
{
    return std::abs(d1 - d2) < tol;
}


bool isRegularPolygon(const std::vector<Vec3>& points, const Vec3& center, const Vec3& normal, double tol)
{
    if (points.size() < 3)
        return false;

    // create an orthonormal basis on the plane
    Vec3 u = (points[0] - center).normalized();
    Vec3 v = normal.cross(u);

    // calculate the angle for each point with respect to the orthonormal basis
    std::vector<double> angles;
    angles.reserve(points.size());
    for (const auto& p : points)
    {
        double cosAngle = (p - center).dot(u);
        double sinAngle = (p - center).dot(v);
        angles.push_back(atan2(sinAngle, cosAngle));
    }

    // sort the angles
    std::sort(angles.begin(), angles.end());

    // check the angles are equally spaced
    double diff = angles[1] - angles[0];
    for (size_t i = 1; i < angles.size() - 1; ++i)
    {
        if (!areEqual(angles[i + 1] - angles[i], diff, tol))
            return false;
    }

    // check the first and last angles are also equally spaced
    return areEqual((2 * M_PI + angles[0]) - angles.back(), diff, tol);
}


// Custom hasher for Vec3 to group similar vectors into the same hash bucket.
struct Vec3Hasher
{
    double tol;
    Vec3Hasher(double tol) : tol(tol) {}

    std::size_t operator()(const Vec3& v) const
    {
        // Create a canonical representation of the vector's direction.
        Vec3 canonical = v; // already normalized

        // Ensure the vector points in a consistent direction.
        // If the first non-zero component is negative, flip the vector (can't be zero vector due to the check in insert()).
        double* firstNonZero = std::find_if(
            canonical.data(), canonical.data() + 3, [this](double val) { return !areEqual(val, 0.0, tol); }
        );
        if (*firstNonZero < -tol)
            canonical = -canonical;

        long ix = static_cast<long>(std::round(canonical.x() / tol));
        long iy = static_cast<long>(std::round(canonical.y() / tol));
        long iz = static_cast<long>(std::round(canonical.z() / tol));

        auto hashCombine = [](std::size_t& seed, const long& v)
        { seed ^= std::hash<long> {}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); };

        size_t seed = 0;
        hashCombine(seed, ix);
        hashCombine(seed, iy);
        hashCombine(seed, iz);
        return seed;
    }
};


// Custom equality comparator for Vec3 to handle floating-point tolerance.
struct Vec3Equal
{
    double tolerance;
    Vec3Equal(double tol) : tolerance(tol) {}
    bool operator()(const Vec3& a, const Vec3& b) const
    {
        return ((a.cross(b)).squaredNorm() < (tolerance * tolerance));
    }
};


class UniqueVec3Set
{
  private:
    std::unordered_set<Vec3, Vec3Hasher, Vec3Equal> set;
    double tol;

  public:
    UniqueVec3Set(size_t initial_buckets, double tol) : set(initial_buckets, Vec3Hasher(tol), Vec3Equal(tol)), tol(tol)
    {
    }

    void insert(const Vec3& v)
    {
        if (v.squaredNorm() > tol * tol)
            set.insert(v.normalized());
    }

    void erase(const Vec3& v) { set.erase(v.normalized()); }

    auto begin() const { return set.begin(); }
    auto end() const { return set.end(); }
    size_t size() const { return set.size(); }
};


// A struct to uniquely identify a plane by its normal and distance from origin.
struct PlaneIdentifier
{
    Vec3 normal;
    double distance;
};


struct PlaneIdentifierEqual
{
    double tol;
    PlaneIdentifierEqual(double tol) : tol(tol) {}

    bool operator()(const PlaneIdentifier& p1, const PlaneIdentifier& p2) const
    {
        long distanceKey_1 = static_cast<long>(std::round(p1.distance / tol));
        long distanceKey_2 = static_cast<long>(std::round(p2.distance / tol));
        if (distanceKey_1 != distanceKey_2)
            return false;
        return Vec3Equal(tol)(p1.normal, p2.normal);
    }
};


// A hasher for the PlaneIdentifier struct.
struct PlaneIdentifierHasher
{
    double tol;
    PlaneIdentifierHasher(double tol) : tol(tol) {}

    std::size_t operator()(const PlaneIdentifier& p) const
    {
        size_t h1 = Vec3Hasher(tol)(p.normal);
        size_t h2 = std::hash<long> {}(p.distance);
        return (h1 ^ (std::hash<long> {}(h2) + 0x9e3779b9 + (h1 << 6) + (h1 >> 2)));
    }
};

} // namespace


namespace Symmetry
{

std::pair<Vec3, Eigen::Matrix3d> diagonalizeInertiaTensor(const std::vector<Atom>& geometry)
{
    auto inertiaTensor = Eigen::MatrixXd(3, 3);
    inertiaTensor.setZero();

    for (const auto& atom : geometry)
    {
        double mass = Utils::atomicMasses.at(atom.atomicNumber);
        Vec3 r      = atom.coords;

        // Calculate the inertia tensor components
        inertiaTensor(0, 0) += mass * (r.y() * r.y() + r.z() * r.z());
        inertiaTensor(1, 1) += mass * (r.x() * r.x() + r.z() * r.z());
        inertiaTensor(2, 2) += mass * (r.x() * r.x() + r.y() * r.y());

        inertiaTensor(0, 1) -= mass * r.x() * r.y();
        inertiaTensor(0, 2) -= mass * r.x() * r.z();
        inertiaTensor(1, 2) -= mass * r.y() * r.z();
    }

    // Since the inertia tensor is symmetric, we can fill in the lower triangle.
    inertiaTensor(1, 0) = inertiaTensor(0, 1);
    inertiaTensor(2, 0) = inertiaTensor(0, 2);
    inertiaTensor(2, 1) = inertiaTensor(1, 2);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertiaTensor);
    Vec3 eigenValues             = solver.eigenvalues();
    Eigen::Matrix3d eigenVectors = solver.eigenvectors();

    return std::pair<Vec3, Eigen::Matrix3d>(eigenValues, eigenVectors);
}


std::vector<Atom> translateCOMToOrigin(const std::vector<Atom>& geometry)
{
    // Find center of mass
    Vec3 centerOfMass = Eigen::Vector3d::Zero();
    double totalMass  = 0.0;
    for (const auto& atom : geometry)
    {
        centerOfMass += atom.coords * Utils::atomicMasses.at(atom.atomicNumber);
        totalMass += Utils::atomicMasses.at(atom.atomicNumber);
    }
    centerOfMass /= totalMass;

    // Translate the center of mass to the origin
    auto newGeometry = geometry;
    for (auto& atom : newGeometry) { atom.coords -= centerOfMass; }

    return newGeometry;
}


void alignCanonically(
    std::vector<Atom>& geometry, PointGroup& pointGroup, const Vec3& principalMoments, const Eigen::Matrix3d& principalAxes, double tol
)
{
    std::vector<SymmetryOperation>& symmetryOperations = pointGroup.operations;

    // sort the principal moments of inertia along with their corresponding axes
    std::vector<std::pair<double, Vec3>> momentAxisPairs = {
        {principalMoments[0], principalAxes.col(0)},
        {principalMoments[1], principalAxes.col(1)},
        {principalMoments[2], principalAxes.col(2)}
    };
    std::sort(
        momentAxisPairs.begin(),
        momentAxisPairs.end(),
        [](const std::pair<double, Vec3>& a, const std::pair<double, Vec3>& b) { return a.first < b.first; }
    );

    Vec3 sortedMoments = {momentAxisPairs[0].first, momentAxisPairs[1].first, momentAxisPairs[2].first};
    Eigen::Matrix3d sortedAxes;
    sortedAxes << momentAxisPairs[0].second.normalized(), momentAxisPairs[1].second.normalized(),
        momentAxisPairs[2].second.normalized();

    // handle linear molecules
    if (areEqual(sortedMoments[0], 0, tol) && areEqual(sortedMoments[1], sortedMoments[2], tol))
    {
        // rotate the geometry to align the molecule with the z-axis
        auto atom = std::find_if(
            geometry.begin(),
            geometry.end(),
            [&tol](const Atom& atom) { return !areEqual(atom.coords, Vec3(0, 0, 0), tol); }
        );
        Vec3 principalAxis   = atom->coords.normalized();
        Vec3 rotationAxis    = principalAxis.cross(Vec3(0, 0, 1));
        double rotationAngle = acos(principalAxis.dot(Vec3(0, 0, 1)));

        Eigen::Matrix3d rotationMatrix = createRotationMatrix(rotationAxis, rotationAngle);
        for (auto& atom : geometry) { atom.coords = rotationMatrix * atom.coords; }
        for (auto& op : symmetryOperations)
        {
            op.element = rotationMatrix * op.element;
            op.matrix  = rotationMatrix * op.matrix * rotationMatrix.transpose();
        }
        return;
    }

    // If not changed, the geometry is aligned with the principal axes of inertia.
    Vec3 xAxis = sortedAxes.col(0).normalized();
    Vec3 yAxis = sortedAxes.col(1).normalized();
    Vec3 zAxis = sortedAxes.col(2).normalized();

    // For C1 and Ci point groups we can simply align the geometry with the principal axes.
    if (pointGroup.name == PointGroup::C1 || pointGroup.name == PointGroup::Ci)
    {
        Eigen::Matrix3d axes;
        axes << xAxis, yAxis, zAxis;
        Eigen::Matrix3d rotationMatrix = axes.transpose();
        for (auto& atom : geometry) { atom.coords = rotationMatrix * atom.coords; }
        for (auto& op : symmetryOperations)
        {
            op.element = rotationMatrix * op.element;
            op.matrix  = rotationMatrix * op.matrix * rotationMatrix.transpose();
        }
        return;
    }

    // For Cs point group, we set the reflection plane to the xy-plane.
    if (pointGroup.name == PointGroup::Cs)
    {
        auto sigmaH = std::find_if(
            symmetryOperations.begin(),
            symmetryOperations.end(),
            [](const SymmetryOperation& op) { return (op.name == SymmetryOperation::sigma_h); }
        );

        zAxis = sigmaH->element.normalized();
    }
    // For non spherical-top molecules, we set the z-axis to the principal axis of rotation.
    else if (pointGroup.name != PointGroup::T && pointGroup.name != PointGroup::Td && pointGroup.name != PointGroup::Th
             && pointGroup.name != PointGroup::O && pointGroup.name != PointGroup::Oh
             && pointGroup.name != PointGroup::I && pointGroup.name != PointGroup::Ih)
    {
        auto principalRotation = std::max_element(
            symmetryOperations.begin(),
            symmetryOperations.end(),
            [](const SymmetryOperation& a, const SymmetryOperation& b)
            {
                if ((a.name != SymmetryOperation::Cn && a.name != SymmetryOperation::Sn) || a.power != 1)
                    return true;
                else if ((b.name != SymmetryOperation::Cn && b.name != SymmetryOperation::Sn) || b.power != 1)
                    return false;
                return (a.n < b.n);
            }
        );

        zAxis = principalRotation->element;
    }

    switch (pointGroup.name)
    {
        case PointGroup::Cnh:
        {
            auto sigmaH = std::find_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                [](const SymmetryOperation& op) { return (op.name == SymmetryOperation::sigma_h); }
            );

            if (zAxis.dot(sigmaH->element) < 0)
                sigmaH->element = -sigmaH->element; // Ensure the reflection plane is in the correct direction.

            [[fallthrough]];
        }
        case PointGroup::Cn: [[fallthrough]];
        case PointGroup::S2n:
        {
            // Get the projection of the geometry on the xy plane (relative to the new z-axis).
            // The x Axis is the eigenvector corresponding to the largest eigenvalue of the inertia tensor

            // rotate the geometry to align the z-axis, then take the projection on the xy plane
            auto projectedGeometry         = geometry;
            Vec3 rotationAxis              = zAxis.cross(Vec3(0, 0, 1));
            double rotationAngle           = acos(zAxis.dot(Vec3(0, 0, 1)));
            Eigen::Matrix3d rotationMatirx = createRotationMatrix(rotationAxis, rotationAngle);
            for (auto& atom : projectedGeometry)
            {
                atom.coords     = rotationMatirx * atom.coords;
                atom.coords.z() = 0.0;
            }

            // Diagonalize the inertia tensor in the xy plane.
            Eigen::Matrix2d inertiaTensor2D = Eigen::Matrix2d::Zero();
            for (const auto& atom : projectedGeometry)
            {
                double mass = Utils::atomicMasses.at(atom.atomicNumber);
                Vec3 r      = atom.coords;

                // Calculate the inertia tensor components in the xy plane
                inertiaTensor2D(0, 0) += mass * (r.y() * r.y());
                inertiaTensor2D(1, 1) += mass * (r.x() * r.x());
                inertiaTensor2D(0, 1) -= mass * r.x() * r.y();
            }
            inertiaTensor2D(1, 0) = inertiaTensor2D(0, 1);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver2D(inertiaTensor2D);
            Eigen::Vector2d eigenvector = solver2D.eigenvectors().col(1);

            xAxis = Vec3(eigenvector.x(), eigenvector.y(), 0.0).normalized();
            xAxis = rotationMatirx.transpose() * xAxis; // Rotate back to the original coordinate system.

            yAxis = zAxis.cross(xAxis); // y axis is determined by the right-hand rule
            break;
        }
        case PointGroup::Cnv:
        {
            // Align a vertical reflection plane with the XZ plane.
            auto sigmaV = std::find_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                [](const SymmetryOperation& op) { return (op.name == SymmetryOperation::sigma_v); }
            );

            yAxis = sigmaV->element.normalized();
            xAxis = yAxis.cross(zAxis);
            break;
        }
        case PointGroup::Dn: [[fallthrough]];
        case PointGroup::Dnd: [[fallthrough]];
        case PointGroup::Dnh:
        {
            // Align the x-axis with the a C2' axis.
            std::vector<SymmetryOperation> perpendicularC2s;
            std::copy_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                std::back_inserter(perpendicularC2s),
                [&zAxis, &tol](const SymmetryOperation& op)
                { return (op.name == SymmetryOperation::Cn && op.n == 2 && areOrthogonal(op.element, zAxis, tol)); }
            );

            // C2' axes are the ones passing through the most atoms, or intersecting the largest number of bonds.
            auto C2Prime = std::max_element(
                perpendicularC2s.begin(),
                perpendicularC2s.end(),
                [&geometry, &tol](const SymmetryOperation& a, const SymmetryOperation& b)
                {
                    size_t atomsIntersected_a = std::count_if(
                        geometry.begin(),
                        geometry.end(),
                        [&a, &tol](const Atom& atom) { return areCollinear(atom.coords, a.element, tol); }
                    );
                    size_t atomsIntersected_b = std::count_if(
                        geometry.begin(),
                        geometry.end(),
                        [&a, &tol](const Atom& atom) { return areCollinear(atom.coords, a.element, tol); }
                    );
                    if (atomsIntersected_a != atomsIntersected_b)
                        return (atomsIntersected_a < atomsIntersected_b);

                    // If the number of atoms intersected is the same, we compare the number of bonds intersected.
                    size_t bondsIntersected_a = 0;
                    size_t bondsIntersected_b = 0;
                    for (size_t i = 0; i < geometry.size(); ++i)
                    {
                        for (size_t j = i + 1; j < geometry.size(); ++j)
                        {
                            if (areCollinear(geometry[i].coords + geometry[j].coords, a.element, tol))
                                bondsIntersected_a++;
                            if (areCollinear(geometry[i].coords + geometry[j].coords, b.element, tol))
                                bondsIntersected_b++;
                        }
                    }
                    return (bondsIntersected_a < bondsIntersected_b);
                }
            );

            xAxis = C2Prime->element;
            yAxis = zAxis.cross(xAxis); // y axis is determined by the right-hand rule
            break;
        }
        case PointGroup::T: [[fallthrough]];
        case PointGroup::Td: [[fallthrough]];
        case PointGroup::Th:
        {
            // Align the three C2 axes with the caretesian Axes.
            std::vector<SymmetryOperation> C2Axes;
            std::copy_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                std::back_inserter(C2Axes),
                [](const SymmetryOperation& op) { return (op.name == SymmetryOperation::Cn && op.n == 2); }
            );

            xAxis = C2Axes[0].element;
            yAxis = C2Axes[1].element;
            zAxis = C2Axes[2].element;
            if (xAxis.cross(yAxis).dot(zAxis) < 0)
            {
                // If the axes are not right-handed, swap the x and y axes.
                std::swap(xAxis, yAxis);
                zAxis = xAxis.cross(yAxis);
            }
            break;
        }
        case PointGroup::O: [[fallthrough]];
        case PointGroup::Oh:
        {
            // Align the three C4 axes with the caretesian Axes.
            std::vector<SymmetryOperation> C4Axes;
            std::copy_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                std::back_inserter(C4Axes),
                [](const SymmetryOperation& op)
                { return (op.name == SymmetryOperation::Cn && op.n == 4 && op.power == 1); }
            );

            zAxis = C4Axes[0].element;
            yAxis = C4Axes[1].element;
            xAxis = C4Axes[2].element;
            break;
        }
        case PointGroup::I: [[fallthrough]];
        case PointGroup::Ih:
        {
            // Align one of the five sets of mutually orthogonal C2 axes
            // with the caretesian Axes, such that the pair of closest C5
            // axes to the z-axis lie on the YZ plane.
            std::vector<SymmetryOperation> C2Axes;
            std::copy_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                std::back_inserter(C2Axes),
                [](const SymmetryOperation& op) { return (op.name == SymmetryOperation::Cn && op.n == 2); }
            );
            zAxis = C2Axes[0].element; // Choose the first C2 axis as the z-axis.

            auto orthogonalC2 = std::find_if(
                C2Axes.begin(),
                C2Axes.end(),
                [&zAxis, &tol](const SymmetryOperation& op) { return areOrthogonal(op.element, zAxis, tol); }
            );
            xAxis = orthogonalC2->element; // Take the first orthogonal C2 axis as the x-axis

            yAxis = zAxis.cross(xAxis); // y axis is determined by the right-hand rule, and will fall on a C2 axis.

            std::vector<SymmetryOperation> C5Axes;
            std::copy_if(
                symmetryOperations.begin(),
                symmetryOperations.end(),
                std::back_inserter(C5Axes),
                [](const SymmetryOperation& op)
                { return (op.name == SymmetryOperation::Cn && op.n == 5 && op.power == 1); }
            );

            // find the pair of C5 axes closest to the z-axis
            auto closestC5_1 = std::min_element(
                C5Axes.begin(),
                C5Axes.end(),
                [&zAxis](const SymmetryOperation& a, const SymmetryOperation& b)
                { return (a.element - zAxis).squaredNorm() < (b.element - zAxis).squaredNorm(); }
            );

            auto closestC5_2 = std::find_if(
                C5Axes.begin(),
                C5Axes.end(),
                [&closestC5_1, &zAxis, &tol](const SymmetryOperation& op)
                {
                    return (
                        !areEqual(op.element, closestC5_1->element, tol)
                        && areEqual((op.element - zAxis).squaredNorm(), (closestC5_1->element - zAxis).squaredNorm(), tol * tol)
                    );
                }
            );

            // If the C5 Axes do not lie on the YZ plane, swap the x and y axes.
            if (!(areEqual(closestC5_1->element.x(), 0, tol) && areEqual(closestC5_2->element.x(), 0, tol)))
            {
                std::swap(xAxis, yAxis);
                zAxis = -zAxis; // maintain the right-hand rule
            }

            break;
        }

        default: break;
    }

    // Apply the rotation
    Eigen::Matrix3d axes;
    axes << xAxis, yAxis, zAxis;
    Eigen::Matrix3d rotationMatrix = axes.transpose();
    for (auto& atom : geometry) { atom.coords = rotationMatrix * atom.coords; }
    for (auto& op : symmetryOperations)
    {
        op.element = rotationMatrix * op.element;
        op.matrix  = rotationMatrix * op.matrix * rotationMatrix.transpose();
    }

    // Swap sign of rotation axes if necessary
    switch (pointGroup.name)
    {
        case PointGroup::T: [[fallthrough]];
        case PointGroup::Td: [[fallthrough]];
        case PointGroup::Th: [[fallthrough]];
        case PointGroup::O: [[fallthrough]];
        case PointGroup::Oh:
        {
            for (auto& op : symmetryOperations)
            {
                if ((op.name == SymmetryOperation::Cn && op.n == 3) || (op.name == SymmetryOperation::Sn && op.n == 6))
                {
                    std::vector<double> sortedComponents = {op.element.x(), op.element.y(), op.element.z()};
                    std::sort(sortedComponents.begin(), sortedComponents.end());
                    if (sortedComponents[2] < -tol || (sortedComponents[1] > tol && sortedComponents[0] < -tol))
                    {
                        op.element = -op.element;
                        op.power   = op.n - op.power;
                    }
                }
            }
            break;
        }
        default: break;
    }
}


std::vector<std::vector<std::reference_wrapper<const Atom>>> findSEAGroups(const std::vector<Atom>& geometry, double tol)
{
    size_t n = geometry.size();

    if (n == 0)
        return {};

    // Use a map to group atoms. The key is the "fingerprint" and the value is the list of atoms sharing that
    // fingerprint. The fingerprint is a vector of pairs, where each pair contains the distance to another atom and its
    // atomic number. Each fingerprint will be sorted to create a canonical form for comparison.
    using Fingerprint = std::vector<std::pair<double, size_t>>;

    auto compareFingerprints = [tol](const Fingerprint& a, const Fingerprint& b)
    {
        for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        {
            if (areEqual(a[i].first, b[i].first, tol) && a[i].second == b[i].second)
                continue;

            if (areEqual(a[i].first, b[i].first, tol))
                return a[i].second < b[i].second;

            return a[i].first < b[i].first;
        }
        return a.size() < b.size();
    };

    std::map<Fingerprint, std::vector<std::reference_wrapper<const Atom>>, decltype(compareFingerprints)> groups(
        compareFingerprints
    );

    for (size_t i = 0; i < n; ++i)
    {
        Fingerprint fingerprint;
        fingerprint.reserve(n);
        for (size_t j = 0; j < n; ++j)
        {
            double distance = (geometry[i].coords - geometry[j].coords).norm();
            fingerprint.emplace_back(distance, geometry[j].atomicNumber);
        }

        std::sort(
            fingerprint.begin(),
            fingerprint.end(),
            [tol](const auto& a, const auto& b)
            {
                if (areEqual(a.first, b.first, tol))
                    return a.second < b.second; // Sort by atomic number if distances are equal

                return a.first < b.first;
            }
        );

        groups[fingerprint].push_back(std::cref(geometry[i]));
    }

    // Collect the final groups from the map's values.
    std::vector<std::vector<std::reference_wrapper<const Atom>>> SEAGroups;
    SEAGroups.reserve(groups.size());
    for (const auto& [fingerprint, atomGroup] : groups) { SEAGroups.push_back(atomGroup); }

    return SEAGroups;
}

std::vector<Vec3> getRotationAxesCandidates(const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups, double tol)
{
    UniqueVec3Set candidateAxes(10, tol);

    // add axes between each pair of symmetrically equivalent atoms
    for (const auto& group : SEAGroups)
    {
        if (group.size() < 2)
            continue;

        // axis going between the two atoms
        for (size_t i = 0; i < group.size(); ++i)
        {
            for (size_t j = i + 1; j < group.size(); ++j)
            {
                candidateAxes.insert((group[i].get().coords + group[j].get().coords));
            }
        }

        if (group.size() < 3)
            continue;

        // find unique normals to planes that contain three atoms
        UniqueVec3Set planeNormals(10, tol); // Use a set to store unique axes
        for (size_t i = 0; i < group.size(); ++i)
        {
            const Vec3& p1 = group[i].get().coords;
            for (size_t j = i + 1; j < group.size(); ++j)
            {
                const Vec3& p2 = group[j].get().coords;
                for (size_t k = j + 1; k < group.size(); ++k)
                {
                    const Vec3& p3 = group[k].get().coords;

                    Vec3 normal = (p2 - p1).cross(p3 - p1);
                    planeNormals.insert(normal);
                }
            }
        }

        // Map planes to the positions of atoms that lie on them.
        std::unordered_map<PlaneIdentifier, std::vector<Vec3>, PlaneIdentifierHasher, PlaneIdentifierEqual> all_planes(
            10, PlaneIdentifierHasher(tol), PlaneIdentifierEqual(tol)
        );
        for (const auto& normal : planeNormals)
        {
            for (const auto& atom : group)
            {
                const Vec3& p   = atom.get().coords;
                double distance = normal.dot(p);

                all_planes[{normal, distance}].push_back(p);
            }
        }

        // Check each planar group for a regular polygon.
        for (const auto& [planeID, pointsOnPlane] : all_planes)
        {
            if (isRegularPolygon(pointsOnPlane, planeID.distance * planeID.normal, planeID.normal, tol))
                candidateAxes.insert(planeID.normal);
        }
    }

    // Add a normal to a plane that contains the center of mass and any two non-collinear atoms.
    // This handles planar Cnh groups.
    Vec3 atom1 = Vec3(0, 0, 0);
    for (const auto& group : SEAGroups)
    {
        auto nonZeroAtom = std::find_if(
            group.begin(),
            group.end(),
            [&tol](const std::reference_wrapper<const Atom>& atom)
            { return !areEqual(atom.get().coords, Vec3(0, 0, 0), tol); }
        );

        if (nonZeroAtom != group.end())
        {
            atom1 = nonZeroAtom->get().coords;
            break;
        }
    }
    for (const auto& group : SEAGroups)
    {
        auto atom2 = std::find_if(
            group.begin(),
            group.end(),
            [&atom1, &tol](const std::reference_wrapper<const Atom>& atom)
            {
                return (!areEqual(atom.get().coords, Vec3(0, 0, 0), tol) && !areCollinear(atom.get().coords, atom1, tol));
            }
        );

        if (atom2 != group.end())
        {
            Vec3 normal = atom1.cross(atom2->get().coords);
            candidateAxes.insert(normal);
            break;
        }
    }


    return {candidateAxes.begin(), candidateAxes.end()};
}

std::vector<SymmetryOperation> findProperRotations(
    const std::vector<Atom>& geometry,
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups,
    const std::vector<Vec3>& candidateAxes,
    double tol
)
{
    // Validate the candidate axes. The proper rotation axis order is limited by the
    // greatest common divisor of the number of off-axis atoms from each symmetrically equivalent group.
    // The order must be a divisor of the GCD of the off-axis counts.
    std::vector<SymmetryOperation> properRotations;
    for (const auto& axis : candidateAxes)
    {
        // get number of off-axis atoms from each symmetrically equivalent group
        std::vector<size_t> offAxisCounts;
        for (const auto& group : SEAGroups)
        {
            size_t count = std::count_if(
                group.begin(),
                group.end(),
                [&axis, &tol](const std::reference_wrapper<const Atom>& atom)
                {
                    return (!areEqual(atom.get().coords, Vec3(0, 0, 0), tol) && !areCollinear(atom.get().coords, axis, tol));
                }
            );
            if (count == 0)
                continue;
            offAxisCounts.push_back(count);
        }

        // get the greatest common divisor of the off-axis counts
        size_t offAxisGCD = std::accumulate(
            offAxisCounts.begin(),
            offAxisCounts.end(),
            0,
            [](size_t currentGCD, size_t count) { return std::gcd(currentGCD, count); }
        );

        for (size_t order = offAxisGCD; order > 1; --order)
        {
            if (offAxisGCD % order != 0)
                continue;

            double angle                   = 2.0 * M_PI / order;
            Eigen::Matrix3d rotationMatrix = createRotationMatrix(axis, angle);

            if (isSymmetryOperation(geometry, rotationMatrix, tol))
            {
                properRotations.emplace_back(SymmetryOperation::Cn, rotationMatrix, axis, order);

                // Add all sub-orders of the rotation (i.e, C5^2, C5^3, etc.)
                for (size_t i = 2; i < order; ++i)
                {
                    int gcd = std::gcd(order, i);

                    rotationMatrix = createRotationMatrix(axis, i * angle);
                    properRotations.emplace_back(SymmetryOperation::Cn, rotationMatrix, axis, order / gcd, i / gcd);
                }

                break;
            }
        }
    }

    return properRotations;
}

std::vector<SymmetryOperation> findImproperRotations(
    const std::vector<Atom>& geometry,
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups,
    const std::vector<Vec3>& candidateAxes,
    double tol
)
{
    // Validate the candidate axes. The improper rotation axis order is limited by the
    // greatest common divisor of the number of off-axis atoms from each symmetrically equivalent group.
    // The order must be a divisor of the GCD of the off-axis counts.
    std::vector<SymmetryOperation> improperRotations;
    for (const auto& axis : candidateAxes)
    {
        // get number of off-axis atoms from each symmetrically equivalent group
        std::vector<size_t> offAxisCounts;
        for (const auto& group : SEAGroups)
        {
            size_t count = std::count_if(
                group.begin(),
                group.end(),
                [&axis, &tol](const std::reference_wrapper<const Atom>& atom)
                {
                    return (!areEqual(atom.get().coords, Vec3(0, 0, 0), tol) && !areCollinear(atom.get().coords, axis, tol));
                }
            );
            if (count == 0)
                continue;
            offAxisCounts.push_back(count);
        }

        // get the greatest common divisor of the off-axis counts
        size_t offAxisGCD = std::accumulate(
            offAxisCounts.begin(),
            offAxisCounts.end(),
            0,
            [](size_t currentGCD, size_t count) { return std::gcd(currentGCD, count); }
        );

        for (size_t order = 3; order <= offAxisGCD; ++order)
        {
            if (offAxisGCD % order != 0)
                continue;

            double angle                           = 2.0 * M_PI / order;
            Eigen::Matrix3d improperRotationMatrix = createImproperRotationMatrix(axis, angle);

            if (isSymmetryOperation(geometry, improperRotationMatrix, tol))
            {
                improperRotations.emplace_back(SymmetryOperation::Sn, improperRotationMatrix, axis, order);
                // Add all sub-orders of the improper rotation.
                for (size_t i = 3; i < (order % 2 == 0 ? order : 2 * order); ++i)
                {
                    if (order % i == 0 || i % 2 == 0)
                        continue;

                    int gcd                = std::gcd(order, i);
                    improperRotationMatrix = createImproperRotationMatrix(axis, i * angle);
                    improperRotations.emplace_back(SymmetryOperation::Sn, improperRotationMatrix, axis, order / gcd, i / gcd);
                }
            }
        }
    }

    return improperRotations;
}


std::vector<Vec3> getReflectionPlaneNormalsCandidates(
    const std::vector<std::vector<std::reference_wrapper<const Atom>>>& SEAGroups, double tol
)
{
    std::vector<Vec3> candidateNormals;

    UniqueVec3Set uniqueNormals(10, tol);

    // add normals to planes bisecting symmetrically equivalent atom pairs
    for (const auto& group : SEAGroups)
    {
        if (group.size() < 2)
            continue;

        for (size_t i = 0; i < group.size(); ++i)
        {
            for (size_t j = i + 1; j < group.size(); ++j)
            {
                uniqueNormals.insert(group[i].get().coords - group[j].get().coords);
            }
        }
    }

    // Add any single plane that containes the center of mass and two non-collinear atoms,
    // this handles planar molcules with no symmetrically equivalent atoms.
    Vec3 atom1 = Vec3(0, 0, 0);
    for (const auto& group : SEAGroups)
    {
        auto nonZeroAtom = std::find_if(
            group.begin(),
            group.end(),
            [&tol](const std::reference_wrapper<const Atom>& atom)
            { return !areEqual(atom.get().coords, Vec3(0, 0, 0), tol); }
        );

        if (nonZeroAtom != group.end())
        {
            atom1 = nonZeroAtom->get().coords;
            break;
        }
    }
    for (const auto& group : SEAGroups)
    {
        auto atom2 = std::find_if(
            group.begin(),
            group.end(),
            [&atom1, &tol](const std::reference_wrapper<const Atom>& atom)
            {
                return (!areEqual(atom.get().coords, Vec3(0, 0, 0), tol) && !areCollinear(atom.get().coords, atom1, tol));
            }
        );

        if (atom2 != group.end())
        {
            Vec3 normal = atom1.cross(atom2->get().coords);
            uniqueNormals.insert(normal);
            break;
        }
    }

    return {uniqueNormals.begin(), uniqueNormals.end()};
}

std::vector<SymmetryOperation> findReflectionPlanes(
    const std::vector<Atom>& geometry, const std::vector<Vec3>& candidateNormals, double tol
)
{
    std::vector<SymmetryOperation> reflections;
    for (const auto& normal : candidateNormals)
    {
        Eigen::Matrix3d reflectionMatrix = createReflectionMatrix(normal);

        if (isSymmetryOperation(geometry, reflectionMatrix, tol))
            reflections.emplace_back(SymmetryOperation::sigma, reflectionMatrix, normal);
    }

    return reflections;
}

void classifyReflectionPlanes(
    const std::vector<Atom>& geometry,
    std::vector<SymmetryOperation>& reflectionPlanes,
    const std::vector<SymmetryOperation>& properRotations,
    const std::vector<SymmetryOperation>& improperRotations,
    double tol
)
{
    if (reflectionPlanes.empty())
        return;

    if (properRotations.empty() && !reflectionPlanes.empty())
    {
        reflectionPlanes[0].name = SymmetryOperation::sigma_h;
        return;
    }

    // find the principal rotation axis
    auto highestN = [](const SymmetryOperation& r1, const SymmetryOperation& r2)
    {
        if (r1.power != 1)
            return true;
        else if (r2.power != 1)
            return false;
        return r1.n < r2.n;
    };

    auto highestCn            = std::max_element(properRotations.begin(), properRotations.end(), highestN);
    auto principalRotation    = highestCn;
    bool useImproperRotations = false;
    if (!improperRotations.empty())
    {
        auto highestSn = std::max_element(improperRotations.begin(), improperRotations.end(), highestN);
        if (highestSn->n > highestCn->n)
        {
            principalRotation    = highestSn;
            useImproperRotations = true;
        }
    }

    size_t principalRotationOrder     = principalRotation->n;
    const Vec3& principalRotationAxis = principalRotation->element;

    // check if multiple principal rotation axes are present
    bool hasMultiplePrincipalAxes = false;

    const std::vector<SymmetryOperation>& rotationsOfInterset = (useImproperRotations) ? improperRotations : properRotations;
    auto secondPrincipalAxis = std::find_if(
        rotationsOfInterset.begin(),
        rotationsOfInterset.end(),
        [&principalRotationAxis, &principalRotationOrder, &tol](const SymmetryOperation& rotation)
        {
            return (
                rotation.n == principalRotationOrder && rotation.power == 1
                && !areCollinear(rotation.element, principalRotationAxis, tol)
            );
        }
    );
    if (secondPrincipalAxis != rotationsOfInterset.end())
        hasMultiplePrincipalAxes = true;

    // Get C2 axes perpendicular to the principal rotation axis.
    std::vector<SymmetryOperation> perpendicularC2s;
    std::copy_if(
        properRotations.begin(),
        properRotations.end(),
        std::back_inserter(perpendicularC2s),
        [&principalRotationAxis, &tol](const SymmetryOperation& r)
        { return (r.n == 2 && r.power == 1 && areOrthogonal(r.element, principalRotationAxis, tol)); }
    );

    // Find the C2' axes, which are the ones that go through the most atoms or bonds.
    std::vector<SymmetryOperation> CPrimesAxes;
    auto C2Prime = std::max_element(
        perpendicularC2s.begin(),
        perpendicularC2s.end(),
        [&geometry, &tol](const SymmetryOperation& a, const SymmetryOperation& b)
        {
            size_t atomsIntersected_a = std::count_if(
                geometry.begin(),
                geometry.end(),
                [&a, &tol](const Atom& atom) { return areCollinear(atom.coords, a.element, tol); }
            );
            size_t atomsIntersected_b = std::count_if(
                geometry.begin(),
                geometry.end(),
                [&a, &tol](const Atom& atom) { return areCollinear(atom.coords, a.element, tol); }
            );
            if (atomsIntersected_a != atomsIntersected_b)
                return (atomsIntersected_a < atomsIntersected_b);

            // If the number of atoms intersected is the same, we compare the number of bonds intersected.
            size_t bondsIntersected_a = 0;
            size_t bondsIntersected_b = 0;
            for (size_t i = 0; i < geometry.size(); ++i)
            {
                for (size_t j = i + 1; j < geometry.size(); ++j)
                {
                    if (areCollinear(geometry[i].coords + geometry[j].coords, a.element, tol))
                        bondsIntersected_a++;
                    if (areCollinear(geometry[i].coords + geometry[j].coords, b.element, tol))
                        bondsIntersected_b++;
                }
            }
            return (bondsIntersected_a < bondsIntersected_b);
        }
    );
    if (C2Prime != perpendicularC2s.end())
    {
        Vec3 C2PrimeGenerator = C2Prime->element;

        std::copy_if(
            perpendicularC2s.begin(),
            perpendicularC2s.end(),
            std::back_inserter(CPrimesAxes),
            [&principalRotationAxis, &principalRotationOrder, &C2PrimeGenerator, &tol](const SymmetryOperation& c2)
            {
                for (size_t i = 0; i < principalRotationOrder; ++i)
                {
                    Eigen::Matrix3d rotationMatrix = createRotationMatrix(
                        principalRotationAxis, i * 2.0 * M_PI / principalRotationOrder
                    );
                    if (areCollinear(rotationMatrix * c2.element, C2PrimeGenerator, tol))
                        return true;
                }
                return false;
            }
        );
    }


    // Special case for Cnv where n > 2 and n is even.
    if (principalRotationOrder > 2 && principalRotationOrder % 2 == 0 && !hasMultiplePrincipalAxes
        && perpendicularC2s.empty() && !hasInversion(geometry, tol))
    {
        // sigma_v planes go through highest number of atoms
        size_t maxAtomsIntersected = 0;
        for (const auto& reflection : reflectionPlanes)
        {
            size_t atomsIntersected = std::count_if(
                geometry.begin(),
                geometry.end(),
                [&reflection, &tol](const Atom& atom) { return areOrthogonal(atom.coords, reflection.element, tol); }
            );
            maxAtomsIntersected = std::max(maxAtomsIntersected, atomsIntersected);
        }
        for (auto& reflection : reflectionPlanes)
        {
            size_t atomsIntersected = std::count_if(
                geometry.begin(),
                geometry.end(),
                [&reflection, &tol](const Atom& atom) { return areOrthogonal(atom.coords, reflection.element, tol); }
            );
            if (atomsIntersected == maxAtomsIntersected)
                reflection.name = SymmetryOperation::sigma_v;
            else
                reflection.name = SymmetryOperation::sigma_d;
        }
        return;
    }

    for (auto& reflection : reflectionPlanes)
    {
        // Check if the reflection plane is a sigma_h plane.
        if (hasMultiplePrincipalAxes)
        {
            if (highestCn->n == 3)
            {
                auto orthogonalC2 = std::find_if(
                    properRotations.begin(),
                    properRotations.end(),
                    [&reflection, &tol](const SymmetryOperation& rotation)
                    { return (rotation.n == 2 && areCollinear(reflection.element, rotation.element, tol)); }
                );
                if (orthogonalC2 != properRotations.end())
                {
                    reflection.name = SymmetryOperation::sigma_h;
                    continue;
                }
            }
            else if (highestCn->n == 4)
            {
                auto orthogonalC4 = std::find_if(
                    properRotations.begin(),
                    properRotations.end(),
                    [&reflection, &tol](const SymmetryOperation& rotation)
                    {
                        return (
                            rotation.n == 4 && rotation.power == 1 && areCollinear(reflection.element, rotation.element, tol)
                        );
                    }
                );
                if (orthogonalC4 != properRotations.end())
                {
                    reflection.name = SymmetryOperation::sigma_h;
                    continue;
                }
            }
            else
            {
                // Ih point group has sigma planes that are not defined as v, d or h planes.
                continue;
            }
        }
        else if (areCollinear(reflection.element, principalRotationAxis, tol))
        {
            reflection.name = SymmetryOperation::sigma_h;
            continue;
        }

        // Check if the reflection plane is a sigma_d or sigma_v plane.
        if (!hasMultiplePrincipalAxes && areOrthogonal(reflection.element, principalRotationAxis, tol))
        {
            auto containedC2Prime = std::find_if(
                CPrimesAxes.begin(),
                CPrimesAxes.end(),
                [&reflection, &tol](const SymmetryOperation& cPrime)
                { return areOrthogonal(reflection.element, cPrime.element, tol); }
            );
            if (containedC2Prime != CPrimesAxes.end()) // If the reflection plane contains a C2' axis it is a sigma_v plane.
            {
                reflection.name = SymmetryOperation::sigma_v;
                continue;
            }

            bool isSigmaD = false;
            for (size_t i = 0; i < perpendicularC2s.size(); ++i)
            {
                for (size_t j = i + 1; j < perpendicularC2s.size(); ++j)
                {
                    Vec3 bisector1 = (perpendicularC2s[i].element + perpendicularC2s[j].element).normalized();
                    Vec3 bisector2 = (perpendicularC2s[i].element - perpendicularC2s[j].element).normalized();
                    if (areOrthogonal(reflection.element, bisector1, tol)
                        || areOrthogonal(reflection.element, bisector2, tol))
                    {
                        reflection.name = SymmetryOperation::sigma_d;
                        isSigmaD        = true;
                        break;
                    }
                }
                if (isSigmaD)
                    break;
            }

            if (!isSigmaD)
                reflection.name = SymmetryOperation::sigma_v;
        }
        else if (hasMultiplePrincipalAxes) // handle spherical top molecules
        {
            bool isSigmaD = false;
            for (size_t i = 0; i < properRotations.size(); ++i)
            {
                if (properRotations[i].n != 2)
                    continue;
                for (size_t j = i + 1; j < properRotations.size(); ++j)
                {
                    if (properRotations[j].n != 2)
                        continue;
                    Vec3 bisector1 = (properRotations[i].element + properRotations[j].element).normalized();
                    Vec3 bisector2 = (properRotations[i].element - properRotations[j].element).normalized();
                    if (areOrthogonal(reflection.element, bisector1, tol)
                        || areOrthogonal(reflection.element, bisector2, tol))
                    {
                        reflection.name = SymmetryOperation::sigma_d;
                        isSigmaD        = true;
                        break;
                    }
                }
                if (isSigmaD)
                    break;
            }
            if (!isSigmaD)
                reflection.name = SymmetryOperation::sigma_v;
        }
    }
}

PointGroup classifyPointGroup(const std::vector<Atom>& geometry, const Vec3& principalMoments, double tol)
{
    auto sortedMoments = principalMoments;
    std::sort(sortedMoments.begin(), sortedMoments.end());

    auto SEAGroups = findSEAGroups(geometry, tol);

    std::vector<Vec3> candidateRotationAxes           = getRotationAxesCandidates(SEAGroups, tol);
    std::vector<Vec3> candidateReflectionPlaneNormals = getReflectionPlaneNormalsCandidates(SEAGroups, tol);

    std::vector<SymmetryOperation> properRotations = findProperRotations(geometry, SEAGroups, candidateRotationAxes, tol);
    std::vector<SymmetryOperation> improperRotations = findImproperRotations(geometry, SEAGroups, candidateRotationAxes, tol);
    std::vector<SymmetryOperation> reflectionPlanes = findReflectionPlanes(geometry, candidateReflectionPlaneNormals, tol);

    classifyReflectionPlanes(geometry, reflectionPlanes, properRotations, improperRotations, tol);

    bool hasInv = hasInversion(geometry, tol);

    std::vector<SymmetryOperation> allOperations;
    allOperations.reserve(properRotations.size() + improperRotations.size() + reflectionPlanes.size() + 2);
    allOperations.push_back(SymmetryOperation::Identity());
    if (hasInversion(geometry, tol))
        allOperations.push_back(SymmetryOperation::Inversion());
    allOperations.insert(allOperations.end(), properRotations.begin(), properRotations.end());
    allOperations.insert(allOperations.end(), improperRotations.begin(), improperRotations.end());
    allOperations.insert(allOperations.end(), reflectionPlanes.begin(), reflectionPlanes.end());

    PointGroup::PointGroupName pointGroupName = PointGroup::C1;

    // Ia = 0, Ib = Ic: linear
    if (std::abs(sortedMoments[0]) < tol * 10 && std::abs(sortedMoments[1] - sortedMoments[2]) < tol * 10)
    {
        if (hasInv)
            pointGroupName = PointGroup::Dinfh;
        else
            pointGroupName = PointGroup::Cinfh;
        return PointGroup(pointGroupName, 1, allOperations);
    }

    bool hasRotation   = properRotations.size() > 0;
    bool hasSigmaPlane = reflectionPlanes.size() > 0;

    if (!hasRotation)
    {
        if (hasSigmaPlane)
            pointGroupName = PointGroup::Cs;
        else if (hasInv)
            pointGroupName = PointGroup::Ci;
        else
            pointGroupName = PointGroup::C1;
        return PointGroup(pointGroupName, 1, allOperations);
    }

    auto principalRotation = std::max_element(
        properRotations.begin(),
        properRotations.end(),
        [](const SymmetryOperation& r1, const SymmetryOperation& r2) { return r1.n < r2.n; }
    );
    const Vec3& principalRotationAxis = principalRotation->element;
    size_t principalRotationOrder     = principalRotation->n;

    auto sigmaHPlane = std::find_if(
        reflectionPlanes.begin(),
        reflectionPlanes.end(),
        [](const SymmetryOperation& r) { return (r.name == SymmetryOperation::sigma_h); }
    );
    bool hasSigmaHPlane = (sigmaHPlane != reflectionPlanes.end());

    auto sigmaVPlane = std::find_if(
        reflectionPlanes.begin(),
        reflectionPlanes.end(),
        [](const SymmetryOperation& r) { return (r.name == SymmetryOperation::sigma_v); }
    );
    bool hasSigmaVPlane = (sigmaVPlane != reflectionPlanes.end());

    size_t numOfSigmaDPlanes = std::accumulate(
        reflectionPlanes.begin(),
        reflectionPlanes.end(),
        0,
        [](size_t acc, const SymmetryOperation& r) { return (r.name == SymmetryOperation::sigma_d ? acc + 1 : acc); }
    );

    // Ia = Ib = Ic: spherical top
    if (areEqual(sortedMoments[0], sortedMoments[1], tol * 10) && areEqual(sortedMoments[1], sortedMoments[2], tol * 10))
    {
        // tetrahedral, octahedral, or icosahedral
        size_t numOfC5Axes = 0;
        size_t numOfC4Axes = 0;
        size_t numOfC3Axes = 0;
        for (const auto& rotation : properRotations)
        {
            if (rotation.n == 5 && rotation.power == 1)
                numOfC5Axes++;
            if (rotation.n == 4 && rotation.power == 1)
                numOfC4Axes++;
            if (rotation.n == 3 && rotation.power == 1)
                numOfC3Axes++;
        }

        if (numOfC5Axes == 6)
        {
            if (hasInv)
                pointGroupName = PointGroup::Ih;
            else
                pointGroupName = PointGroup::I;
        }
        else if (numOfC4Axes == 3)
        {
            if (hasInv)
                pointGroupName = PointGroup::Oh;
            else
                pointGroupName = PointGroup::O;
        }
        else if (numOfC3Axes == 4)
        {
            if (hasInv)
                pointGroupName = PointGroup::Th;
            else if (numOfSigmaDPlanes == 6)
                pointGroupName = PointGroup::Td;
            else
                pointGroupName = PointGroup::T;
        }
        return PointGroup(pointGroupName, principalRotationOrder, allOperations);
    }

    // Ia = Ib != Ic
    else if (areEqual(sortedMoments[0], sortedMoments[1], tol * 10)
             || areEqual(sortedMoments[1], sortedMoments[2], tol * 10))
    {
        size_t numOfPerpendicularC2Axes = std::count_if(
            properRotations.begin(),
            properRotations.end(),
            [&principalRotationAxis, &tol](const SymmetryOperation& rotation)
            { return (rotation.n == 2 && areOrthogonal(rotation.element, principalRotationAxis, tol)); }
        );

        // symmetric top
        if (numOfPerpendicularC2Axes == principalRotationOrder)
        {
            if (hasSigmaHPlane)
                pointGroupName = PointGroup::Dnh;
            else if (numOfSigmaDPlanes == principalRotationOrder)
                pointGroupName = PointGroup::Dnd;
            else
                pointGroupName = PointGroup::Dn;
        }
        else
        {
            if (hasSigmaHPlane)
                pointGroupName = PointGroup::Cnh;
            else if (hasSigmaVPlane)
                pointGroupName = PointGroup::Cnv;
            else if (!improperRotations.empty())
                pointGroupName = PointGroup::S2n;
            else
                pointGroupName = PointGroup::Cn;
        }

        return PointGroup(pointGroupName, principalRotationOrder, allOperations);
    }

    // Ia != Ib != Ic: asymmetric top
    size_t numOfC2Axes = std::count_if(
        properRotations.begin(),
        properRotations.end(),
        [](const SymmetryOperation& rotation) { return (rotation.n == 2); }
    );

    if (numOfC2Axes == 3)
    {
        if (!reflectionPlanes.empty())
            pointGroupName = PointGroup::Dnh;
        else
            pointGroupName = PointGroup::Dn;
    }
    else
    {
        if (hasSigmaHPlane)
            pointGroupName = PointGroup::Cnh;
        else if (hasSigmaVPlane)
            pointGroupName = PointGroup::Cnv;
        else
            pointGroupName = PointGroup::Cn;
        // We alrady handled C1, Cs, Ci and S2n
    }
    return PointGroup(pointGroupName, 2, allOperations);
}


std::vector<std::vector<SymmetryOperation>> getClasses(const std::vector<SymmetryOperation>& operations, double tol)
{
    std::vector<std::vector<SymmetryOperation>> classes;
    std::vector<bool> isClassified(operations.size(), false);

    for (size_t i = 0; i < operations.size(); ++i)
    {
        if (isClassified[i])
            continue;

        std::vector<SymmetryOperation> operationsClass;
        operationsClass.push_back(operations[i]);
        isClassified[i] = true;

        auto A = operations[i].matrix;

        for (size_t j = 0; j < operations.size(); ++j)
        {
            auto X = operations[j].matrix;
            auto B = X.transpose() * A * X;
            for (size_t k = i + 1; k < operations.size(); ++k)
            {
                if (!isClassified[k] && areEqual(B, operations[k].matrix, tol))
                {
                    operationsClass.push_back(operations[k]);
                    isClassified[k] = true;
                    break;
                }
            }
        }

        classes.push_back(operationsClass);
    }

    return classes;
}


std::vector<std::string> nameClasses(
    const std::vector<std::vector<SymmetryOperation>>& classes, const PointGroup& pointGroup, double tol
)
{
    std::vector<std::string> classNames;

    for (const auto& operationsClass : classes)
    {
        size_t classSize                 = operationsClass.size();
        SymmetryOperation representative = operationsClass[0];
        for (const auto& op : operationsClass)
        {
            if (op.n < representative.n || (op.n == representative.n && op.power < representative.power))
                representative = op;
        }
        std::string className = (classSize > 1 ? std::to_string(classSize) : "") + representative.getName();

        // handle some special cases
        if (pointGroup.name == PointGroup::Cnv && pointGroup.n == 2 && className == "sigmav") // C2v
        {
            if (areCollinear(representative.element, Vec3(1, 0, 0), tol))
                className = "sigmayz";
            else if (areCollinear(representative.element, Vec3(0, 1, 0), tol))
                className = "sigmaxz";
            else if (areCollinear(representative.element, Vec3(0, 0, 1), tol))
                className = "sigmaxy";
        }
        else if (pointGroup.name == PointGroup::Dnh && pointGroup.n == 2 && (className == "sigma")) // D2h
        {
            if (areCollinear(representative.element, Vec3(0, 0, 1), tol))
                className = "sigmaxy";
            else if (areCollinear(representative.element, Vec3(0, 1, 0), tol))
                className = "sigmaxz";
            else if (areCollinear(representative.element, Vec3(1, 0, 0), tol))
                className = "sigmayz";
        }

        if ((pointGroup.name == PointGroup::Dn || pointGroup.name == PointGroup::Dnh || pointGroup.name == PointGroup::Dnd)
            && representative.name == SymmetryOperation::Cn && representative.n == 2) // Dn
        {
            if (pointGroup.n == 2 && pointGroup.name != PointGroup::Dnd) // D2 or D2h
            {
                if (areCollinear(representative.element, Vec3(1, 0, 0), tol))
                    className = "C2x";
                else if (areCollinear(representative.element, Vec3(0, 1, 0), tol))
                    className = "C2y";
                else if (areCollinear(representative.element, Vec3(0, 0, 1), tol))
                    className = "C2z";
            }
            else // Dn or D2h, n != 2
            {
                bool isC2Prime = false;
                for (const auto& c2 : operationsClass)
                {
                    if (areEqual(c2.element, Vec3(1, 0, 0), tol))
                    {
                        className += "'";
                        isC2Prime = true;
                        break;
                    }
                }

                if (!isC2Prime && !areEqual(representative.element, Vec3(0, 0, 1), tol))
                    className += "''";
            }
        }

        if (pointGroup.name == PointGroup::O && className == "6C2") // O
            className = "6C2'";

        // if (pointGroup.name == PointGroup::Ih && className == "15sigmad") // Ih
        //     className = "15sigma";


        classNames.push_back(className);
    }

    return classNames;
}


bool areGeometriesSuperimposable(const std::vector<Atom>& geometry1, const std::vector<Atom>& geometry2, double tol)
{
    if (geometry1.size() != geometry2.size())
        return false;

    std::vector<bool> matched(geometry2.size(), false);

    for (const auto& atom1 : geometry1)
    {
        bool foundMatch = false;

        for (size_t j = 0; j < geometry2.size(); ++j)
        {
            if (matched[j])
                continue;

            const auto& atom2 = geometry2[j];

            if (atom1.atomicNumber != atom2.atomicNumber)
                continue;

            if (areEqual(atom1.coords, atom2.coords, tol))
            {
                matched[j] = true;
                foundMatch = true;
                break;
            }
        }

        if (!foundMatch)
            return false;
    }

    return true;
}


Eigen::MatrixXd getAOBasisRep(const SymmetryOperation& op, const std::vector<AtomicOrbital>& aos, double tol)
{
    size_t n_aos             = aos.size();
    Eigen::MatrixXd D        = Eigen::MatrixXd::Zero(n_aos, n_aos);
    const Eigen::Matrix3d& R = op.matrix;

    for (size_t i = 0; i < n_aos; ++i)
    {
        const AtomicOrbital& ao                = aos[i];
        const Vec3& center                     = ao.center;
        const Eigen::Vector3i& angularMomentum = ao.angularMomentum;

        Vec3 newCenter = R * center;

        // get indicies of atomic orbitals that are on the new center after the symmetry operation
        std::vector<size_t> AOsOnNewCenter;
        for (size_t j = 0; j < n_aos; ++j)
        {
            if (areEqual(aos[j].center, newCenter, tol))
                AOsOnNewCenter.push_back(j);
        }

        if (AOsOnNewCenter.empty())
            continue;

        int L = angularMomentum.sum();

        std::vector<Vec3> components;
        if (L > 0)
        {
            components.reserve(L);
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < angularMomentum[k]; ++l) { components.emplace_back(R * Eigen::Vector3d::Unit(k)); }
            }
        }

        for (size_t j : AOsOnNewCenter)
        {
            if (!AtomicOrbital::sameSubshell(aos[j], ao))
                continue;

            const Eigen::Vector3i& targetAM = aos[j].angularMomentum;

            if (L == 0)
                D(j, i) = 1.0;
            else
            {
                // recursive function to calculate the coefficient
                std::function<double(size_t, const Eigen::Vector3i&, double)> calculateCartesianCoefficient =
                    [&](size_t currentIndex, const Eigen::Vector3i& currentAMCounts, double currentProduct) -> double
                {
                    // Base case: All compnents have been processed
                    if (currentIndex == components.size())
                    {
                        // If the accumulated angular momentum matches the target, return the product.
                        if (currentAMCounts == targetAM)
                            return currentProduct;
                        else
                            return 0.0; // Otherwise, this path doesn't contribute.
                    }

                    double sumOfContributions = 0.0;

                    const Vec3& comp = components[currentIndex];

                    // Try to assign the current compnent's x-coordinate
                    if (currentAMCounts.x() < targetAM.x())
                    {
                        sumOfContributions += calculateCartesianCoefficient(
                            currentIndex + 1, currentAMCounts + Eigen::Vector3i(1, 0, 0), currentProduct * comp.x()
                        );
                    }

                    // Try to assign the current component's y-coordinate
                    if (currentAMCounts.y() < targetAM.y())
                    {
                        sumOfContributions += calculateCartesianCoefficient(
                            currentIndex + 1, currentAMCounts + Eigen::Vector3i(0, 1, 0), currentProduct * comp.y()
                        );
                    }

                    // Try to assign the current component's z-coordinate
                    if (currentAMCounts.z() < targetAM.z())
                    {
                        sumOfContributions += calculateCartesianCoefficient(
                            currentIndex + 1, currentAMCounts + Eigen::Vector3i(0, 0, 1), currentProduct * comp.z()
                        );
                    }

                    return sumOfContributions;
                };

                D(j, i) = calculateCartesianCoefficient(0, Eigen::Vector3i::Zero(), 1.0);
            }
        }
    }
    return D;
}


bool isSymmetryOperation(const std::vector<Atom>& geometry, const Eigen::Matrix3d& operation, double tol)
{
    std::vector<Atom> transformedGeometry;
    transformedGeometry.reserve(geometry.size());

    for (const auto& atom : geometry) { transformedGeometry.emplace_back(atom.atomicNumber, operation * atom.coords); }

    return areGeometriesSuperimposable(geometry, transformedGeometry, tol);
}


Eigen::Matrix3d createRotationMatrix(const Vec3& axis, double angle)
{
    Eigen::Matrix3d R;
    double s            = sin(angle);
    double c            = cos(angle);
    double t            = 1.0 - c;
    auto normalizedAxis = axis.normalized();
    double x            = normalizedAxis.x();
    double y            = normalizedAxis.y();
    double z            = normalizedAxis.z();

    R(0, 0) = c + x * x * t;
    R(0, 1) = x * y * t - z * s;
    R(0, 2) = y * s + x * z * t;
    R(1, 0) = z * s + x * y * t;
    R(1, 1) = c + y * y * t;
    R(1, 2) = -x * s + y * z * t;
    R(2, 0) = -y * s + x * z * t;
    R(2, 1) = x * s + z * y * t;
    R(2, 2) = c + z * z * t;

    return R;
}


Eigen::Matrix3d createReflectionMatrix(const Vec3& normal)
{
    return Eigen::Matrix3d::Identity() - 2.0 * normal.normalized() * normal.normalized().transpose();
}


Eigen::Matrix3d createImproperRotationMatrix(const Vec3& axis, double angle)
{
    Eigen::Matrix3d rotation   = createRotationMatrix(axis, angle);
    Eigen::Matrix3d reflection = createReflectionMatrix(axis);
    return reflection * rotation; // Combine rotation and inversion
}


bool hasInversion(const std::vector<Atom>& geometry, double tol)
{
    Eigen::Matrix3d inversionMatrix = -Eigen::Matrix3d::Identity();
    return isSymmetryOperation(geometry, inversionMatrix, tol);
}

} // namespace Symmetry
