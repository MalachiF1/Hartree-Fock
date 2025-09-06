#include "Symmetry/Symmetry.hpp"

#include "Eigen/Core"
#include "Symmetry/AlignGeometry.hpp"
#include "Symmetry/Helpers.hpp"

#include <algorithm>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <map>

namespace Symmetry
{

std::vector<SEA> findSEAs(const std::vector<Atom>& geometry, double tol)
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

    std::map<Fingerprint, std::vector<size_t>, decltype(compareFingerprints)> groups(compareFingerprints);

    for (size_t i = 0; i < n; ++i)
    {
        Fingerprint fingerprint;
        fingerprint.reserve(n);
        for (size_t j = 0; j < n; ++j)
        {
            double distance = (geometry[i].coords - geometry[j].coords).norm();
            fingerprint.emplace_back(distance, geometry[j].atomicNumber);
        }

        std::ranges::sort(
            fingerprint,
            [&tol](const auto& a, const auto& b)
            {
                if (areEqual(a.first, b.first, tol))
                    return a.second < b.second; // Sort by atomic number if distances are equal
                return a.first < b.first;
            }
        );

        groups[fingerprint].push_back(i);
    }

    // Collect the final groups from the map's values.
    std::vector<SEA> SEAs;
    for (const auto& [fingerprint, indices] : groups)
    {
        std::vector<Atom> SEAAtoms;
        SEAAtoms.reserve(indices.size());
        for (const auto& index : indices) { SEAAtoms.push_back(geometry[index]); }

        // Calculate center of the SEAs
        Vec3 center = Vec3::Zero();
        for (const auto& atom : SEAAtoms) { center += atom.coords; }
        center /= static_cast<double>(SEAAtoms.size());

        // Calculate inertia tensor of the SEA
        translateOrigin(SEAAtoms, center);
        InertiaTensor inertiaTensor = diagonalizeInertiaTensor(SEAAtoms);

        SEAs.push_back(SEA {.indices = indices, .center = center, .inertiaTensor = inertiaTensor});
    }

    return SEAs;
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


bool isSymmetryOperation(const std::vector<Atom>& geometry, const Eigen::Matrix3d& operation, double tol)
{
    std::vector<Atom> transformedGeometry;
    transformedGeometry.reserve(geometry.size());

    for (const auto& atom : geometry) { transformedGeometry.emplace_back(atom.atomicNumber, operation * atom.coords); }

    return areGeometriesSuperimposable(geometry, transformedGeometry, tol);
}


std::vector<std::pair<Vec3, size_t>> findRotations(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol)
{
    // Each inner vector contains rotation ax3s within a SEA, along with the order of rotation (0 = any order).
    std::vector<std::vector<std::pair<Vec3, size_t>>> candidateAxesPerSEA;

    for (const auto& sea : SEAs)
    {
        size_t k = sea.indices.size();

        auto [tensor, moments, axes] = sea.inertiaTensor;

        double Ia = moments[0], Ib = moments[1], Ic = moments[2];
        Vec3 Iav = axes.col(0), Ibv = axes.col(1), Icv = axes.col(2);

        if ((areEqual(Ia, 0.0, tol) || areEqual(Ia, Ib, tol)) && areEqual(Ib, Ic, tol)) // linear/spherical
            continue;

        if (k == 1 && sea.center.squaredNorm() > tol * tol)
        {
            candidateAxesPerSEA.push_back({{sea.center, 0}});
            continue;
        }

        if (areEqual(Ia + Ib, Ic, tol)) // polygonal
        {
            bool regular = areEqual(Ia, Ib, tol);
            std::vector<std::pair<Vec3, size_t>> axes;
            std::vector<size_t> divisors = findDivisors(k);

            Vec3 axis = Icv.normalized();
            if (regular)
                axes.emplace_back(axis.normalized(), k);

            for (size_t n : divisors) { axes.emplace_back(axis, n); };

            candidateAxesPerSEA.push_back(axes);
        }
        else // polyhedral
        {
            if (!areEqual(Ia, Ib, tol) && !areEqual(Ib, Ic, tol)) // asymmetric top
            {
                candidateAxesPerSEA.push_back({{Iav.normalized(), 2}, {Ibv.normalized(), 2}, {Icv.normalized(), 2}});
            }
            else
            {
                Vec3 axis = (areEqual(Ia, Ib, tol)) ? Icv.normalized() : Iav.normalized(); // oblate or prolate symmetric top

                std::vector<std::pair<Vec3, size_t>> axes;
                axes.emplace_back(axis, k / 2);

                std::vector<size_t> divisors = findDivisors(k / 2);
                for (size_t n : divisors) { axes.emplace_back(axis, n); };

                candidateAxesPerSEA.push_back(axes);
            }
        }
    }

    if (candidateAxesPerSEA.empty())
        return {};

    // Find common axes among all SEAs
    std::vector<std::pair<Vec3, size_t>> commonAxes = candidateAxesPerSEA.front();
    for (size_t i = 1; i < candidateAxesPerSEA.size(); ++i)
    {
        std::vector<std::pair<Vec3, size_t>> newCommonAxes;
        for (const auto& [axis1, n1] : commonAxes)
        {
            for (const auto& [axis2, n2] : candidateAxesPerSEA[i])
            {
                if (areCollinear(axis1, axis2, tol))
                {
                    if (n1 == n2)
                    {
                        newCommonAxes.emplace_back(axis1, n1);
                    }
                    else if (n1 == 0 && n2 == 0)
                    {
                        newCommonAxes.emplace_back(axis1, 0);
                    }
                    else if (n2 == 0 || n1 == 0)
                    {
                        newCommonAxes.emplace_back(axis1, n1);
                        newCommonAxes.emplace_back(axis2, n2);
                    }
                }
            }
        }
        commonAxes = newCommonAxes;
        if (commonAxes.empty())
            break;
    }

    // Remove n = 0
    auto it = std::ranges::remove_if(commonAxes, [](const auto& pair) { return pair.second == 0; });
    commonAxes.erase(it.begin(), it.end());


    std::vector<std::pair<Vec3, size_t>> rotations;
    for (const auto& [axis, n] : commonAxes)
    {
        if (isSymmetryOperation(geometry, Cn(axis, n), tol))
            rotations.emplace_back(axis, n);
    }

    return rotations;
}


std::vector<Vec3> findC2s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol)
{
    UniqueVec3Set singleAtomAxes(10, tol);
    UniqueVec3Set linearSEAaxes(5, tol);
    for (const auto& sea : SEAs)
    {
        if (sea.indices.size() == 1)
        {
            singleAtomAxes.insert(sea.center); // handles zero vector and duplicates
            continue;
        }

        auto [tensor, moments, axes] = sea.inertiaTensor;
        double Ia                    = moments[0];
        double Ib                    = moments[1];
        double Ic                    = moments[2];

        if (areEqual(Ia, 0.0, tol) && areEqual(Ib, Ic, tol)) // linear
            linearSEAaxes.insert(axes.col(0));
    }

    // If more than 2 linear SEAs, they all must be coplanar for a C2 to exist.
    if (linearSEAaxes.size() > 2)
    {
        for (const auto& axis1 : linearSEAaxes)
        {
            for (const auto& axis2 : linearSEAaxes)
            {
                if (areCollinear(axis1, axis2, tol))
                    continue;

                Vec3 normal = axis1.cross(axis2);
                for (const auto& axis3 : linearSEAaxes)
                {
                    if (areCollinear(axis3, axis1, tol) || areCollinear(axis3, axis2, tol))
                        continue;

                    if (!areOrthogonal(normal, axis3, tol))
                        return {}; // Not coplanar, no C2s.
                }
            }
        }
    }

    // If there are at least two linear SEAs (and they're all coplanar), there can be only one C2 axis perpendicular to both.
    if (linearSEAaxes.size() > 1)
    {
        auto it    = linearSEAaxes.begin();
        Vec3 axis1 = *it;
        ++it;
        Vec3 axis2 = *it;

        Vec3 C2axis = axis1.cross(axis2).normalized();

        if (isSymmetryOperation(geometry, Cn(C2axis, 2), tol))
            return {C2axis};
        else
            return {};
    }

    // Collect candidate C2 axes
    UniqueVec3Set C2Candidates(10, tol);

    // Axes through atoms
    for (const auto& atom : geometry)
    {
        if (atom.coords.squaredNorm() < tol * tol)
            continue;

        // If there is a linear SEA, any C2 axis must be perpendicular to it.
        if (linearSEAaxes.size() == 1 && !areOrthogonal(atom.coords, *linearSEAaxes.begin(), tol))
            continue;

        // If the axis does not go through all single atom SEAs, it cannot be a C2 axis.
        bool passesSingleAtomSEAs = std::ranges::all_of(
            singleAtomAxes, [&atom, &tol](const Vec3& seaAxis) { return areCollinear(atom.coords, seaAxis, tol); }
        );
        if (!passesSingleAtomSEAs)
            continue;

        // hndles zero vector and duplicates
        C2Candidates.insert(atom.coords);
    }

    // Axes through midpoints of atom pairs in each SEA
    for (const auto& sea : SEAs)
    {
        if (sea.indices.size() < 2)
            continue;

        for (size_t i = 0; i < sea.indices.size(); ++i)
        {
            for (size_t j = i + 1; j < sea.indices.size(); ++j)
            {
                Vec3 axis = geometry[sea.indices[i]].coords + geometry[sea.indices[j]].coords;

                if (linearSEAaxes.size() == 1 && !areOrthogonal(axis, *linearSEAaxes.begin(), tol))
                    continue;

                // If the axis does not go through all single atom SEAs, it cannot be a C2 axis.
                bool passesSingleAtomSEAs = std::ranges::all_of(
                    singleAtomAxes, [&axis, &tol](const Vec3& seaAxis) { return areCollinear(axis, seaAxis, tol); }
                );
                if (!passesSingleAtomSEAs)
                    continue;

                // handles zero vector and duplicates
                C2Candidates.insert(axis);
            }
        }
    }

    // Check each candidate C2 axis
    std::vector<Vec3> C2Axes;
    for (const auto& axis : C2Candidates)
    {
        if (isSymmetryOperation(geometry, Cn(axis, 2), tol))
            C2Axes.push_back(axis);
    }

    return C2Axes;
}


std::vector<Vec3> findIcosahedralC3s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol)
{
    UniqueVec3Set C3Candidates(10, tol);

    auto isEqualateralTriangle = [](const Vec3& v1, const Vec3& v2, const Vec3& v3, double tol) -> bool
    {
        double d12 = (v1 - v2).norm();
        double d23 = (v2 - v3).norm();
        double d31 = (v3 - v1).norm();
        return areEqual(d12, d23, tol) && areEqual(d23, d31, tol);
    };

    for (const auto& sea : SEAs)
    {
        if (sea.indices.size() < 3)
            continue;

        for (size_t i = 0; i < sea.indices.size(); ++i)
        {
            for (size_t j = i + 1; j < sea.indices.size(); ++j)
            {
                for (size_t k = j + 1; k < sea.indices.size(); ++k)
                {
                    Vec3 v1 = geometry[sea.indices[i]].coords;
                    Vec3 v2 = geometry[sea.indices[j]].coords;
                    Vec3 v3 = geometry[sea.indices[k]].coords;

                    if (!isEqualateralTriangle(v1, v2, v3, tol))
                        continue;

                    C3Candidates.insert(v1 + v2 + v3);
                }
            }
        }
    }

    // Check each candidate C3 axis
    std::vector<Vec3> C3Axes;
    for (const auto& axis : C3Candidates)
    {
        if (isSymmetryOperation(geometry, Cn(axis, 3), tol))
            C3Axes.push_back(axis);
    }

    return C3Axes;
}


std::vector<Vec3> findOctahedralC4s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol)
{

    UniqueVec3Set C4Candidates(10, tol);

    for (const auto& sea : SEAs)
    {
        size_t size = sea.indices.size();
        if (size < 4) // looking for squares
            continue;

        const auto& idx = sea.indices;

        // find unique normals to planes that at least contain four atoms
        UniqueVec3Set planeNormals(10, tol);
        for (size_t i = 0; i < size; ++i)
        {
            const Vec3& p1 = geometry[idx[i]].coords;
            for (size_t j = i + 1; j < size; ++j)
            {
                const Vec3& p2 = geometry[idx[j]].coords;
                for (size_t k = j + 1; k < size; ++k)
                {
                    const Vec3& p3 = geometry[idx[k]].coords;

                    Vec3 normal = (p2 - p1).cross(p3 - p1);
                    for (size_t l = k + 1; l < size; ++l)
                    {
                        const Vec3& p4 = geometry[idx[l]].coords;
                        if (areOrthogonal(p4 - p1, normal, tol))
                        {
                            planeNormals.insert(normal);
                            break;
                        }
                    }
                }
            }
        }

        // Map planes to the positions of atoms that lie on them.
        std::unordered_map<PlaneIdentifier, std::vector<Vec3>, PlaneIdentifierHasher, PlaneIdentifierEqual> all_planes(
            10, PlaneIdentifierHasher(tol), PlaneIdentifierEqual(tol)
        );
        for (const auto& normal : planeNormals)
        {
            for (const auto& i : idx)
            {
                const Vec3& p   = geometry[i].coords;
                double distance = normal.dot(p);

                all_planes[{normal, distance}].push_back(p);
            }
        }

        // Check each planar group for a square.
        for (const auto& [planeID, pointsOnPlane] : all_planes)
        {
            if (pointsOnPlane.size() != 4)
                continue;

            if (isRegularPolygon(pointsOnPlane, planeID.distance * planeID.normal, planeID.normal, tol))
                C4Candidates.insert(planeID.normal);
        }
    }

    // Check each candidate C4 axis
    std::vector<Vec3> C4Axes;
    for (const auto& axis : C4Candidates)
    {
        if (isSymmetryOperation(geometry, Cn(axis, 4), tol))
            C4Axes.push_back(axis);
    }

    return C4Axes;
}

std::vector<Vec3> findIcosahedralC5s(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, double tol)
{

    UniqueVec3Set C5Candidates(10, tol);

    for (const auto& sea : SEAs)
    {
        size_t size = sea.indices.size();
        if (size < 5) // looking for pentagons
            continue;

        const auto& idx = sea.indices;

        // find unique normals to planes that at least contain five atoms
        UniqueVec3Set planeNormals(10, tol);
        for (size_t i = 0; i < size; ++i)
        {
            const Vec3& p1 = geometry[idx[i]].coords;
            for (size_t j = i + 1; j < size; ++j)
            {
                const Vec3& p2 = geometry[idx[j]].coords;
                for (size_t k = j + 1; k < size; ++k)
                {
                    const Vec3& p3 = geometry[idx[k]].coords;

                    Vec3 normal = (p2 - p1).cross(p3 - p1);
                    for (size_t l = k + 1; l < size; ++l)
                    {
                        const Vec3& p4 = geometry[idx[l]].coords;
                        if (areOrthogonal(p4 - p1, normal, tol))
                        {
                            for (size_t m = l + 1; m < size; ++m)
                            {
                                const Vec3& p5 = geometry[idx[m]].coords;
                                if (areOrthogonal(p5 - p1, normal, tol))
                                {
                                    planeNormals.insert(normal);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Map planes to the positions of atoms that lie on them.
        std::unordered_map<PlaneIdentifier, std::vector<Vec3>, PlaneIdentifierHasher, PlaneIdentifierEqual> all_planes(
            10, PlaneIdentifierHasher(tol), PlaneIdentifierEqual(tol)
        );
        for (const auto& normal : planeNormals)
        {
            for (const auto& i : idx)
            {
                const Vec3& p   = geometry[i].coords;
                double distance = normal.dot(p);

                all_planes[{normal, distance}].push_back(p);
            }
        }

        // Check each planar group for a regular pentagon.
        for (const auto& [planeID, pointsOnPlane] : all_planes)
        {
            if (pointsOnPlane.size() != 5)
                continue;

            if (isRegularPolygon(pointsOnPlane, planeID.distance * planeID.normal, planeID.normal, tol))
                C5Candidates.insert(planeID.normal);
        }
    }

    // Check each candidate C5 axis
    std::vector<Vec3> C5Axes;
    for (const auto& axis : C5Candidates)
    {
        if (isSymmetryOperation(geometry, Cn(axis, 5), tol))
            C5Axes.push_back(axis);
    }

    return C5Axes;
}


Vec3 findC2Prime(const std::vector<Atom>& geometry, const std::vector<Vec3>& C2s, const Vec3& pAxis, double tol)
{
    std::vector<Vec3> orthoC2s;
    std::ranges::copy_if(
        C2s, std::back_inserter(orthoC2s), [&pAxis, &tol](const Vec3& axis) { return areOrthogonal(axis, pAxis, tol); }
    );

    if (orthoC2s.empty())
        return Vec3::Zero();

    // Choose the C2' that passes through the most atoms, and if tied, the most bonds.
    auto c2prime = std::ranges::max_element(
        orthoC2s,
        [&geometry, &tol](const Vec3& a, const Vec3& b)
        {
            size_t atomsA = 0;
            size_t atomsB = 0;
            for (const auto& atom : geometry)
            {
                if (areCollinear(atom.coords, a, tol))
                    atomsA++;
                if (areCollinear(atom.coords, b, tol))
                    atomsB++;
            }
            if (atomsA != atomsB)
                return atomsA < atomsB;

            size_t bondsA = 0;
            size_t bondsB = 0;
            for (size_t i = 0; i < geometry.size(); ++i)
            {
                for (size_t j = i + 1; j < geometry.size(); ++j)
                {
                    Vec3 bondCenter = geometry[i].coords + geometry[j].coords;
                    if (areCollinear(bondCenter, a, tol))
                        bondsA++;
                    if (areCollinear(bondCenter, b, tol))
                        bondsB++;
                }
            }
            return bondsA < bondsB;
        }
    );

    return *c2prime;
}

std::vector<Vec3> findSigmaVs(const std::vector<Atom>& geometry, const std::vector<SEA>& SEAs, const Vec3& paxis, double tol)
{
    UniqueVec3Set sigmaVCandidates(10, tol);

    for (const auto& sea : SEAs)
    {
        if (sea.indices.size() < 2)
            continue;

        for (size_t i = 0; i < sea.indices.size(); ++i)
        {
            const Vec3& p1 = geometry[sea.indices[i]].coords;

            for (size_t j = i + 1; j < sea.indices.size(); ++j)
            {
                const Vec3& p2 = geometry[sea.indices[j]].coords;

                Vec3 n = (p2 - p1).normalized();

                if (!areOrthogonal(n, paxis, tol))
                    continue;

                if (isSymmetryOperation(geometry, createReflectionMatrix(n), tol))
                    sigmaVCandidates.insert(n);
            }
        }
    }

    // Handle planar molecules with no SEAs.
    auto [tensor, moments, axes] = diagonalizeInertiaTensor(geometry);
    double Ia = moments[0], Ib = moments[1], Ic = moments[2];
    if (areEqual(Ia + Ib, Ic, tol))
        sigmaVCandidates.insert(axes.col(2));

    std::vector<Vec3> sigmaVs;
    for (const auto& axis : sigmaVCandidates)
    {
        if (isSymmetryOperation(geometry, createReflectionMatrix(axis), tol))
            sigmaVs.push_back(axis);
    }

    return sigmaVs;
}


Vec3 findSigmaVPrime(const std::vector<Atom>& geometry, const std::vector<Vec3>& sigmaVs, double tol)
{
    if (sigmaVs.empty())
        return Vec3::Zero();

    // Choose the sigmaV' that contains the most atoms, and if tied, passes through the most bonds.
    auto sigmavprime = std::ranges::max_element(
        sigmaVs,
        [&geometry, &tol](const Vec3& a, const Vec3& b)
        {
            size_t atomsA = 0;
            size_t atomsB = 0;
            for (const auto& atom : geometry)
            {
                if (areOrthogonal(atom.coords, a, tol))
                    atomsA++;
                if (areOrthogonal(atom.coords, b, tol))
                    atomsB++;
            }
            if (atomsA != atomsB)
                return atomsA < atomsB;

            size_t bondsA = 0;
            size_t bondsB = 0;
            for (size_t i = 0; i < geometry.size(); ++i)
            {
                for (size_t j = i + 1; j < geometry.size(); ++j)
                {
                    Vec3 bondCenter = geometry[i].coords + geometry[j].coords;
                    if (areOrthogonal(bondCenter, a, tol))
                        bondsA++;
                    if (areOrthogonal(bondCenter, b, tol))
                        bondsB++;
                }
            }
            return bondsA < bondsB;
        }
    );

    return *sigmavprime;
}


bool hasSigmaH(const std::vector<Atom>& geometry, const Vec3& paxis, double tol)
{
    Vec3 n = paxis.normalized();
    return isSymmetryOperation(geometry, createReflectionMatrix(n), tol);
}


MolecularPG findPointGroup(const std::vector<Atom>& geometry, double tol)
{
    auto [tensor, moments, axes] = diagonalizeInertiaTensor(geometry);
    double Ia = moments[0], Ib = moments[1], Ic = moments[2];
    Vec3 Iav = axes.col(0), Ibv = axes.col(1), Icv = axes.col(2);

    std::string symbol;
    Vec3 primaryAxis   = Vec3::Zero();
    Vec3 secondaryAxis = Vec3::Zero();

    if (geometry.size() == 1)
        return {.symbol = "Kh", .primaryAxis = Vec3(0, 0, 1), .secondaryAxis = Vec3(1, 0, 0)};

    if (Ia < tol && areEqual(Ib, Ic, tol)) // Linear
    {
        if (isSymmetryOperation(geometry, inversionMatrix(), tol))
            symbol = "Dinfh";
        else
            symbol = "Cinfv";
        return {.symbol = symbol, .primaryAxis = Iav, .secondaryAxis = Ibv};
    }

    std::vector<SEA> SEAs = findSEAs(geometry, tol);
    if (areEqual(Ia, Ib, tol) && areEqual(Ib, Ic, tol)) // Spherical top
    {
        std::vector<Vec3> c2Axes = findC2s(geometry, SEAs, tol);
        if (c2Axes.size() == 15)
        {
            std::vector<Vec3> c3Axes = findIcosahedralC3s(geometry, SEAs, tol);
            std::vector<Vec3> c5Axes = findIcosahedralC5s(geometry, SEAs, tol);

            primaryAxis = c5Axes[0];

            // Secondary axis a perpendicular C2 axis, such that a counter-clockwise rotation about the C2 brings the C5
            // axis to the C3 axis closest to the C5 axis which is also perpendicular to the C2 axis.
            std::vector<Vec3> pC2s; // C2s that are perpindicular to primary axis
            std::ranges::copy_if(
                c2Axes,
                std::back_inserter(pC2s),
                [&primaryAxis, &tol](const Vec3& axis) { return areOrthogonal(axis, primaryAxis, tol); }
            );

            secondaryAxis = pC2s[0];

            std::vector<Vec3> pC3s;
            std::ranges::copy_if(
                c3Axes,
                std::back_inserter(pC3s),
                [&secondaryAxis, &tol](const Vec3& axis) { return areOrthogonal(axis, secondaryAxis, tol); }
            );

            auto angle = [](const Vec3& a, const Vec3& b)
            {
                double cosTheta = a.dot(b);
                return std::acos(cosTheta);
            };

            Vec3 closestC3 = angle(primaryAxis, pC3s[0]) < angle(primaryAxis, -pC3s[0]) ? pC3s[0] : -pC3s[0];
            double C3Angle = angle(primaryAxis, closestC3);
            for (size_t i = 1; i < pC3s.size(); ++i)
            {
                double tempAngle = angle(primaryAxis, pC3s[i]);
                if (tempAngle < C3Angle)
                {
                    closestC3 = pC3s[i];
                    C3Angle   = tempAngle;
                }

                tempAngle = angle(primaryAxis, -pC3s[i]);
                if (tempAngle < C3Angle)
                {
                    closestC3 = -pC3s[i];
                    C3Angle   = tempAngle;
                }
            }

            Eigen::Matrix3d R = createRotationMatrix(secondaryAxis, 0.6523581392114948);
            if (!areEqual(R * primaryAxis, closestC3, tol))
                secondaryAxis = -secondaryAxis;

            if (isSymmetryOperation(geometry, inversionMatrix(), tol))
                symbol = "Ih";
            else
                symbol = "I";
        }
        else if (c2Axes.size() == 9)
        {
            std::vector<Vec3> c4Axes = findOctahedralC4s(geometry, SEAs, tol);
            primaryAxis              = c4Axes[0];
            secondaryAxis            = c4Axes[1];
            if (isSymmetryOperation(geometry, inversionMatrix(), tol))
                symbol = "Oh";
            else
                symbol = "O";
        }
        else if (c2Axes.size() == 3)
        {
            // The cartesian axes are the three C2 axes.
            primaryAxis   = c2Axes[0];
            secondaryAxis = c2Axes[1];

            if (isSymmetryOperation(geometry, inversionMatrix(), tol))
                symbol = "Th";
            else
            {
                std::vector<Vec3> sigmaVs = findSigmaVs(geometry, SEAs, primaryAxis, tol);
                if (sigmaVs.empty())
                    symbol = "T";
                else
                    symbol = "Td";
            }
        }
    }
    else
    {
        std::vector<std::pair<Vec3, size_t>> rotations = findRotations(geometry, SEAs, tol);
        std::vector<Vec3> c2axes                       = findC2s(geometry, SEAs, tol);

        size_t order = 1;
        if (rotations.empty())
        {
            if (!c2axes.empty())
            {
                order       = 2;
                primaryAxis = c2axes[0];
            }
        }
        else
        {
            auto pRotation = std::ranges::max_element(
                rotations, [](const auto& a, const auto& b) { return a.second < b.second; }
            );
            order       = pRotation->second;
            primaryAxis = pRotation->first;
        }

        if (order == 1)
        {
            if (isSymmetryOperation(geometry, inversionMatrix(), tol))
                symbol = "Ci";
            else
            {
                std::vector<Vec3> sigmaVs = findSigmaVs(geometry, SEAs, Vec3::Zero(), tol);
                if (sigmaVs.empty())
                    symbol = "C1";
                else
                {
                    symbol      = "Cs";
                    primaryAxis = findSigmaVPrime(geometry, sigmaVs, tol);
                }
            }
        }
        else
        {
            Vec3 c2prime              = findC2Prime(geometry, c2axes, primaryAxis, tol);
            std::vector<Vec3> sigmaVs = findSigmaVs(geometry, SEAs, Vec3::Zero(), tol);
            if (c2prime != Vec3::Zero())
            {
                secondaryAxis = c2prime;

                if (hasSigmaH(geometry, primaryAxis, tol))
                    symbol = fmt::format("D{}h", order);
                else if (!sigmaVs.empty())
                    symbol = fmt::format("D{}d", order);
                else
                    symbol = fmt::format("D{}", order);
            }
            else if (hasSigmaH(geometry, primaryAxis, tol))
            {
                symbol = fmt::format("C{}h", order);
            }
            else if (!sigmaVs.empty())
            {
                symbol = fmt::format("C{}v", order);
                if (areEqual(Ia + Ib, Ic, tol)) // planar
                {
                    secondaryAxis = axes.col(2); // normal to plane
                }
                else
                {
                    Vec3 sigmaVPrime = findSigmaVPrime(geometry, sigmaVs, tol);
                    secondaryAxis    = primaryAxis.cross(sigmaVPrime).normalized();
                }
            }
            else
            {
                Eigen::Matrix3d S2n = Sn(primaryAxis, order * 2);
                if (isSymmetryOperation(geometry, S2n, tol))
                    symbol = fmt::format("S{}", 2 * order);
                else
                    symbol = fmt::format("C{}", order);
            }
        }
    }

    // If axes are arbitrary (C1, Ci) set the axes to the axes of inertia.
    if (primaryAxis == Vec3::Zero())
    {
        primaryAxis   = Icv;
        secondaryAxis = Iav;
    }
    else if (secondaryAxis == Vec3::Zero())
    {
        // If secondary axis is arbitrary, take the projection on the plane orthogonal to the primary axis, and set the
        // secondary axis to the highest moment of inertia axis in that plane.

        // Rotate the geometry to align the z-axis, then take the projection on the xy plane.
        auto projectedGeometry         = geometry;
        Vec3 rotationAxis              = primaryAxis.cross(Vec3(0, 0, 1));
        double rotationAngle           = acos(primaryAxis.dot(Vec3(0, 0, 1)));
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

        secondaryAxis = rotationMatirx.transpose() * Vec3(eigenvector.x(), eigenvector.y(), 0.0).normalized();
    }

    return {.symbol = symbol, .primaryAxis = primaryAxis, .secondaryAxis = secondaryAxis};
}

} // namespace Symmetry
