#include "Symmetry.hpp"

#include <sstream>

SymmetryOperation::SymmetryOperation(SymmetryType type, size_t order, const Vec3& axis) :
    type(type), order(order), axis(axis)
{
}

std::string SymmetryOperation::toString() const
{
    std::string typeStr;
    switch (type)
    {
        case Cn: typeStr = "Cn"; break;
        case sigma: typeStr = "sigma"; break;
        case i: typeStr = "i"; break;
        case Sn: typeStr = "Sn"; break;
    }

    std::stringstream ss;
    ss << typeStr << " (order: " << order << ", axis: [" << axis.x() << ", " << axis.y() << ", " << axis.z() << "])";
    return ss.str();
}

PointGroup::PointGroup(PointGroupClass pointGroupClass, size_t order) : pointGroupClass(pointGroupClass), order(order)
{
}

std::string PointGroup::toString() const
{
    std::string n = std::to_string(order);
    switch (pointGroupClass)
    {
        case Dinfh: return "D∞h";
        case Cinfh: return "C∞h";
        case I: return "I";
        case Ih: return "Ih";
        case O: return "O";
        case Oh: return "Oh";
        case T: return "T";
        case Td: return "Td";
        case Th: return "Th";
        case Dn: return "D" + n;
        case Dnh: return "D" + n + "h";
        case Dnd: return "D" + n + "d";
        case Cn: return "C" + n;
        case Cnh: return "C" + n + "h";
        case Cnv: return "C" + n + "v";
        case Cs: return "Cs";
        case Ci: return "Ci";
        case S2n: return "S" + std::to_string(this->order * 2);
        case C1: return "C1";
        default: return "Unknown point group";
    }
}

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


std::vector<Atom> translateAndRotate(const std::vector<Atom>& geometry, const Eigen::Matrix3d& principalAxes)
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
    auto newGeometry = geometry; // Create a copy to avoid modifying the original geometry
    for (auto& atom : newGeometry) { atom.coords -= centerOfMass; }

    // Rotate the geometry to align with the principal axes
    auto R = principalAxes.transpose(); // rotation matrix
    for (auto& atom : newGeometry) { atom.coords = R * atom.coords; }

    return newGeometry;
}

std::vector<std::vector<size_t>> findSEAGroups(const std::vector<Atom>& geometry, double tol)
{

    // struct to hold pairs of atomic numbers to match with the distance matrix
    struct AtomTypesMatrix
    {
        std::vector<std::vector<std::pair<unsigned, unsigned>>> atomTypes;

        AtomTypesMatrix() = default;
        AtomTypesMatrix(size_t rowNumber, size_t colNumber) :
            atomTypes(rowNumber, std::vector<std::pair<unsigned, unsigned>>(colNumber))
        {
        }

        std::pair<unsigned, unsigned>& operator()(size_t row, size_t col) { return atomTypes[row][col]; }

        std::vector<std::pair<unsigned, unsigned>> col(size_t col) const
        {
            std::vector<std::pair<unsigned, unsigned>> column(atomTypes.size());
            for (size_t i = 0; i < atomTypes.size(); ++i) { column[i] = atomTypes[i][col]; }
            return column;
        }

        bool colsEqual(size_t col1, size_t col2) const
        {
            return std::all_of(
                atomTypes.begin(),
                atomTypes.end(),
                [col1, col2](const std::vector<std::pair<unsigned, unsigned>>& row) { return row[col1] == row[col2]; }
            );
        }

        std::string toString() const
        {
            std::stringstream ss;
            for (const auto& row : atomTypes)
            {
                for (const auto& pair : row) { ss << "(" << pair.first << ", " << pair.second << ") "; }
                ss << "\n";
            }
            return ss.str();
        }
    };

    // map atomic numbers to indices of atoms with that atomic number in the geometry
    std::unordered_map<int, std::vector<size_t>> atomGroups;
    for (size_t i = 0; i < geometry.size(); ++i) { atomGroups[geometry[i].atomicNumber].push_back(i); }

    // calculate the interatomic distance matrix
    Eigen::MatrixXd distanceMatrix(geometry.size(), geometry.size());
    AtomTypesMatrix distanceMatrixAtomTypes(geometry.size(), geometry.size()); // keep track of atom types in the distance matrix
    for (size_t i = 0; i < geometry.size(); ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            distanceMatrix(i, j) = distanceMatrix(j, i) = (geometry[i].coords - geometry[j].coords).norm();

            auto [atomicNumber1, atomicNumber2] = std::minmax(geometry[i].atomicNumber, geometry[j].atomicNumber);
            distanceMatrixAtomTypes(i, j) = distanceMatrixAtomTypes(j, i) = std::make_pair(atomicNumber1, atomicNumber2);
        }
    }

    // sort the distance matrix columns to create a canonical form
    Eigen::MatrixXd canonicalDistanceMatrix(geometry.size(), geometry.size());
    AtomTypesMatrix canonicalDistanceMatrixAtomTypes(geometry.size(), geometry.size());
    for (int i = 0; i < distanceMatrix.cols(); ++i)
    {
        Eigen::VectorXd col                                     = distanceMatrix.col(i);
        std::vector<std::pair<unsigned, unsigned>> colAtomTypes = distanceMatrixAtomTypes.col(i);

        // sort the atomic types toegther with the distances so they stay aligned
        std::vector<std::pair<double, std::pair<int, int>>> pairs;
        for (int j = 0; j < col.size(); ++j) { pairs.emplace_back(col(j), colAtomTypes[j]); }
        std::sort(
            pairs.begin(),
            pairs.end(),
            [tol](
                const std::pair<double, std::pair<unsigned, unsigned>>& a,
                const std::pair<double, std::pair<unsigned, unsigned>>& b
            )
            {
                if (std::abs(a.first - b.first) < tol) // if distances are equal, sort by atom types
                    return a.second.first + a.second.second < b.second.first + b.second.second;
                return a.first < b.first;
            }
        );

        // Fill the sorted column back into the matrix
        for (int j = 0; j < col.size(); ++j)
        {
            canonicalDistanceMatrix(j, i)          = pairs[j].first;
            canonicalDistanceMatrixAtomTypes(j, i) = pairs[j].second;
        }
    }

    // Find indices of Symmetrically Equivalent Atoms (SEAs)
    std::vector<std::vector<size_t>> SEAIndeces;
    for (const auto& [atomicNumber, indices] : atomGroups)
    {
        std::vector<bool> inGroup(indices.size(), false);
        for (size_t i = 0; i < indices.size(); ++i)
        {
            if (inGroup[i])
                continue;

            std::vector<size_t> group;
            group.push_back(indices[i]);
            inGroup[i] = true;

            // Compare the distance vectors of the current atom with all other atoms of the same type
            for (size_t j = i + 1; j < indices.size(); ++j)
            {
                // make sure that the disstances are to the same atom type as well
                if (canonicalDistanceMatrixAtomTypes.colsEqual(indices[i], indices[j])
                    && (canonicalDistanceMatrix.col(indices[i]) - canonicalDistanceMatrix.col(indices[j])).norm() < tol)
                {
                    group.push_back(indices[j]);
                    inGroup[j] = true;
                }
            }

            // add the group to the SEAIndeces
            SEAIndeces.push_back(group);
        }
    }

    return SEAIndeces;
}


bool areMoleculesSuperimposable(const std::vector<Atom>& geometry1, const std::vector<Atom>& geometry2, double tol)
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
                continue; // Skip already matched atoms

            const auto& atom2 = geometry2[j];

            if (atom1.atomicNumber != atom2.atomicNumber)
                continue; // Atomic numbers must match

            if ((atom1.coords - atom2.coords).norm() < tol)
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

    for (const auto& atom : geometry)
    {
        Atom newAtom;
        newAtom.atomicNumber = atom.atomicNumber;
        newAtom.coords       = operation * atom.coords;
        transformedGeometry.push_back(newAtom);
    }

    return areMoleculesSuperimposable(geometry, transformedGeometry, tol);
}


bool hasInvertion(const std::vector<Atom>& geometry, double tol)
{
    Eigen::Matrix3d inversionMatrix = -Eigen::Matrix3d::Identity();
    return isSymmetryOperation(geometry, inversionMatrix, tol);
}

// anonymous namespace for helper functions
namespace
{

bool areColinear(const Vec3& v1, const Vec3& v2, double tol)
{
    return (v1.normalized().cross(v2.normalized())).norm() < tol;
}

void addUniqueVector(std::vector<Vec3>& vecList, const Vec3& vec, double tol)
{
    if (vec.norm() < tol)
        return;

    for (const auto& existingVec : vecList)
    {
        if (areColinear(existingVec, vec, tol))
            return;
    }

    vecList.push_back(vec.normalized());
}

} // namespace

std::vector<SymmetryOperation> findProperRotations(
    const std::vector<Atom>& geometry, const std::vector<std::vector<size_t>>& SEAIndecies, double tol
)
{
    std::vector<Vec3> candidateAxes;

    // add principle axes as candidates, geometry is assumed to be aligned with principal axes.
    addUniqueVector(candidateAxes, {1.0, 0.0, 0.0}, tol);
    addUniqueVector(candidateAxes, {0.0, 1.0, 0.0}, tol);
    addUniqueVector(candidateAxes, {0.0, 0.0, 1.0}, tol);

    // add axes from origin to each atom
    for (const auto& atom : geometry) { addUniqueVector(candidateAxes, atom.coords, tol); }

    // add axes between each pair of symmetrically equivalent atoms
    for (const auto& group : SEAIndecies)
    {
        if (group.size() < 2)
            continue;

        for (size_t i = 0; i < group.size(); ++i)
        {
            for (size_t j = i + 1; j < group.size(); ++j)
            {
                // axis going between the two atoms
                addUniqueVector(candidateAxes, geometry[group[i]].coords + geometry[group[j]].coords, tol);
                // axis from cross product of the two atoms
                addUniqueVector(candidateAxes, geometry[group[i]].coords.cross(geometry[group[j]].coords), tol);
            }
        }
    }

    // validate the candidate axes
    std::vector<SymmetryOperation> symmetryOperations;
    for (const auto& axis : candidateAxes)
    {
        for (size_t order = 2; order <= geometry.size(); ++order)
        {
            double angle                   = 2.0 * M_PI / order;
            Eigen::Matrix3d rotationMatrix = createRotationMatrix(axis, angle);

            if (isSymmetryOperation(geometry, rotationMatrix, tol))
            {
                symmetryOperations.emplace_back(SymmetryOperation::Cn, order, axis);
            }
        }
    }

    return symmetryOperations;
}

std::vector<SymmetryOperation> findReflectionPlanes(
    const std::vector<Atom>& geometry, const std::vector<std::vector<size_t>>& SEAIndecies, double tol
)
{
    std::vector<Vec3> candidateNormals;

    // add cartesian planes, cartesian axes are assumed to be aligned with principal axes.
    addUniqueVector(candidateNormals, {1.0, 0.0, 0.0}, tol);
    addUniqueVector(candidateNormals, {0.0, 1.0, 0.0}, tol);
    addUniqueVector(candidateNormals, {0.0, 0.0, 1.0}, tol);

    // add normals to planes bisecting symmetrically equivalent atom pairs
    for (const auto& group : SEAIndecies)
    {
        if (group.size() < 2)
            continue;

        for (size_t i = 0; i < group.size(); ++i)
        {
            for (size_t j = i + 1; j < group.size(); ++j)
            {
                addUniqueVector(candidateNormals, geometry[group[i]].coords - geometry[group[j]].coords, tol);
            }
        }
    }

    // add normals to planes bisecting (not necessarly identical) atom triplets
    if (geometry.size() >= 3)
    {
        for (size_t i = 0; i < geometry.size(); ++i)
        {
            for (size_t j = i + 1; j < geometry.size(); ++j)
            {
                for (size_t k = j + 1; k < geometry.size(); ++k)
                {
                    Vec3 v1 = geometry[j].coords - geometry[i].coords;
                    Vec3 v2 = geometry[k].coords - geometry[i].coords;
                    addUniqueVector(candidateNormals, v1.cross(v2), tol);
                }
            }
        }
    }

    // validate the candidate axes
    std::vector<SymmetryOperation> symmetryOperations;
    for (const auto& normal : candidateNormals)
    {
        // Check for proper rotations
        Eigen::Matrix3d rotationMatrix = createReflectionMatrix(normal);

        if (isSymmetryOperation(geometry, rotationMatrix, tol))
        {
            symmetryOperations.emplace_back(SymmetryOperation::sigma, 0, normal);
        }
    }

    return symmetryOperations;
}


std::vector<SymmetryOperation> findImproperRotations(
    const std::vector<Atom>& geometry,
    const std::vector<SymmetryOperation>& properRotations,
    const std::vector<SymmetryOperation>& reflectionPlanes,
    double tol
)
{
    if (properRotations.empty() || reflectionPlanes.empty())
        return {};

    std::vector<SymmetryOperation> symmetryOperations;
    for (const auto& rotation : properRotations)
    {
        const auto& axis = rotation.axis;
        size_t order     = rotation.order;
        if (order <= 2)
            continue; // We don't return sigma (S1) of i (S2) operations

        double angle                           = 2.0 * M_PI / order;
        Eigen::Matrix3d rotationMatrix         = createRotationMatrix(axis, angle);
        Eigen::Matrix3d reflectionMatrix       = createReflectionMatrix(axis);
        Eigen::Matrix3d improperRotationMatrix = reflectionMatrix * rotationMatrix;

        if (isSymmetryOperation(geometry, improperRotationMatrix, tol))
        {
            // Find the corresponding reflection plane
            for (const auto& reflection : reflectionPlanes)
            {
                if (areColinear(reflection.axis, axis, tol))
                {
                    symmetryOperations.emplace_back(SymmetryOperation::Sn, order, reflection.axis);
                }
            }
        }
    }

    return symmetryOperations;
}


PointGroup classifyPointGroup(const std::vector<SymmetryOperation>& symmetryOperations, const Vec3& principalMoments, double tol)
{
    auto sortedMoments = principalMoments;
    std::sort(sortedMoments.begin(), sortedMoments.end());

    bool hasInversion = false;
    for (SymmetryOperation op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::i)
        {
            hasInversion = true;
            break;
        }
    }

    // Ia = 0, Ib = Ic: linear
    if (std::abs(sortedMoments[0]) < tol && std::abs(sortedMoments[1] - sortedMoments[2]) < tol)
    {
        if (hasInversion)
            return PointGroup(PointGroup::Dinfh, 0);
        else
            return PointGroup(PointGroup::Cinfh, 0);
    }

    bool hasRotation = false;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::Cn && op.order > 1)
        {
            hasRotation = true;
            break;
        }
    }

    bool hasSigmaPlane = false;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::sigma)
        {
            hasSigmaPlane = true;
            break;
        }
    }

    if (!hasRotation)
    {
        if (hasSigmaPlane)
            return PointGroup(PointGroup::Cs, 0);
        if (hasInversion)
            return PointGroup(PointGroup::Ci, 0);
        else
            return PointGroup(PointGroup::C1, 0);
    }


    size_t principalRotationIndex = 0;
    for (size_t i = 0; i < symmetryOperations.size(); ++i)
    {
        if (symmetryOperations[i].type == SymmetryOperation::Cn
            && symmetryOperations[i].order > symmetryOperations[principalRotationIndex].order)
        {
            principalRotationIndex = i;
        }
    }
    SymmetryOperation principalRotation = symmetryOperations[principalRotationIndex];

    std::vector<SymmetryOperation> perpendicularC2Axes;
    // for (size_t i = 0; i < symmetryOperations.size(); ++i)
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::Cn && op.order == 2)
        {
            if (op.axis.dot(principalRotation.axis) < tol)

            {
                perpendicularC2Axes.push_back(op);
            }
        }
    }

    bool hasSigmaHPlane = false;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::sigma
            && areColinear(op.axis, symmetryOperations[principalRotationIndex].axis, tol))
        {
            hasSigmaHPlane = true;
            break;
        }
    }

    size_t numOfSigmaDPlanes = 0;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::sigma && op.axis.dot(principalRotation.axis) < tol)
        {
            for (size_t i = 0; i < perpendicularC2Axes.size(); ++i)
            {
                for (size_t j = i + 1; j < perpendicularC2Axes.size(); ++j)
                {
                    Vec3 bisector1 = (perpendicularC2Axes[i].axis + perpendicularC2Axes[j].axis).normalized();
                    Vec3 bisector2 = (perpendicularC2Axes[i].axis - perpendicularC2Axes[j].axis).normalized();
                    if (op.axis.dot(bisector1) < tol || op.axis.dot(bisector2) < tol)
                    {
                        numOfSigmaDPlanes++;
                    }
                }
            }
        }
    }

    bool hasSigmaVPlane = false;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::sigma && op.axis.dot(symmetryOperations[principalRotationIndex].axis) < tol)
        {
            hasSigmaVPlane = true;
            break;
        }
    }


    // Ia = Ib = Ic: spherical top
    if (std::abs(sortedMoments[0] - sortedMoments[1]) < tol && std::abs(sortedMoments[1] - sortedMoments[2]) < tol)
    {
        // tetrahedral, octahedral, or icosahedral
        size_t numOfC5Axes = 0;
        for (const auto& op : symmetryOperations)
        {
            if (op.type == SymmetryOperation::Cn && op.order == 5)
            {
                numOfC5Axes++;
            }
        }

        if (numOfC5Axes == 6)
        {
            if (hasInversion)
                return PointGroup(PointGroup::Ih, 0);
            else
                return PointGroup(PointGroup::I, 0);
        }

        size_t numOfC4Axes = 0;
        for (const auto& op : symmetryOperations)
        {
            if (op.type == SymmetryOperation::Cn && op.order == 4)
            {
                numOfC4Axes++;
            }
        }

        if (numOfC4Axes == 3)
        {
            if (hasInversion)
                return PointGroup(PointGroup::Oh, 0);
            else
                return PointGroup(PointGroup::O, 0);
        }

        size_t numOfC3Axes = 0;
        for (const auto& op : symmetryOperations)
        {
            if (op.type == SymmetryOperation::Cn && op.order == 3)
            {
                numOfC3Axes++;
            }
        }

        if (numOfC3Axes == 4)
        {
            if (hasInversion)
                return PointGroup(PointGroup::Th, 0);
            else if (numOfSigmaDPlanes == 6)
                return PointGroup(PointGroup::Td, 0);
            else
                return PointGroup(PointGroup::T, 0);
        }
    }

    // Ia = Ib != Ic
    else if (std::abs(sortedMoments[0] - sortedMoments[1]) < tol || std::abs(sortedMoments[1] - sortedMoments[2]) < tol)
    {
        // symmetric top
        if (perpendicularC2Axes.size() == principalRotation.order)
        {
            if (hasSigmaHPlane)
            {
                return PointGroup(PointGroup::Dnh, principalRotation.order);
            }
            if (numOfSigmaDPlanes == principalRotation.order)
                return PointGroup(PointGroup::Dnd, principalRotation.order);
            return PointGroup(PointGroup::Dn, principalRotation.order);
        }
        else
        {
            if (hasSigmaHPlane)
                return PointGroup(PointGroup::Cnh, principalRotation.order);
            if (hasSigmaVPlane)
                return PointGroup(PointGroup::Cnv, principalRotation.order);
            bool hasImproperRotation = false;
            for (const auto& op : symmetryOperations)
            {
                if (op.type == SymmetryOperation::Sn && op.order > 2)
                {
                    hasImproperRotation = true;
                    break;
                }
            }
            if (hasImproperRotation)
                return PointGroup(PointGroup::S2n, principalRotation.order); // Improper rotation
            return PointGroup(PointGroup::Cn, principalRotation.order);
        }
    }

    // Ia != Ib != Ic: asymmetric top
    size_t numOfC2Axes = 0;
    for (const auto& op : symmetryOperations)
    {
        if (op.type == SymmetryOperation::Cn && op.order == 2)
        {
            numOfC2Axes++;
        }
    }

    if (numOfC2Axes == 3)
    {
        if (hasSigmaHPlane)
            return PointGroup(PointGroup::Dnh, 2);
        else
            return PointGroup(PointGroup::Dn, 2);
    }
    else
    {
        if (hasSigmaHPlane)
            return PointGroup(PointGroup::Cnh, 2);
        if (hasSigmaVPlane)
            return PointGroup(PointGroup::Cnv, 2);
        return PointGroup(PointGroup::Cn, 2);
        // We alrady handled C1, Cs, Ci and S2n
    }
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
} // namespace Symmetry
