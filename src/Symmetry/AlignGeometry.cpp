#include "Symmetry/AlignGeometry.hpp"


namespace Symmetry
{

InertiaTensor diagonalizeInertiaTensor(const std::vector<Atom>& geometry)
{
    Eigen::Matrix3d inertiaTensor = Eigen::Matrix3d::Zero();

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

    return InertiaTensor {.tensor = inertiaTensor, .moments = eigenValues, .axes = eigenVectors};
}

Vec3 findCOM(const std::vector<Atom>& geometry)
{
    Vec3 com         = Eigen::Vector3d::Zero();
    double totalMass = 0.0;
    for (const auto& atom : geometry)
    {
        com += atom.coords * Utils::atomicMasses.at(atom.atomicNumber);
        totalMass += Utils::atomicMasses.at(atom.atomicNumber);
    }
    com /= totalMass;
    return com;
}


void translateOrigin(std::vector<Atom>& geometry, const Vec3& point)
{
    // Translate the center of mass to the origin
    for (auto& atom : geometry) { atom.coords -= point; }
}


void rotateToNewAxes(std::vector<Atom>& geometry, Vec3 newX, Vec3 newY, Vec3 newZ)
{
    Eigen::Matrix3d axes;
    axes << newX.normalized(), newY.normalized(), newZ.normalized();
    Eigen::Matrix3d rotationMatrix = axes.transpose();
    for (auto& atom : geometry) { atom.coords = rotationMatrix * atom.coords; }
}

} // namespace Symmetry
