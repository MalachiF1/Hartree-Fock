#pragma once
#include "Utils.hpp"

namespace Symmetry
{

struct InertiaTensor
{
    Eigen::Matrix3d tensor;
    Eigen::Vector3d moments;
    Eigen::Matrix3d axes;
};

InertiaTensor diagonalizeInertiaTensor(const std::vector<Atom>& geometry);

Vec3 findCOM(const std::vector<Atom>& geometry);

void translateOrigin(std::vector<Atom>& geometry, const Vec3& point);

void rotateToNewAxes(std::vector<Atom>& geometry, Vec3 newX, Vec3 newY, Vec3 newZ);


} // namespace Symmetry
