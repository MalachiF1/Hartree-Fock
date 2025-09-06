#pragma once
#include "Utils.hpp"

#include <unordered_set>

namespace Symmetry
{

bool areCollinear(const Vec3& v1, const Vec3& v2, double tol);

bool areOrthogonal(const Vec3& v1, const Vec3& v2, double tol);

bool areEqual(const Vec3& v1, const Vec3& v2, double tol);

bool areEqual(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2, double tol);

bool areEqual(double d1, double d2, double tol);

// Custom hasher for Vec3 to group similar vectors into the same hash bucket.
struct Vec3Hasher
{
    double tol;
    Vec3Hasher(double tol);
    std::size_t operator()(const Vec3& v) const;
};


// Custom equality comparator for Vec3 to handle floating-point tolerance.
struct Vec3Equal
{
    double tol;
    Vec3Equal(double tol);
    bool operator()(const Vec3& a, const Vec3& b) const;
};


class UniqueVec3Set
{
  private:
    std::unordered_set<Vec3, Vec3Hasher, Vec3Equal> set;
    double tol;

  public:
    UniqueVec3Set(size_t initial_buckets, double tol);

    void insert(const Vec3& v);

    void erase(const Vec3& v);

    auto begin() const { return set.begin(); }
    auto end() const { return set.end(); }
    size_t size() const { return set.size(); }
    bool empty() const { return set.empty(); }
};


struct PlaneIdentifier
{
    Vec3 normal;
    double distance;
};


struct PlaneIdentifierEqual
{
    double tol;
    PlaneIdentifierEqual(double tol);

    bool operator()(const PlaneIdentifier& p1, const PlaneIdentifier& p2) const;
};


// A hasher for the PlaneIdentifier struct.
struct PlaneIdentifierHasher
{
    double tol;
    PlaneIdentifierHasher(double tol);

    std::size_t operator()(const PlaneIdentifier& p) const;
};


Eigen::Matrix3d createRotationMatrix(const Vec3& axis, double angle);

Eigen::Matrix3d createImproperRotationMatrix(const Vec3 &axis, double angle);

Eigen::Matrix3d createReflectionMatrix(const Vec3& normal);

Eigen::Matrix3d Cn(const Vec3& axis, size_t n);

Eigen::Matrix3d Sn(const Vec3& axis, size_t n);

Eigen::Matrix3d inversionMatrix();

std::vector<size_t> findDivisors(size_t n);

bool isRegularPolygon(const std::vector<Vec3>& points, const Vec3& center, const Vec3& normal, double tol);

} // namespace Symmetry
