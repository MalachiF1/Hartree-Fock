#include "Symmetry/Helpers.hpp"

namespace Symmetry
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


Vec3Hasher::Vec3Hasher(double tol) : tol(tol) {}


std::size_t Vec3Hasher::operator()(const Vec3& v) const
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


Vec3Equal::Vec3Equal(double tol) : tol(tol) {}
bool Vec3Equal::operator()(const Vec3& a, const Vec3& b) const
{
    return ((a.cross(b)).squaredNorm() < (this->tol * this->tol));
}


UniqueVec3Set::UniqueVec3Set(size_t initial_buckets, double tol) :
    set(initial_buckets, Vec3Hasher(tol), Vec3Equal(tol)), tol(tol)
{
}


void UniqueVec3Set::insert(const Vec3& v)
{
    if (v.squaredNorm() > tol * tol)
        set.insert(v.normalized());
}


void UniqueVec3Set::erase(const Vec3& v)
{
    set.erase(v.normalized());
}


PlaneIdentifierEqual::PlaneIdentifierEqual(double tol) : tol(tol) {}


bool PlaneIdentifierEqual::operator()(const PlaneIdentifier& p1, const PlaneIdentifier& p2) const
{
    long distanceKey_1 = static_cast<long>(std::round(p1.distance / tol));
    long distanceKey_2 = static_cast<long>(std::round(p2.distance / tol));
    if (distanceKey_1 != distanceKey_2)
        return false;
    return Vec3Equal(tol)(p1.normal, p2.normal);
}


PlaneIdentifierHasher::PlaneIdentifierHasher(double tol) : tol(tol) {}

std::size_t PlaneIdentifierHasher::operator()(const PlaneIdentifier& p) const
{
    size_t h1 = Vec3Hasher(tol)(p.normal);
    size_t h2 = std::hash<long> {}(p.distance);
    return (h1 ^ (std::hash<long> {}(h2) + 0x9e3779b9 + (h1 << 6) + (h1 >> 2)));
}


Eigen::Matrix3d createRotationMatrix(const Vec3& axis, double angle)
{
    Eigen::Matrix3d R;
    double s = sin(angle);
    double c = cos(angle);
    double t = 1.0 - c;
    double x = axis.x(), y = axis.y(), z = axis.z();

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


Eigen::Matrix3d createImproperRotationMatrix(const Vec3& axis, double angle)
{
    Eigen::Matrix3d rotation   = createRotationMatrix(axis, angle);
    Eigen::Matrix3d reflection = createReflectionMatrix(axis);
    return reflection * rotation;
}

Eigen::Matrix3d createReflectionMatrix(const Vec3& normal)
{
    return Eigen::Matrix3d::Identity() - 2.0 * normal.normalized() * normal.normalized().transpose();
}


Eigen::Matrix3d Cn(const Vec3& axis, size_t n)
{
    double angle = 2.0 * M_PI / static_cast<double>(n);
    return createRotationMatrix(axis, angle);
}

Eigen::Matrix3d Sn(const Vec3& axis, size_t n)
{
    double angle = 2.0 * M_PI / static_cast<double>(n);
    return createImproperRotationMatrix(axis, angle);
}

Eigen::Matrix3d inversionMatrix()
{
    return -Eigen::Matrix3d::Identity();
}


std::vector<size_t> findDivisors(size_t n)
{
    // Perhaps this isn't the most efficient way to find divisors, but it's good enough for small n.
    // 1 and n are excluded.
    std::vector<size_t> divisors;
    for (size_t i = 2; i <= sqrt(n); ++i)
    {
        if (n % i == 0)
        {
            divisors.push_back(i);

            if (i * i != n)
                divisors.push_back(n / i);
        }
    }

    return divisors;
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
    std::ranges::sort(angles);

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

} // namespace Symmetry
