#include "AtomicOrbital.hpp"

#include <cmath>
#include <map>
#include <tuple>

// Memoization caches for the recursive integral helper functions.
namespace
{
std::map<std::tuple<int, int, int, double, double, double>, double> E_cache;
std::map<std::tuple<int, int, int, int, double, double, double, double, double>, double> R_cache;
} // namespace


AtomicOrbital::AtomicOrbital(const Vec3& center, const std::vector<PrimitiveGaussian>& primitives) :
    center(center), primitives(primitives)
{
}

std::string AtomicOrbital::toString() const
{
    std::stringstream ss;
    ss << "AtomicOrbital with " << primitives.size() << " primitive Gaussians at center (" << center.x() << ", "
       << center.y() << ", " << center.z() << "):\n{\n";

    for (const auto& p : primitives)
    {
        // Indent the primitive gaussian string for readability
        std::string prim_str = p.toString();
        std::stringstream prim_ss(prim_str);
        std::string line;
        while (std::getline(prim_ss, line)) { ss << "  " << line << "\n"; }
    }
    ss << "}";
    return ss.str();
}

double AtomicOrbital::overlap(const AtomicOrbital& ao1, const AtomicOrbital& ao2)
{
    double total_overlap = 0.0;
    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            double p = p1.exponent + p2.exponent;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;

            double prefactor = std::pow(M_PI / p, 1.5) * std::exp(-q * (p1.coords - p2.coords).squaredNorm());

            double Sx = E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent);
            double Sy = E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent);
            double Sz = E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent);

            total_overlap += p1.normConst * p2.normConst * p1.coeff * p2.coeff * prefactor * Sx * Sy * Sz;
        }
    }
    return total_overlap;
}

double AtomicOrbital::kinetic(const AtomicOrbital& ao1, const AtomicOrbital& ao2)
{
    double total_kinetic = 0.0;
    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            double p = p1.exponent + p2.exponent;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;

            double term1 = p2.exponent * (2 * (p2.l + p2.m + p2.n) + 3) * E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                         * E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                         * E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent);

            double term2 = -2.0 * std::pow(p2.exponent, 2)
                         * (E(0, p1.l, p2.l + 2, Q.x(), p1.exponent, p2.exponent)
                                * E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                                * E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent)
                            + E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, p1.m, p2.m + 2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent)
                            + E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, p1.n, p2.n + 2, Q.z(), p1.exponent, p2.exponent));

            double term3 = -0.5
                         * (p2.l * (p2.l - 1) * E(0, p1.l, p2.l - 2, Q.x(), p1.exponent, p2.exponent)
                                * E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                                * E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent)
                            + p2.m * (p2.m - 1) * E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, p1.m, p2.m - 2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent)
                            + p2.n * (p2.n - 1) * E(0, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, p1.n, p2.n - 2, Q.z(), p1.exponent, p2.exponent));

            double prefactor = std::pow(M_PI / p, 1.5) * std::exp(-q * Q.squaredNorm());
            total_kinetic += p1.normConst * p2.normConst * p1.coeff * p2.coeff * prefactor * (term1 + term2 + term3);
        }
    }
    return total_kinetic;
}

double AtomicOrbital::nuclearAttraction(const AtomicOrbital& ao1, const AtomicOrbital& ao2, const Vec3& nucleusCoords)
{
    double total_attraction = 0.0;
    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            double p = p1.exponent + p2.exponent;
            Vec3 P   = (p1.exponent * p1.coords + p2.exponent * p2.coords) / p;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;
            Vec3 PC  = P - nucleusCoords;

            double prefactor = -1.0 * (2.0 * M_PI / p) * std::exp(-q * (p1.coords - p2.coords).squaredNorm());

            double sum = 0.0;
            for (int t = 0; t <= p1.l + p2.l; ++t)
            {
                for (int u = 0; u <= p1.m + p2.m; ++u)
                {
                    for (int v = 0; v <= p1.n + p2.n; ++v)
                    {
                        sum += E(t, p1.l, p2.l, Q.x(), p1.exponent, p2.exponent)
                             * E(u, p1.m, p2.m, Q.y(), p1.exponent, p2.exponent)
                             * E(v, p1.n, p2.n, Q.z(), p1.exponent, p2.exponent)
                             * R(t, u, v, 0, p, PC, p * PC.squaredNorm());
                    }
                }
            }
            total_attraction += p1.normConst * p2.normConst * p1.coeff * p2.coeff * prefactor * sum;
        }
    }

    // We multiply by the nuclear charge when calculating the nuclear attraction matrix in the Molecule object
    return total_attraction;
}


double AtomicOrbital::electronRepulsion(
    const AtomicOrbital& ao1, const AtomicOrbital& ao2, const AtomicOrbital& ao3, const AtomicOrbital& ao4
)
{
    double total_repulsion = 0.0;
    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            for (const auto& p3 : ao3.primitives)
            {
                for (const auto& p4 : ao4.primitives)
                {
                    double p_ab = p1.exponent + p2.exponent;
                    double p_cd = p3.exponent + p4.exponent;

                    Vec3 P_ab = (p1.exponent * p1.coords + p2.exponent * p2.coords) / p_ab;
                    Vec3 P_cd = (p3.exponent * p3.coords + p4.exponent * p4.coords) / p_cd;

                    double q_ab = (p1.exponent * p2.exponent) / p_ab;
                    double q_cd = (p3.exponent * p4.exponent) / p_cd;

                    Vec3 Q_ab = p1.coords - p2.coords;
                    Vec3 Q_cd = p3.coords - p4.coords;

                    double delta_eri = (p_ab * p_cd) / (p_ab + p_cd);
                    double T         = delta_eri * (P_ab - P_cd).squaredNorm();

                    double prefactor = (2.0 * std::pow(M_PI, 2.5) / (p_ab * p_cd * std::sqrt(p_ab + p_cd)))
                                     * std::exp(-q_ab * (p1.coords - p2.coords).squaredNorm())
                                     * std::exp(-q_cd * (p3.coords - p4.coords).squaredNorm());

                    double sum = 0.0;
                    for (int t1 = 0; t1 <= p1.l + p2.l; ++t1)
                    {
                        for (int u1 = 0; u1 <= p1.m + p2.m; ++u1)
                        {
                            for (int v1 = 0; v1 <= p1.n + p2.n; ++v1)
                            {
                                for (int t2 = 0; t2 <= p3.l + p4.l; ++t2)
                                {
                                    for (int u2 = 0; u2 <= p3.m + p4.m; ++u2)
                                    {
                                        for (int v2 = 0; v2 <= p3.n + p4.n; ++v2)
                                        {
                                            double Ex1 = E(t1, p1.l, p2.l, Q_ab.x(), p1.exponent, p2.exponent);
                                            double Ey1 = E(u1, p1.m, p2.m, Q_ab.y(), p1.exponent, p2.exponent);
                                            double Ez1 = E(v1, p1.n, p2.n, Q_ab.z(), p1.exponent, p2.exponent);

                                            double Ex2 = E(t2, p3.l, p4.l, Q_cd.x(), p3.exponent, p4.exponent);
                                            double Ey2 = E(u2, p3.m, p4.m, Q_cd.y(), p3.exponent, p4.exponent);
                                            double Ez2 = E(v2, p3.n, p4.n, Q_cd.z(), p3.exponent, p4.exponent);

                                            double R_val = R(t1 + t2, u1 + u2, v1 + v2, 0, delta_eri, P_ab - P_cd, T);

                                            sum += Ex1 * Ey1 * Ez1 * Ex2 * Ey2 * Ez2 * std::pow(-1.0, t2 + u2 + v2) * R_val;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    total_repulsion += p1.normConst * p2.normConst * p3.normConst * p4.normConst * p1.coeff * p2.coeff
                                     * p3.coeff * p4.coeff * prefactor * sum;
                }
            }
        }
    }
    return total_repulsion;
}

double AtomicOrbital::E(int i, int l1, int l2, double Q, double exponentA, double exponentB)
{
    // Use a tuple as the key for the cache
    auto key = std::make_tuple(i, l1, l2, Q, exponentA, exponentB);
    // Check if the result is already cached
    if (E_cache.count(key))
    {
        return E_cache[key];
    }

    double p = exponentA + exponentB;
    double q = (exponentA * exponentB) / p;

    double result = 0.0;
    if (i < 0 || i > (l1 + l2))
    {
        // out of bounds
        return 0.0;
    }
    else if (i == 0 && l1 == 0 && l2 == 0)
    {
        // base case
        return 1.0;
    }
    else if (l2 == 0)
    {
        // decrement l1
        result = ((1.0 / (2.0 * p)) * E(i - 1, l1 - 1, l2, Q, exponentA, exponentB))
               - ((q * Q / exponentA) * E(i, l1 - 1, l2, Q, exponentA, exponentB))
               + ((i + 1) * E(i + 1, l1 - 1, l2, Q, exponentA, exponentB));
    }
    else
    {
        // decrement l2
        result = ((1.0 / (2.0 * p)) * E(i - 1, l1, l2 - 1, Q, exponentA, exponentB))
               + ((q * Q / exponentB) * E(i, l1, l2 - 1, Q, exponentA, exponentB))
               + ((i + 1) * E(i + 1, l1, l2 - 1, Q, exponentA, exponentB));
    }

    // Store the result in the cache
    E_cache[key] = result;
    return result;
}

double AtomicOrbital::R(int t, int u, int v, int n, double p, const Vec3& PC, double T)
{
    // Use a tuple as the key for the cache
    auto key = std::make_tuple(t, u, v, n, p, PC.x(), PC.y(), PC.z(), T);
    // Check if the result is already cached
    if (R_cache.count(key))
    {
        return R_cache[key];
    }

    if (t < 0 || u < 0 || v < 0)
    {
        // out of bounds
        return 0.0;
    }
    else if (t == 0 && u == 0 && v == 0)
    {
        // base case
        return std::pow(-2.0 * p, n) * boys(n, T);
    }

    double result = 0.0;
    if (t > 0)
    {
        result = (t - 1) * R(t - 2, u, v, n + 1, p, PC, T) + PC.x() * R(t - 1, u, v, n + 1, p, PC, T);
    }
    else if (u > 0)
    {
        result = (u - 1) * R(t, u - 2, v, n + 1, p, PC, T) + PC.y() * R(t, u - 1, v, n + 1, p, PC, T);
    }
    else
    { // v > 0
        result = (v - 1) * R(t, u, v - 2, n + 1, p, PC, T) + PC.z() * R(t, u, v - 1, n + 1, p, PC, T);
    }

    R_cache[key] = result;
    return result;
}
