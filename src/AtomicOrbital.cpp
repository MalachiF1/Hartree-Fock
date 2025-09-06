#include "AtomicOrbital.hpp"

#include "Utils.hpp"

#include <cmath>
#include <fmt/core.h>
#include <iostream>
#include <ranges>

AtomicOrbital::AtomicOrbital(
    const Vec3& center, const Eigen::Vector3i& angularMomentum, const std::vector<PrimitiveGaussian>& primitives
) :
    center(center), angularMomentum(angularMomentum), primitives(primitives)
{
}

std::string AtomicOrbital::toString() const
{
    return fmt::format(
        "AtomicOrbital with {} primitive Gaussians at center ({:.8f}, {:.8f}, {:.8f}):\n{{\n",
        primitives.size(),
        center.x(),
        center.y(),
        center.z()
    );
}


bool AtomicOrbital::operator==(const AtomicOrbital& other) const
{
    if (this->center != other.center || this->angularMomentum != other.angularMomentum
        || this->primitives.size() != other.primitives.size())
        return false;

    std::vector<PrimitiveGaussian> primitives1 = this->primitives;
    std::vector<PrimitiveGaussian> primitives2 = other.primitives;
    auto comparator                            = [](const PrimitiveGaussian& a, const PrimitiveGaussian& b)
    {
        return std::tie(a.exponent, a.coeff, a.angularMomentum.x(), a.angularMomentum.y(), a.angularMomentum.z(), a.normConst)
             < std::tie(b.exponent, b.coeff, b.angularMomentum.x(), b.angularMomentum.y(), b.angularMomentum.z(), b.normConst);
    };
    std::ranges::sort(primitives1, comparator);
    std::ranges::sort(primitives2, comparator);
    return primitives1 == primitives2;
}


bool AtomicOrbital::sameSubshell(const AtomicOrbital& ao1, const AtomicOrbital& ao2)
{
    if (ao1.angularMomentum.sum() != ao2.angularMomentum.sum())
        return false;

    std::vector<std::pair<double, double>> exponentCoeffPairs1;
    std::vector<std::pair<double, double>> exponentCoeffPairs2;

    for (const auto& primitive : ao1.getPrimitives())
        exponentCoeffPairs1.emplace_back(primitive.exponent, primitive.coeff);
    for (const auto& primitive : ao2.getPrimitives())
        exponentCoeffPairs2.emplace_back(primitive.exponent, primitive.coeff);

    std::ranges::sort(exponentCoeffPairs1);
    std::ranges::sort(exponentCoeffPairs2);


    return exponentCoeffPairs1 == exponentCoeffPairs2;
}


double AtomicOrbital::overlap(const AtomicOrbital& ao1, const AtomicOrbital& ao2)
{
    int l1 = ao1.angularMomentum.x();
    int m1 = ao1.angularMomentum.y();
    int n1 = ao1.angularMomentum.z();
    int l2 = ao2.angularMomentum.x();
    int m2 = ao2.angularMomentum.y();
    int n2 = ao2.angularMomentum.z();

    double total_overlap = 0.0;
    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            double p = p1.exponent + p2.exponent;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;

            double prefactor = std::pow(M_PI / p, 1.5) * std::exp(-q * (p1.coords - p2.coords).squaredNorm());

            double Sx = E(0, l1, l2, Q.x(), p1.exponent, p2.exponent);
            double Sy = E(0, m1, m2, Q.y(), p1.exponent, p2.exponent);
            double Sz = E(0, n1, n2, Q.z(), p1.exponent, p2.exponent);

            total_overlap += p1.normConst * p2.normConst * p1.coeff * p2.coeff * prefactor * Sx * Sy * Sz;
        }
    }
    return total_overlap;
}

double AtomicOrbital::kinetic(const AtomicOrbital& ao1, const AtomicOrbital& ao2)
{
    double total_kinetic = 0.0;

    int l1 = ao1.angularMomentum.x();
    int m1 = ao1.angularMomentum.y();
    int n1 = ao1.angularMomentum.z();
    int l2 = ao2.angularMomentum.x();
    int m2 = ao2.angularMomentum.y();
    int n2 = ao2.angularMomentum.z();

    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            double p = p1.exponent + p2.exponent;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;

            double term1 = p2.exponent * (2 * (l2 + m2 + n2) + 3) * E(0, l1, l2, Q.x(), p1.exponent, p2.exponent)
                         * E(0, m1, m2, Q.y(), p1.exponent, p2.exponent) * E(0, n1, n2, Q.z(), p1.exponent, p2.exponent);

            double term2 = -2.0 * std::pow(p2.exponent, 2)
                         * (E(0, l1, l2 + 2, Q.x(), p1.exponent, p2.exponent) * E(0, m1, m2, Q.y(), p1.exponent, p2.exponent)
                                * E(0, n1, n2, Q.z(), p1.exponent, p2.exponent)
                            + E(0, l1, l2, Q.x(), p1.exponent, p2.exponent) * E(0, m1, m2 + 2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, n1, n2, Q.z(), p1.exponent, p2.exponent)
                            + E(0, l1, l2, Q.x(), p1.exponent, p2.exponent) * E(0, m1, m2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, n1, n2 + 2, Q.z(), p1.exponent, p2.exponent));

            double term3 = -0.5
                         * (l2 * (l2 - 1) * E(0, l1, l2 - 2, Q.x(), p1.exponent, p2.exponent)
                                * E(0, m1, m2, Q.y(), p1.exponent, p2.exponent) * E(0, n1, n2, Q.z(), p1.exponent, p2.exponent)
                            + m2 * (m2 - 1) * E(0, l1, l2, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, m1, m2 - 2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, n1, n2, Q.z(), p1.exponent, p2.exponent)
                            + n2 * (n2 - 1) * E(0, l1, l2, Q.x(), p1.exponent, p2.exponent)
                                  * E(0, m1, m2, Q.y(), p1.exponent, p2.exponent)
                                  * E(0, n1, n2 - 2, Q.z(), p1.exponent, p2.exponent));

            double prefactor = std::pow(M_PI / p, 1.5) * std::exp(-q * Q.squaredNorm());
            total_kinetic += p1.normConst * p2.normConst * p1.coeff * p2.coeff * prefactor * (term1 + term2 + term3);
        }
    }
    return total_kinetic;
}

double AtomicOrbital::nuclearAttraction(const AtomicOrbital& ao1, const AtomicOrbital& ao2, const Vec3& nucleusCoords)
{
    double total_attraction = 0.0;

    int l1 = ao1.angularMomentum.x();
    int m1 = ao1.angularMomentum.y();
    int n1 = ao1.angularMomentum.z();
    int l2 = ao2.angularMomentum.x();
    int m2 = ao2.angularMomentum.y();
    int n2 = ao2.angularMomentum.z();

    int max_t = l1 + l2;
    int max_u = m1 + m2;
    int max_v = n1 + n2;
    RCache RCache(max_t, max_u, max_v);

    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {
            RCache.clear();

            double p = p1.exponent + p2.exponent;
            Vec3 P   = (p1.exponent * p1.coords + p2.exponent * p2.coords) / p;
            double q = (p1.exponent * p2.exponent) / p;
            Vec3 Q   = p1.coords - p2.coords;
            Vec3 PC  = P - nucleusCoords;

            double prefactor = -1.0 * (2.0 * M_PI / p) * std::exp(-q * (p1.coords - p2.coords).squaredNorm());

            double sum = 0.0;
            for (int t = 0; t <= l1 + l2; ++t)
            {
                for (int u = 0; u <= m1 + m2; ++u)
                {
                    for (int v = 0; v <= n1 + n2; ++v)
                    {
                        sum += E(t, l1, l2, Q.x(), p1.exponent, p2.exponent)
                             * E(u, m1, m2, Q.y(), p1.exponent, p2.exponent) * E(v, n1, n2, Q.z(), p1.exponent, p2.exponent)
                             * R(t, u, v, 0, p, PC, p * PC.squaredNorm(), RCache);
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

    int l1 = ao1.angularMomentum.x();
    int m1 = ao1.angularMomentum.y();
    int n1 = ao1.angularMomentum.z();
    int l2 = ao2.angularMomentum.x();
    int m2 = ao2.angularMomentum.y();
    int n2 = ao2.angularMomentum.z();
    int l3 = ao3.angularMomentum.x();
    int m3 = ao3.angularMomentum.y();
    int n3 = ao3.angularMomentum.z();
    int l4 = ao4.angularMomentum.x();
    int m4 = ao4.angularMomentum.y();
    int n4 = ao4.angularMomentum.z();

    int max_t = l1 + l2 + l3 + l4;
    int max_u = m1 + m2 + m3 + m4;
    int max_v = n1 + n2 + n3 + n4;
    RCache RCache(max_t, max_u, max_v);

    for (const auto& p1 : ao1.primitives)
    {
        for (const auto& p2 : ao2.primitives)
        {

            double p_ab = p1.exponent + p2.exponent;
            Vec3 P_ab   = (p1.exponent * p1.coords + p2.exponent * p2.coords) / p_ab;
            double q_ab = (p1.exponent * p2.exponent) / p_ab;
            Vec3 Q_ab   = p1.coords - p2.coords;

            std::vector<double> Ex1_values(l1 + l2 + 1);
            std::vector<double> Ey1_values(m1 + m2 + 1);
            std::vector<double> Ez1_values(n1 + n2 + 1);

            // precompute E values so we don't recompute them in the nested loops
            for (int t1 = 0; t1 <= l1 + l2; ++t1) Ex1_values[t1] = E(t1, l1, l2, Q_ab.x(), p1.exponent, p2.exponent);
            for (int u1 = 0; u1 <= m1 + m2; ++u1) Ey1_values[u1] = E(u1, m1, m2, Q_ab.y(), p1.exponent, p2.exponent);
            for (int v1 = 0; v1 <= n1 + n2; ++v1) Ez1_values[v1] = E(v1, n1, n2, Q_ab.z(), p1.exponent, p2.exponent);

            for (const auto& p3 : ao3.primitives)
            {
                for (const auto& p4 : ao4.primitives)
                {
                    RCache.clear();

                    double p_cd = p3.exponent + p4.exponent;
                    Vec3 P_cd   = (p3.exponent * p3.coords + p4.exponent * p4.coords) / p_cd;
                    double q_cd = (p3.exponent * p4.exponent) / p_cd;
                    Vec3 Q_cd   = p3.coords - p4.coords;

                    std::vector<double> Ex2_values(l3 + l4 + 1);
                    std::vector<double> Ey2_values(m3 + m4 + 1);
                    std::vector<double> Ez2_values(n3 + n4 + 1);

                    for (int t2 = 0; t2 <= l3 + l4; ++t2)
                        Ex2_values[t2] = E(t2, l3, l4, Q_cd.x(), p3.exponent, p4.exponent);
                    for (int u2 = 0; u2 <= m3 + m4; ++u2)
                        Ey2_values[u2] = E(u2, m3, m4, Q_cd.y(), p3.exponent, p4.exponent);
                    for (int v2 = 0; v2 <= n3 + n4; ++v2)
                        Ez2_values[v2] = E(v2, n3, n4, Q_cd.z(), p3.exponent, p4.exponent);

                    double delta_eri = (p_ab * p_cd) / (p_ab + p_cd);
                    double T         = delta_eri * (P_ab - P_cd).squaredNorm();

                    double prefactor = (2.0 * std::pow(M_PI, 2.5) / (p_ab * p_cd * std::sqrt(p_ab + p_cd)))
                                     * std::exp(-q_ab * (p1.coords - p2.coords).squaredNorm())
                                     * std::exp(-q_cd * (p3.coords - p4.coords).squaredNorm());

                    double sum = 0.0;
                    for (int t1 = 0; t1 <= l1 + l2; ++t1)
                    {
                        double Ex1 = Ex1_values[t1];
                        for (int u1 = 0; u1 <= m1 + m2; ++u1)
                        {
                            double Ey1 = Ey1_values[u1];
                            for (int v1 = 0; v1 <= n1 + n2; ++v1)
                            {
                                double Ez1 = Ez1_values[v1];
                                for (int t2 = 0; t2 <= l3 + l4; ++t2)
                                {
                                    double Ex2 = Ex2_values[t2];
                                    for (int u2 = 0; u2 <= m3 + m4; ++u2)
                                    {
                                        double Ey2 = Ey2_values[u2];
                                        for (int v2 = 0; v2 <= n3 + n4; ++v2)
                                        {
                                            double Ez2 = Ez2_values[v2];
                                            double R_val = R(t1 + t2, u1 + u2, v1 + v2, 0, delta_eri, P_ab - P_cd, T, RCache);
                                            int signFactor = (t2 + u2 + v2) % 2 == 0 ? 1 : -1;
                                            sum += Ex1 * Ey1 * Ez1 * Ex2 * Ey2 * Ez2 * signFactor * R_val;
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
    if (i < 0 || i > (l1 + l2) || l1 < 0 || l2 < 0)
    {
        // out out bounds
        return 0.0;
    }
    else if (i == 0 && l1 == 0 && l2 == 0)
    {
        // base case
        return 1.0;
    }

    double p = exponentA + exponentB;
    double q = (exponentA * exponentB) / p;

    double result = 0.0;
    if (l2 == 0)
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

    return result;
}

double AtomicOrbital::R(int t, int u, int v, int n, double p, const Vec3& PC, double T, RCache& cache)
{
    if (t < 0 || u < 0 || v < 0)
    {
        return 0.0;
    }

    double cachedValue = cache(t, u, v, n);
    if (cachedValue > -0.5) // Use -1.0 as the sentinel value
        return cachedValue;

    double result = 0.0;
    if (t == 0 && u == 0 && v == 0)
    {
        // base case
        result = std::pow(-2.0 * p, n) * boys(n, T);
    }
    else if (t > 0)
    {
        result = (t - 1) * R(t - 2, u, v, n + 1, p, PC, T, cache) + PC.x() * R(t - 1, u, v, n + 1, p, PC, T, cache);
    }
    else if (u > 0)
    {
        result = (u - 1) * R(t, u - 2, v, n + 1, p, PC, T, cache) + PC.y() * R(t, u - 1, v, n + 1, p, PC, T, cache);
    }
    else // v > 0
    {
        result = (v - 1) * R(t, u, v - 2, n + 1, p, PC, T, cache) + PC.z() * R(t, u, v - 1, n + 1, p, PC, T, cache);
    }

    return cache(t, u, v, n) = result;
}
