#pragma once

#include "PrimitiveGaussian.hpp"

#include <vector>

/**
 * @brief Represents an atomic orbital composed of cartesian primitive Gaussian functions.
 *
 * This class represents an atomic orbital, which is a linear combination of primitive Gaussian functions.
 * An atomic orbital is defined by its center and a set of primitive Gaussians, each with its own exponent,
 * contraction coefficient, and angular momentum numbers. This class contains methods to calculate
 * various integrals involving atomic orbitals, such as overlap, kinetic energy, nuclear attraction,
 * and electron-electron repulsion integrals.
 */
class AtomicOrbital
{
  public:
    const Vec3 center; // center of the atomic orbital (the nucleus position on which the orbital is centered)
    Eigen::Vector3i angularMomentum; // angular momentum quantum numbers for the orbital (l, m, n)

    AtomicOrbital(const Vec3& center, const std::vector<PrimitiveGaussian>& primitives);

    /**
     * @return A string containing the center and the primitive gaussians of the orbital.
     */
    std::string toString() const;

    bool operator==(const AtomicOrbital& other) const;

    /**
     * Get the primitive gaussians that make up the atomic orbital.
     * @return the reference to the vector of primitive gaussians that make up the atomic orbital.
     */
    const std::vector<PrimitiveGaussian>& getPrimitives() const { return primitives; }

    /**
     * Calculates the overlap integral between two atomic orbitals.
     *
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @return The value of the overlap integral.
     */
    static double overlap(const AtomicOrbital& ao1, const AtomicOrbital& ao2);

    /**
     * Calculates the kinetic energy integral between two atomic orbitals.
     *
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @return The value of the kinetic energy integral.
     */
    static double kinetic(const AtomicOrbital& ao1, const AtomicOrbital& ao2);

    /**
     * Calculates the kinetic energy integral between two atomic orbitals.
     *
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @param nucleusCoords The coordinates of the nucleus.
     * @return The value of the kinetic energy integral.
     */
    static double nuclearAttraction(const AtomicOrbital& ao1, const AtomicOrbital& ao2, const Vec3& nucleusCoords);

    /**
     * Calculates the electron-electron repulsion integral between four atomic orbitals.
     *
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @param ao3 The third atomic orbital.
     * @param ao4 The fourth atomic orbital.
     * @return The value of the electron-electron repulsion integral.
     */
    static double electronRepulsion(
        const AtomicOrbital& ao1, const AtomicOrbital& ao2, const AtomicOrbital& ao3, const AtomicOrbital& ao4
    );

    static bool sameSubshell(const AtomicOrbital& ao1, const AtomicOrbital& ao2);


  private:
    /**
     * Helper function for calculating Hermite gaussians.
     * Used for the electron-nuclear attraction and electron-electron repulsion integrals.
     *
     * @param i The index of the Hermite polynomial.
     * @param l1 The angular momentum quantum number of the first orbital.
     * @param l2 The angular momentum quantum number of the second orbital.
     * @param Q The distance between the centers of the two orbitals along the axis.
     * @return The value of the Hermite expansion coefficient.
     */
    static double E(int i, int l1, int l2, double Q, double exponentA, double exponentB);

    /**
     * Helper function for calculating the auxiliary integrals R_tuv.
     * Used for the electron-nuclear attraction and electron-electron repulsion integrals.
     *
     * @param t The combined angular momentum index for the x-dimension.
     * @param u The combined angular momentum index for the y-dimension.
     * @param v The combined angular momentum index for the z-dimension.
     * @param n The order of the Boys function.
     * @param p The sum of the exponent constants of the two Gaussians.
     * @param PC The vector of the distance between the centers of the two Gaussians along each axis.
     * @param T The argument for Boys function, F_m(T).
     * @return The value of the auxiliary integral R_tuv.
     */
    static double R(int t, int u, int v, int n, double p, const Vec3& PC, double T);

    const std::vector<PrimitiveGaussian> primitives;
};
