#pragma once

#include "PrimitiveGaussian.hpp"

#include <vector>

class AtomicOrbital
{
  public:
    const Vec3 center;
    const std::vector<PrimitiveGaussian> primitives;

    AtomicOrbital(const Vec3& center, const std::vector<PrimitiveGaussian>& primitives);


    /**
     * Returns a string representation of the atomic orbital.
     * @return A string containing the center and the primitive gaussians of the orbital.
     */
    std::string toString() const;

    /**
     * Calculates the overlap integral between two atomic orbitals.
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @return The value of the overlap integral.
     */
    static double overlap(const AtomicOrbital& ao1, const AtomicOrbital& ao2);

    /**
     * Calculates the kinetic energy integral between two atomic orbitals.
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @return The value of the kinetic energy integral.
     */
    static double kinetic(const AtomicOrbital& ao1, const AtomicOrbital& ao2);

    /**
     * Calculates the kinetic energy integral between two atomic orbitals.
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @param nucleusCoords The coordinates of the nucleus.
     * @return The value of the kinetic energy integral.
     */
    static double nuclearAttraction(const AtomicOrbital& ao1, const AtomicOrbital& ao2, const Vec3& nucleusCoords);

    /**
     * Calculates the electron-electron repulsion integral between four atomic orbitals.
     * @param ao1 The first atomic orbital.
     * @param ao2 The second atomic orbital.
     * @param ao3 The third atomic orbital.
     * @param ao4 The fourth atomic orbital.
     * @return The value of the electron-electron repulsion integral.
     */
    static double electronRepulsion(
        const AtomicOrbital& ao1, const AtomicOrbital& ao2, const AtomicOrbital& ao3, const AtomicOrbital& ao4
    );

  private:
    /**
     * Helper function for calculating Hermite expansion coefficients.
     * Used for the electron-nuclear attraction and electron-electron repulsion integrals.
     * This function is memoized to improve performance.
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
     * This function is memoized to improve performance.
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
};
