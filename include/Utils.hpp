#pragma once

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

#include <Eigen/Dense>

/*
 * In the std namespace, we specialize the `hash` struct for tuples.
 *This is used for hashing results of E and R auxiliary functions
 * in the AtomicOrbital class, in an unordered_map.
 */
namespace std
{
template<typename... Args>
struct hash<tuple<Args...>>
{
    size_t operator()(const tuple<Args...>& t) const
    {
        size_t seed = 0;       // Start with a seed of 0
        hash_combine(seed, t); // Combine the hash of each element into the seed
        return seed;
    }
  private:
    // This is the core logic for combining hashes
    template<typename T, size_t I>
    void hash_combine_impl(size_t& seed, const T& t) const
    {
        // 1. Get the hash of the current tuple element (get<I>(t))
        // 2. Add a "magic number" (derived from the golden ratio) and bit-shifted versions ofthe seed
        // 3. XOR the result back into the seed
        seed ^= hash<typename tuple_element<I, T>::type>()(get<I>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    // These helpers apply the above logic to each element of the tuple
    template<typename T, size_t... I>
    void hash_combine(size_t& seed, const T& t, index_sequence<I...>) const
    {
        (hash_combine_impl<T, I>(seed, t), ...);
    }
    template<typename... Ts>
    void hash_combine(size_t& seed, const tuple<Ts...>& t) const
    {
        hash_combine(seed, t, make_index_sequence<sizeof...(Ts)>());
    }
};
}; // namespace std


using Vec3 = Eigen::Vector3d;

/**
 * Represents an atom with its atomic number and coordinates.
 */
struct Atom
{
    int atomicNumber;
    Vec3 coords;
};

/**
 * Struct to hold shell information from a basis set
 */
struct Shell
{
    int angularMomentum;
    std::vector<double> exponents;
    std::vector<double> coefficients;
};

/**
 * Computes the double factorial of an integer n.
 * If n is -1, returns 1.
 *
 * @param n The integer for which to compute the double factorial.
 * @return The double factorial of n, or 1 if n is -1.
 */
int dfact(int n);

/**
 * Computes the Boys function F_m(T) using the hypergeometric function.
 *
 * @param m The order of the Boys function.
 * @param T The argument of the Boys function.
 * @return The value of the Boys function F_m(T).
 */
double boys(int m, double T);

/*
 * Compute the inverse square root of a symmetric positive definite matrix S.
 * This is used for symmetric orthogonalization (Löwdin orthogonalization).
 *
 * @param S The symmetric positive definite matrix.
 * @return The inverse square root matrix S^(-1/2).
 */
Eigen::MatrixXd inverseSqrtMatrix(const Eigen::MatrixXd& S);
