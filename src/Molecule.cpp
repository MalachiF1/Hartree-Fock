#include "Molecule.hpp"

#include "Basis.hpp"
#include "Eigen/Core"
#include "IntegralEngine.hpp"
#include "Symmetry/Symmetry.hpp"
#include "Utils.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <sstream>

Molecule::Molecule(
    int charge, int multiplicity, const std::string& basisName, const std::vector<Atom>& geometry, bool detectSymmetry, double symmetryTolerance
) :
    charge(charge), multiplicity(multiplicity), symmetryTolerance(symmetryTolerance)
{

    if (detectSymmetry)
    {
        this->geometry = geometry;
        Symmetry::translateOrigin(this->geometry, Symmetry::findCOM(this->geometry));

        auto SEAs                   = Symmetry::findSEAs(this->geometry, symmetryTolerance);
        auto [symbol, paxis, saxis] = Symmetry::findPointGroup(this->geometry, symmetryTolerance);
        Symmetry::rotateToNewAxes(this->geometry, saxis, paxis.cross(saxis), paxis);
    }
    else
    {
        this->geometry = geometry;
    }

    this->electronCount = countElectrons();
    buildBasis(basisName);
    this->basis              = Basis(basisName, this->geometry);
    this->basisFunctionCount = this->basis.getAOCount();
    fmt::println("Number of basis functions (original): {}", this->atomicOrbitals.size());
    fmt::println("Number of basis functions (new): {}", this->basisFunctionCount);
    fmt::println("Number of shells: {}", this->basis.getShellCount());
    for (const auto& shell : this->basis.getShells()) { fmt::println("  {}", shell.toString()); }
    fmt::print("Exponents ({}): ", this->basis.getExponents().size());
    for (const auto& exp : this->basis.getExponents()) { fmt::print("{:.4f} ", exp); }
    fmt::println("\nAngular momenta ({}): ", this->basis.getLx().size());
    for (size_t i = 0; i < this->basis.getAOCount(); ++i)
    {
        fmt::print("({}, {}, {}) ", this->basis.getLx()[i], this->basis.getLy()[i], this->basis.getLz()[i]);
    }
    fmt::println("\nCoefficients ({}): ", this->basis.getCoefficients().size());
    for (const auto& coeff : this->basis.getCoefficients()) { fmt::print("{:.4f} ", coeff); }
    fmt::println("");
}

std::string Molecule::toString() const
{
    return fmt::format(
        "Molecule with charge {} and multiplicity {}:\nNumber of electrons: {}\nNumber of basis functions: "
        "{}\nGeometry (atom, x, y, z):\n{}",
        charge,
        multiplicity,
        electronCount,
        basisFunctionCount,
        fmt::join(
            geometry
                | std::views::transform(
                    [](const Atom& atom)
                    {
                        return fmt::format(
                            "\t{:<2} {:>3} {:>13.8f} {:>13.8f} {:>13.8f}\n",
                            Utils::atomicNumberToName.at(atom.atomicNumber),
                            atom.atomicNumber,
                            atom.coords.x(),
                            atom.coords.y(),
                            atom.coords.z()
                        );
                    }
                ),
            ""
        )
    );
}


std::vector<std::string> Molecule::getAOLabels() const
{
    std::vector<std::string> aoLabels;
    aoLabels.reserve(atomicOrbitals.size());
    for (size_t i = 0; i < atomicOrbitals.size(); ++i)
    {
        const AtomicOrbital& ao = atomicOrbitals[i];
        std::stringstream label;
        auto atom = std::ranges::find_if(
            geometry, [&ao](const Atom& a) { return std::abs((a.coords - ao.center).squaredNorm()) < 1e-12; }
        );
        // Get the number of atoms of same atomic number before this atom.
        size_t atomCount = std::ranges::count_if(
            geometry.begin(), atom, [&atom](const Atom& a) { return a.atomicNumber == atom->atomicNumber; }
        );

        label << fmt::format("{:<3} {:<2} {:<2} ", i + 1, Utils::atomicNumberToName.at(atom->atomicNumber), atomCount + 1);

        switch (ao.angularMomentum.sum())
        {
            case 0: label << "s"; break; // s-orbital
            case 1:
            {
                std::string sublabel;
                if (ao.angularMomentum.x() == 1)
                    sublabel = "x";
                else if (ao.angularMomentum.y() == 1)
                    sublabel = "y";
                else if (ao.angularMomentum.z() == 1)
                    sublabel = "z";
                label << fmt::format("p{:<4}", sublabel);
                break;
            }
            case 2:
            {
                std::string sublabel;
                if (ao.angularMomentum.x() == 2)
                    sublabel = "x2";
                else if (ao.angularMomentum.y() == 2)
                    sublabel = "y2";
                else if (ao.angularMomentum.z() == 2)
                    sublabel = "z2";
                else if (ao.angularMomentum.x() == 1 && ao.angularMomentum.y() == 1)
                    sublabel = "xy";
                else if (ao.angularMomentum.x() == 1 && ao.angularMomentum.z() == 1)
                    sublabel = "xz";
                else if (ao.angularMomentum.y() == 1 && ao.angularMomentum.z() == 1)
                    sublabel = "yz";
                label << fmt::format("d{:<4}", sublabel);
                break;
            }
            case 3:
            {
                std::string sublabel;
                if (ao.angularMomentum.x() == 3)
                    sublabel = "x3";
                else if (ao.angularMomentum.y() == 3)
                    sublabel = "y3";
                else if (ao.angularMomentum.z() == 3)
                    sublabel = "z3";
                else if (ao.angularMomentum.x() == 2 && ao.angularMomentum.y() == 1)
                    label << "x2y";
                else if (ao.angularMomentum.x() == 2 && ao.angularMomentum.z() == 1)
                    label << "x2z";
                else if (ao.angularMomentum.y() == 2 && ao.angularMomentum.x() == 1)
                    label << "y2x";
                else if (ao.angularMomentum.y() == 2 && ao.angularMomentum.z() == 1)
                    label << "y2z";
                else if (ao.angularMomentum.z() == 2 && ao.angularMomentum.x() == 1)
                    label << "z2x";
                else if (ao.angularMomentum.z() == 2 && ao.angularMomentum.y() == 1)
                    label << "z2y";
                else if (ao.angularMomentum.x() == 1 && ao.angularMomentum.y() == 1 && ao.angularMomentum.z() == 1)
                    sublabel = "xyz";
                label << fmt::format("f{:<4}", sublabel);
                break;
            }
            case 4:
            {
                std::string sublabel;
                if (ao.angularMomentum.x() == 4)
                    sublabel = "x4";
                else if (ao.angularMomentum.y() == 4)
                    sublabel = "y4";
                else if (ao.angularMomentum.z() == 4)
                    sublabel = "z4";
                else if (ao.angularMomentum.x() == 3 && ao.angularMomentum.y() == 1)
                    sublabel = "x3y";
                else if (ao.angularMomentum.x() == 3 && ao.angularMomentum.z() == 1)
                    sublabel = "x3z";
                else if (ao.angularMomentum.y() == 3 && ao.angularMomentum.x() == 1)
                    sublabel = "y3x";
                else if (ao.angularMomentum.y() == 3 && ao.angularMomentum.z() == 1)
                    sublabel = "y3z";
                else if (ao.angularMomentum.z() == 3 && ao.angularMomentum.x() == 1)
                    sublabel = "z3x";
                else if (ao.angularMomentum.z() == 3 && ao.angularMomentum.y() == 1)
                    sublabel = "z3y";
                else if (ao.angularMomentum.x() == 2 && ao.angularMomentum.y() == 2)
                    sublabel = "x2y2";
                else if (ao.angularMomentum.x() == 2 && ao.angularMomentum.z() == 2)
                    sublabel = "x2z2";
                else if (ao.angularMomentum.y() == 2 && ao.angularMomentum.z() == 2)
                    sublabel = "y2z2";
                else if (ao.angularMomentum.x() == 1 && ao.angularMomentum.y() == 1 && ao.angularMomentum.z() == 2)
                    sublabel = "xyz2";
                else if (ao.angularMomentum.x() == 1 && ao.angularMomentum.z() == 1 && ao.angularMomentum.y() == 2)
                    sublabel = "xy2z";
                else if (ao.angularMomentum.y() == 1 && ao.angularMomentum.z() == 1 && ao.angularMomentum.x() == 2)
                    sublabel = "x2yz";
                label << fmt::format("g{:<4}", sublabel);
                break;
            }
            default: break;
        }

        aoLabels.emplace_back(label.str());
    }

    return aoLabels;
}

size_t Molecule::countElectrons() const
{
    size_t totalElectrons = 0;
    for (const auto& atom : geometry) { totalElectrons += atom.atomicNumber; }
    return totalElectrons - charge;
}


void Molecule::buildBasis(const std::string& basisName)
{
    std::vector<int> elements;
    for (const auto& atom : geometry)
    {
        // Ensure that no duplicate atomic numbers are added to the elements vector.
        if (std::ranges::find(elements, atom.atomicNumber) == elements.end())
            elements.push_back(atom.atomicNumber);
    }

    auto basisData = Basis::readBasis(basisName, elements);

    for (const auto& atom : geometry)
    {
        const auto& shellsForAtom = basisData.at(atom.atomicNumber);
        for (const auto& shell : shellsForAtom)
        {
            // For a given angular momentum L, generate all components (e.g., L=1 -> px, py, pz)
            // This loop handles the generation of Cartesian Gaussians.
            for (int i = shell.angularMomentum; i >= 0; --i)
            {
                for (int j = shell.angularMomentum - i; j >= 0; --j)
                {
                    int k = shell.angularMomentum - i - j;

                    std::vector<PrimitiveGaussian> primitives;
                    for (size_t p = 0; p < shell.exponents.size(); ++p)
                    {
                        primitives.emplace_back(
                            shell.exponents[p], shell.coefficients[p], atom.coords, Eigen::Vector3i(i, j, k)
                        );
                    }
                    atomicOrbitals.emplace_back(atom.coords, Eigen::Vector3i(i, j, k), primitives);
                }
            }
        }
    }
}


double Molecule::nuclearRepulsion() const
{
    // The nuclear repulsion energy is calculated using the formula:
    // E_NN = sum_{i < j} (Z_i * Z_j) / r_ij
    double energy = 0.0;
    for (size_t i = 0; i < geometry.size(); ++i)
    {
        for (size_t j = i + 1; j < geometry.size(); ++j)
        {
            double dist = (geometry[i].coords - geometry[j].coords).norm();
            energy += (geometry[i].atomicNumber * geometry[j].atomicNumber) / dist;
        }
    }
    return energy;
}


Eigen::MatrixXd Molecule::overlapMatrix() const
{
    Eigen::MatrixXd S(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as S is symmetric
        {
            S(i, j) = S(j, i) = AtomicOrbital::overlap(atomicOrbitals[i], atomicOrbitals[j]);
        }
    }
    // std::cout << "Overlap matrix S:\n" << std::fixed << std::setprecision(6) << S << std::endl;

    Eigen::MatrixXd S_test = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::overlap(basis, basis.getShells()[i], basis.getShells()[j], S_test);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { S_test(j, i) = S_test(i, j); }
    }
    // std::cout << "\nOverlap matrix S_test:\n" << std::fixed << std::setprecision(6) << S_test << "\n" << std::endl;

    fmt::println("Max abs diff S - S_test: {:.6e}", (S - S_test).cwiseAbs().maxCoeff());
    return S_test;
}


Eigen::MatrixXd Molecule::kineticMatrix() const
{
    Eigen::MatrixXd T(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as T is symmetric
        {
            T(i, j) = T(j, i) = AtomicOrbital::kinetic(atomicOrbitals[i], atomicOrbitals[j]);
        }
    }
    // std::cout << "Kinetic matrix T:\n" << std::fixed << std::setprecision(6) << T << std::endl;

    Eigen::MatrixXd T_test = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::kinetic(basis, basis.getShells()[i], basis.getShells()[j], T_test);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { T_test(j, i) = T_test(i, j); }
    }
    // std::cout << "\nKinetic matrix T_test:\n" << std::fixed << std::setprecision(6) << T_test << "\n" << std::endl;

    fmt::println("Max abs diff T - T_test: {:.6e}", (T - T_test).cwiseAbs().maxCoeff());

    return T_test;
}


Eigen::MatrixXd Molecule::nuclearAttractionMatrix() const
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as V is symmetric
        {
            double v_ij = 0.0;
            for (const auto& atom : geometry)
            {
                v_ij += atom.atomicNumber
                      * AtomicOrbital::nuclearAttraction(atomicOrbitals[i], atomicOrbitals[j], atom.coords);
            }
            V(i, j) = V(j, i) = v_ij;
        }
    }
    // std::cout << "Nuclear attraction matrix V:\n" << std::fixed << std::setprecision(6) << V << "\n" << std::endl;

    Eigen::MatrixXd V_test = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::nuclearAttraction(geometry, basis, basis.getShells()[i], basis.getShells()[j], V_test);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { V_test(j, i) = V_test(i, j); }
    }
    // std::cout << "\nNuclear attraction matrix V_test:\n"
    //           << std::fixed << std::setprecision(6) << V_test << "\n"
    //           << std::endl;

    fmt::println("Max abs diff V - V_test: {:.6e}", (V - V_test).cwiseAbs().maxCoeff());

    return V_test;
}


ElectronRepulsionTensor Molecule::electronRepulsionTensor(double threshold) const
{
    size_t N_ao = basisFunctionCount;
    ElectronRepulsionTensor Vee(N_ao);

    // Pre-calculate Schwartz screening matrix
    Eigen::MatrixXd Q(N_ao, N_ao);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double integral = AtomicOrbital::electronRepulsion(
                atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[i], atomicOrbitals[j]
            );
            Q(i, j) = Q(j, i) = integral;

            // We can set these elements in the tensor instead of calculating them again in the next loop.
            // The ElectronRepulsionTensor class handles the symmetry
            Vee(i, j, i, j) = integral;
        }
    }

#pragma omp parallel for collapse(4) schedule(dynamic, 64)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            for (size_t k = 0; k < N_ao; ++k)
            {
                for (size_t l = 0; l < N_ao; ++l)
                {
                    // no need to recalculate identical elements of the tensor (enforce quartet symmetry)
                    if (j > i || l > k || ((i * (i + 1) / 2) + j) < ((k * (k + 1) / 2) + l))
                        continue;

                    // make sure we did not already calculate this integral in Schwartz screening loop
                    if ((i == k && j == l) || (i == l && j == k))
                        continue;

                    // apply Schwartz screening at the AO level
                    if (Q(i, j) * Q(k, l) < threshold * threshold)
                        continue;

                    double integral = AtomicOrbital::electronRepulsion(
                        atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[k], atomicOrbitals[l]
                    );
                    // Store with 8-fold symmetry
                    Vee(i, j, k, l) = integral;
                }
            }
        }
    }

    size_t sc = basis.getShellCount();
    Eigen::MatrixXd Q_test(sc, sc);
    ElectronRepulsionTensor Vee_test(N_ao);
    const std::vector<Shell>& shells = basis.getShells();

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < sc; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[i], shells[j], Vee_test);

            double maxVal = 0;
            for (size_t a = 0; a < shells[i].nao; ++a)
            {
                for (size_t b = 0; b < shells[j].nao; ++b)
                {
                    double val = Vee_test(
                        shells[i].aoOffset + a, shells[j].aoOffset + b, shells[i].aoOffset + a, shells[j].aoOffset + b
                    );
                    maxVal = std::max(maxVal, val);
                }
            }

            Q_test(i, j) = Q_test(j, i) = maxVal;
        }
    }

#pragma omp parallel for collapse(4) schedule(dynamic, 64)
    for (size_t i = 0; i < sc; ++i)
    {
        for (size_t j = 0; j < sc; ++j)
        {
            for (size_t k = 0; k < sc; ++k)
            {
                for (size_t l = 0; l < sc; ++l)
                {
                    // no need to recalculate identical elements of the tensor (enforce quartet symmetry)
                    if (j > i || l > k || (i * (i + 1) / 2 + j) < (k * (k + 1) / 2 + l))
                        continue;

                    // make sure we did not already calculate this integral in Schwartz screening loop
                    if ((i == k && j == l) || (i == l && j == k))
                        continue;

                    // apply Schwartz screening at the shell level
                    if (Q_test(i, j) * Q_test(k, l) < threshold * threshold)
                        continue;

                    IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[k], shells[l], Vee_test);
                }
            }
        }
    }

    fmt::println("Electron repulsion tensor diff:");
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            for (size_t k = 0; k < N_ao; ++k)
            {
                for (size_t l = 0; l < N_ao; ++l)
                {
                    if (j > i || l > k || ((i * (i + 1) / 2) + j) < ((k * (k + 1) / 2) + l))
                        continue;
                    double diff = std::abs(Vee(i, j, k, l) - Vee_test(i, j, k, l));
                    if (diff > 1e-8)
                        fmt::println(
                            "({}, {}, {}, {}): {:.8f} vs {:.8f} (diff: {:.8f})",
                            i,
                            j,
                            k,
                            l,
                            Vee(i, j, k, l),
                            Vee_test(i, j, k, l),
                            diff
                        );
                }
            }
        }
    }


    return Vee_test;
}

Eigen::MatrixXd Molecule::schwartzScreeningMatrix() const
{
    // size_t N_ao = basisFunctionCount;

    // Pre-calculate Schwartz screening matrix
    //     Eigen::MatrixXd Q(N_ao, N_ao);
    // #pragma omp parallel for collapse(2)
    //     for (size_t i = 0; i < N_ao; ++i)
    //     {
    //         for (size_t j = 0; j <= i; ++j)
    //         {
    //             double integral = AtomicOrbital::electronRepulsion(
    //                 atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[i], atomicOrbitals[j]
    //             );
    //             Q(i, j) = Q(j, i) = std::sqrt(integral);
    //         }
    //     }
    //     return Q;

    size_t sc  = basis.getShellCount();
    size_t Nao = basis.getAOCount();
    Eigen::MatrixXd Q_test(sc, sc);
    ElectronRepulsionTensor Vee_test(Nao);
    const std::vector<Shell>& shells = basis.getShells();

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < sc; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[i], shells[j], Vee_test);

            double maxVal = 0;
            for (size_t a = 0; a < shells[i].nao; ++a)
            {
                for (size_t b = 0; b < shells[j].nao; ++b)
                {
                    double val = Vee_test(
                        shells[i].aoOffset + a, shells[j].aoOffset + b, shells[i].aoOffset + a, shells[j].aoOffset + b
                    );
                    maxVal = std::max(maxVal, val);
                }
            }

            Q_test(i, j) = Q_test(j, i) = std::sqrt(maxVal);
        }
    }
    return Q_test;
}
