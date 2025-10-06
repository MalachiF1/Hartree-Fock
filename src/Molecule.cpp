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

    this->electronCount      = countElectrons();
    this->basis              = Basis(basisName, this->geometry);
    this->basisFunctionCount = this->basis.getAOCount();
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
    aoLabels.reserve(basisFunctionCount);

    const auto& shells = basis.getShells();
    const auto& lx     = basis.getLx();
    const auto& ly     = basis.getLy();
    const auto& lz     = basis.getLz();
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    for (const auto& shell : shells)
    {
        const double x = cx[shell.primOffset];
        const double y = cy[shell.primOffset];
        const double z = cz[shell.primOffset];
        auto atom      = std::ranges::find_if(
            geometry,
            [&x, &y, &z](const Atom& a) { return std::abs((a.coords - Eigen::Vector3d(x, y, z)).squaredNorm()) < 1e-12; }
        );

        // Get the number of atoms of same atomic number before this atom.
        size_t atomCount = std::ranges::count_if(
            geometry.begin(), atom, [&atom](const Atom& a) { return a.atomicNumber == atom->atomicNumber; }
        );

        for (size_t i = 0; i < shell.nao; ++i)
        {
            std::stringstream label;
            label << fmt::format(
                "{:<3} {:<2} {:<2} ", aoLabels.size() + 1, Utils::atomicNumberToName.at(atom->atomicNumber), atomCount + 1
            );

            const unsigned l = lx[shell.aoOffset + i];
            const unsigned m = ly[shell.aoOffset + i];
            const unsigned n = lz[shell.aoOffset + i];

            switch (shell.l)
            {
                case 0: label << "s"; break; // s-orbital
                case 1:
                {
                    std::string sublabel;
                    if (l == 1)
                        sublabel = "x";
                    else if (m == 1)
                        sublabel = "y";
                    else if (n == 1)
                        sublabel = "z";
                    label << fmt::format("p{:<4}", sublabel);
                    break;
                }
                case 2:
                {
                    std::string sublabel;
                    if (l == 2)
                        sublabel = "x2";
                    else if (m == 2)
                        sublabel = "y2";
                    else if (n == 2)
                        sublabel = "z2";
                    else if (l == 1 && m == 1)
                        sublabel = "xy";
                    else if (l == 1 && n == 1)
                        sublabel = "xz";
                    else if (m == 1 && n == 1)
                        sublabel = "yz";
                    label << fmt::format("d{:<4}", sublabel);
                    break;
                }
                case 3:
                {
                    std::string sublabel;
                    if (l == 3)
                        sublabel = "x3";
                    else if (m == 3)
                        sublabel = "y3";
                    else if (n == 3)
                        sublabel = "z3";
                    else if (l == 2 && m == 1)
                        label << "x2y";
                    else if (l == 2 && n == 1)
                        label << "x2z";
                    else if (m == 2 && l == 1)
                        label << "y2x";
                    else if (m == 2 && n == 1)
                        label << "y2z";
                    else if (n == 2 && l == 1)
                        label << "z2x";
                    else if (n == 2 && m == 1)
                        label << "z2y";
                    else if (l == 1 && m == 1 && n == 1)
                        sublabel = "xyz";
                    label << fmt::format("f{:<4}", sublabel);
                    break;
                }
                case 4:
                {
                    std::string sublabel;
                    if (l == 4)
                        sublabel = "x4";
                    else if (m == 4)
                        sublabel = "y4";
                    else if (n == 4)
                        sublabel = "z4";
                    else if (l == 3 && m == 1)
                        sublabel = "x3y";
                    else if (l == 3 && n == 1)
                        sublabel = "x3z";
                    else if (m == 3 && l == 1)
                        sublabel = "y3x";
                    else if (m == 3 && n == 1)
                        sublabel = "y3z";
                    else if (n == 3 && l == 1)
                        sublabel = "z3x";
                    else if (n == 3 && m == 1)
                        sublabel = "z3y";
                    else if (l == 2 && m == 2)
                        sublabel = "x2y2";
                    else if (l == 2 && n == 2)
                        sublabel = "x2z2";
                    else if (m == 2 && n == 2)
                        sublabel = "y2z2";
                    else if (l == 1 && m == 1 && n == 2)
                        sublabel = "xyz2";
                    else if (l == 1 && n == 1 && m == 2)
                        sublabel = "xy2z";
                    else if (m == 1 && n == 1 && l == 2)
                        sublabel = "x2yz";
                    label << fmt::format("g{:<4}", sublabel);
                    break;
                }
                default:
                {
                    label << fmt::format("{:<4}", fmt::format("({}{}{})", l, m, n));
                    break;
                }
            }

            aoLabels.emplace_back(label.str());
        }
    }
    return aoLabels;
}


size_t Molecule::countElectrons() const
{
    size_t totalElectrons = 0;
    for (const auto& atom : geometry) { totalElectrons += atom.atomicNumber; }
    return totalElectrons - charge;
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
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::overlap(basis, basis.getShells()[i], basis.getShells()[j], S);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { S(j, i) = S(i, j); }
    }

    return S;
}


Eigen::MatrixXd Molecule::kineticMatrix() const
{

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::kinetic(basis, basis.getShells()[i], basis.getShells()[j], T);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { T(j, i) = T(i, j); }
    }

    return T;
}


Eigen::MatrixXd Molecule::nuclearAttractionMatrix() const
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basis.getShellCount(); ++i)
    {
        for (size_t j = i; j < basis.getShellCount(); ++j) // only compute upper triangle as S is symmetric
        {
            IntegralEngine::nuclearAttraction(geometry, basis, basis.getShells()[i], basis.getShells()[j], V);
        }
    }
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) { V(j, i) = V(i, j); }
    }

    return V;
}


ElectronRepulsionTensor Molecule::electronRepulsionTensor(double threshold) const
{
    size_t sc = basis.getShellCount();
    Eigen::MatrixXd Q(sc, sc);
    ElectronRepulsionTensor Vee(basisFunctionCount);
    const auto& shells = basis.getShells();

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < sc; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[i], shells[j], Vee);

            double maxVal = 0;
            for (size_t a = 0; a < shells[i].nao; ++a)
            {
                for (size_t b = 0; b < shells[j].nao; ++b)
                {
                    double val = Vee(
                        shells[i].aoOffset + a, shells[j].aoOffset + b, shells[i].aoOffset + a, shells[j].aoOffset + b
                    );
                    maxVal = std::max(maxVal, val);
                }
            }

            Q(j, i) = maxVal;
        }
    }

    Q = Q.selfadjointView<Eigen::Upper>();
    Q = Q.cwiseSqrt();

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
                    if (Q(i, j) * Q(k, l) < threshold)
                        continue;

                    IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[k], shells[l], Vee);
                }
            }
        }
    }

    return Vee;
}

Eigen::MatrixXd Molecule::schwartzScreeningMatrix() const
{
    size_t sc = basis.getShellCount();
    Eigen::MatrixXd Q(sc, sc);
    const auto& shells = basis.getShells();

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < sc; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            std::vector<double> ERIs(shells[i].nao * shells[j].nao * shells[i].nao * shells[j].nao, 0.0);
            IntegralEngine::electronRepulsion(basis, shells[i], shells[j], shells[i], shells[j], ERIs);

            const size_t aOffset = shells[j].nao * shells[i].nao * shells[j].nao + shells[j].nao;
            const size_t bOffset = shells[i].nao * shells[j].nao + 1;

            double maxVal = 0;
            for (size_t a = 0; a < shells[i].nao; ++a)
            {
                for (size_t b = 0; b < shells[j].nao; ++b)
                {
                    double val = ERIs[a * aOffset + b * bOffset];
                    maxVal     = std::max(maxVal, val);
                }
            }

            Q(j, i) = maxVal;
        }
    }

    Q = Q.selfadjointView<Eigen::Upper>();
    Q = Q.cwiseSqrt();

    return Q;
}
