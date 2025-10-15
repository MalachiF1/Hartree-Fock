#include "Basis.hpp"

#include "Utils.hpp"

#include <charconv>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <json.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

BasisSet::BasisSet(const std::string& name, const std::vector<Atom>& geometry) : name(name)
{
    std::map<int, int> elementCount;
    for (const auto& atom : geometry) { elementCount[atom.atomicNumber]++; }

    std::vector<int> elements;
    for (const auto& [element, count] : elementCount) { elements.push_back(element); }
    std::map<int, std::vector<RawShell>> basisData = readBasis(name, elements);

    nPrimitives  = 0;
    nAOs         = 0;
    nNormFactors = 0;
    for (size_t i = 0; i < geometry.size(); ++i)
    {
        const auto& shellsForAtom = basisData.at(geometry[i].atomicNumber);
        for (const auto& rawShell : shellsForAtom)
        {
            Shell shell;
            shell.l            = rawShell.angularMomentum;
            shell.nPrimitives  = rawShell.exponents.size();
            shell.nAOs         = (rawShell.angularMomentum + 1) * (rawShell.angularMomentum + 2) / 2;
            shell.nNormFactors = shell.nPrimitives * shell.nAOs;
            shell.atomIndex    = i;
            shell.center       = geometry[i].coords;

            shell.primitives.reserve(shell.nPrimitives);
            for (size_t j = 0; j < rawShell.exponents.size(); ++j)
            {
                shell.primitives.emplace_back(rawShell.exponents[j], rawShell.coefficients[j]);
            }

            shell.normalizationFactors.resize(shell.nPrimitives * shell.nAOs);
            shell.angularMomentum.resize(3, shell.nAOs);
            size_t aoIndex = 0;
            for (int i = rawShell.angularMomentum; i >= 0; --i)
            {
                for (int j = rawShell.angularMomentum - i; j >= 0; --j)
                {
                    int k = rawShell.angularMomentum - i - j;

                    shell.angularMomentum(0, aoIndex) = i;
                    shell.angularMomentum(1, aoIndex) = j;
                    shell.angularMomentum(2, aoIndex) = k;
                    aoIndex++;
                }
            }

            for (size_t l = 0; l < rawShell.exponents.size(); ++l)
            {
                aoIndex = 0;
                for (int i = rawShell.angularMomentum; i >= 0; --i)
                {
                    for (int j = rawShell.angularMomentum - i; j >= 0; --j)
                    {
                        int k              = rawShell.angularMomentum - i - j;
                        double prefactor1  = std::pow(2.0 * rawShell.exponents[l] / M_PI, 0.75);
                        double prefactor2  = std::pow(4.0 * rawShell.exponents[l], (i + j + k) / 2.0);
                        double denominator = std::sqrt(dfact((2 * i) - 1) * dfact((2 * j) - 1) * dfact((2 * k) - 1));
                        double normFactor  = prefactor1 * prefactor2 / denominator;
                        shell.normalizationFactors[l * shell.nAOs + aoIndex] = normFactor;
                        aoIndex++;
                    }
                }
            }

            shells.push_back(shell);
            nPrimitives += shell.nPrimitives;
            nAOs += shell.nAOs;
            nNormFactors += shell.nPrimitives * shell.nAOs;
        }
    }
    nShells = shells.size();
}


std::map<int, std::vector<RawShell>> BasisSet::readBasis(const std::string& name, const std::vector<int>& elements)
{
    // Convert name to lowercase for case-insensitive comparison.
    auto view          = name | std::ranges::views::transform(::tolower);
    auto lowercaseName = std::string(view.begin(), view.end());

    // Replace '*' with '_st_'.
    size_t startPos = 0;
    while ((startPos = lowercaseName.find("*", startPos)) != std::string::npos)
    {
        lowercaseName.replace(startPos, 1, "_st_");
        startPos += 4;
    }

    std::string basisPath = "basis_sets/" + lowercaseName + ".json";

    // Check if the basis set exists.
    if (!std::filesystem::exists(basisPath))
    {
        throw std::runtime_error("Basis set '" + name + "' not found.");
    }

    // Load the basis set from the JSON file.
    std::ifstream basisSetFile(basisPath);
    std::stringstream buffer;
    buffer << basisSetFile.rdbuf();
    auto basisJson = nlohmann::json::parse(buffer.str());

    // Store the basis set data in a map where the key is the atomic number and the value is a vector of Shell objects.
    std::map<int, std::vector<RawShell>> basisResult;
    for (const auto& element : elements)
    {
        if (basisJson["elements"].contains(std::to_string(element)))
        {
            for (const auto& shell : basisJson["elements"][std::to_string(element)]["electron_shells"])
            {
                for (size_t i = 0; i < shell["angular_momentum"].get<std::vector<int>>().size(); ++i)
                {
                    int angularMomentum = shell["angular_momentum"].get<std::vector<int>>()[i];

                    // parse exponents
                    std::vector<std::string> exponentsStrings = shell["exponents"].get<std::vector<std::string>>();
                    std::vector<double> exponents;
                    exponents.reserve(exponentsStrings.size());
                    for (const auto& expStr : exponentsStrings)
                    {
                        double value;
                        auto [ptr, ec] = std::from_chars(expStr.data(), expStr.data() + expStr.size(), value);
                        if (ec != std::errc())
                        {
                            throw std::runtime_error("Invalid exponent value: " + expStr);
                        }
                        exponents.push_back(value);
                    }

                    // parse coefficients
                    std::vector<std::string> coefficientsStrings = shell["coefficients"]
                                                                       .get<std::vector<std::vector<std::string>>>()[i];
                    std::vector<double> coefficients;
                    coefficients.reserve(coefficientsStrings.size());
                    for (const auto& coeffStr : coefficientsStrings)
                    {
                        double value;
                        auto [ptr, ec] = std::from_chars(coeffStr.data(), coeffStr.data() + coeffStr.size(), value);
                        if (ec != std::errc())
                        {
                            throw std::runtime_error("Invalid coefficient value: " + coeffStr);
                        }
                        coefficients.push_back(value);
                    }
                    basisResult[element].emplace_back(angularMomentum, exponents, coefficients);
                }
            }
        }
        else
        {
            throw std::runtime_error(
                fmt::format("{} not found in basis set '{}'.", Utils::atomicNumberToName.at(element), name)
            );
        }
    }
    return basisResult;
}
