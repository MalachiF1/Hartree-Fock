#include "Basis.hpp"

#include <filesystem>
#include <fstream>
#include <json.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

std::map<int, std::vector<Shell>> Basis::getBasis(const std::string& name, const std::vector<int>& elements)
{
    // convert name to lowercase for case-insensitive comparison
    std::string lowercaseName = name;
    std::transform(
        lowercaseName.begin(),
        lowercaseName.end(),
        lowercaseName.begin(),
        [](unsigned char c) { return std::tolower(c); }
    );

    std::string basisPath = "basis_sets/" + lowercaseName + ".json";

    // Check if the basis set exists
    if (!std::filesystem::exists(basisPath))
    {
        throw std::runtime_error("Basis set '" + name + "' not found.");
    }

    // Load the basis set from the JSON file
    std::ifstream basisSetFile(basisPath);
    std::stringstream buffer;
    buffer << basisSetFile.rdbuf();
    auto basisJson = nlohmann::json::parse(buffer.str());

    std::map<int, std::vector<Shell>> basisResult;

    for (const auto& element : elements)
    {
        if (basisJson["elements"].contains(std::to_string(element)))
        {
            for (const auto& shell : basisJson["elements"][std::to_string(element)]["electron_shells"])
            {
                for (size_t i = 0; i < shell["angular_momentum"].get<std::vector<int>>().size(); ++i)
                {
                    {
                        // Create a new Shell object for each angular momentum
                        Shell newShell;

                        // parse anguar momentum
                        newShell.angularMomentum = shell["angular_momentum"].get<std::vector<int>>()[i];

                        // parse exponents
                        std::vector<std::string> exponentsStrings = shell["exponents"].get<std::vector<std::string>>();
                        std::vector<double> exponents;
                        for (const auto& expStr : exponentsStrings)
                        {
                            try
                            {
                                exponents.push_back(std::stod(expStr));
                            }
                            catch (const std::invalid_argument& e)
                            {
                                throw std::runtime_error("Invalid exponent value: " + expStr);
                            }
                        }
                        newShell.exponents = exponents;

                        // parse coefficients
                        std::vector<std::string>
                            coefficientsStrings = shell["coefficients"].get<std::vector<std::vector<std::string>>>()[i];
                        std::vector<double> coefficients;
                        for (const auto& coeffStr : coefficientsStrings)
                        {
                            try
                            {
                                coefficients.push_back(std::stod(coeffStr));
                            }
                            catch (const std::invalid_argument& e)
                            {
                                throw std::runtime_error("Invalid coefficient value: " + coeffStr);
                            }
                        }
                        newShell.coefficients = coefficients;
                        basisResult[element].push_back(newShell);
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error(
                "Element with atomic number " + std::to_string(element) + " not found in basis set '" + name + "'."
            );
        }
    }
    return basisResult;
}
