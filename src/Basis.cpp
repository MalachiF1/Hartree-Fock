#include "Basis.hpp"

#include <charconv>
#include <filesystem>
#include <fstream>
#include <json.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

std::map<int, std::vector<Shell>> Basis::getBasis(const std::string& name, const std::vector<int>& elements)
{
    // convert name to lowercase for case-insensitive comparison
    auto toLowerString = [](std::string_view sv) -> std::string
    {
        auto view = sv | std::ranges::views::transform(::tolower);
        return std::string(view.begin(), view.end());
    };
    std::string lowercaseName = toLowerString(name);

    // replace '*' with '_st_'
    size_t startPos = 0;
    while ((startPos = lowercaseName.find("*", startPos)) != std::string::npos)
    {
        lowercaseName.replace(startPos, 1, "_st_");
        startPos += 4;
    }

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

    // store the basis set data in a map where the key is the atomic number and the value is a vector of Shell objects
    std::map<int, std::vector<Shell>> basisResult;
    for (const auto& element : elements)
    {
        if (basisJson["elements"].contains(std::to_string(element)))
        {
            for (const auto& shell : basisJson["elements"][std::to_string(element)]["electron_shells"])
            {
                // Only support cartesian GTOs for now
                if (shell["function_type"] != "gto" && shell["function_type"] != "gto_cartesian")
                {
                    throw std::runtime_error(
                        "Unsupported function type '" + shell["function_type"].get<std::string>()
                        + "' for element with atomic number " + std::to_string(element) + " in basis set '" + name + "'."
                    );
                }

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
                "Element with atomic number " + std::to_string(element) + " not found in basis set '" + name + "'."
            );
        }
    }
    return basisResult;
}
