#include "functions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>


// Function to read data from a text file and return it as a vector of doubles
std::vector<std::vector<double>> readData(const std::string& filePath) {
    std::ifstream file(filePath);
    std::vector<std::vector<double>> data;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<double> currentLineData;
            double value;
            while (iss >> value) {
                currentLineData.push_back(value);
            }
            data.push_back(currentLineData);
        }
        file.close();
    } else {
        std::cerr << "Failed to open file: " << filePath << std::endl;
    }

    return data;
}

// Function to read scale factors from a text file and return them as a vector of strings
std::vector<std::string> readScaleFactors(const std::string& filePath) {
    std::ifstream file(filePath);
    std::vector<std::string> scaleFactors;
    std::string line;

    if (file.is_open()) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::string factor;
        while (std::getline(iss, factor, ',')) {
            scaleFactors.push_back(factor);
        }
        file.close();
    } else {
        std::cerr << "Failed to open file: " << filePath << std::endl;
    }

    return scaleFactors;
}

// Function to rescale and translate data using scale and translation factors
std::vector<std::vector<double>> rescaleAndTranslateData(const std::vector<std::vector<double>>& data, const std::vector<std::string>& scaleFactors) {
    if (scaleFactors.size() != 6) {
        std::cerr << "Invalid number of scale factors." << std::endl;
        return data;
    }

    double scaleX = std::stod(scaleFactors[0]);
    double scaleZ = std::stod(scaleFactors[2]);
    double translateX = std::stod(scaleFactors[3]);
    double translateZ = std::stod(scaleFactors[5]);

    std::vector<std::vector<double>> result = data;

    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = 0; j < result[i].size(); j += 2) {
            result[i][j] = (result[i][j] - translateX) / scaleX;     // Scale x-values
            result[i][j + 1] = (result[i][j + 1] - translateZ) / scaleZ; // Scale z-values
        }
    }

    return result;
}


