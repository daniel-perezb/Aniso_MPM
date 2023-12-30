#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>

struct TV {
    float x, y, z;
};

double computeMSE(const std::vector<TV>& set1, const std::vector<TV>& set2);
std::pair<TV, int> findClosestPoint(const std::vector<double> &interestPoint, const std::vector<TV> &XP_data);
void readF(const std::string& filename, std::vector<TV>& Xp_data);
std::vector<std::vector<double>> readData(const std::string& filePath);
std::vector<std::string> readScaleFactors(const std::string& filePath);
std::vector<std::vector<double>> rescaleAndTranslateData(const std::vector<std::vector<double>>& data, const std::vector<std::string>& scaleFactors);


#endif
