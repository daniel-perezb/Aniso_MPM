#include "functions.h"
#include <stdexcept>

double computeMSE(const std::vector<TV>& set1, const std::vector<TV>& set2) {
    
    if (set1.size() != set2.size()) {
        throw std::invalid_argument("The two sets must be of the same size.");
    }

    double mse = 0.0;
    for (size_t i = 0; i < set1.size(); ++i) {
        double dx = set1[i].x - set2[i].x;
        double dy = set1[i].y - set2[i].y;
        double dz = set1[i].z - set2[i].z;

        mse += (dx * dx + dy * dy + dz * dz);
    }

    return mse / set1.size();


}