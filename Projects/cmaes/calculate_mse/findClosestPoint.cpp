#include "functions.h"
#include <cmath>
#include <limits>

std::pair<TV, int> findClosestPoint(const std::vector<double> &interestPoint, const std::vector<TV> &XP_data) {
    double minDistance = std::numeric_limits<double>::max();
    TV closestPoint;
    int index = -1;

    for (size_t i = 0; i < XP_data.size(); i++)
    {
        double dx = interestPoint[0] - XP_data[i].x;
        double dz = interestPoint[1] - XP_data[i].z;
        double distance = dx * dx + dz * dz;

        if (distance < minDistance)
        {
            minDistance = distance;
            closestPoint = XP_data[i];
            index = i;
        }
    }

    return {closestPoint, index};
}
