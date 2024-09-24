#include "functions.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>

// Parameters for the image
const int image_size = 1200; // Image will be 1200x1200 pixels
const int point_radius = 1;  // Radius of the points in the image

// Enum for selecting the viewpoint
enum Viewpoint {
    XY_VIEW, // Top-down view (Z-axis projection)
    XZ_VIEW, // View from Y-axis
    YZ_VIEW  // View from X-axis
};

void createBinaryImage(const std::vector<TV>& Xp_data, Viewpoint viewpoint) {
    // Create a binary image with all zeros (black)
    cv::Mat image = cv::Mat::zeros(image_size, image_size, CV_8UC1);

    // Find min/max values to normalize the points to the image size
    float minX = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float minY = std::numeric_limits<float>::max();
    float maxY = std::numeric_limits<float>::lowest();

    for (const auto& point : Xp_data) {
        float x, y;
        
        // Select coordinates based on the viewpoint
        switch (viewpoint) {
            case XY_VIEW:
                x = point.x;
                y = point.y;
                break;
            case XZ_VIEW:
                x = point.x;
                y = point.z;
                break;
            case YZ_VIEW:
                x = point.y;
                y = point.z;
                break;
            default:
                x = point.x;
                y = point.y;
        }

        // Update min/max based on the selected axes
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    // Scale and translate points to fit them into the image
    for (const auto& point : Xp_data) {
        float x, y;
        
        switch (viewpoint) {
            case XY_VIEW:
                x = point.x;
                y = point.y;
                break;
            case XZ_VIEW:
                x = point.x;
                y = point.z;
                break;
            case YZ_VIEW:
                x = point.y;
                y = point.z;
                break;
            default:
                x = point.x;
                y = point.y;
        }

        // Normalize the points to the image size
        int imgX = static_cast<int>(((x - minX) / (maxX - minX)) * (image_size - 1));
        int imgY = static_cast<int>(((y - minY) / (maxY - minY)) * (image_size - 1));

        // Draw the point as a small circle in the image
        cv::circle(image, cv::Point(imgX, imgY), point_radius, cv::Scalar(255), -1);
    }

    // Save the binary image to a file
    cv::imwrite("../binary_image.png", image);
}

int main() {
    std::vector<TV> Xp_data;
    readF("../data_7.dat", Xp_data);
    
    // Choose the viewpoint: XY_VIEW (top-down), XZ_VIEW (from side Y), YZ_VIEW (from side X)
    createBinaryImage(Xp_data, XY_VIEW); 
    
    return 0;
}
