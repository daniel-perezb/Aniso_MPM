#include "functions.h"
#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>

// Parameters for the image
const int image_size =
    1200; // Image will be 1200x1200 pixels
          // Rescaling to original image will be applied later

const int point_radius =
    5; // Radius of the points in the image
       // A value of 5 to have a the gaps filled in
       // Value can be reduced if more than 50k points are used for mesh

void createBinaryImage(const std::vector<TV> &Xp_data, Viewpoint viewpoint,
                       const std::string &image_name) {
  // Create a binary image with all zeros (black)
  cv::Mat image = cv::Mat::zeros(image_size, image_size, CV_8UC1);

  // Find min/max values to normalize the points to the image size
  float minX = std::numeric_limits<float>::max();
  float maxX = std::numeric_limits<float>::lowest();
  float minY = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::lowest();

  for (const auto &point : Xp_data) {
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
    if (x < minX)
      minX = x;
    if (x > maxX)
      maxX = x;
    if (y < minY)
      minY = y;
    if (y > maxY)
      maxY = y;
  }

  // Scale and translate points to fit them into the image
  for (const auto &point : Xp_data) {
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
    int imgX =
        static_cast<int>(((x - minX) / (maxX - minX)) * (image_size - 1));
    int imgY =
        static_cast<int>(((y - minY) / (maxY - minY)) * (image_size - 1));

    // Draw the point as a small circle in the image
    cv::circle(image, cv::Point(imgX, imgY), point_radius, cv::Scalar(255), -1);
  }

  // Save the binary image to a file
  cv::imwrite(image_name, image);
}

int main() {

  // Create binary imagw from the first frame of animation
  std::vector<TV> Xp_data1;
  readF("output/RIG_fleece_fibres/data_0.dat", Xp_data1);
  createBinaryImage(Xp_data1, XZ_VIEW,
                    "../../Data/TetMesh/fleece_files/output_binary_image1.png");

  // Find last file in output
  int num_files = findLastFileNumber();
  std::string filename = "output/RIG_fleece_fibres/data_" + std::to_string(num_files) + ".dat";
  
  // Create binary imagw from the last frame of animation
  std::vector<TV> Xp_data2;
  readF(filename, Xp_data2);
  createBinaryImage(Xp_data2, XZ_VIEW,
                    "../../Data/TetMesh/fleece_files/output_binary_image2.png");

  std::cout << "Images created succesfully" << std::endl;
  return 0;
}
