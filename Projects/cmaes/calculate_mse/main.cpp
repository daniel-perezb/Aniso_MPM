#include "functions.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

int main() {
  // Using the scale module
  std::string dataFilePath =
      "../../Data/TetMesh/mesh_files/initial_final_points.txt";
  std::string scaleFilePath = "../../Data/TetMesh/mesh_files/scale_values.txt";
  std::vector<std::vector<double>> data = readData(dataFilePath);
  std::vector<std::string> scaleFactors = readScaleFactors(scaleFilePath);
  std::vector<std::vector<double>> rescaledData =
      rescaleAndTranslateData(data, scaleFactors);

  std::vector<TV> Xp_data;
  std::string filename = "../anisofracture/output/RIG_fleece_fibres/data_0.dat";
  readF(filename, Xp_data);

  // Find the closest point to each point in rescaledData
  std::vector<int> indexOfClosestPoints;
  // Array to save files
  std::vector<TV> all_data;

  for (const auto &row : rescaledData) {
    if (row.size() >= 2) { // Ensure there are at least two values in the row
      std::vector<double> interestPoint = {row[0], row[1]};
      auto result = findClosestPoint(interestPoint, Xp_data);
      TV closestPoint = result.first;
      int indexOfClosestPoint = result.second;

      // Store the indexOfClosestPoint value
      indexOfClosestPoints.push_back(indexOfClosestPoint);

      // Append the closest point to the vector
      all_data.push_back(closestPoint);
    }
  }

  // Value below is obtained using the number of frames desired for the output simulation
  int num_files = 7;

  std::vector<TV> new_data;
  std::vector<TV> final_values;
  // Using the readF module

  std::vector<TV> currentXp_data;
  std::string final_filename =
      "../anisofracture/output/RIG_fleece_fibres/data_" + std::to_string(num_files) + ".dat";
  readF(final_filename, currentXp_data);
  // Use the stored indexOfClosestPoint values to extract corresponding values

  for (int index : indexOfClosestPoints) {
    TV closestPoint = currentXp_data[index];
    // Append the closest point to the vector
    new_data.push_back(closestPoint);
  }

  int row_index = 2 * (2 - 1);
  for (const auto &row : rescaledData) {
    if (row_index + 1 < row.size()) {
      TV video_point = {row[row_index], 1.84, row[row_index + 1]};
      final_values.push_back(video_point);
    }
  }

  double mse_final = computeMSE(new_data, final_values);
  std::cout << "MSE error: " << mse_final << std::endl;

  // Save to a txt file
  std::ofstream outFile("mse.txt");

  if (outFile.is_open()) {
    outFile << mse_final;
    outFile.close();
    std::cout << "MSE written to txt" << std::endl;

  } else {
    std::cout << "Unable to open file" << std::endl;
  }

  return 0;
}
