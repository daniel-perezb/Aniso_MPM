#include "functions.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

// Function to find the last file number instead of hyperparameter
int findLastFileNumber() {
  int lastFileNumber = 0;
  for (int i = 0; i < 1000; i++) {
    std::string filename =
        "output/RIG_fleece_fibres/data_" + std::to_string(i) + ".dat";
    std::ifstream file(filename);
    if (!file) {
      lastFileNumber = i - 1;
      break;
    }
  }
  return lastFileNumber;
}

int main() {

  // Points below are obtained using TAPIR
  std::string dataFilePath = "../../Data/TetMesh/fleece_files/data.txt";
  std::string visiblesFilePath = "../../Data/TetMesh/fleece_files/visibles.txt";
  
  // Read the data and remove points that are not visible at some point in the
  // video
  std::vector<std::vector<double>> data =
      readData(dataFilePath, visiblesFilePath);

  // Crop frame rate based on animation frames and TAPIR frames
  int animation_frame_rate = 1;
  int TAPIR_frame_rate = 10;
  std::vector<std::vector<double>> croppedData =
      cropData(data, animation_frame_rate, TAPIR_frame_rate);

  // Rescale the data using ooutpput of py file to create stl
  std::vector<std::vector<double>> rescaledData = rescaleData(croppedData, "../../Data/TetMesh/fleece_files/camera_params.txt", "../../Data/TetMesh/fleece_files/initial_mask.png");

  // Read the data from the output of the simulation
  std::vector<TV> Xp_data;
  std::string filename = "output/RIG_fleece_fibres/data_0.dat";
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

  // Value below is obtained using the number of frames desired for the output
  // simulation
  int num_files = findLastFileNumber();

  std::vector<TV> new_data;
  std::vector<TV> final_values;
  std::vector<TV> currentXp_data;

  std::string final_filename =
      "output/RIG_fleece_fibres/data_" + std::to_string(num_files) + ".dat";
  readF(final_filename, currentXp_data);

  // Use the stored indexOfClosestPoint values to extract corresponding values
  for (int index : indexOfClosestPoints) {
    TV closestPoint = currentXp_data[index];
    // Append the closest point to the vector
    new_data.push_back(closestPoint);
  }

  int row_index = 2;
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
