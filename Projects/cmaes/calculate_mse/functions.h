#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <tuple>
#include <vector>

struct TV {
  float x, y, z;
};

struct CameraParameters {
  double camera_distance;
  double focal_length_mm;
  double sensor_width_mm;
  double sensor_height_mm;
};

// Enum for selecting the viewpoint
enum Viewpoint {
  XY_VIEW, // Top-down view (Z-axis projection)
  XZ_VIEW, // View from Y-axis
  YZ_VIEW  // View from X-axis
};

/**
Compute MSE between two set of points
@return MSE error
*/
double computeMSE(const std::vector<TV> &set1, const std::vector<TV> &set2);

/**
Find closest point between TAPIR and simulation frames
@return index of closest point in simulation
*/
std::pair<TV, int> findClosestPoint(const std::vector<double> &interestPoint,
                                    const std::vector<TV> &XP_data);

/**
Read file from a txt
@return data from txt
*/
std::vector<std::vector<double>> readFile(const std::string &filename);

/**
Read pointcloud data from the simulation
@return data from simulation files
*/
void readF(const std::string &filename, std::vector<TV> &Xp_data);

/**
Read the data from TAPIR files
@return data from TAPIR
*/
std::vector<std::vector<double>> readData(const std::string &dataFile,
                                          const std::string &visiblesFile, int num_files, 
                                          int &animation_frame_rate, int &TAPIR_frame_rate);

/**
Crop the data from TAPIR based on the frame rate from simulation and TAPIR
@return cropped data from TAPIR
*/
std::vector<std::vector<double>> cropData(std::vector<std::vector<double>> data,
                                          int &animation_frame_rate,
                                          int &TAPIR_frame_rate);

/**
Camera Params used for rescaling the data
@return camera params
*/
bool readCameraParameters(const std::string &filename,
                          CameraParameters &params);


/**
Rescale data from TAPIR to match the simulation parameters
@return rescaled data
*/
std::vector<std::vector<double>>
rescaleData(const std::vector<std::vector<double>> &data, const std::string &filename, const std::string &camera_filename);

/**
Find the last file in the output directory
@return last file num
*/
int findLastFileNumber(); 

#endif
