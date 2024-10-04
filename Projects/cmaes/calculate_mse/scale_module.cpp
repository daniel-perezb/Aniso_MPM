#include "functions.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> 
#include <opencv2/opencv.hpp>

using namespace std;

// Function to read camera parameters from a file
bool readCameraParameters(const std::string& filename, CameraParameters& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open camera parameters file: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool inCameraParametersSection = false;

    while (std::getline(file, line)) {
        // Trim whitespace from the beginning and end of the line
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), line.end());

        // Skip empty lines or comments
        if (line.empty() || line[0] == ';' || line[0] == '#') {
            continue;
        }

        // Check for section headers
        if (line[0] == '[' && line[line.size() - 1] == ']') {
            std::string section = line.substr(1, line.size() - 2);
            if (section == "CameraParameters") {
                inCameraParametersSection = true;
            } else {
                inCameraParametersSection = false;
            }
            continue;
        }

        if (!inCameraParametersSection) {
            continue;
        }

        // Split line into key and value
        size_t pos = line.find('=');
        if (pos == std::string::npos) {
            std::cerr << "Invalid line in camera parameters file: " << line << std::endl;
            continue;
        }

        // Extract key and value, and trim whitespace
        std::string key = line.substr(0, pos);
        std::string value_str = line.substr(pos + 1);

        key.erase(key.begin(), std::find_if(key.begin(), key.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        key.erase(std::find_if(key.rbegin(), key.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), key.end());

        value_str.erase(value_str.begin(), std::find_if(value_str.begin(), value_str.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        value_str.erase(std::find_if(value_str.rbegin(), value_str.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), value_str.end());

        // Convert value to double
        try {
            double value = std::stod(value_str);

            // Assign value to the corresponding parameter
            if (key == "camera_distance") {
                params.camera_distance = value;
            } else if (key == "focal_length_mm") {
                params.focal_length_mm = value;
            } else if (key == "sensor_width_mm") {
                params.sensor_width_mm = value;
            } else if (key == "sensor_height_mm") {
                params.sensor_height_mm = value;
            } else {
                std::cerr << "Unknown parameter: " << key << std::endl;
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid value for parameter '" << key << "': " << value_str << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Value out of range for parameter '" << key << "': " << value_str << std::endl;
        }
    }

    file.close();
    return true;
}

// Function to read data from a file and store it in a vector of vectors (2D vector)
vector<vector<double>> readFile(const string &filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
    vector<double> row;
    stringstream ss(line);
    double value;
    while (ss >> value) {
        row.push_back(value);
    }
    data.push_back(row);
    }
    return data;
}

// Function to filter the data from 'data.txt' based on 'visibles.txt' and return the filtered data
std::vector<std::vector<double>> readData(const std::string &dataFile,
                                          const std::string &visiblesFile, int num_files, 
                                          int &animation_frame_rate, int &TAPIR_frame_rate) {
  
    // Read both files into vectors of vectors
    vector<vector<double>> data = readFile(dataFile);
    vector<vector<double>> visibles = readFile(visiblesFile);

    // Vector to store the filtered data
    vector<vector<double>> filteredData;

    // Calculate the number of frames to check in visibles
    int total_visibles_frames = visibles[0].size();
    int frames_to_check = num_files * (TAPIR_frame_rate / animation_frame_rate);
    frames_to_check = std::min(frames_to_check, total_visibles_frames);
    
    // Iterate over each row in visibles
    for (size_t i = 0; i < visibles.size(); ++i) {
        bool keepRow = true;
        // Check the first frames_to_check values in visibles[i]
        for (int j = 0; j < frames_to_check; ++j) {
            if (visibles[i][j] == 0) {
                keepRow = false;
                break;
            }
        }
        // If the point is always visible in the first frames_to_check frames, keep the row
        if (keepRow) {
            filteredData.push_back(data[i]);
        }
    }

    return filteredData;
}
// Function to crop the data based on the frame rate of both videos
vector<vector<double>> cropData(vector<vector<double>> data,
                                int &animation_frame_rate,
                                int &TAPIR_frame_rate) {
    // Vector to store cropped data
    vector<vector<double>> croppedData;

    // Iterate over each row of data
    for (const auto &row : data) {
    vector<double> croppedRow;

    // Extract points based on frame rates
    for (size_t i = 0; i < row.size();
            i += 2 * TAPIR_frame_rate) { 
        if (i + 1 < row.size()) {      
        croppedRow.push_back(row[i]); 
        croppedRow.push_back(row[i + 1]); 
        }
    }

    croppedData.push_back(croppedRow);
    }

    return croppedData;
}


// Function to rescale data
std::vector<std::vector<double>> rescaleData(
    const std::vector<std::vector<double>> &data,
    const std::string &filename, const std::string &camera_filename) {

    std::vector<std::vector<double>> rescaledData;

    // Read camera parameters
    CameraParameters camParams;
    if (!readCameraParameters(filename, camParams)) {
        std::cerr << "Error reading camera parameters." << std::endl;
        return rescaledData;
    }

    // Read image to get dimensions
    cv::Mat image = cv::imread(camera_filename, cv::IMREAD_GRAYSCALE);
    if (image.empty()) {
        std::cerr << "Failed to read image." << std::endl;
        return rescaledData;
    }

    int image_width = image.cols;
    int image_height = image.rows;

    // Define target centroid coordinates
    double target_centroid_x = 2.1;
    double target_centroid_y = 1.832;
    double target_centroid_z = 2.0;

    // Compute scaling factors
    double focal_length_m = camParams.focal_length_mm / 1000.0;
    double real_width = (camParams.camera_distance * camParams.sensor_width_mm) / camParams.focal_length_mm;
    double real_height = (camParams.camera_distance * camParams.sensor_height_mm) / camParams.focal_length_mm;

    double voxel_size_x = real_width / image_width;
    double voxel_size_y = real_height / image_height;

    // Thickness and voxel size z
    double thickness = 0.05; // Desired thickness in meters
    double extrusion_height = 50.0; // As the Python code
    double voxel_size_z = thickness / extrusion_height;

    // Define rotation angles in degrees
    double rotation_deg_x = -90.0;
    double rotation_deg_y = 180.0;
    double rotation_deg_z = 180.0;

    // Convert degrees to radians
    double rotation_rad_x = rotation_deg_x * M_PI / 180.0;
    double rotation_rad_y = rotation_deg_y * M_PI / 180.0;
    double rotation_rad_z = rotation_deg_z * M_PI / 180.0;

    // Create rotation matrices
    cv::Matx33d Rx(1, 0, 0,
                    0, cos(rotation_rad_x), -sin(rotation_rad_x),
                    0, sin(rotation_rad_x), cos(rotation_rad_x));

    cv::Matx33d Ry(cos(rotation_rad_y), 0, sin(rotation_rad_y),
                    0, 1, 0,
                    -sin(rotation_rad_y), 0, cos(rotation_rad_y));

    cv::Matx33d Rz(cos(rotation_rad_z), -sin(rotation_rad_z), 0,
                    sin(rotation_rad_z), cos(rotation_rad_z), 0,
                    0, 0, 1);

    // Combined rotation matrix
    cv::Matx33d R = Rz * Ry * Rx;

    // Variables for centroid computation
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    int total_points = 0;

    std::vector<cv::Vec3d> all_scaled_points;

    for (const auto& row : data) {
        for (size_t i = 0; i + 1 < row.size(); i += 2) {
            double x = row[i];
            double y = row[i + 1];
            double z = extrusion_height / 2.0; // Midpoint of extrusion height

            // Apply scaling
            double rescaled_x = x * voxel_size_x;
            double rescaled_y = y * voxel_size_y;
            double rescaled_z = z * voxel_size_z;

            // Collect scaled points
            all_scaled_points.emplace_back(rescaled_x, rescaled_y, rescaled_z);

            // Accumulate for centroid computation
            sum_x += rescaled_x;
            sum_y += rescaled_y;
            sum_z += rescaled_z;
            total_points++;
        }
    }

    // Compute centroid
    double centroid_x = sum_x / total_points;
    double centroid_y = sum_y / total_points;
    double centroid_z = sum_z / total_points;

    // Process each row individually
    int point_index = 0; // Index into all_scaled_points

    for (const auto& row : data) {
        std::vector<double> rescaledRow;
        for (size_t i = 0; i + 1 < row.size(); i += 2) {
            // Get the scaled point
            cv::Vec3d scaled_pt = all_scaled_points[point_index];

            // Center the point
            cv::Vec3d centered_pt = scaled_pt - cv::Vec3d(centroid_x, centroid_y, centroid_z);

            // Apply rotation
            cv::Vec3d rotated_pt = R * centered_pt;

            // Apply translation
            double final_x = rotated_pt[0] + target_centroid_x;
            double final_y = rotated_pt[1] + target_centroid_y;
            double final_z = rotated_pt[2] + target_centroid_z;

            // Store the final x and z values
            rescaledRow.push_back(final_x);
            rescaledRow.push_back(final_z);

            point_index++;
        }
        rescaledData.push_back(rescaledRow);
    }

    return rescaledData;
}
