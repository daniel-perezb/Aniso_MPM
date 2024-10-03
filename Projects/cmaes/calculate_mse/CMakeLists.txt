cmake_minimum_required(VERSION 3.10)

# Project name
project(OpenFile)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)

# Find OpenCV
find_package(OpenCV REQUIRED)

# Include OpenCV headers
include_directories(${OpenCV_INCLUDE_DIRS})

# MSE executable
add_executable(mse main.cpp readf_module.cpp scale_module.cpp computeMSE.cpp findClosestPoint.cpp)
target_link_libraries(mse ${OpenCV_LIBS})

# Create_image executable
add_executable(create_image main_image.cpp readf_module.cpp)
target_link_libraries(create_image ${OpenCV_LIBS})
