cmake_minimum_required(VERSION 3.10)

# Project name
project(OpenFile)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)

# MSE executable
add_executable(open main.cpp scale_module.cpp readf_module.cpp computeMSE.cpp findClosestPoint.cpp)

# Create_image executabler
add_executable(create_image main_image.cpp readf_module.cpp)

# Additional libraries 
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# Link OpenCV to create_image:
target_link_libraries(create_image ${OpenCV_LIBS})
