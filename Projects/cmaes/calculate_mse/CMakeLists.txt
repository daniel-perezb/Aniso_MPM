cmake_minimum_required(VERSION 3.10)
project(OpenFile)

# Set C++ version
set(CMAKE_CXX_STANDARD 14)

# Add the main executable that uses the two modules
add_executable(open main.cpp scale_module.cpp readf_module.cpp computeMSE.cpp findClosestPoint.cpp)


