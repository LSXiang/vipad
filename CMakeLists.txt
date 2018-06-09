##
# @File     CMakeLists.txt
# @Version  v1.0.0
# @Date     May,28 2018
# @Author   akatsuki (jacob.lsx)
# @Brief    Visual-Inertial Positioning with AprilTag for Drones
#
# Copyright (c) 2018, The Akatsuki. All rights reserved.
#

# set the minimum required version of cmake for the project
cmake_minimum_required(VERSION 2.8.3)

# project name 
project(vipad)

# set property
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -O3")
set(CMAKE_BUILD_TYPE "Debug")

# add path 
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# find OpenCV library
find_package(OpenCV REQUIRED)

# include header file path 
include_directories(
    ${OpenCV_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/apriltag/inc
    ${PROJECT_SOURCE_DIR}/src/ippe
    )

# add subdirectories
add_subdirectory(src)
add_subdirectory(app)

