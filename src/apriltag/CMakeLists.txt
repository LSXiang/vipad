##
# @File     src/apriltag/CMakeLists.txt
# @Version  v1.0.0
# @Date     May,22 2018
# @Author   akatsuki (jacob.lin)
# @Brief    
#
# Copyright (c) 2018, The Akatsuki. All rights reserved.
#

# add header file path
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

# list apriltag directory all source file 
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src apriltag_src)
# message(STATUS "${apriltag_src}")

# create ArUco library 
add_library(AprilTag_lib SHARED STATIC ${apriltag_src})
