##
# @File     src/CMakeLists.txt
# @Version  v1.0.0
# @Date     May,22 2018
# @Author   akatsuki (jacob.lin)
# @Brief    
#
# Copyright (c) 2018, The Akatsuki. All rights reserved.
#

add_subdirectory(apriltag)
add_subdirectory(ippe)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
add_library(vipad_lib ${source_files})
target_link_libraries(vipad_lib AprilTag_lib ippe_lib)
