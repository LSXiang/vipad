##
# @File     app/CMakeLists.txt
# @Version  v1.0.0
# @Date     May,22 2018
# @Author   akatsuki (jacob.lin)
# @Brief    
#
# Copyright (c) 2018, The Akatsuki. All rights reserved.
#

# create executable file 
add_executable(locationEstimation main.cpp)

target_link_libraries(locationEstimation 
    ${OpenCV_LIBS}
    AprilTag_lib
    ippe_lib
    vipad_lib
    )

# create executable file 
add_executable(lpe_test test.cpp)

target_link_libraries(lpe_test 
    ${OpenCV_LIBS}
    AprilTag_lib
    ippe_lib
    vipad_lib
    )
    
