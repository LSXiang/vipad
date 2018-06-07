#file : /src/ippe/CMakeLists.txt

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ippe_src)
add_library(ippe_lib ${ippe_src})
