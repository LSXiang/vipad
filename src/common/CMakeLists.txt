# file: /src/common/CMakeLists.txt

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} common_src)
add_library(common_lib ${common_src})
