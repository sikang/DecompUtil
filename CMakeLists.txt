cmake_minimum_required(VERSION 2.8.3)
project(decomp_util)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wno-deprecated-declarations")

find_package(Eigen3 REQUIRED)
include_directories(./include ${EIGEN3_INCLUDE_DIRS})

add_library(geometry_utils src/geometry_utils.cpp)

add_library(ellipse_decomp
	src/line_segment.cpp
	src/ellipse_decomp.cpp)
target_link_libraries(ellipse_decomp maxdet geometry_utils)
add_library(iterative_decomp src/iterative_decomp.cpp)
target_link_libraries(iterative_decomp ellipse_decomp geometry_utils)

add_library(ellipsoid_utils src/ellipsoid_utils.cpp)
target_link_libraries(ellipsoid_utils geometry_utils)
