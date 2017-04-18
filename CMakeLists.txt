cmake_minimum_required(VERSION 2.8.3)
project(decomp_util)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wno-deprecated-declarations")

set(USE_ROS TRUE)


IF(USE_ROS)
  message(STATUS "Compile ROS utils.")
  # roslib cmake_modules
  FIND_PACKAGE(catkin REQUIRED COMPONENTS roscpp)
  catkin_package(
    DEPENDS Eigen3 
    CATKIN_DEPENDS 
    INCLUDE_DIRS include
    LIBRARIES ellipse_decomp ellipsoid_utils iterative_decomp
    )
  INCLUDE_DIRECTORIES(${catkin_INCLUDE_DIRS})
  LIST(APPEND LINK_LIBS ${catkin_LIBRARIES})
ENDIF()

find_package(Eigen3 REQUIRED)
include_directories(./include ${EIGEN3_INCLUDE_DIRS})

add_library(geometry_utils src/geometry_utils.cpp)

add_library(ellipse_decomp
	src/line_segment.cpp
	src/ellipse_decomp.cpp)
target_link_libraries(ellipse_decomp geometry_utils)
add_library(iterative_decomp src/iterative_decomp.cpp)
target_link_libraries(iterative_decomp ellipse_decomp geometry_utils)

add_library(ellipsoid_utils src/ellipsoid_utils.cpp)
target_link_libraries(ellipsoid_utils geometry_utils)
