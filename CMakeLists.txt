cmake_minimum_required(VERSION 2.8.3)
project(decomp_util)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations")

set(BUILD_SHARED_LIBS ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS} include)

#add_library(ellipsoid_decomp
#  src/line_segment.cpp
#  src/ellipse_decomp.cpp)

#add_library(iterative_decomp src/iterative_decomp.cpp)
#target_link_libraries(iterative_decomp ellipse_decomp)

#add_library(ellipsoid_utils src/ellipsoid_utils.cpp)
#target_link_libraries(ellipsoid_utils geometry_utils)

include(CTest)

add_executable(test_seed_decomp test/test_seed_decomp.cpp)
target_link_libraries(test_seed_decomp ${Boost_LIBRARIES})
add_test(test_seed_decomp test_seed_decomp)

install(FILES "${PROJECT_NAME}Config.cmake" "${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "share/${PROJECT_NAME}/cmake")

install(TARGETS
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY include/decomp_basis include/decomp_geomtry include/decomp_util
  DESTINATION include)
