cmake_minimum_required(VERSION 3.0)

project(gnr1d)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)

# Look for OpenMP package
find_package(OpenMP REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -O3 ${OpenMP_CXX_FLAGS} -pg")

set(SOURCES
  ./Model_wall.hpp 
  ./Model_wall.cpp 
  ./Sys_Tools.hpp 
  ./Sys_Tools.cpp 
  ./Time_solver.hpp 
  ./Time_solver.cpp
)

include_directories(.)

add_library(gnr1d_lib ${SOURCES})

add_executable(gnr1d main.cpp)

target_link_libraries(gnr1d gnr1d_lib ${OpenMP_CXX_LIBRARIES})

set_target_properties(gnr1d PROPERTIES LINK_FLAGS "-pg")
