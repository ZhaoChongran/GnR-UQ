PROJECT(gnr1d_cpp)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

SET( CMAKE_CXX_COMPILER "g++")

SET( CMAKE_CXX_FLAGS "-g -Wall -O3")

SET( file ./Model_wall.cpp 
  ./Time_solver.cpp ./Sys_Tools.cpp )

INCLUDE_DIRECTORIES(.)

ADD_LIBRARY( mat ${file} )

ADD_EXECUTABLE( gnr1d main.cpp )

#ADD_EXECUTABLE( test test.cpp )

TARGET_LINK_LIBRARIES( gnr1d mat )

#TARGET_LINK_LIBRARIES( test mat )

# EOF
