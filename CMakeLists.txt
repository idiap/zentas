# Copyright (c) 2016 Idiap Research Institute
# Written by James Newling <jnewling@idiap.ch>

cmake_minimum_required(VERSION 3.0)
project(zentas)

set(${PROJECT_NAME}_VERSION 0.1)
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)

# this flag should not really be set in CMakeLists.txt.
SET(CMAKE_BUILD_TYPE Release)

find_package (Threads REQUIRED)

option(BUILD_SHARED_LIBS "Build shared library" ON)

# optional, comment out if clang is not available
SET (CMAKE_CXX_COMPILER "clang++")

add_definitions(-DPATH_DATADIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")

add_subdirectory(zentas)
add_subdirectory(testsexamples)

# python code. It requires cython (and only builds with certain versions). You might want to comment this out.
add_subdirectory(python)

#R code. You might want to comment this out.
#add_subdirectory(R)
