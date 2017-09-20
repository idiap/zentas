####################################################
#   Copyright (c) 2016 Idiap Research Institute    #
#   Written by James Newling <jnewling@idiap.ch>   #
####################################################

cmake_minimum_required(VERSION 3.0)


#preferred compiler. TODO, make this optional
# (see https://cmake.org/Wiki/CMake_FAQ#How_do_I_use_a_different_compiler.3F).

#SET (CMAKE_CXX_COMPILER "clang++")
#SET (CMAKE_CXX_COMPILER "g++")

project(zentas)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(${PROJECT_NAME}_VERSION 0.1)
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)

option(BUILD_PYTHON_LIB "Build Python library" ON)
option(BUILD_R_LIB "Build R library (NOT yet functioning)" OFF)
option(BUILD_SHARED_LIBS "Build shared library" ON)

set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." )

#Where to find FindCython.cmake
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )



find_package (Threads REQUIRED)

add_definitions(-DPATH_DATADIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
add_subdirectory(zentas)
add_subdirectory(testsexamples)

if (BUILD_PYTHON_LIB)
  find_package(Cython REQUIRED)
  add_subdirectory(python)
endif()

if (BUILD_R_LIB)
  add_subdirectory(R)
endif()
