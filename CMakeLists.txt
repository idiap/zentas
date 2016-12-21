# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.


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
