# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.


cmake_minimum_required(VERSION 3.0)
project(zentas)

set(${PROJECT_NAME}_VERSION 0.1)
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)

find_package (Threads REQUIRED)

option(BUILD_SHARED_LIBS "Build shared library" ON)

SET (CMAKE_CXX_COMPILER "clang++")

add_definitions(-DPATH_DATADIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")

add_subdirectory(zentas)
add_subdirectory(testsexamples)



# python code. It requires cython (and only builds with certain versions). You probably want to comment this out.
add_subdirectory(python)

#This flag should not really be set in CMakeLists.txt.
SET(CMAKE_BUILD_TYPE Release)



#cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++11 ")

##The name of the project is not important, it will only appear in fancy IDEs
#project(ProjectTemplate)



#find_package (Threads REQUIRED)

##This flag should not really be set in CMakeLists.txt.
#SET(CMAKE_BUILD_TYPE Release)
## The bin and lib dirs where all executables and libraries will be put will be `build'/cleaninstall
#SET(JINSTALL_DIR cleaninstall)
#SET(JBINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
#SET(WITHBLAS "1101" CACHE STRING "noob")
#SET(CMAKE_INSTALL_PREFIX  ${JBINARY_DIR}/${JINSTALL_DIR})
#message("Installation directory : ${CMAKE_INSTALL_PREFIX}")

#MACRO(SUBDIRLIST result curdir)
  #FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  #SET(dirlist "")
  #FOREACH(child ${children})
    #IF(IS_DIRECTORY ${curdir}/${child})
      #LIST(APPEND dirlist ${child})
    #ENDIF()
  #ENDFOREACH()
  #SET(${result} ${dirlist})
#ENDMACRO()

#SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)

#LIST(LENGTH SUBDIRS NSUBDIRS)

#if(NSUBDIRS EQUAL 1)
  #SET(JLIB1_NAME ${SUBDIRS})
  #SET(JEXE1_NAME x${JLIB1_NAME})
#else()
  #error("More than 1 subdirectory of source, cannot determine the correct name for the library")
#endif()
#file(GLOB PYX_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/python/*pyx)
#get_filename_component(JLIBPY_NAME ${PYX_FILENAME} NAME_WE)
#SET(JLIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
#SET(JBIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
#SET(JSRCDIR_NAME ${JLIB1_NAME})


#include_directories(include)
#add_subdirectory(include)
#add_subdirectory(src)

##comment out this line if you don't want the python version:
#add_subdirectory(python)

