file(GLOB source_files *cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++14")

#create the makefile for making the library
add_library(${JLIB1_NAME} SHARED ${source_files})
target_link_libraries (${JLIB1_NAME} ${CMAKE_THREAD_LIBS_INIT}) #pthread)

#how to install
install(TARGETS ${JLIB1_NAME} DESTINATION ${JLIB_INSTALL_DIR})
