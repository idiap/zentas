include(EnableCompilerWarnings)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11  -Wall -Wextra -pedantic " )

include_directories("include")

file(GLOB_RECURSE source_files src/*cpp)#include/*)

add_library("zentas" SHARED ${source_files})



# Inspiration from https://cmake.org/examples/
target_include_directories (zentas PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # ${OPENCL_INCLUDE_DIRS})
target_link_libraries(zentas  ${CMAKE_THREAD_LIBS_INIT})
