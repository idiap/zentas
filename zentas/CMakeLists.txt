#include(EnableCompilerWarnings)
add_definitions(-std=c++11)

include_directories("include")
file(GLOB_RECURSE source_files src/*cpp)
add_library("zentas" SHARED ${source_files})

# Inspiration from https://cmake.org/examples/
target_include_directories (zentas PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(zentas  ${CMAKE_THREAD_LIBS_INIT})
