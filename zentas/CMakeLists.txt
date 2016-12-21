SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11  -Wall -Wextra -pedantic " )

include_directories("include")

file(GLOB_RECURSE source_files src/* include/*)

add_library("zentas" SHARED ${source_files})



# Inspiration from https://cmake.org/examples/
target_include_directories (zentas PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # ${OPENCL_INCLUDE_DIRS})

target_link_libraries(zentas  ${CMAKE_THREAD_LIBS_INIT})




# ${OPENCL_LIBRARIES})

#target_link_libraries (${JLIB1_NAME} ${CMAKE_THREAD_LIBS_INIT}) #pthread)

#python/* opencl/*

# ${OpenCL_INCLUDE_DIR})
#add_definitions(-DPATH_MAKE_KERNEL_CMDL_PY="${CMAKE_CURRENT_SOURCE_DIR}/python/make_kernel_cmdl.py")
#add_definitions(-DPATH_BETAC_KERNEL_DIR="${CMAKE_CURRENT_SOURCE_DIR}/opencl/")

