add_subdirectory(${JSRCDIR_NAME})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++14")

SET(JPROJECT_LIBS ${JLIB1_NAME})
# We need this if bin != lib. Only works on linux & mac...
SET(CMAKE_INSTALL_RPATH "${JLIB_INSTALL_DIR}")
message(" the RPATH : ${CMAKE_INSTALL_RPATH} ")

add_executable(${JEXE1_NAME} main.cpp)

TARGET_LINK_LIBRARIES(${JEXE1_NAME}  ${OPENBLAS_LINK_DIR} ${JPROJECT_LIBS}) 
install(TARGETS ${JEXE1_NAME} DESTINATION ${JBIN_INSTALL_DIR})

