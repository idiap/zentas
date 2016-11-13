SET(CYTHONGEN_CPP ${JBINARY_DIR}/cythongen/${JLIBPY_NAME}.cpp)
SET(CYTHON_PYX ${CMAKE_CURRENT_SOURCE_DIR}/${JLIBPY_NAME}.pyx)

file(MAKE_DIRECTORY ${JBINARY_DIR}/cythongen)

SET(CMAKE_INSTALL_RPATH "${JLIB_INSTALL_DIR}")
message(" the RPATH : ${CMAKE_INSTALL_RPATH} ")

add_library(${JLIBPY_NAME} SHARED ${CYTHONGEN_CPP})

SET_TARGET_PROPERTIES(${JLIBPY_NAME} PROPERTIES LINKER_LANGUAGE CXX PREFIX "")

#How and when to cythonise
add_custom_command(
  OUTPUT ${CYTHONGEN_CPP}
  DEPENDS ${CYTHON_PYX}
  COMMAND cython --cplus ${CYTHON_PYX} -o ${CYTHONGEN_CPP}
  )

#The python library relies on the main library
TARGET_LINK_LIBRARIES(${JLIBPY_NAME} ${JLIB1_NAME})

#This is linux specific! Getting the include and link files for making a python library.
execute_process(COMMAND python-config --cflags OUTPUT_VARIABLE CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE "-Wall" "" CFLAGS ${CFLAGS})

SET(CMAKE_CXX_FLAGS  " ${CFLAGS} --std=c++11 ")

message("(In python dir) CMAKE_CXX_FLAGS set to : ${CFLAGS}")

include_directories(include)

install(TARGETS ${JLIBPY_NAME} DESTINATION lib)
