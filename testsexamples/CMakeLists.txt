SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pedantic" )

include_directories("../zentas/include")

add_executable(exdense exdense.cpp)
target_link_libraries(exdense LINK_PUBLIC zentas)


add_executable(exsparse exsparse.cpp)
target_link_libraries(exsparse LINK_PUBLIC zentas)


add_executable(exwords exwords.cpp)
target_link_libraries(exwords LINK_PUBLIC zentas)

add_executable(exvdimap exvdimap.cpp)
target_link_libraries(exvdimap LINK_PUBLIC zentas)


#add_executable(
#deepbench deepbench.cpp)

#add_executable(
#backconvwrw backconvwrw.cpp)

#add_executable(
#smallgeometrytests smallgeometrytests.cpp)

#add_executable(
#redirectionexample redirectionexample.cpp)

#add_executable(
#devtest devtest.cpp)



#target_link_libraries(
#devtest LINK_PUBLIC tinygemm ${OpenCL_LIBRARY})


#target_link_libraries(
#deepbench LINK_PUBLIC tinygemm ${OpenCL_LIBRARY})

#target_link_libraries(
#backconvwrw LINK_PUBLIC tinygemm ${OpenCL_LIBRARY})

#target_link_libraries(
#redirectionexample LINK_PUBLIC tinygemm ${OpenCL_LIBRARY})

#target_link_libraries(
#smallgeometrytests LINK_PUBLIC tinygemm ${OpenCL_LIBRARY})
