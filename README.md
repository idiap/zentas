WHAT:
-------------------------------------------
Accelerated CLARANS as described at https://arxiv.org/abs/1609.04723
(1) Linked library with C++ headers
(2) Python library

REQUIREMENTS:
-------------------------------------------
Minimal installation requirements:
-- cmake version 3.0 or greater

In addition, for the python version:
-- Python and Cython


CONFIGURATION:
-------------------------------------------
If you do NOT want the python version, comment out the final line in the topmost CMakeLists.txt file, so that it reads
#add_subdirectory(python)


BUILDING:
-------------------------------------------
mkdir build
cd ./build
cmake ..
make
The shared library is then in ./build/zentas, and the python library is in ./build/python. 

USING: 
-------------------------------------------
See python/examples and testsexamples (with executables in build/testsexamples).


DOESN'T WORK?
-------------------------------------------
Please contact me at jnewling@idiap.ch or raise an issue
