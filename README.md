what
-------------------------------------------
Fast and versatile CLARANS, as described at https://arxiv.org/abs/1609.04723
(1) Shared library with C++ headers
(2) Python library

requirements
-------------------------------------------
Minimal installation requirements,
-- Cmake version 3.0 or greater

In addition, for the python version,
-- Python and Cython


configuring
-------------------------------------------
If you do NOT want the python version, comment out the final line in the topmost CMakeLists.txt file, so that it reads
`#add_subdirectory(python)`


building
-------------------------------------------
`mkdir build`
`cd ./build`
`cmake ..`
`make`
The shared library is then in ./build/zentas, and the python library (if built) is in ./build/python. These can be moved around.

using
-------------------------------------------
See testsexamples (with executables in build/testsexamples) python/examples and 


doesn't work?
-------------------------------------------
Please contact me at jnewling@idiap.ch or raise an issue.
