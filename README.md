What
-------------------------------------------
Fast and versatile CLARANS, as described at https://arxiv.org/abs/1609.04723

(1) Shared library with C++ headers

(2) Python library

Requirements
-------------------------------------------
Minimal installation requirements:

-- Cmake version 3.0 or greater


In addition, for the python library:

-- Python and Cython


Configuring
-------------------------------------------
If you do NOT want the Python library, comment out the final line in the topmost CMakeLists.txt file, so that it reads
`#add_subdirectory(python)`


Building
-------------------------------------------

```cmake 
mkdir build
cd build
cmake ..
make 
``` 


The shared library is in ./build/zentas (libzentas.so in Linux) and the Python shared library (pyzentas.so in Linux) is in ./build/python. These can be moved manually. 


Using
-------------------------------------------
Example use cases of the library in C++ are in testsexamples, with the corresponding executables in build/testsexamples. There are examples of clustering dense vectors, sparse vectors, and strings using the Levenshtein distance.

Example use cases of the library in Python are in python/examples.py. To use the Python library, make sure pyzentas.so is on PYTHONPATH, for example you can use `sys.path.append(/path/to/pyzentas.so)`



Doesn't work or missing a feature?
-------------------------------------------
Please raise an issue in the zentas repository
