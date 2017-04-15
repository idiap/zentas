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


The shared library is in ./build/zentas (libzentas.so in Linux) and the Python shared library is in ./build/python (pyzentas.so in Linux). These can be moved/copied elsewhere manually. 


Using
-------------------------------------------
Example use cases of the C++ library and headers are in testsexamples, with the corresponding executables in build/testsexamples. There is an example of clustering dense vectors (exdense.cpp), sparse vectors (exsparse.cpp), and strings (exwords.cpp).

To use the Python library, make sure pyzentas.so is on PYTHONPATH, for example you can use `sys.path.append(/path/to/pyzentas.so)`. There is only one clustering function for all data types, metrics etc. To use it pyzentas interactively (using iPython for example), try

```import pyzentas
pyzentas.pyzentas?
```

At the bottom of the function string are example uses. More examples are in python/examples.py


Doesn't work or missing a feature?
-------------------------------------------
Please raise an issue in the zentas repository
