# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

"""
Further using pyzentas.
"""

import sys
import random
import numpy as np
import numpy.random as npr

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../../build/python")
import pyzentas

from IPython.core.debugger import Tracer 


def mnist_example():
  """
  Showing how do_vdimap can help. 
  """
  import mnist
  reload(mnist)
  
  ndata = int(1e4)
  X = mnist.read_MNIST(dataset = "original", ndata = ndata)
  dimension = X[0].shape[-1]
  npr.seed(1011)
  
  z = pyzentas.pyzen(K = 1e3, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_time = 10000, max_rounds = 4, seed = 1011, patient = True, nthreads = 1)
  do_vdimap = False
  tangerine =  z.den(X, do_vdimap)
