# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

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
  
  ndata = int(2e3)
  X = mnist.read_MNIST(dataset = "original", ndata = ndata)/1000.
  dimension = X[0].shape[-1]
  npr.seed(1000)
  
  
  z = pyzentas.pyzen(K = 200, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_time = 0, max_rounds = 5000, seed = 1011, patient = True, nthreads = 1, init = "uniform", with_tests = False)
  do_vdimap = False
  tangerine =  z.den(X, do_vdimap, do_refinement = True)


def rna_example():
  """
  pass
  """
