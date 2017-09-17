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
  pass

import time
import rna
reload(rna)
K = 200
X = rna.get_rna()[0:10000, :]
from sklearn.cluster import KMeans
t0 = time.time()

print "running skl..."
sklc = KMeans(n_clusters = K, init = "k-means++", max_iter = 7, tol = 1e-9, algorithm = "full")
sklc.fit(X)
print sklc.score(X)
print "done."

t1 = time.time()
print "sklearn took : " , t1 - t0


print "running zentas..."
z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_time = 0, max_rounds = 5000, seed = 1011, patient = True, nthreads = 3, init = "kmeans++-5", with_tests = False, capture_output = True)
tangerine =  z.den(X, do_vdimap = False, do_refinement = True)
print "done."

t2 = time.time()

print "zentas took : " , t2 - t1

