# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

"""
Some simple test cases
"""

import sys
import random
import numpy as np
import numpy.random as npr

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../build/python")
import pyzentas

def from_file_example():
  """
  read lists of words from files and cluster using Levenshtein
  """
  seed =  1012
  random.seed(seed)
  npr.seed(seed)
  max_proposals = 100000
  maxtime = 100.0
  K = 3
  maxrounds = 10
  
    
  print "string settings (1)"
  root = "../data/"
  cl1 = pyzentas.pyzentas(filenames_list = [root + "words1.txt", root + "words2.txt"], outfilename = root + "output1.txt", costfilename = root + "costs.txt", K = K, algorithm = "clarans", level = 1, max_proposals = 100, capture_output = False, maxtime = 10, metric = "levenshtein", nthreads = 1, maxrounds = 10, patient = False, energy = "identity", rooted = True, seed = seed)
  





def generated_sequence_example():
  seed =  1011 #npr.randint(1000)
  print "seed :", seed
  random.seed(seed)
  npr.seed(seed)
  max_proposals = 100000
  maxtime = 100.0
  ndata = 10000
  K = 120
  maxrounds = 10
  sizes = np.array(npr.randint(2, 8, size = ndata), dtype = np.uint64)
  data = []
  for i in range(sizes.sum()):
    if npr.rand() < 0.25:
      data.append('A')
    elif npr.rand() < 0.5:
      data.append('C')
    elif npr.rand() < 0.75:
      data.append('G')
    else:
      data.append('T')
  data = np.array(data, dtype = 'c')
  #data.dtype = np.int8 #or keep it as |s1 (c). 
  
  
  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()
  
    
  print "string settings (1)"
  cl1 = pyzentas.pyzentas(ndata, dimension = None, sizes = sizes, X = data, K = K, indices_init = indices_init, algorithm = "clarans", level = 1, max_proposals = max_proposals, capture_output = False, seed = 1011, maxtime = maxtime, nthreads = 1, maxrounds = maxrounds, patient = False, metric = "normalised levenshtein", rooted = False, energy = 'identity', with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None)
  
  


def dense_vector_example():
  seed =  npr.randint(10000)
  random.seed(seed)
  npr.seed(seed)
  max_proposals = 1000#000000
  maxtime = 60.0

  ndata = 50000
  K = 500
  maxrounds = 2
  dimension = 2

  basedata = 1*npr.randn(ndata, dimension)
  data = np.array(basedata, dtype = np.float64)
  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()
  
  print " settings (1) "

  pyzentas.pyzentas(ndata = ndata, dimension = dimension, sizes = None, X = data, K = K, indices_init = indices_init, algorithm = "clarans", level = 0, max_proposals = max_proposals, capture_output = False, seed = seed, maxtime = maxtime, nthreads = 3, maxrounds = maxrounds, patient = True, metric = "l2", rooted = False, energy = 'quadratic') #, energy = 'squarepotential', critical_radius = 5.) #, energy = 'exp', exponent_coeff = 0.1)#

  


def sparse_data_example():
  seed =  1011 #npr.randint(1000)
  print "seed :", seed
  random.seed(seed)
  npr.seed(seed)
  max_proposals = 1000
  maxtime = 100.0
  maxrounds = 100
  ndata = 5
  K = 2
  sizes = np.array([2,3,2,1,2], dtype = np.uint64)
  indices_s = np.array([1,2,
  1,2,3,
  8,9,
  8,
  8,9000], dtype = np.uint64)
  data = np.ones(sizes.sum()) #npr.randn(sizes.size)  
  indices_init = np.array([0,1], dtype = np.uint64)
  indices_init.sort()
  
    
  #cl1 = pyzentas.pyzentas(ndata = ndata, sizes = sizes, X = data, K = K, indices_init = indices_init, algorithm = "clarans", level = 3, max_proposals = max_proposals, capture_output = False, seed = seed, maxtime = maxtime, metric = "l0", nthreads = 1, maxrounds = maxrounds, patient = False, energy = "identity", rooted = False, indices_s = indices_s)

  z = pyzentas.pyzen(K)
  return z.spa(sizes, indices_s, data)  
  
  #return cl1
  
  
  
  
  
def function_string_example():
  """
  Thus is the example in the function string
  """
  import numpy as np
  import numpy.random as npr
  ndata = 100000
  dimension = 10
  X = npr.randn(ndata, dimension)
  K = 500
  indices_init = range(K)
  results = pyzentas.pyzentas(X = X, K = K, indices_init = indices_init, maxrounds = 200000, maxtime = 180, capture_output = True)

  times = []
  mEs = []
  for l in results['output'].split("\n"):
    if "mE" in l and "nprops" in l:
      mE = l.split("mE: ")[1].split()[0]
      ttime = l.split("ttime: ")[1].split()[0]
      mEs.append(mE)
      times.append(ttime)
  
  import matplotlib.pyplot as pl
  pl.plot(times, mEs)
  
  return times, mEs, results 
  
def kmeanspp(X, K):
  d2_nearest = 1e55*np.ones(X.shape[0])

  Cs = [X[0]]
  
  for i in range(1, K):
    d2_new = np.sum((X - Cs[-1])**2, axis = 1)
    new_is_nearer = d2_new < d2_nearest
    d2_nearest = d2_nearest*(1 - new_is_nearer) + d2_new*new_is_nearer
    
    d2_cum = np.cumsum(d2_nearest)
    rind = npr.rand()*d2_cum[-1]
    
    #Cs.append(X[i])
    #print np.where(rind < d2_cum)[0][0]
    Cs.append(X[np.where(rind < d2_cum)[0][0]])
  
  print np.mean(d2_nearest)
    
  #print Cs
    
  
  
  


