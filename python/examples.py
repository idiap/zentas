# Copyright (c) 2016 Idiap Research Institute
# Written by James Newling <jnewling@idiap.ch>

"""
Some examples using pyzentas.
"""

import sys
import random
import numpy as np
import numpy.random as npr

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../build/python")
import pyzentas

from IPython.core.debugger import Tracer 


def easy_kmeans_example():
  output = pyzentas.kmeans(X = npr.randn(10000, 3), K = 100, maxtime = 10, capture_output = False)
  return output 
  
def tests():
  """
  run some tests to confirm the algorithms are working correctly. 
  """
  
  #Test 1
  seed = 107
  npr.seed(seed)
  data = 10*npr.rand(25, 2)
  K = 7
  z = pyzentas.pyzen(K = K, metric = 'l2', algorithm = "clarans", level = 3, energy = 'quadratic', exponent_coeff = 0,  max_time = 2, max_rounds = 30, seed = 1011, nthreads = 1, with_tests = True, patient = False, init = "uniform", do_balance_labels = False, rooted = False)
  tangerine =  z.den(data, False, False)

  #Test 2
  z = pyzentas.pyzen(init = np.arange(K), K = K, algorithm = "voronoi", level = 0,  metric = 'l2', energy = 'identity', exponent_coeff = 0,  max_time = 1, max_rounds = 40, seed = 1011, nthreads = 1, with_tests = True, patient = False)
  tangerine =  z.den(data, True)


def dense_data_example():
  """
  cluster dense data points.
  """

  import random 
  npr.seed(1011)
  random.seed(1011)
  K = int(1e2)
  ndata = int(7e3)
  dimension = 5
  data = np.array(1 + npr.randn(ndata, dimension), dtype = np.float64)
  seed = 1011
  z = pyzentas.pyzen(init = "kmeans++~1", K = K, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_rounds = 20, max_time = 1., max_itok = 3.0, seed = seed, nthreads = 1, patient = False, with_tests = False, algorithm = "clarans", level = 3)
  do_vdimap = False
  do_refinement = True
  refinement_algorithm = "yinyang"
  rf_max_rounds = 3000;
  rf_max_time = 10000.;
  tangerine =  z.den(data, do_vdimap, do_refinement, refinement_algorithm, rf_max_rounds, rf_max_time)

  run_eakmeans = False
  if run_eakmeans:
    sys.path.append(datapaths.datapaths["eaklibdir"])
    import kmeans
    indices = random.sample(xrange(ndata), K)
    indices = np.array(indices, dtype = np.uint64)
    indices.sort()
    X = kmeans.get_clustering(X = data, n_clusters = K, init = indices, verbose = 1, n_threads = 1, algorithm = "yin-sn")

def sparse_data_example():
  """ 
  cluster sparse data points.
  """
  ndata = 20000
  sizes = np.array(npr.randint(5,9, size = (ndata,)), dtype = np.uint64)
  indices_s = []
  for i in range(ndata):
    new_indices = random.sample(xrange(10), sizes[i])
    new_indices.sort()
    indices_s.extend(new_indices)
    
  indices_s = np.array(indices_s, dtype = np.uint64)
  data = npr.randn(sizes.sum())
  z = pyzentas.pyzen(K = 100, max_rounds = 10000, max_itok = 15.001, seed = npr.randint(1011), with_tests = False)
  do_refinement = True
  refinement_algorithm = "yinyang"
  rf_max_rounds = 3000;
  rf_max_time = 10000.;
  
  X = z.spa(sizes, indices_s, data, do_refinement, refinement_algorithm, rf_max_rounds, rf_max_time)

  


def generated_sequence_example():
  """
  generate random sequences of chars/ints and cluster using levenshtein
  """
  ndata = 2000
  
  sizes = np.array(npr.randint(10, 30, size = ndata), dtype = np.uint64)
  
  #the values of the sequences
  data = []
  usechars = False
  for i in range(sizes.sum()):
    if npr.rand() < 0.25:
      data.append('A' if usechars else 0)
    elif npr.rand() < 0.5:
      data.append('C' if usechars else 1)      
    elif npr.rand() < 0.75:
      data.append('G' if usechars else 2)
    else:
      data.append('T' if usechars else 3)
  
  data = np.array(data, dtype = 'c' if usechars else np.int32)
  
  
  #The cost of mutating chars/ints (a 4x4 matrix)
  cost_switch = np.array(
  [[ 0.,   10,   8, 11],
  [ 10,     0.,   10, 9],
  [ 8,   10,   0,   10],
  [ 11, 9, 10,   0. ]], dtype = np.float64)
  
  
  #The cost of inserting or deleting a char/int
  cost_indel = np.array([10, 11, 10, 9], dtype = np.float64)
  
  z = pyzentas.pyzen(K = 400, metric = 'levenshtein', max_proposals = 100000, max_rounds = 2, energy = 'quadratic', seed = npr.randint(1000), nthreads = 1, init = "kmeans++-2")
  tangerine =  z.seq(sizes = sizes, values = data, cost_indel = cost_indel, cost_switch = cost_switch) 


def from_file_example():
  """
  read lists of words from files 
  and cluster using Normalised Levenshtein.
  """
  
  root = "../data/"  
  filenames_list = [root + "words1.txt", root + "words2.txt"]
  outfilename = root + "output1.txt"  
  #the costs of indels and switches
  costfilename = root + "costs.txt"
  
  z = pyzentas.pyzen(K = 5, metric = 'normalised levenshtein', max_proposals = 1000, energy = 'cubic', seed = npr.randint(1000))
  tangerine =  z.txt_seq(filenames_list, outfilename + "N", costfilename)


import matplotlib.pyplot as pl  
def capture_example():
  """
  extracting statistics from runs, with plots
  """
  
  import matplotlib.pyplot as pl  
  K = int(1e3)
  ndata = int(1e4)
  dimension = 3
  #so that distance to nearest ~ 1
  centers = K**(1./dimension)*npr.rand(K, dimension)
  indices = npr.randint(K, size = (ndata,))
  data = centers[indices] + 0.5*npr.randn(ndata,dimension)
  
  pl.clf()
  pl.ion()
  
  for max_itok in [x/2. for x in range(10)]:
    for init in ["kmeans++-5"]: #["kmeans++-10", "uniform"]:# "kmeans++", "afk-mc2-10", "afk-mc2-100", ]:
    
      print "with ", init, max_itok, "."
      z = pyzentas.pyzen(K = K, capture_output = True, init = init, max_itok  = max_itok, seed = 1011, level = 3)
      tangerine = z.den(data, do_vdimap = False, do_refinement = True)        
      oo = pyzentas.get_processed_output(tangerine['output'])
      pl.plot(oo['tT'], oo['mE'],  label = init + "[" + str(max_itok) + "]", linestyle = 'none', marker = '.')
    
  #pl.xscale('log')
  pl.xlabel("time [ms]")
  pl.ylabel("mean energy")
  #pl.legend()
  

  
  

def afk_mc2_subopt_mode():
  """
  We present an example where afk-mc^2 requires extremely
  long chain lengths to match k-means++. There are K centers,
  K-1 are uniformly distributed in [0,1] x [0,1], and the
  K'th center is at (10,10). Data points are drawn as follows:
  (1) select one of the K centers at random, 
  (2) add N(0,sigma) to it, where sigma is small (see below).
  
  Essentially, afk-mc^2 reduces to k-mc^2 
  (see Bachem 2016 for algorithms).
  """
  
  #################
  # Generate data #
  #################
  K = 200
  ndata = K*50
  dimension = 2
  centers = npr.rand(K,2)
  centers[-1] = [10, 10]

  sigma = 1e-5  
  labels = npr.randint(K, size = (ndata,))
  data = centers[labels] + sigma*npr.randn(ndata,2)
  data/=sigma
  
  
  ################
  # Run kmeans++ #
  ################
  E_kmeanspp = []
  print "kmeans++ energies: ", 
  for i in range(5):
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', seed = npr.randint(1000), max_rounds = 0, init = "kmeans++-1", capture_output = True)
    tangerine =  z.den(data)
    E_kmeanspp.append(tangerine['output'].split("R=0")[1].split("mE=")[1].split()[0])
    print E_kmeanspp[-1], " ",
  
  print "\nwith afk-mc2:"

  ###############
  # Run afk-mc2 #
  ###############
  E_afkmc2 = []
  l_afkmc2 = []
  for i in range(30):
    chain_length = i**2
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', seed = npr.randint(1000), max_rounds = 0, init = "afk-mc2-%d"%(chain_length,), capture_output = True)
    tangerine =  z.den(data)
    E_afkmc2.append(tangerine['output'].split("R=0")[1].split("mE=")[1].split()[0])
    l_afkmc2.append(chain_length)
    print "chain length: ", l_afkmc2[-1], "\tenergy: ", E_afkmc2[-1]
    
  return {"E_kmeanspp":E_kmeanspp, "E_afkmc2":E_afkmc2, "l_afkmc2":l_afkmc2}
  
  

