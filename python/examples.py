# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

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


def dense_data_example():
  """
  cluster dense data ndata points in dimension `dimension'.
  """
  
  import numpy as np
  import mnist
  reload(mnist)
  
  ndata = int(1e4)
  X = mnist.read_MNIST(dataset = "original", ndata = ndata)
  dimension = X[0].shape[-1]
  npr.seed(1011)
  #data = np.array(npr.randn(ndata, dimension), dtype = np.float32)
  #data[:, 5:-5] *= 0.01
  
  z = pyzentas.pyzen(K = 1e3, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_time = 10000, max_rounds = 4, seed = 1011, patient = True, nthreads = 1)
  
  do_vdimap = True
  tangerine =  z.den(X, do_vdimap)

def sparse_data_example():

  ndata = 5
  sizes = np.array([2,3,2,1,2], dtype = np.uint64)
  indices_s = np.array([1,2,
  1,2,3,
  8,9,
  8,
  8,9000], dtype = np.uint64)
  data = np.ones(sizes.sum())

  z = pyzentas.pyzen(K = 2, max_proposals = 1e3, seed = npr.randint(1011))
  return z.spa(sizes, indices_s, data)

  


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


def capture_example():
  """
  extracting statistics from runs, with plots
  """

  import matplotlib.pyplot as pl  
  K = 1e2
  ndata = 1e4
  centers = np.sqrt(K)*npr.rand(K, 2)
  indices = npr.randint(K, size = (ndata,))
  data = centers[indices] + 0.5*npr.randn(ndata,2)

  pl.clf()
  pl.ion()
  
  for init in ["kmeans++", "kmeans++-10", "afk-mc2-10", "afk-mc2-100", "uniform"]:
    z = pyzentas.pyzen(K = K, capture_output = True, max_time = 0.2, init = init)
    tangerine = z.den(data)
    energies = [float(x.split()[0]) for x in tangerine["output"].split("mE=")[2::]]
    times = [float(x.split()[0]) for x in tangerine["output"].split("Tt=")[2::]]
    pl.plot(times, energies,  label = init, linestyle = ':', marker = '.')
    
  pl.xscale('log')
  pl.xlabel("time [ms]")
  pl.ylabel("mean energy")
  pl.legend()
  
  

def afk_mc2_failure_mode():
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
  
  

