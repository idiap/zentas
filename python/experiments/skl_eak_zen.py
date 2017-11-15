# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

"""
experiments comparing: 
 -   scikit-learn, 
 -   eakmeans and 
 -   zentas.  
"""

import matplotlib.pyplot as pl
import sys
import random
import numpy as np
import numpy.random as npr

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../../build/python")
import pyzentas
from IPython.core.debugger import Tracer 
import time
import rna
reload(rna)


def go(X, K, withskl, witheak, withzen):
  """
  X : data
  K : number of clusters
  withskl, witheak, withzen : bools indicating whether to run with em.
  """
  indices_init = np.arange(K, dtype = np.uint64)
  C_init = X[indices_init]

  results = {}
  if withskl == True:
    results["skl"] = {}
    from sklearn.cluster import KMeans
    # run until convergence, initialise with scikit-learn's special version of k-means++ (see zentas wiki entry for discussion). 
    sklc = KMeans(n_clusters = K, init = "k-means++", max_iter = 100000000, tol = 1e-20, verbose = 0, n_init = 1)
    tsk0 = time.time()
    sklc.fit(X)
    tsk1 = time.time()
    sklacc = np.sum(np.min(np.sum((np.expand_dims(X, axis = 1) - np.expand_dims(sklc.cluster_centers_, axis = 0))**2, axis = 2), axis = 1)) / X.shape[0]
    results["skl"]["t"] = tsk1 - tsk0
    results["skl"]["mse"] = sklacc


  if witheak:
    results["eak"] = {}
    sys.path.append(datapaths.datapath["eaklibdir"])
    import kmeans
    teak0 = time.time()
    eak = kmeans.get_clustering(X, K, verbose = 1, init = "kmeans++", n_threads = 4)
    teak1 = time.time()
    results["eak"]["t"] = teak1 - teak0
    results["eak"]['mse'] = eak["mse"]
  

  if withzen:
    results["zen"] = {}
    # run with zentas. pipeline here is (1) kmeans++ (2) clarans (3) lloyd.
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 10.0, max_time = 5.0, max_proposals = K**2, seed = npr.randint(1000), patient = True, nthreads = 4, init = "kmeans++-4", with_tests = False, capture_output = True, rooted = False)
    tzen0 = time.time()
    tangerine =  z.den(X, do_vdimap = True, do_refinement = True, rf_max_rounds = 10000000)
    tzen1 = time.time()
    results["zen"]["t"] = tzen0 - tzen1
    results["zen"]["out"] = pyzentas.get_processed_output(tangerine['output'])
    results["zen"]['mse'] = results["zen"]["out"]["mE"][-1]
    
  return results


def experiment1():


  #set the experiment settings here : 
  withskl = True
  witheak = False
  withzen = True
  K = 200
  ndata = 30*K
  seed = 107
  nruns = 10
  npr.seed(seed)

  dataset = "rna" #"random" # or "rna" or "mnist".

  if dataset == "mnist":
    import mnist
    reload(mnist)
    X = mnist.read_MNIST(dataset = "original", ndata = ndata)/1000.

  
  elif dataset == "rna":
    X = rna.get_rna()[0:ndata, 2::]
    X += 0.001*npr.randn(X.shape[0], X.shape[1])
    npr.seed()

  elif dataset == "random":
    X = npr.randn(ndata, 2)
    
  else:
    raise RuntimeError("unrecognised dataset string")

  for i in range(nruns):
    npr.seed()
    results = go(X, K, withskl, witheak, withzen)
    
    
    if withskl:
      label = "scikit-learn" if i == 0 else None
      pl.plot(results["skl"]["t"], results["skl"]["mse"], marker = "o", color = 'k', markersize = 15, label = label)
  
    if witheak:
      label = "eakmeans" if i == 0 else None
      pl.plot(results["eak"]["t"], results["eak"]["mse"], marker = "x", color = 'k', markersize = 15, label = label)
  
    
    if withzen:
      pl.plot(results["zen"]["out"]["Tt"]/1000., results["zen"]["out"]["mE"], color = 'k', marker = "+", markersize = 2, linestyle = "none")
      label = "zentas" if i == 0 else None
      pl.plot(results["zen"]["out"]["Tt"][-1]/1000., results["zen"]["out"]["mE"][-1], color = 'k', marker = "+", markersize = 15, linestyle = "none", label = label)

  pl.ion()
  pl.figure(1)
  pl.legend()
  pl.show()

def sklearn_elkan():
  """
  a scikit-learn inconsisitency due to numerical rounding
  """
  import numpy.random as npr
  from sklearn.cluster import KMeans
  seed = 51220
  npr.seed(seed)
  N = 1200
  K = 100
  X = npr.randn(N, 2)**7
  indices_init = np.arange(K, dtype = np.uint64)
  C_init = X[indices_init]
  for alg in ["elkan", "full"]:
    sklc = KMeans(n_clusters = K, init = C_init, max_iter = int(1e6), tol = 1e-20, verbose = 0, n_init = 1, algorithm = alg)
    sklc.fit(X)
    print "final E with algorithm ", alg, "\t : \t", np.sum(np.min(np.sum((np.expand_dims(X, axis = 1) - np.expand_dims(sklc.cluster_centers_, axis = 0))**2, axis = 2), axis = 1)) / X.shape[0]

def mnist_example():
  """
  Showing how do_vdimap can help. 
  """
  
  import mnist
  reload(mnist)
  
  ndata = int(1e3)
  X = mnist.read_MNIST(dataset = "original", ndata = ndata)/1000.
  dimension = X[0].shape[-1]
  npr.seed(1000)
  
  
  z = pyzentas.pyzen(K = 100, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_itok = 10000, max_rounds = 20, seed = 1011, nthreads = 1, init = "kmeans++-5", with_tests = False, patient = False)
  do_vdimap = True
  tangerine =  z.den(X, do_vdimap, do_refinement = True)
