# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

import matplotlib.pyplot as pl
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

import time

#def get_sklearn_mses(X, K, max_iter):
  #print "time taken : ", time.time() - t0
  
  #Tracer()()


##def a_random_example():
#X = npr.randn(10000, 3)
#K = 500

import rna
reload(rna)
K = 100
X = rna.get_rna()[10000:20000, :]


indices_init = np.arange(K, dtype = np.uint64)
C_init = X[indices_init]

#my_score = np.sum((np.expand_dims(X, axis = 1) - np.expand_dims(sklc.cluster_centers_, axis = 0))**2, axis = 1).shape
#print np.sum(np.min(np.sum((np.expand_dims(X, axis = 1) - np.expand_dims(sklc.cluster_centers_, axis = 0))**2, axis = 2), axis = 1)) / X.shape[0]
#get_sklearn_mses(X, K, 100000)

from sklearn.cluster import KMeans
t0 = time.time()
sklc = KMeans(n_clusters = K, init = C_init, max_iter = 100000000, tol = 1e-9, verbose = 0, n_init = 1, algorithm = "full")#, algorithm = "elkan") #k-means++
sklc.fit(X)
t1 = time.time()


#"random"


z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0., max_time = 1000.0, max_rounds = 50000, seed = npr.randint(1000), patient = True, nthreads = 1, init = indices_init, with_tests = False, capture_output = True)
tangerine =  z.den(X, do_vdimap = False, do_refinement = True, rf_max_rounds = 1000000000, rf_alg = "exponion")

t2 = time.time()

sys.path.append("/home/james/clustering/idiap/eakmeans/lib")
import kmeans
bla = kmeans.get_clustering(X, K, verbose = 1, init = indices_init)

t3 = time.time()

print "skl : ", t1 - t0, "   :   ", -sklc.score(X)/X.shape[0]
print "zen : ", t2 - t1, "   :   ", pyzentas.get_processed_output(tangerine['output'])["mE"][-1]
print "eak : ", t3 - t2, "   :   ", "hmm"





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


def the_rna_example():
  """
  pass
  """
  

  import time
  import rna
  reload(rna)
  K = 300
  X = rna.get_rna()[0:6000, :]

  with_skl = True
  
  if with_skl:

    t0 = time.time()
    print "running skl..."
    sklc = KMeans(n_clusters = K, init = "k-means++", max_iter = 20, tol = 1e-9, algorithm = "full") #k-means++
    sklc.fit(X)
    print sklc.score(X)
    print "done."    
    t1 = time.time()
    print "sklearn took : " , t1 - t0
  
  
  #resultsprint "running zentas..."
  results = {}
  for init in ["kmeans++-5"]:      #, "uniform"]:
    results[init] = {}
    for max_itok in [0, 0.3, 0.6, 0.9, 1.2]:
      print init, max_itok
      z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_itok = max_itok, max_time = 1000, max_rounds = 1000, seed = 1011, patient = True, nthreads = 1, init = init, with_tests = False, capture_output = True) #kmeans++-5
      
      tangerine = z.den(X, do_vdimap = False, do_refinement = True, rf_max_rounds = 10000, rf_alg = "exponion")
      
      results[init][max_itok] = pyzentas.get_processed_output(tangerine['output'])
  
    #print "done."
  
  #t2 = time.time()
  
  pl.figure(3)
  pl.clf()
  pl.ion()
  
  for k in results.keys():
    for k2 in results[k].keys():
      pl.plot(results[k][k2]["Tt"], results[k][k2]["mE"], label = "%s %s"%(k, k2), linestyle = "none", marker = ".")
  pl.legend()
  
  #print "zentas took : " , t2 - t1
  
  
  
  
  
