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





#def get_sklearn_mses(X, K, max_iter):
  #print "time taken : ", time.time() - t0  
  #Tracer()()

##def a_random_example():
#X = npr.randn(10000, 3)
#K = 500


#X = 0*rna.get_rna()[0:40000, 2::]


#noise = 1e-1*npr.randn(X.shape[0], X.shape[1])
#X += noise


#npr.seed(1000)
#X = npr.randn(40000, 8)

#npr.randn(20000, 8) #


def go(X, K):
  

  indices_init = np.arange(K, dtype = np.uint64)
  C_init = X[indices_init]

  
  withskl = False
  witheak = False
  withzen = True

  results = {}
  if withskl == True:
    results["skl"] = {}
    from sklearn.cluster import KMeans
    sklc = KMeans(n_clusters = K, init = "k-means++", max_iter = 100000000, tol = 1e-20, verbose = 0, n_init = 1)
    tsk0 = time.time()
    sklc.fit(X)
    tsk1 = time.time()
    sklacc = np.sum(np.min(np.sum((np.expand_dims(X, axis = 1) - np.expand_dims(sklc.cluster_centers_, axis = 0))**2, axis = 2), axis = 1)) / X.shape[0]
    results["skl"]["t"] = tsk1 - tsk0
    results["skl"]["mse"] = sklacc

  

  if witheak:
    results["eak"] = {}
    sys.path.append("/home/james/clustering/idiap/eakmeans/lib")
    import kmeans
    teak0 = time.time()
    eak = kmeans.get_clustering(X, K, verbose = 1, init = "kmeans++", n_threads = 4)
    teak1 = time.time()
    results["eak"]["t"] = teak1 - teak0
    results["eak"]['mse'] = eak["mse"]
  

  if withzen:
    results["zen"] = {}
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 15.0, max_time = 10.0, max_rounds = 10000, seed = npr.randint(1000), patient = True, nthreads = 4, init = "kmeans++-5", with_tests = False, capture_output = False, rooted = False)
    tzen0 = time.time()
    tangerine =  z.den(X, do_vdimap = False, do_refinement = True, rf_max_rounds = 10000000)
    tzen1 = time.time()
    results["zen"]["t"] = tzen0 - tzen1
    results["zen"]["out"] = pyzentas.get_processed_output(tangerine['output'])
    results["zen"]['mse'] = results["zen"]["out"]["mE"][-1]
    

  return results
  #if withskl:
    
    #print "skl (time) ", tsk1 - tsk0, "     (accuracy)   ", 
  
  #if witheak:
    #print "eak (time) ", teak1 - teak0, "   (accuracy)   ", bla['mse']

  ##for x in tangerine['output'].split("\n")[-10::]:
    ##print x
  
  #print "zen (time) ", tzen1 - tzen0, "   (accuracy)   ", zenout["mE"][-1]
  


K = 2000
seed = 107 #npr.randint(1000)
npr.seed(seed)
#X = npr.randn(K*16, 3)

X = rna.get_rna() #[0:K*16, 2::]
X += 0.001*npr.randn(X.shape[0], X.shape[1])


#X = npr.rand(K*320, 2)**1


for i in range(1):
  npr.seed()
  results = go(X, K)
  
  pl.ion()
  
  #label = "scikit-learn" if i == 0 else None
  #pl.plot(results["skl"]["t"], results["skl"]["mse"], marker = "o", color = 'k', markersize = 15, label = label)
  #label = "eakmeans" if i == 0 else None
  #pl.plot(results["eak"]["t"], results["eak"]["mse"], marker = "x", color = 'k', markersize = 15, label = label)

  
  pl.plot(results["zen"]["out"]["Tt"]/1000., results["zen"]["out"]["mE"], color = 'k', marker = "+", markersize = 2, linestyle = "none")
  label = "zentas" if i == 0 else None
  pl.plot(results["zen"]["out"]["Tt"][-1]/1000., results["zen"]["out"]["mE"][-1], color = 'k', marker = "+", markersize = 15, linestyle = "none", label = label)
  
  #print results["skl"]["mse"], results["eak"]["mse"], results["zen"]["out"]["mE"][-1]
  
  pl.show()

pl.legend()

def sklearn_elkan():
  """
  scikit learn inconsisitency due to numerical rounding. 
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



