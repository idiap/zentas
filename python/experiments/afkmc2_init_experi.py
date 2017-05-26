import hardpaths
reload(hardpaths)
import sys
sys.path.append(hardpaths.zentas_lib_dir)
sys.path.append("../../build/python")
import pyzentas
import kmedoids
import kmeans
import numpy as np
import numpy.random as npr
import time
import random
import copy
import os
import cPickle
import kmc2 
import copy
import hardpaths
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib 
from IPython.core.debugger import Tracer
pl.ion()


def load_all_rna():
  data = {}
  for fn in ["cod-rna",  "cod-rna.r",  "cod-rna.t"]:
    filly = open("/idiap/user/jnewling/bachemetaldata/" + fn, "r")
    lines = filly.readlines()
    filly.close()
    alldata = []
    for l in lines:
      subl = []
      for fr in l.split()[1::]:
        subl.append(float(fr.split(":")[1]))
      alldata.append(subl)
    data[fn] = np.array(alldata)
    print fn

  all_data = np.vstack(data.values())
  return all_data


rna = load_all_rna()
K = 2000

#ndata = 100
#dimension = 8
#rna = 0.01*npr.randn(ndata,dimension)
#psi = 1

#for i in random.sample(xrange(ndata), psi):
  #rna[i] = 10000.*npr.randn(8,)
  #print i,
#print "."
    


ndata, dimension = rna.shape

#afk-mc2-25

pyzentas.pyzentas(X = rna, K = K, initialisation_method = "kmeans++", indices_init = None, energy = "quadratic", metric = "l2", maxrounds = 3, seed = npr.randint(1000), maxtime = 100, exponent_coeff = 0)

#pyzentas.pyzentas(X = rna, K = K, indices_init = random.sample(xrange(ndata), K), initialisation_method = "from_indices_init", energy = "quadratic", metric = "l2", maxrounds = 20, seed = npr.randint(1000))

#kmeans_out_kmeanspp = kmeans.get_clustering(X = rna, n_clusters = K, algorithm = "exp-sn", init = 'kmeans++', verbose = 2, capture_verbose = False, seed = None, n_threads = 3, max_iter = 1)
