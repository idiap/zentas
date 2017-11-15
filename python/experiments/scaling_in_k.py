# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

"""
experiment determining scaling in k and n 
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


firstrun = True
plotresults = True

#obtain results
if firstrun:
  dim = 3
  N = 20000
  X = npr.randn(N, dim)
  results = {}
  
  for K in [int(x) for x in 2**np.linspace(np.log2(40), np.log2(2500), 40)]:
  
    results[K] = {}
    max_itoks = [0, 1, 2, 4]
    for max_itok in max_itoks:
      results[K][max_itok] = {}
      z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = max_itok, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "kmeans++-4", with_tests = False, capture_output = True, rooted = False)
    
      tzen0 = time.time()
      tangerine =  z.den(X, do_vdimap = True, do_refinement = True, rf_max_rounds = 10000000)
      tzen1 = time.time()    
      results[K][max_itok]["t"] = tzen0 - tzen1
      results[K][max_itok]["out"] = pyzentas.get_processed_output(tangerine['output'])
      results[K][max_itok]['mse'] = results[K][max_itok]["out"]["mE"][-1]
    
    print K, results[K][max_itoks[0]]['mse']/results[K][max_itoks[1]]['mse'], results[K][max_itoks[0]]['mse']/results[K][max_itoks[2]]['mse'], results[K][max_itoks[0]]['mse']/results[K][max_itoks[3]]['mse']

  Ks = np.array(results.keys())
  Ks.sort()
  itoks = np.array(results[Ks[0]].keys())
  itoks.sort()
  
  
if plotresults:  
  pl.figure(figsize = (8,4))
  pl.ion()
  for itok in itoks[1::]:
    vals = []
    for K in Ks:
      vals.append(results[K][itok]['mse'] / results[K][itoks[0]]['mse'])
    pl.plot(Ks, vals, label = "itok = %d"%(itok,), linestyle = ":", marker = "+")
  
  pl.legend(loc = "lower left")
  pl.xlabel("$K$")
  pl.ylabel("(mse with clarans) / (mse without)")
  pl.subplots_adjust(bottom = 0.2, top = 0.9)
  pl.savefig(datapaths.datapaths["kscalingfigpath"])

