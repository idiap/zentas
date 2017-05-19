# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

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
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.weight'] = 'light'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif', serif='cm10')

from IPython.core.debugger import Tracer
pl.ion()


def get_default_zentas_dict(X, K, alg, n_threads, seed, init_time_limit = None):

  switch_type =  alg.split('_')[-1][0]
  ndata, dimension = X.shape
  
  if switch_type == 't':
    if init_time_limit != None:
      raise RuntimeError("switch_type is t, init_time_limit should be None")
    maxtime = kmeanspp_runtime*float(alg.replace("-mc2", "").split('_')[-1][1::])
    maxrounds = 1000000
    max_proposals = 1000000
  
  elif switch_type == 's':
    maxtime = init_time_limit if init_time_limit != None else 10**9
    maxrounds = int(K*float(alg.replace("-mc2", "").split('_')[-1][1::]))
    max_proposals = 1000 + int(K**2 / 2.)
  
  elif switch_type == 'u':
    maxtime = init_time_limit if init_time_limit != None else 10**9
    maxrounds = 100*K
    max_proposals = K**2
    
  else:
    raise RuntimeError("Unrecognised switch_type")
    
  #case 1.1 : start with kmc2, then do clarans.
  indices_init = None
  if "-mc2" in alg:
    initcenters = kmc2.kmc2(X, K, chain_length = 100, afkmc2 = True)
    indices_init = []
    for ikik , initcent in enumerate(initcenters):
      indices_init.append(np.where(np.sum(initcent != X, axis = 1) == 0)[0][0])
    indices_init = np.array(indices_init, dtype = np.uint64)
    
  #case 1.2 : do clarans with random initialisation
  else:
    indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)

  indices_init.sort()
      
  #run clarans.
  #TODO : patient should be True
  return pyzentas.pyzentas(ndata = ndata, dimension = dimension, X = X, K = K, indices_init = indices_init, algorithm = "clarans", level = 3, max_proposals = max_proposals, capture_output = True, seed = seed, maxtime = maxtime, nthreads = n_threads, maxrounds = maxrounds, patient = False, metric = 'l2', rooted = False, energy = 'quadratic')
  
  
def get_mses(X, K, alg, n_threads, seed, init_time_limit = None, holdout_error = True, ntest = 5000, ntrain = 5000):
  """
  run alg + lloyd
  """

  cloutput = "NONO"
  
  ndata, dimension = X.shape
  random.seed(seed)
  npr.seed(seed)
 
  if holdout_error == False:
    X_train = X
    X_test = X
    
  else:
    if ndata < ntest + ntrain:
      raise RuntimeError("To run with holdout_error == True, user should provide data array with more than %d elements (as per get_mses call, ntest + ntrain)"%(ntest + ntrain,))  

    #data has been shuffled already
    X_test = X[0:ntest]
    train_indices = random.sample(xrange(ntrain), ntrain)
    X_train = X[ntest::][train_indices]

  #case 0 :    
  if alg in ["pp", "BF", "BFS"]:
    #will just be a direct call to k-means.

    if alg == "pp":
      init = "kmeans++"
   
    elif alg in ["BF", "BFS"]:
      init = alg

    if holdout_error == True:
      #a first call to k-means on X_train to get initialisation centers.
      #As BF initialisation does not return initialisation centers, we raise an error.  
      if init == "BF":
        #TODO : raise an error rather.
        indices_init = np.array(random.sample(xrange(ntrain), K), dtype = np.uint64)
        indices_init.sort()
 
      else:
        indices_init = kmeans.get_clustering(X = X_train, n_clusters = K, algorithm = 'syin-ns', init = init, seed = seed, verbose = 2, capture_verbose = True, n_threads = n_threads, max_iter = 0)['I']
 
      init = X_train[indices_init]

    #will just use string "kmeans++" etc.
    else:
      pass


  #get **indices_init** on X_train. 
  #If holdout_error == False, init is them just indices_init.
  #If holdout_error == True, init is X_train [indices_init]
  else:
    
    #case 1 : uses clarans for initialisation
    if 'cl_' in alg:
      zentas_dict = get_default_zentas_dict(X_train, K, alg, n_threads, seed, init_time_limit = init_time_limit)
      cloutput = zentas_dict['output']
      indices_init = zentas_dict['indices_final']

    
    #case 2 : medlloyd initialisation, case 3 : uniform 
    elif "medlloyd" == alg or alg == "un":
      indices_init = np.array(random.sample(xrange(ntrain), K), dtype = np.uint64)
      indices_init.sort()
      
      if "medlloyd" == alg:
        maxrounds = 1000000
        maxtime = 100000000.
        indices_init = pyzentas.pyzentas(ndata = ntrain, dimension = dimension, X = X_train, K = K, indices_init = indices_init, algorithm = "voronoi", level = 0, capture_output = True, seed = seed, maxtime = maxtime, nthreads = n_threads, maxrounds = maxrounds, metric = 'l2', rooted = False, energy = 'quadratic')['indice_final']
      
    #case 4 : just do mc2
    elif "mc2" == alg:
      initcenters = kmc2.kmc2(X_train, K, chain_length = 100, afkmc2 = True)
      indices_init = []
      for ikik , initcent in enumerate(initcenters):
        indices_init.append(np.where(np.sum(initcent != X_train, axis = 1) == 0)[0][0])
      indices_init = np.array(indices_init, dtype = np.uint64)
    
  
    else:
      raise RuntimeError("unrecognised alg:" + alg)
    
    indices_init.sort()
    if holdout_error == True:
      init = X_train[indices_init]
    else:
      init = indices_init
  
 
  outp = kmeans.get_clustering(X = X_test, n_clusters = K, algorithm = 'auto', init = init, verbose = 2, capture_verbose = True, seed = seed, n_threads = n_threads)

  start_mse = float(outp['output'].split("valmse")[-1].split("\n")[1].split()[-2])
  mse = outp['mse']
  
  return {'start_mse': start_mse, 'mse':mse, 'I':outp['I'], 'cloutput':cloutput}

  
def get_kmeans_results(X, K, algs = None, n_threads = 3, init_time_fraction_of_total = 0.3, holdout_error = True):
  """
  (1) run k means ++ and lloyd with K centers on X. Set TL to be NRUNS_KMEANSPP times this.
  (2) for alg in algs : run series of alg + lloyd while TL not exceeded.
  (3) return the results
  """
  
  

  ndata, dimension = X.shape

  
  ntest, ntrain = ndata, ndata
  if holdout_error == True:
    ntest = 20000
    ntrain = 10000
    if ntest + ntrain > ndata:
      raise RuntimeError("ntest + ntrain > ndata in get_kmeans_results")
    
  seed = npr.randint(100000)
  random.seed(seed)
  npr.seed(seed)

  
  if holdout_error == True:
    timesetsize = ntest + ntrain
    
  else:
    timesetsize = ndata

  #see how long a run with kmeans++ takes. 
  print "getting time of complete run with kmeans++ on %d data points ... "%(timesetsize,),
    
  tstart = time.time()
  
  out_kmeanspp = kmeans.get_clustering(X = X[0:timesetsize], n_clusters = K, algorithm = 'auto', init = "kmeans++", verbose = 2, capture_verbose = True, n_threads = n_threads)
  
  tstop = time.time()
  t_elapsed = tstop - tstart
  print t_elapsed
  kmeanspp_runtime = t_elapsed
  time_per_alg = kmeanspp_runtime*NRUNS_KMEANSPP

  #cl_x : x is 
  # (if txxx) time spent in clarans, fraction of time taken for 1 run of out_kmeanspp.
  # (if kxxx) number of swaps, 
  results = dict.fromkeys(algs)
  
  subkeys = ['times', 'mses', "start_mses"]
  for alg in algs:
    results[alg] = dict.fromkeys(subkeys)
    for sk in subkeys:
      results[alg][sk] = []    

  for alg in algs: 
    t_start = time.time()
    
    print "\n%s \t "%(alg,),
    while time.time() - t_start < time_per_alg:
      print ".",
      both_mses = get_mses(X, K, alg, n_threads, seed = npr.randint(100000), init_time_limit = time_per_alg*init_time_fraction_of_total
      , holdout_error = holdout_error, ntest = ntest, ntrain = ntrain)
      results[alg]["mses"].append(both_mses['mse'])
      results[alg]["start_mses"].append(both_mses['start_mse'])
      results[alg]["times"].append(time.time() - t_start)
    
    results[alg]["I"] = both_mses["I"]
    for sk in subkeys:
      results[alg][sk] = np.array(results[alg][sk])
    
    print "%.6f %.6f [%.6f %.6f]"%(results[alg]["start_mses"].mean(), results[alg]['mses'].mean(), results[alg]["start_mses"].min(), results[alg]['mses'].min())

  results['t_elapsed_1_kmeanspp'] = t_elapsed
  return results
  
            ###########################
            ### 2D grid experiments ###
            ###########################


#TODO : set to all 4
grid_algorithms = ['un','cl_s3.0','pp', 'medlloyd']
#grid_algorithms = ['cl_s3.0']

grid_lsigmas = np.linspace(-10, 0, 401)
#grid_lsigmas = [-6] #-4]

grid_alg_plot_parms = {
  'un':{'color':"#dab595", 'marker':'.', 'markersize' : 2}, 
  'pp':{'color':"#008b50",'marker':'.', 'markersize' : 2}, 
  'medlloyd':{'color':"#48464a",'marker':'.', 'markersize' : 2}, 
  'cl_s3.0':{'color':"#cc252a",'marker':'.', 'markersize' : 2}}
grid_results_dir = "/idiap/user/jnewling/zentasoutput/icml/gridkmeans"
grid_kmeans_results_pickle_fn = os.path.join(grid_results_dir, "scores.pkl")
grid_K = 400
grid_nperc = 100

def get_grid_data(lsigma, K = grid_K, nperc = grid_nperc):
  """
  return data which is on a grid of unit length 1/lsigma, sigma = 1.
  """  
  rK = int(np.sqrt(K))
  ndata = nperc*K
  X_base = np.empty(shape = (ndata, 2), dtype = np.float64)
  X_base[:,0] = np.repeat(np.arange(rK), nperc*rK)
  X_base[:,1] = np.tile(np.arange(rK), nperc*rK)  
  return X_base/(2.**lsigma) + npr.randn(ndata, 2)


def get_energy_time(output_string):
  """
  extract energy and time from output string.
  
  1190    mE: 0.6364063   itime: 39       ctime: 6641     utime: 85       rtime: 8        ttime: 6789     lg2 nc(c): 25.532       lg2 nc: 25.554   pc: 1.0795      nprops : 12205


  """
  E = []
  ctime = []
  ncalcs = []
  ttime = []
  nprops = []
  #lg2_nc = []
  
  for l in output_string.split("\n")[0:-1]:
    if "mE:" in l:
      E.append(float(l.split("E: ")[1].split("\t")[0]))
      ctime.append(int(l.split("ctime: ")[1].split("\t")[0]))
      ttime.append(int(l.split("ttime: ")[1].split("\t")[0]))
      ncalcs.append(float(l.split("lg2 nc: ")[1].split("\t")[0]))
      nprops.append(int(l.split("nprops : ")[1].split("\t")[0]))
    #lg2_nc.append(float(l.split("lg2 nc: ")[1].split("\t")[0]))
  
  return {'ttime' : np.array(ttime), 'ctime' : np.array(ctime), 'E': np.array(E), 'ncalcs' :np.array(ncalcs), 'nprops':np.array(nprops)} #, 'lg2_nc':lg2_nc}
  
def get_grid_runtime_clarans():
  X = get_grid_data(lsigma = -4, K = grid_K, nperc = grid_nperc)
  seed = npr.randint(1000000)
  K = grid_K
  ndata, dimension = X.shape

  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()
  
  zentas_dicts = {}
  for level in [0,1,2]:#,1,2]:#[0,1,2,3]:
    print level
    zentas_dicts[level] = pyzentas.pyzentas(ndata = ndata, dimension = dimension, X = X, K = K, indices_init = indices_init, algorithm = "clarans", level = level, max_proposals = 100000, capture_output = True, seed = seed, maxtime = 10000000, nthreads = 1, maxrounds = 3*grid_K, patient = False, metric = 'l2', rooted = False, energy = 'quadratic') #K**2/4.
    
  
  return zentas_dicts # get_energy_time(zentas_dict['output'])

def get_grid_speed_table(zentas_dicts):
  
  en_tim = {}
  
  for level in [0,1,2]:
    en_tim[level] = get_energy_time(zentas_dicts[level]['output'])
    print "level : ", level, "    ",
    print "ncalcs", en_tim[level]["ncalcs"][-2], " \t ",
    print "ttime", en_tim[level]["ttime"][-2], " \t ",    
    print "ctime", en_tim[level]["ctime"][-2], " \t ",   
    print "nprops", en_tim[level]["nprops"].sum(), " \t ",   
    print "nswaps", en_tim[level]["nprops"].size
    
 
    
  return en_tim 
  
  

if False:
  zentas_dict = pyzentas.pyzentas(X = X, K = 400, indices_init = range(400), algorithm = "clarans", level = 3, max_proposals = 400*400, capture_output = True, seed = 1011, maxtime = 2000, nthreads = 3, maxrounds = 1000000, patient = False, metric = 'l2', rooted = False, energy = 'quadratic')
  
def plot_runtime_clarans(zentas_dict, fnsave = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/eval_vs_impl.pdf"):

  pl.figure(figsize = (10,1.2), num = 1011)
  pl.clf()
  
  pl.subplots_adjust(top = 0.85, bottom = 0.4, left = 0.32)
  energies_times = {}

  en_tim = get_energy_time(zentas_dict['output'])
  npoints = en_tim['nprops'][1::].size + 100
  pl.plot(en_tim['nprops'][1::], marker = '.', linestyle = 'none', color = 'k', markersize = 1)
  
  #pl.ylabel(, horizontalalignment = 'center')
  ax = pl.gca()
  pl.text(x = -0.04, y = -0.3, s ="evaluations", rotation = 90, transform=ax.transAxes, horizontalalignment = 'right', verticalalignment = 'bottom' )
  pl.xlabel("accepted swaps (implementations)", labelpad = 0)
  pl.yscale('log', basey = 2)
  lymax = np.log2(400.**2)
  ytocks = np.array([0, 10, lymax])


  pl.yticks((2.**ytocks).tolist(), ["$2^{%d}$"%(e,) for e in ytocks[0:-1]] + ["$N_r$"])

  pl.ylim(ymin = 2**(-0.8), ymax = 2**lymax*4)
  pl.plot([0, npoints], 2*[2**lymax], linestyle = '-', color = 'r')#(0.5, 0.5, 0.5))
  
  
  #xtocks = np.array([0, 200, 400, 600, 800, 1000, 1200])
  #pl.xticks((xtocks).tolist(), ["$%d$"%(e,) for e in xtocks[0:-1]] + ["$N_s$"])
  pl.xlim(xmin = 0, xmax = npoints)
  #pl.plot([1200, 1200], [2**-5, 2**20], linestyle = '-', color = (0.5, 0.5, 0.5))
  
  #Tracer()()
  if fnsave:
    pl.savefig(fnsave)
    import commands 
    commands.getstatusoutput("pdfcrop %s %s"%(fnsave, fnsave))
    
  


  
def get_grid_kmeans_results():
  """
  100 centers (0,0) ... (0,rK) ... (rK,0) ... (rK,rK). points 100 per cluster, center + N(0,I_2)
  """

  scores = {}
  for lsigma in grid_lsigmas:
    scores[lsigma] = dict.fromkeys(grid_algorithms)

  for lsigma in grid_lsigmas:
    for alg in grid_algorithms:
      scores[lsigma][alg] = dict.fromkeys(['mse', 'start_mse'])
  
  for lsigma in grid_lsigmas:
    X = get_grid_data(lsigma) 

    print "lsigma : ", lsigma, " \t "
    seed = npr.randint(100000)
    for alg in grid_algorithms:
      both_mses = get_mses(X, grid_K, alg, n_threads = 3, seed = seed, holdout_error = False)
      mse = both_mses['mse']*2.**(2.*lsigma)
      start_mse = both_mses['start_mse']*2.**(2.*lsigma)
      
      scores[lsigma][alg]['mse'] = mse
      scores[lsigma][alg]['start_mse'] = start_mse
      scores[lsigma][alg]['init'] = X[both_mses['I']]
      scores[lsigma][alg]['cloutput'] = both_mses['cloutput']
      
      print alg, ":", scores[lsigma][alg]['start_mse'], '->',  scores[lsigma][alg]['mse'], " \t ",

    print " "

  return scores


def write_grid_kmeans_results():
  print "Obtaining grid kmeans results ..."
  scores = get_grid_kmeans_results() 
  print "...results obtained..."
  filly = open(grid_kmeans_results_pickle_fn, "w")
  #cloutput = scores.pop('cloutput', None)
  cPickle.dump(scores, filly)
  filly.close()
  
  
def get_label(alg):
  label = r"\texttt{%s}"%(alg.replace("_", "-").replace("-s","").replace(".0","").replace("un", "uni").replace("pp", "++"))
  if "3" in label or "cl" in label:
    label = r"\texttt{clarans}"
  return label


def plot_grid_data(fnsave = "/idiap/home/jnewling/colftex/gridplots.pdf"):
  pl.figure(figsize = (7.1,2.6), num = 1443)
  pl.clf()
  for ai, lsigma in enumerate([-6, -4, -2]):
    data = get_grid_data(lsigma, nperc = 10)*2.**lsigma
    pl.subplot(1,3,ai+1)
    pl.plot(data[:,0], data[:,1], marker = '.', color = 'k', linestyle = 'none', markersize = 1.2)
    pl.xlim([-2, 21])
    pl.ylim([-2, 21])
    pl.xticks([])
    pl.yticks([])
    if ai == 0:
      pl.yticks([0,19])
      #pl.gca().xaxis.tick_top()
      #pl.xticks([0,19])
      
    pl.xlabel("$\sigma = 2^{%d}$"%(lsigma,), fontsize = 17)
  
  pl.subplots_adjust(bottom = 0.25, top = 0.89)
  if fnsave:
    pl.savefig(fnsave)
    import commands 
    commands.getstatusoutput("pdfcrop %s %s"%(fnsave, fnsave))


def initial_grid_locations(fnsave = "/idiap/home/jnewling/colftex/gridinits.pdf"):
  
  import copy
  pl.figure(figsize = (7,6), num = 31209)
  pl.subplots_adjust(left = 0.25, top = 0.7, bottom = 0.25, right = 0.75)
  
  pl.clf()
  
  filly = open(grid_kmeans_results_pickle_fn, "r")
  scores = cPickle.load(filly)
  filly.close()



  lsigma = -2.0 #grid_lsigmas[16]

  #pl.suptitle("$\sigma = 2^{-4}$")
  
  lsigmas = [-6.0, -4.0, -2.0]
   
  #l_grid_lsigmas = grid_lsigmas.tolist()[13:21]
  for lsigmai, lsigma in enumerate(lsigmas):
    orders = [0,3,2,1]
    for algi, alg in enumerate([grid_algorithms[i] for i in orders]):
      pl.subplot(3, 4, 4*lsigmai+ orders[algi] + 1)  #len(l_grid_lsigmas) lsigmai*4
      data = scores[lsigma][alg]['init'].T*2**lsigma
      baseargs = copy.copy(grid_alg_plot_parms[alg])
      #baseargs['markersize'] = 1.
      
      pl.plot(data[0], data[1], linestyle = 'none',  **baseargs)  
  
      if orders[algi] == 0:
        pl.ylabel("$\sigma = 2^{%d}$"%(lsigma,))
      
      if lsigmai == 2:
        pl.xlabel(get_label(alg))
  
      pl.xlim([-2, 21])
      pl.ylim([-2, 21])
      
      
      pl.xticks([])
      pl.yticks([])


  import commands
  if fnsave:
    pl.savefig(fnsave)
    commands.getstatusoutput('pdfcrop %s %s'%(fnsave, fnsave))
  
  
def plot_grid_init_final(fnsave = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/gridexperi.pdf"):

  orders = [0,3,2,1]  
  # data 

  gs1 = gridspec.GridSpec(4,6)
  gs1.update(left=0.1, right=0.99, top = 0.9, bottom = 0.23, wspace=0.05, hspace = 0.05)
  pl.figure(figsize = (7, 7), num = 1443)
  pl.clf()
  for ai, lsigma in enumerate([-6, -4, -2]):
    data = get_grid_data(lsigma, nperc = 10)*2.**lsigma
    print ai
    ax1 = pl.subplot(gs1[ai:ai+1, 0:1])
    pl.plot(data[:,0], data[:,1], marker = '.', color = 'k', linestyle = 'none', markersize = 1.2)
    pl.xlim([-2, 21])
    pl.ylim([-2, 21])
    pl.xticks([])
    pl.yticks([])
    if ai == 2:
      pl.xticks([0,19])
    #if ai == 2:
      #pl.xlabel('data')
    pl.ylabel("$\sigma = 2^{%d}$"%(lsigma,), fontsize = 15)

  #initialisations
  import copy
  filly = open(grid_kmeans_results_pickle_fn, "r")
  scores = cPickle.load(filly)
  filly.close()
  lsigmas = [-6.0, -4.0, -2.0]
  for lsigmai, lsigma in enumerate(lsigmas):
    for algi, alg in enumerate([grid_algorithms[i] for i in orders]):
      ax1 = pl.subplot(gs1[lsigmai:lsigmai+1, algi+1:algi+2])
      data = scores[lsigma][alg]['init'].T*2**lsigma
      baseargs = copy.copy(grid_alg_plot_parms[alg])
      pl.plot(data[0], data[1], linestyle = 'none',  **baseargs)  
      
      if lsigmai == 2:
        pl.xlabel(get_label(alg), size = 'large')
  
      pl.xlim([-2, 21])
      pl.ylim([-2, 21])
      
      
      pl.xticks([])
      pl.yticks([])
      
  import commands
  if fnsave:

    pl.savefig(fnsave)
    commands.getstatusoutput('pdfcrop %s %s'%(fnsave, fnsave))
  


  
  

def plot_grid_kmeans_results(fnsave = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/gridres.pdf"):
  
  filly = open(grid_kmeans_results_pickle_fn, "r")
  scores = cPickle.load(filly)
  filly.close()


  #pl.figure(figsize = (5,2.5), num = 312)
  pl.figure(figsize = (7,1.9), num = 312)
  pl.clf()
  lines = []
  for algi, alg in enumerate(grid_algorithms):
    algseries = {'mse':[], 'start_mse':[]}
    for lsigma in grid_lsigmas:
      for mse_type in ['mse', 'start_mse']:
        algseries[mse_type].append(0.5*scores[lsigma][alg][mse_type]/2.**(2*lsigma))
 
    
    baseargs = copy.copy(grid_alg_plot_parms[alg])
    baseargs['markersize'] = 1.8

    pl.subplot(1,2,1)
    #pl.title("initialisation", fontsize = 'small')
    pl.plot(2.**grid_lsigmas, algseries['start_mse'],  linestyle = "None",  alpha = 0.6, **baseargs) #[0]
    
    pl.subplot(1,2,2)
    #pl.title(r"after  $\texttt{lloyd}$", fontsize = 'small')
    lines.append(pl.plot(2.**grid_lsigmas, algseries['mse'], label = get_label(alg),  linestyle = "none", alpha = 0.6, **baseargs)[0])
    

  #pl.legend(handles = [lines[i] for i in [0,3,2,1]], numpoints = 1, frameon = False, labelspacing = 0.1)  
  pl.legend(handles = [lines[i] for i in [3,0,2,1]], numpoints = 1, frameon = False, labelspacing = 0.00, handletextpad = -0.45, fontsize = 'small', loc = 'upper right')# loc = (0.35 , 0.48))
  pl.subplots_adjust(bottom = 0.25, top = 0.75, left = 0.12, wspace = 0.32, hspace = 0.32)
  
  for spi in [1,2]:
    pl.subplot(1,2,spi)
    pl.xscale('log', basex = 2)  
    pl.yscale('log', basey = 2) 
    
    pl.xlim(xmax = 2.**-1.0)
    pl.ylim(ymax = 2**19.5)
    ytocks = np.array([-4, 0, 4, 8, 12,16])
    pl.yticks((2.**ytocks).tolist(), ["$2^{%d}$"%(e,) for e in ytocks])
    
  pl.subplot(1,2,1)
  #pl.xticks([])
  pl.xlabel("$\sigma$", labelpad = -3)
  ax = pl.gca()
  #pl.text(x = 0, y = 0, s = r"initial $MSE/\sigma^2$", rotation = 90, transform=ax.transAxes, verticalalignment = 'left')
  pl.ylabel(r"init $MSE/\sigma^2$")#, horizontalalignment = 'center') 
  
  pl.subplot(1,2,2)
  pl.xlabel("$\sigma$", labelpad = -3)
  pl.ylabel("final $MSE/\sigma^2$")
  #
    
  


  #pl.subplot(2,1,2)
  #pl.yticks([])

  import commands
  if fnsave:
    pl.savefig(fnsave)
    commands.getstatusoutput('pdfcrop %s %s'%(fnsave, fnsave))


            ###########################
            ### joensuu experiments ###
            ###########################    


import load_joensuu_data
reload(load_joensuu_data)

NRUNS_KMEANSPP = 80
ALGS =  ['pp', "cl_u", "un", "BF", "mc2"]#,  'cl_s1.0', 'cl_s3.0', 'cl_s1.0-mc2', 'cl_s3.0-mc2'] ["un", "BF"]#



#= []

#bloop = [

# from joensuu website. This is ALL the data.
joensuu_datasets = [
'birch1',
'birch2',
'birch3',
's1',
's2',
's3',
's4',
'a1',
'a2',
'a3',
'dim032',
'dim064',
'dim1024',
'KDDCUP04Bio', 
'ConfLongDemo_JSI_164860',
'MopsiLocationsUntil2012-Finland',
'europediff', 
'housec8', 
'mnist',
'yeast'
]
joensuu_datasets.sort(key = str.lower)

joensuu_roman = {}
for k, ds in enumerate(joensuu_datasets):
  joensuu_roman[ds] = '%d'%(k,)


bachem_datasets = [
"susy",
"song",
"rna"
]

joensuu_datasets.sort(key = str.lower)

holdout_datasets = ["susy", "song"]  #'KDDCUP04Bio'] #'MopsiLocationsUntil2012-Finland']





# The K-s to run with
trueks = {}
trueks['ConfLongDemo_JSI_164860'] =  22
trueks['birch1'] = 200
trueks['birch2'] = 200
trueks['birch3'] = 200
trueks['s1']  = 30
trueks['s2']  = 30
trueks['s3']  = 30
trueks['s4']  = 30
trueks['a1']  = 40
trueks['a2']  = 70
trueks['a3']  = 100
trueks['dim032']  = 32
trueks['dim064']  = 32
trueks['dim1024'] = 32
trueks['housec8'] = 400 #-1
trueks['MopsiLocationsUntil2012-Finland'] =  100 #100#-1
trueks['mnist'] = 300 #10
trueks['europediff'] = 1000 #-1
trueks['KDDCUP04Bio'] = 200 #2000
trueks['yeast'] =  40

trueks["rna"] = 200
trueks["song"] = 500
trueks["susy"] = 500


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
    
def load_bachem_data(dataset):
  data = None
  if dataset == "rna":
    filly = open("/idiap/user/jnewling/bachemetaldata/rna.txt", "r")
    lines = filly.readlines()
    filly.close()
    alldata = []
    for l in lines:
      subl = []
      for fr in l.split()[1::]:
        subl.append(float(fr.split(":")[1]))
      alldata.append(subl)
    data = np.array(alldata)
  
  elif dataset == "song":
    data = np.load("/idiap/user/jnewling/bachemetaldata/YearPredictionMSD.npy")[:,1::]

  elif dataset == "susy":
    data = np.load("/idiap/user/jnewling/bachemetaldata/SUSY.npy")[:,1::]
  

    
  #ndata, dimension = X.shape
  #indices = random.sample(xrange(ndata), ndata)
  #X = X[indices][0:-250000, 1::]
  

  return data



def joensuu_bachem_experiment(dataset = "MopsiLocationsUntil2012-Finland", K = 15, writedata = True, holdout_error = True):
  """
  (1) Read data
  (2) Get results via initialisation_test
  (3) Write results
  """
  if dataset in joensuu_datasets:
    X = load_joensuu_data.load_data(dataset)
 
  elif dataset in bachem_datasets:
    X = load_bachem_data(dataset)
  
  else:
    raise RuntimeError("Unrecognised dataset")
  
  indices = random.sample(xrange(X.shape[0]), X.shape[0])
  X = X[indices] 
    
  scores = get_kmeans_results(X, K, algs = ALGS, holdout_error = dataset in holdout_datasets)
  
  if writedata : 
    resultsdir = "/idiap/user/jnewling/zentasoutput/nips"
    elapsed_dir = os.path.join(resultsdir, "elapsed_times_pp")
    pickles_dir = os.path.join(resultsdir, "pickles")
    
    if not os.path.exists(elapsed_dir):
      os.makedirs(elapsed_dir)
  
    if not os.path.exists(pickles_dir):
      os.makedirs(pickles_dir)
    
    filly = open(os.path.join(elapsed_dir, "%s.txt"%(dataset,)), "w")
    value = scores['t_elapsed_1_kmeanspp']
    filly.write("%.5f"%(value))
    filly.close()
    
    filly = open(os.path.join(pickles_dir, "%s.pkl"%(dataset,)), "w")
    cPickle.dump(scores, filly)
    filly.close()
    
    
  
def all_joensuu_bachem_experiments(writedata = True):  
  """
  All rank vs energy plots.
  """
  
  for ik, k in enumerate(joensuu_datasets):
    print "\n joensuu ---------------------", ik, "  :  ", k, "-----------------------"
    joensuu_bachem_experiment(k, trueks[k], writedata = writedata)
  
  for ik, k in enumerate(bachem_datasets):
    print "\n bachem ---------------------", ik, "  :  ", k, "-----------------------"
    joensuu_bachem_experiment(k, trueks[k], writedata = writedata)
  
  
def print_joensuu_dataset_table():
  """
  Write the statistics : dataset, K, N, dim, TL. 
  """
  
  def get_abridged(dataset):
    text = dataset
    if 'Mopsi' in text:
      text = 'Mopsi'
    elif 'ConfLong' in text:
      text = 'ConfLong'
    elif 'KDDCUP04' in text:
      text = 'KDDCUP04'
    elif 'europe' in text:
      text = 'europe'
    return text
    
    
  resultsdir = "/idiap/user/jnewling/zentasoutput/icml"
  elapsed_dir = os.path.join(resultsdir, "elapsed_times_pp")
  datasets = joensuu_datasets 
  datasets.sort(key = str.lower)
  taboo = r"""
\begin{table}
\centering
\begin{tabular}{cccccc}
\hline
dataset & $\#$ & N & dim & K & TL $[s]$\\
\hline
"""
  for ik, k in enumerate(datasets):
    X = load_joensuu_data.load_data(k)
    ndata, dimension = X.shape
    ndata += (dimension == 2)
    K = trueks[k]
    
    filly2 = open(os.path.join(elapsed_dir, "%s.txt"%(k,)), "r")
    line = filly2.readline()
    filly2.close()
    elapsed_time = float(line.strip())
    allocated_time = elapsed_time*NRUNS_KMEANSPP
      
    taboo +=  r""" %s & %d &  %d &  %d &  %d &   %.2f \\
"""%(get_abridged(k), ik, ndata, dimension, K, allocated_time)
  
  taboo += r"""
\hline
\end{tabular}
\caption{datasets}
\label{tab:datasets}
\end{table}
"""
  print taboo




def load_results(datasets):
  resultsdir = "/idiap/user/jnewling/zentasoutput/icml2"
  elapsed_dir = os.path.join(resultsdir, "elapsed_times_pp")
  pickles_dir = os.path.join(resultsdir, "pickles")
  datasets.sort(key = str.lower)
  results = {}  
  for ikop, dataset in enumerate(datasets):
    pklfn = os.path.join(pickles_dir, "%s.pkl"%(dataset,))
    print pklfn
    if os.path.exists(pklfn):
      filly = open(pklfn, "r")
      scores = cPickle.load(filly)
      filly.close()
      results[dataset] = scores
  return results


all_datasets = joensuu_datasets + bachem_datasets

def load_all_results():
  return  load_results(all_datasets)




colors = ["#2faced" , '#fd1f58', 'g']
def experiplot(sfn = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/bars.pdf"):
  results = load_all_results()


  pl.figure(figsize = (9,4.7))
  pl.clf()
  ds = 's3'

  ndatasets = len(results.keys())

  ylabs = ["initialisation MSE", "final MSE"]


  gs1 = gridspec.GridSpec(4,1)
  gs1.update(left=0.1, right=0.99, top = 0.9, bottom = 0.23, wspace=0.05, hspace = 0.05)


  for mse_key_i, mse_key in enumerate(['start_mses', 'mses']):
    
    if mse_key_i == 0:
      ax1 = pl.subplot(gs1[0:2, 0:1])
  
    else:
      ax1 = pl.subplot(gs1[2:4, 0:1])
  
    #pl.subplot(2,1,1 + mse_key_i)
    pl.ylabel(ylabs[mse_key_i], size = 'large')

    pl.xlim([-1.6, 23])
    
    for dsi, ds in enumerate(all_datasets):
      bar_width = 0.33
          
      valuesx = results[ds]["pp"]['start_mses']
      nfactor = valuesx.mean() 

    
    
    
      algs = ["pp","cl_u"]
      
      for algi, alg in enumerate(algs):
  

        values = results[ds][alg][mse_key].copy()
        x0 = dsi + bar_width*(algi) # + 2.5*mse_key_i)
        values /= nfactor
        
        mean = values.mean()
        std = values.std()
        maxval = values.max()
        minval = values.min()
        
        if mse_key_i == 0 and dsi == 0 and algi == 0:
          pl.text(-0.36, mean, '(2)', horizontalalignment = 'right', verticalalignment = 'center')
          pl.text(-0.36, minval, '(1)', horizontalalignment = 'right', verticalalignment = 'center')
          pl.text(-0.36, mean + std, '(3)', horizontalalignment = 'right', verticalalignment = 'center')
        
        #pl.plot([x0 - bar_width, x0], [mean, mean], color = colors[algi], linewidth = 2)
        
        if (dsi == 0 and mse_key_i == 1):
          label = get_label(alg)
          if "+" in label:
            label = r"$\texttt{km++}$"
        else:
          label = None
        
        pl.plot([x0 - bar_width, x0], [minval, minval], color = colors[algi], linewidth = 1, label = label)
        pl.fill_between([x0 - bar_width, x0], [0, 0], [minval, minval], edgecolor = colors[algi], facecolor = colors[algi],  alpha = 0.5)
        
        if algi == 0:
          alpha = 0.25
        else:
          alpha = 0.1
          
        pl.fill_between([x0 - bar_width, x0], [minval, minval], [mean, mean], edgecolor = colors[algi], facecolor = colors[algi],  alpha = alpha)
        
        pl.plot([x0 - bar_width/2., x0 - bar_width/2.], [mean + std, mean], color = colors[algi], linewidth = 1, linestyle = '-')       
        pl.plot([x0 - bar_width, x0], [mean, mean], color = colors[algi], linewidth = 1, linestyle = '-')
        pl.plot([x0 - bar_width, x0], [mean + std, mean + std], color = colors[algi], linewidth = 1, linestyle = '-')
      
    pl.xticks(range(23), [x + 1  for x in  range(23)])
    

      
    if mse_key_i == 1:
      pl.legend(loc = 'upper left', labelspacing = 0.06)
      pl.ylim([0.5, 1.16])
      #pl.yticks([0.5, 0.6, 0.7, 0.8, 0.9])
      
    else:
      pl.ylim([0.5, 1.16])
      pl.xticks([])
        

        #pl.plot([x0 - bar_width/2., x0 + bar_width/2.], [mean - std, mean - std], color = colors[algi], linewidth = 1, linestyle = '-')
        
        #pl.plot([x0, x0], [minval, minval], markerfacecolor = colors[algi], markeredgecolor = colors[algi], linewidth = 3, marker = "o", markersize = 2)
    
  pl.subplots_adjust(wspace = 0.06, hspace = 0.06)
  pl.xlabel('dataset', size = 'large')
  
  
  pl.savefig(sfn)
  import commands
  commands.getstatusoutput('pdfcrop %s %s'%(sfn, sfn))

    #pl.subplot(ndatasets, 2, 2*dsi + 1)
    #pl.boxplot([results[ds][alg]['start_mses'] for alg in algs])

    #pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #pl.xticks([])
    #yto = pl.yticks()[0]
    #pl.yticks([yto[1], yto[-2]])


    #pl.subplot(ndatasets, 2, 2*dsi + 2)
    
    #pl.text(0.95,0.9, ds, horizontalalignment = 'right', verticalalignment = 'top', transform=pl.gca().transAxes)
    #for algi, alg in enumerate(algs):
      #pl.text((algi + 0.92)/(len(algs) + 0.), 0.0, '%d'%(results[ds][alg]['mses'].size), transform=pl.gca().transAxes, verticalalignment = 'bottom', horizontalalignment = 'right', size = 'small')
    
    #pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #pl.xticks([])
    #yto = pl.yticks()[0]
    #pl.yticks([yto[1], yto[-2]])

    
    #pl.subplots_adjust(left = 0.3, wspace = 0.3, hspace = 0.4, bottom = 0.3)
  
    
  

def print_joensuu_result_tables():

  TABALGS =  ALGS #["pp", "un", "BF", "mc2",  'cl_s1.0', 'cl_s1.0-mc2', 'cl_s3.0']

  shortalgname = {
  "un":r'\begin{sideways} \texttt{uni}   \end{sideways}',
  "pp":r'\begin{sideways} ++  \end{sideways}',
  "BF":r'\begin{sideways} \texttt{bf}    \end{sideways}',
  "cl_s1.0":r'\begin{sideways}   \texttt{cl1}    \end{sideways}',
  "cl_s3.0":r'\begin{sideways}   \texttt{cl3}    \end{sideways}',
  "mc2":r'\begin{sideways}       \texttt{mc2}     \end{sideways}',
  "cl_s1.0-mc2":r" \begin{sideways} \hspace{-2.0pt} \texttt{mc2} \end{sideways}  \begin{sideways} \texttt{+cl1} \end{sideways}",
  "cl_s3.0-mc2":r'\texttt{mc2cl3}',
  "cl_u":r'\texttt{cla}'
  }
  
  results = load_joensuu_results()
  
  for table in ["min", "mean", "sigma", "count"]:
    
    
    print r"\begin{table}"
    print "\centering"
    print r"\begin{tabular}{cccccccc}" #ccccccccc
    print "\hline"
    print " ",
    for alg in TABALGS:
      print " &  %s "%(shortalgname[alg])  ,
    print r"\\"
    print "\hline",
  
    for ikop, dataset in enumerate(joensuu_datasets):
      print "\n", joensuu_roman[dataset], " & ",
      
    
      
      if table in ["min", "mean"]:
        
        if table == "min":
          foffer = lambda alg : results[dataset][alg]['mses'].min()
      
        if table == "mean":
          foffer = lambda alg : results[dataset][alg]['mses'].mean()
        
        nfactor = foffer('pp')
        minval = min([foffer(alg)/nfactor for alg in TABALGS])
        formatprint = lambda v :  r"\textbf{%.3g}"%(v,) if v == minval and v<1. else r'%.3g'%(v)
        
        for  alg in TABALGS:
          print formatprint(foffer(alg)/nfactor,), " &  " if alg != TABALGS[-1] else r"\\",
  
      elif table == "sigma":
        for  alg in TABALGS:
          print r'%.3g'%( np.abs(results[dataset]['pp']['mses'].mean() - results[dataset][alg]['mses'].mean()) / results[dataset][alg]['mses'].std(),), " & " if alg != TABALGS[-1] else r"\\",
                
      else:
        for  alg in TABALGS:
          print results[dataset][alg]['mses'].size, " & " if alg != TABALGS[-1] else r"\\",
  
    print ""
    print "\hline"
    print "\end{tabular}"
    print r"\caption{%s}"%(table,)
    print r"\label{tab::%s}"%(table,)
    print "\end{table}"
    print ""
  


            ##########################
            ### bachem experiments ###
            ##########################    


def bachem_experiment(X, K):

  ndata, dimension = X.shape  
  print X.shape

  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()
  
  
  
  if True:
    
    methods = ["kmeans++", "uniform"]
    
    results = dict.fromkeys(methods)
    for m in methods:
      results[m] = {"init":[], "final":[]}
      for i in range(10):
        bla = kmeans.get_clustering(X = X, algorithm = 'auto', n_clusters = 200, init = m, n_threads = 3, capture_verbose = True, seed = npr.randint(1000000))
        results[m]["init"].append(float(bla['output'].split("initial mse : ")[-1].split()[0])*ndata)
        results[m]["final"].append(bla['mse']*ndata)
        print m, "%.3e"%(results[m]["init"][-1]), "%.3e"%(results[m]["final"][-1])
  

  #bla = pyzentas.pyzentas(dimension = dimension, ndata = ndata, X = X.ravel(), K = K, metric = "l2", energy = "quadratic", algorithm = "clarans", level = 3, capture_output = False, indices_init = indices_init, max_proposals = 10000000, maxtime = 100000, maxrounds = 200000, patient = True, rooted = False)

  
  return results



def get_kmeans_energy_time(kmo):
  
  energies = []
  times = []
  
  good_zone = False
  for l in kmo.split("\n"):
    if "roundchanges" in l:
      good_zone = True
      continue
      
    
    if good_zone:
      frags = l.split()
      if len(frags) != 7:
        good_zone = False
        break

      times.append(float(frags[4]))
      energies.append(float(frags[5]))
  
  return {"E":np.array(energies), "T":np.array(times)/1000.}

if False:
  all_data = load_all_rna()

K = 2000
n_runs = 5
maxtime = 18
  
def big_experiment():
  
  results = []
  for runi in range(n_runs):

    print "run ", runi    
    
    #indices_init = kmeans.get_clustering(X = all_data, n_clusters = K, init = "kmeans++", seed = seed, verbose = 2, capture_verbose = True, n_threads = n_threads, max_iter = 0)['I']
    
    
    alg = 'ham' #'selk-ns'

    print "in kmeans ++..."    
    kmeans_out_kmeanspp = kmeans.get_clustering(X = all_data, n_clusters = K, algorithm = alg, init = 'kmeans++', verbose = 2, capture_verbose = True, seed = None, n_threads = 3)
    kmok = get_kmeans_energy_time(kmeans_out_kmeanspp['output'])

    print "in zentas..."
    zentas_dict = pyzentas.pyzentas(X = all_data, indices_init = kmeans_out_kmeanspp["I"], K = K, algorithm = "clarans", level = 3, max_proposals = K*K, capture_output = True, seed = 1011, maxtime = maxtime, nthreads = 3, maxrounds = 1000000, patient = False, metric = 'l2', rooted = False, energy = 'quadratic')      
    zentas_series = get_energy_time(zentas_dict['output'])

    print "in kmeans z..."    
    kmeans_out_zentas = kmeans.get_clustering(X = all_data, n_clusters = K, algorithm = alg, init = zentas_dict['indices_final'], verbose = 2, capture_verbose = True, seed = None, n_threads = 3)    
    kmoz = get_kmeans_energy_time(kmeans_out_zentas['output'])


    #print "in kmeans u..."    
    #kmeans_out_uniform = kmeans.get_clustering(X = all_data, n_clusters = K, algorithm = alg, init = 'uniform', verbose = 2, capture_verbose = True, seed = None, n_threads = 3)
    #kmou = get_kmeans_energy_time(kmeans_out_uniform['output'])

  
  
    results.append({'zentas_series':zentas_series, 'kmok':kmok, 'kmoz':kmoz}) #, 'kmou':kmou})
    
  return results
  

def plot_big_experiment(be_all, fn = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/precede.pdf"):
  

  pl.figure(figsize = (5.2,3))
  pl.subplots_adjust(bottom = 0.3, left = 0.3)

  mio = 0
  for bei, be in enumerate(be_all):
  
    label = None
    label2 = None
    if bei == 0:
      label = r"\texttt{km++}+\texttt{clarans}+\texttt{lloyd}"
      label2 = r"\texttt{km++}+\texttt{lloyd}"
    
    pl.plot((be['zentas_series']['ttime'] - be['zentas_series']['ttime'][0])/1000., be['zentas_series']['E'], color = colors[1], linestyle = '-', alpha = 1.0, label = label)    
    pl.plot(be['kmoz']['T'] - be['kmoz']['T'][0] + (be['zentas_series']['ttime'][-1] - be['zentas_series']['ttime'][0])/1000., be['kmoz']['E'], color = colors[1], label = None)

    #pl.plot((be['kmoz']['T'] + be['zentas_series']['ttime'][-1]/1000.), be['kmoz']['E'] - mio, color = 'b')
    pl.plot(be['kmok']['T'] - be['kmok']['T'][0], be['kmok']['E'], color = colors[0], label = label2)
    
    #pl.plot(be['kmou']['T'], be['kmou']['E'] - mio, color = 'g')
    
    print be['kmoz']['E'].min(),
    print be['kmok']['E'].min()
    #print be['kmou']['E'].min()
    
  pl.xlabel("time [s]", size = 'large')
  pl.ylabel("$MSE$")
  pl.xlim(xmin = -0.4)
  pl.legend(borderpad = 0.5, labelspacing = 0.09, handlelength = 0.07, borderaxespad = 0.2)
  pl.savefig(fn)
  import commands
  commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))

  #p1 = pl.Line2D([1,1], [1,1])
  #p2 = pl.Line2D([1,1], [1,1])
  #pl.legend([p2, p1], ["line 2", "line 1"])

    
  #pl.xscale('log', xscale = 2)
  #pl.yscale('log', yscale = 2)
  #pl.ylim(ymin = 1)
  
#def expkdd04():  
  #"""
  #done by joensuu
  #"""
  #print "loading data..."
  #X = np.loadtxt(fname = "/idiap/user/jnewling/bachemetaldata/kdd04/bio_train.dat")[:,3::]
  
    
  ##X = np.loadtxt(fname = "/idiap/user/jnewling/joensuudata/txts/KDDCUP04Bio.txt")
  #print "loaded."
  #bachem_experiment(X, K = 200)
  

#def exprna():

  
  #bachem_experiment(data, K = 200)
  #

#def expsusy():
  #X = np.load("/idiap/user/jnewling/bachemetaldata/SUSY.npy")[:,1::]
  
  ##problem with bachem's numbers.
  ##todo : k-means and uniform. Then, kmedoids on 500,000 datapoints VALIDATED on all data.
  
  #for i in range(100):
    #index = npr.randint(5*10**6)
    #d2s = np.sum((X - X[index])*(X - X[index]), axis = 1)
    #d2s.sort()
    #print d2s[int(5*10**6 / 2000)]
  
  #indices = random.sample(xrange(5*10**6), 5*10**6)
  #X = X[indices]
  #X.reshap
  #bachem_experiment(X, K = 2000)

#def expsong():

  #X = np.load("/idiap/user/jnewling/bachemetaldata/YearPredictionMSD.npy")#[:,1::]
  #ndata, dimension = X.shape
  #indices = random.sample(xrange(ndata), ndata)
  #X = X[indices][0:-250000, 1::]
    
  #print "entering experiment"
  #bachem_experiment(X = X, K = 2000)
  
#for k in results.keys():
  #for k2 in ["init", "final"]:
    #print k, "-", k2, ' \t mean sse over 10 runs : %.2e'%(np.mean(results[k][k2])) #, '\t std : %.2e'%(np.std(results[k][k2]))
    




def all_levels():
  seed =  1012
  random.seed(seed)
  npr.seed(seed)
  max_proposals = 1000000#000000
  maxtime = 200.0
  ndata = 500000
  K = 500
  maxrounds = 10000000
  dimension = 4
  basedata = 1*npr.randn(ndata, dimension)
  data = np.array(basedata, dtype = np.float64)
  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()
  
  results = {}
  for level in [0,1,2,3]:
    results[level] = pyzentas.pyzentas(ndata = ndata, dimension = dimension, sizes = None, X = data, K = K, indices_init = indices_init, algorithm = "clarans", level = level, max_proposals = max_proposals, capture_output = True, seed = seed, maxtime = maxtime, nthreads = 3, maxrounds = maxrounds, patient = True, metric = "l2", rooted = False, energy = 'quadratic') #, energy = 'squarepotential', critical_radius = 5.) #, energy = 'exp', exponent_coeff = 0.1)#

  return results


def series_plots_app(all_level_results):
  
  pl.figure(figsize = (9,3))
  pl.subplots_adjust(bottom = 0.2, left = 0.2, wspace = 0.3, hspace = 0.3)
  fn = "/idiap/home/jnewling/kmlocal/kmlocal-1.7.2/bin/out.txt"
  filly = open(fn)
  lines = filly.readlines()
  mses = []
  times = []
  for l in lines:
    if "<stage:" in l:
      mses.append(float(l.split("best: ")[1].split()[0]))
      times.append(float(l.split("time: ")[1].split()[0]))
  
  times = np.array(times)
  times += 0.4

  for k in range(2):
    pl.subplot(1,2,k+1)
    pl.plot(times, mses, label = "kmlocal")
    
    
  
    for l in [0,1,2,3]:
      X = get_energy_time(all_level_results[l]['output'])
      pl.plot(X['ttime']/1000., X['E'], label = 'level %s'%(l,))
      
    pl.xlabel("time [s]")
    pl.ylabel("MSE")
    pl.ylim(ymax = 0.4)
    pl.xlim(xmax = 180)
      #Tracer()()
      #pl.plot
    if k == 1:
      pl.xscale('log', basex = 10)
      pl.xlim(0.2, 180)
    
    else:
      pl.legend(labelspacing = 0.05)

  fn = "/idiap/home/jnewling/ftex/papers/arxiv/clarans/v_nips/comparem.pdf"
  pl.savefig(fn)
  import commands
  commands.getstatusoutput("pdfcrop %s %s"%(fn,fn))
