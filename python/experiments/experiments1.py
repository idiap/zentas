# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>


import matplotlib.pyplot as pl
import sys
import os
import random
import numpy as np
import numpy.random as npr
import cPickle

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../../build/python")
import pyzentas
from IPython.core.debugger import Tracer 
import time

import matplotlib 
import datapaths


# Font solution, based on stack-overflow page
# https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
matplotlib.rcParams['text.latex.preamble']= r"\usepackage{cmbright} \renewcommand{\familydefault}{fos} \renewcommand{\seriesdefault}{l} \renewcommand{\bfdefault}{sb}"
matplotlib.rcParams['axes.linewidth'] = 0.4

#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : "sans-serif",
          'text.latex.unicode': True,
}
matplotlib.rcParams.update(params) 


import matplotlib.gridspec as gridspec
pl.ion()


#for loading datasets.
import rna
reload(rna)
import load_csv
reload(load_csv)
import joensuu
reload(joensuu)
import mnist
reload(mnist)



import datapaths

joensuu_datasets = ["dim032", "MopsiLocationsUntil2012-Finland", "europediff", "yeast", "housec8", "KDDCUP04Bio", "a1", "a2", "a3", "birch2"]
all_datasets = joensuu_datasets + ["htru", "yearpredictionmsd", "mopac", "mnist", "dna"]


#all_datasets = ["dim032", "yeast", "housec8", "a1", "a2", "a3", "birch2"]

print_names = {
"KDDCUP04Bio" : "kdd04", 
"a1" : "a1", 
"a2" : "a2", 
"a3" : "a3", 
"birch2" : "birch2", 
"dim032" : "dim032", 
"MopsiLocationsUntil2012-Finland": "mopsi", 
"europediff":"europe", 
"yeast":"yeast", 
"housec8":"housec8", 
"dna":"dna", 
"htru":"htru", 
"yearpredictionmsd":"song", 
"mopac" : "motion", 
"mnist": "mnist-10"
}

alg_colors = {"kmpps" : "r", "clas" : "#000066"}
alg_labels = {"kmpps" : "$K$-means++ - $K$-means", "clas" : "$K$-means++ - CLARANS - $K$-means"}

def go(X, K, alg, seed, kmpp_greedy):
  """
  One run.
  X : data
  K : number of clusters
  """
  results = {}
  
  if "kmpp-cla" in alg:
    max_itok = float(alg.split("-")[-1])
    
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = max_itok, max_time = 500.0, max_proposals = K**2, seed = seed, patient = True, nthreads = 4, init = "kmeans++-3", with_tests = False, capture_output = True, rooted = False)

  elif alg == "kmpp":
    
    if (kmpp_greedy > 1):
      init = "kmeans++~%d"%(kmpp_greedy)
    
    else:
      init = "kmeans++-1"
    
     
    z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0.0, max_time = 5.0, max_proposals = K**2, seed = seed, patient = True, nthreads = 4, init = init, with_tests = False, capture_output = True, rooted = False)
    
  else:
    raise RuntimeError("unrecognised alg in go")
    
  tzen0 = time.time()
  if X.shape[1] > 15:
    rf_alg = "yinyang"
    do_vdimap = True
  else:
    rf_alg = "exponion"
    do_vdimap = False
    
  tangerine =  z.den(X, do_vdimap = do_vdimap, do_refinement = True, rf_alg = rf_alg, rf_max_rounds = 100) 
  tzen1 = time.time()
  results["t"] = tzen1 - tzen0
  results["out"] = pyzentas.get_processed_output(tangerine['output'])
  results['mse'] = results["out"]["mE"][-1]

  return results


def get_dimensions(dataset, ignore_cache = False):
  
  cache_fn = os.path.join(datapaths.datapaths["ds_dimcache"], "%s.dims"%(dataset,))
  if ignore_cache == False and os.path.exists(cache_fn):
    print "loading dimensions of dataset ", dataset , " from ", cache_fn
    filly = open(cache_fn, "r")
    l = filly.readline()
    filly.close() 
    N, d = [int(x) for x in l.split()]
  
  else:
    print "loading dataset ", dataset, " extracting N,d, saving and returning" 
    X = get_dataset(dataset)
    N, d = X.shape
    filly = open(cache_fn, "w")
    filly.write("%d %d"%(N, d))
    filly.close()
  
  return {"N":N, "d":d}
    
def get_dataset(dataset):
  """
  return numpy array, given string
  """
  
  if dataset in joensuu_datasets:
    X = joensuu.get_joensuu(dataset)
    if dataset == "KDDCUP04Bio":
      X = X[0:10000, :]
    
  elif dataset is "dna":
    X = rna.get_rna()#[0:160000, :]

  elif dataset is "htru":
    X = load_csv.get_htru2()
 
  elif dataset is "yearpredictionmsd":
    X = load_csv.get_yearpredictionmsd()[0:30000, :]

  elif dataset is "mopac":
    X = load_csv.get_mopac()
    
  elif dataset is "mnist":
    X = mnist.read_MNIST(dataset = "projected", ndata = 10000, dimension = 10)/1000.  
  
  else:
    raise RuntimeError("Unrecognised dataset ", dataset)


  ndata = X.shape[0] #min(X.shape[0], 50000)
  
  return X[0:ndata, :]

  
def get_pkl_fn(dataset, N, K, n_kmpps, kmpp_greedy):
  if kmpp_greedy == 0:
    return "%s_N%s_K%s_Nkmpp%d.pkl"%(dataset, N, K, n_kmpps)
  else:
    return "%s_N%s_K%s_Nkmpp%d_greed%d.pkl"%(dataset, N, K, n_kmpps, kmpp_greedy)

def get_nips_results(ignore_cache = False, dataset = "dna", n_kmpps = 10, sched = None, N = "full", K = "sqrt", kmpp_greedy = 0):  
  """
  ignore_cache : if True even if there is a pickle matching the control parameters (dataset, etc), run afresh.
  n_kmpps : number of runs of $K$-means++
  sched : one of fixed_schedule1 and time_schedule1, control how many clarans to run, and what max_itok is. 
  N : either an integer, or "full" (all the data)
  K : either an integer, or "sqrt" (sqrt(N)) or divP for K = int(N/P) or like "0.1sqrt" or "10sqrt"
  """
  
  # check if the results in stored, in which case just load them up.
  pkl_path = os.path.join(datapaths.datapaths["pkl_results_dir"], get_pkl_fn(dataset, N, K, n_kmpps, kmpp_greedy))

  if ignore_cache == False and os.path.exists(pkl_path):
    print "loading cPickled results from", pkl_path
    run_from_scratch = False
    kmpps = []
    clas = []
    filly = open(pkl_path, "r")
    n_kmpps = cPickle.load(filly)
    n_clas = cPickle.load(filly)    

    for i in range(n_kmpps):
      kmpps.append(cPickle.load(filly))

    for i in range(n_clas):
      clas.append(cPickle.load(filly))
                
  else:
    print "no pickle being loaded"
    run_from_scratch = True


  # if the results need to be generated from scratch,
  if run_from_scratch == True:
    X = get_dataset(dataset)
    if N != "full":
      N = int(N)
      if N > X.shape[0]:
        raise RuntimeError("N is larger than dataset size")
      else:
        X = X[0:N, :]

    if isinstance(K, str) and "sqrt" in K:
      if K == "sqrt":
        K = int(np.sqrt(X.shape[0]))
      else:
        K = int(float(K.split("sqrt")[0])*np.sqrt(X.shape[0]))
    elif isinstance(K, str) and "div" in K:
      factor = int(K.split("div")[-1])
      K = int(X.shape[0]/factor)
    else:
      K = int(K)
      
    
    print dataset, " ndata : ", X.shape[0], "  K : ", K  

    seed = npr.randint(10000)      
    # K-means++ --> K-means
    t0 = time.time()
    kmpps = []
    
    print "$K$-means++ [ ",
    for i in range(n_kmpps):
      print i, 
      
      aDict = go(X, K, "kmpp", seed, kmpp_greedy)
      kmpps.append(np.array([aDict["out"]["Tt"], aDict["out"]["mE"]]))
      seed = npr.randint(10000)
    
    print "]"
    time_kmpp = time.time() - t0
    
    # K-means++ --> CLARANS --> K-means
    clas = sched.get(time_kmpp, X, K, 0)
        
    print "\nopening file to cPickle...",
    filly = open(pkl_path, "w")

    cPickle.dump(len(kmpps), filly)
    cPickle.dump(len(clas), filly)
        
    for x in kmpps:
      cPickle.dump(x, filly)

    for x in clas:
      cPickle.dump(x, filly)
      
    filly.close()
  
    print "done."
  return {"kmpps":kmpps, "clas":clas}


  
def get_nips_scaled_results(ignore_cache, dataset, n_kmpps, sched, N, K, kmpp_greedy):
  """
  scale to have min E = 1 and for time to be in seconds.
  """
  results = get_nips_results(ignore_cache, dataset, n_kmpps, sched, N, K, kmpp_greedy)

  max_E, min_E = 0, 10**17
  for x in results["kmpps"]:
    max_E, min_E = max(max_E, x[1].max()), min(min_E, x[1].min())

  for x in results["clas"]:
    max_E, min_E = max(max_E, x[1].max()), min(min_E, x[1].min())

  scaled_results = {}
  for alg in results.keys():
    scaled_results[alg] = []
    for x in results[alg]:
      scaled_results[alg].append([x[0]/1000., x[1]/min_E])

  return scaled_results
  


class fixed_schedule1():
  """
  Run clarans with a fixed series of max_itoks
  """
  def get(self, time_limit, X, K, kmpp_greedy):
    seed = npr.randint(10000)
    clas = []
    max_itoks = [0.5,1,2,3,4,5,6,7,8]
    print "clarans [ "
    for i in range(len(max_itoks)):    
      print i, 
      aDict = go(X, K, "kmpp-cla-%.2f"%(max_itoks[i],), seed, kmpp_greedy)
      clas.append(np.array([aDict["out"]["Tt"], aDict["out"]["mE"]]))
      seed = npr.randint(10000)
    
    print "]"
    return clas

fix_sch = fixed_schedule1()


class time_schedule1():
  """
  Run clarans looping, through itoks while the time_limit permits.
  """
  
  def __init__(self, itoks = [2, 5, 1, 10, 16]):
    self.itoks = itoks
    
  def get(self, time_limit, X, K, kmpp_greedy):    
    seed = npr.randint(10000)
    clas = []
    t0 = time.time()
    n_runs = 0
      
    print "clarans [ ",
    while time.time() - t0 < time_limit:
      print n_runs,
      aDict = go(X, K, "kmpp-cla-%d"%(self.itoks[n_runs%len(self.itoks)]), seed, kmpp_greedy)
      clas.append(np.array([aDict["out"]["Tt"], aDict["out"]["mE"]]))
      n_runs += 1
      seed = npr.randint(10000)
    
    print " ] ", time.time() - t0, "[s]"
    return clas

tim_sch = time_schedule1()


def get_extrema(results):

  min_mse = {}
  max_mse = {}
  max_time = 0

  for alg in ["kmpps", "clas"]:
    min_mse[alg] = 10**17
    max_mse[alg] = 0
    for i in range(len(results[alg])):
      max_time = max(max_time, results[alg][i][0][-1])
      min_mse[alg] = min(min_mse[alg], results[alg][i][1][-1])
      max_mse[alg] = max(max_mse[alg], results[alg][i][1][-1])

  return {"min_mse":min_mse, "max_mse":max_mse, "max_time":max_time}

def nips_base_plot(results, allpoints = True, trace = False):
  
  def plot_curve(Time, Energy, alg, allpoints = False, trace = False, withlabel = False):
    
    if allpoints:
      pl.plot(Time, Energy, markersize = 0.4, marker = ".", linestyle = "none", alpha = 0.5, color = alg_colors[alg], label = None)
    
    if trace:
      pl.plot(Time, Energy, color = "0.4", linestyle = "-", linewidth = 0.1, alpha = 0.1, label = None)  

    label = None if not withlabel else alg_labels[alg]
    pl.plot([Time[-1]], [Energy[-1]], linestyle = "none", marker = "x", markersize = 4, alpha = 0.7, color = alg_colors[alg], label = label)

  
  for alg in ["kmpps", "clas"]:
    for i in range(len(results[alg])):
      plot_curve(results[alg][i][0], results[alg][i][1], alg, allpoints = allpoints, trace = trace, withlabel = (i == 0))
  
  
  


    
def nips_plot1(dataset = "dna", savefig = True):

  sched = time_schedule1(itoks = [3])
  
  #results = get_nips_scaled_results(ignore_cache = False, dataset = dataset, n_kmpps = 47, sched = sched, N = "full", K = 400, kmpp_greedy = 0)


#"KDDCUP04Bio" : "kdd04", 
#"a1" : "a1", 
#"a2" : "a2", 
#"a3" : "a3", 
#"birch2" : "birch2", 
#"dim032" : "dim032", 
#"MopsiLocationsUntil2012-Finland": "mopsi", 
#"europediff":"europe", 
#"yeast":"yeast", 
#"housec8":"housec8", 
#"dna":"dna", 
#"htru":"htru", 
#"yearpredictionmsd":"song", 
#"mopac" : "motion", 

  
  results = get_nips_scaled_results(ignore_cache = False, dataset = "yeast", n_kmpps = 50, sched = sched, N = "full", K = "sqrt", kmpp_greedy = 5)

  pl.figure(1, figsize = (6,1.9))
  pl.ion()
  pl.clf()
    
  extrema = get_extrema(results)
  nips_base_plot(results)

  pl.ylim(ymin = 0.95)
  xmax = 1.02*extrema["max_time"]
  pl.xlim(xmin = 0, xmax = xmax)
  pl.plot([0,xmax], [extrema["min_mse"]["kmpps"], extrema["min_mse"]["kmpps"]], color = alg_colors["kmpps"], alpha = 0.2, linewidth = 0.4)
  pl.plot([0,xmax], [extrema["min_mse"]["clas"], extrema["min_mse"]["clas"]], color = alg_colors["clas"], alpha = 0.2, linewidth = 0.4)



  pl.xlabel("time [s]")#, fontproperties = fontprop_normal)
  pl.ylabel("$E$")#, fontproperties = fontprop_italic) #, fontproperties = fontpropx)
  leg = pl.legend(fontsize = "medium", borderpad = 0.2, numpoints = 1, labelspacing = 0.3, columnspacing = 0.2)
  frame = leg.get_frame()
  frame.set_linewidth(0.1)
  frame.set_edgecolor("none")

  pl.subplots_adjust(bottom = 0.25, left = 0.2)
  
  if savefig:
    fn = datapaths.datapaths["nips_plot1"]
    pl.savefig(fn)
    print "saving as", fn
    import commands
    bla = commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))


def nips_plot2():

  sched = time_schedule1()
  savefig = False

  results = {}
  for ds in all_datasets:
    print "In nips_plot2, dataset", ds    
    results[ds] = get_nips_scaled_results(ignore_cache = False, dataset = ds, n_kmpps = 50, sched = sched, N = "full", K = "10sqrt")

  pl.figure(1)
  gs = gridspec.GridSpec(3,5)
  
  for dsi, ds in enumerate(all_datasets):

    xloc = dsi%3
    yloc = dsi/3
    ax = pl.subplot(gs[xloc:xloc+1, yloc:yloc+1])
    nips_base_plot(results[ds], False, False)

    with_descriptors = False
    if with_descriptors:
      extrema = get_extrema(results[ds])
      min_mses = extrema["min_mse"]
      max_mses = extrema["max_mse"]
  
  
      glo_max = max(max_mses["clas"], max_mses["kmpps"])
      delta_E = glo_max - 1
      pl.xlim(xmin = 0)
      pl.ylim(ymin = 1 - 0.08*delta_E, ymax = 1 + delta_E*1.33)
      ytick2 = max(1.01, float("%.2f"%(1 + 0.5*delta_E,)))
      
      xmax = 1.5*extrema["max_time"]
      pl.plot([0,xmax], [min_mses["kmpps"], min_mses["kmpps"]], color = alg_colors["kmpps"], alpha = 0.2, linewidth = 0.5)
      pl.plot([0,xmax], [1, 1], color = alg_colors["clas"], alpha = 0.2, linewidth = 0.5)    
      pl.yticks([min_mses["kmpps"]], ["%.3f"%(min_mses["kmpps"],)])    
      xtickas = pl.xticks()
      pl.xticks([0, xtickas[0][-3]], ["0", xtickas[0][-3]])
      pl.subplots_adjust(hspace = 0.3, wspace = 0.25)
    
    else:
      pl.xticks([])
      pl.yticks([])
      pl.subplots_adjust(hspace = 0.01, wspace = 0.01)
  
  if savefig:
    fn = datapaths.datapaths["nips_plot2"]
    pl.savefig(fn)
    import commands
    bla = commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))
    

def plot3_finishing(savekey):
  dimensions = {}
  for dsi, ds in enumerate(all_datasets):
    #can ignore cache here, there's barely any speed advantage in not ignoring it.
    dimensions[ds] = get_dimensions(ds, ignore_cache = False)
  
  for dsi, ds in enumerate(all_datasets):
    print dimensions[ds], ds
  
  pl.xlabel("($E$ without CLARANS) / ($E$ with CLARANS)", weight = "light") #fontproperties=fontprop)
  pl.ylabel("fraction of datasets", weight = 'light')
  pl.xticks([1, 1.05, 1.15, 1.25], ["1", "1.05","1.15","1.25"])
  pl.subplots_adjust(bottom = 0.24, left = 0.24)
  pl.xlim(xmin = 0.995)
  leg = pl.legend(fontsize = "medium", labelspacing = 0.3)
  frame = leg.get_frame()
  frame.set_linewidth(0.3)
  

  if savekey != None and savekey != "none":
    fn = datapaths.datapaths[savekey]
    pl.savefig(fn)
    import commands
    bla = commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))


def nips_plot3x(kmpp_greedy, fignum = 5, in_gray = False):    
    
  sched = time_schedule1()  

  colors = {"10sqrt":"#800909", "sqrt":"#f27d0c", "0.1sqrt":"#fdcf58"}
  labels = {"10sqrt":"$K=10\sqrt{N}$", "sqrt":"$K=\sqrt{N}$", "0.1sqrt":"$K=0.1\sqrt{N}$"}
  alphas = {}
  for k in colors.keys():
    alphas[k] = 1
    
  if in_gray:
    for k in colors.keys():
      labels[k] = None
      alphas[k] = 0.1
      
  results = {}
  ratios = {}

  pl.figure(fignum, figsize = (4.3,2.4))
  for Ki, K in enumerate(colors.keys()):
      
    results[K] = {}    
    for ds in all_datasets:
      print "In nips_plot2, dataset", ds    
      
      if ds == "dna" and kmpp_greedy != 0 and K == "10sqrt" :
        N = 10000
      else:
        N = "full"
      results[K][ds] = get_nips_scaled_results(ignore_cache = False, dataset = ds, n_kmpps = 50, sched = sched, N = N, K = K, kmpp_greedy = kmpp_greedy)
  
    ratios[K] = []
    for dsi, ds in enumerate(all_datasets):
      extrema = get_extrema(results[K][ds])  
      min_mses = extrema["min_mse"]
      max_mses = extrema["max_mse"]
      ratios[K].append(min_mses["kmpps"]/min_mses["clas"])
  
    ratios[K] = np.array(ratios[K])
    ratios[K].sort()
    
    Y = np.linspace(1, 0, len(ratios[K]))
    pl.plot(ratios[K], Y, marker = '.', markersize = 10, color = colors[K], label = labels[K], alpha = alphas[K])
  
  
    

def nips_plot3():
  nips_plot3x(kmpp_greedy = 0)
  plot3_finishing(savekey = "nips_plot3")

def nips_plot3_greedy():
  """
  version where the k-means++ -> k-means pipeline uses k-means++~5, that is 5 attempts
  for every addition of a sample to the seeding set, the one (of five) which reduced the
  energy the most is actually added. 
  """
  nips_plot3x(kmpp_greedy = 0, in_gray = True)
  nips_plot3x(kmpp_greedy = 5)
  plot3_finishing(savekey = "nips_plot3_greedy")


