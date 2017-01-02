# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

import hardpaths
reload(hardpaths)

import sys
sys.path.append(hardpaths.zentas_lib_dir)

sys.path.append("../../build/python")
import pyzentas


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


import copy


import hardpaths

import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec

import matplotlib 
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

from IPython.core.debugger import Tracer
pl.ion()



import load_joensuu_data
reload(load_joensuu_data)


NRUNS_KMEANSPP = 40
ALGS =  ["un", 'pp', "BF", 'cl_s1.0', 'cl_s3.0']


colors ={"un":"#dab595",'pp':"#48464a","BF":"#008b50",'cl_s1.0':"#f19e34",'cl_s3.0':"#cc252a", "none":'w', 'cl_s3.0-mc2':'k', 'cl_s1.0-mc2':'k'}

def get_alg_color(alg):
  if alg == 'none':
    return 'w'
  else:
    return colors[ALGS.index(alg)]


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

def get_label(alg):
  label = alg
  k = alg
  if k == 'pp':
    label = r'\texttt{k-means}++' #$\mathbf{k}$-means++'
  elif k == 'un':
    label = r'\texttt{uni}'
  elif k == 'cl_s1.0':
    label = r'\texttt{cl-1}'
  elif k == 'cl_s3.0':
    label = r'\texttt{cl-3.0}'
  elif k == 'BF':
    label = r'\texttt{bf}'
  elif k == 'cl_s1.0-mc2':
    label = r'\texttt{cl-1.0-mc2}'
  elif k == 'cl_s3.0-mc2':
    label = r'\texttt{cl-3.0-mc2}'



  return label

  


linestyles = {
"un":"-",
'pp':"-",
"BF":"-",
'cl_s1.0':"-", 
'cl_s3.0':"-"}


linewidths = {
"un":0.5,
"pp":0.5,
"BF":0.5,
"cl_s1.0":0.5,
"cl_s3.0":0.5
}

markers = {
"un":"+",
"pp":"+",
"BF":"+",
"cl_s1.0":"x",
"cl_s3.0":"x",
"cl_s3.0-mc2":"+",
"cl_s1.0-mc2":"+"
}


markersizes = {
"un":5,
"pp":5,
"BF":5,
"cl_s1.0":5,
"cl_s3.0":5,
"cl_s1.0-mc2":5,
"cl_s3.0-mc2":5
}


fullalgname = {
"un":r'\texttt{uni}',
"pp":r'\texttt{k-means}++',
"BF":r'\texttt{bf}',
"cl_s1.0":r'\texttt{cl-1.0}',
"cl_s3.0":r'\texttt{cl-3.0}',
"cl_s1.0-mc2":r'\texttt{cl-1.0-mc2}',
"cl_s3.0-mc2":r'\texttt{cl-3.0-mc2}'
}


# The K-s to run with
trueks = {}
trueks['ConfLongDemo_JSI_164860'] =  11
trueks['birch1'] = 100
trueks['birch2'] = 100
trueks['birch3'] = 100
trueks['s1']  = 15
trueks['s2']  = 15
trueks['s3']  = 15
trueks['s4']  = 15
trueks['a1']  = 20
trueks['a2']  = 35
trueks['a3']  = 50
trueks['dim032']  = 16
trueks['dim064']  = 16
trueks['dim1024'] = 16
#trueks['yeast'] =  20
trueks['housec8'] = 200 #-1
trueks['MopsiLocationsUntil2012-Finland'] =  50 #100#-1
trueks['mnist'] = 150 #10
trueks['europediff'] = 600 #-1
trueks['KDDCUP04Bio'] = 300 #2000

for a in trueks.keys():
  trueks[a]*=2
  trueks[a] = int(trueks[a])


  
  

def initialisation_test(X, K, algs = ALGS, dataset = None, n_threads = 3):
  """
  (1) run k means ++ and lloyd with K centers on X. Set TL to be NRUNS_KMEANSPP times this.
  (2) for alg in algs : run series of alg + lloyd while TL not exceeded.
  (3) return the results
  """
  
  ndata, dimension = X.shape
  seed = npr.randint(100000)
  random.seed(seed)
  npr.seed(seed)
  print "seed : " , seed

  #see how long a run with kmeans++ takes. 
  print "getting time of complete run with kmeans++"
  
  tstart = time.time()
  
  out_kmeanspp = kmeans.get_clustering(X = X, n_clusters = K, algorithm = 'auto', init = "kmeans++", verbose = 2, capture_verbose = True, n_threads = n_threads)
  
  
  tstop = time.time()
  
  t_elapsed = tstop - tstart
  print "t_elapsed : ", t_elapsed
  kmeanspp_runtime = t_elapsed 
  
  time_per_alg = kmeanspp_runtime*NRUNS_KMEANSPP

  #cl_x : x is 
  # (if txxx) time spent in clarans, fraction of time taken for 1 run of out_kmeanspp.
  # (if kxxx) number of swaps, 
  results = dict.fromkeys(algs)
  for k in algs:
    times = []
    mses = []
    t_start = time.time()
    
    while time.time() - t_start < time_per_alg:
      print k, time.time() - t_start, " / ", time_per_alg
      seed = npr.randint(100000)
     
      if 'cl_' in k:
        switch_type =  k.split('_')[-1][0]
       
        if switch_type == 't':
          maxtime = kmeanspp_runtime*float(k.split('_')[-1][1::])
          maxrounds = 1000000
          max_proposals = 1000000
        elif switch_type == 's':
          maxtime = 100000000.
          maxrounds = int(K*float(k.split('_')[-1][1::]))
          max_proposals = 10000
          
        indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
        indices_init.sort()
        
        
    
        
        bla3 = pyzentas.pyzentas(ndata = ndata, dimension = dimension, X = X, K = K, indices_init = indices_init, algorithm = "clarans", level = 3, max_proposals = max_proposals, capture_output = True, seed = seed, maxtime = maxtime, nthreads = n_threads, maxrounds = maxrounds, patient = True, metric = 'l2', rooted = False, energy = 'quadratic')
        bla3['indices_final'].sort()
        out_clarans_3 = kmeans.get_clustering(X = X, n_clusters = K, algorithm = 'auto', init = bla3['indices_final'], verbose = 2, capture_verbose = True, n_threads = n_threads)
        mses.append(float(out_clarans_3['output'].split("\n")[-3].split(" : ")[1].split()[0]))
        times.append(time.time() - t_start)
    
        
        
      else:
        if k == 'un':
          init = "uniform"
        elif k == "pp":
          init = "kmeans++"
        elif k == "BF":
          init = "BF"
          
          
        else:
          raise RuntimeError("unrecognised cl : " + cl)
        
        out = kmeans.get_clustering(X = X, n_clusters = K, algorithm = 'auto', init = init, verbose = 2, capture_verbose = True, seed = seed, n_threads = n_threads)
        mses.append(float(out['output'].split("\n")[-3].split(" : ")[1].split()[0]))
        times.append(time.time() - t_start)
    
    results[k] = dict.fromkeys(['times', 'mses'])
    results[k]['times'] = np.array(times)
    results[k]['mses'] = np.array(mses)
  
  
  if dataset != None:
    runtime_write_fn = os.path.join(hardpaths.initialisation_experiments_runtimes_dir, dataset + ".txt")
    
    filly = open(runtime_write_fn, "w")
    filly.write(t_elapsed)
    filly.close()

  results['t_elapsed_1_kmeanspp'] = t_elapsed
  return results
  


def joensuu_experiment(dataset = "MopsiLocationsUntil2012-Finland", K = 15, writedata = True):
  """
  (1) Read data
  (2) Get results via initialisation_test
  (3) Write results
  """
  X = load_joensuu_data.load_data(dataset)
  scores = initialisation_test(X, K)
  resultsdir = hardpaths.elapsed_time_kmeanspp
  
  if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)
    
  filly = open(os.path.join(resultsdir, "%s.txt"%(dataset,)), "w")
  value = scores['t_elapsed_1_kmeanspp']
  #print value
  filly.write("%.5f"%(value))
  filly.close()
  
  if writedata:
    resultsdir = hardpaths.initialisation_result_pickles        
    filly = open(os.path.join(resultsdir, "%s.pkl"%(dataset,)), "w")
    cPickle.dump(scores, filly)
    filly.close()
  
  

def all_experiments(writedata = True):  
  """
  All rank vs energy plots.
  """
  alldatasets = trueks.keys()
  alldatasets.sort(key = str.lower)
  for ik, k in enumerate(alldatasets):
    print "\n"
    print "********************    ", ik, "  :  ", k, "    ***********************"
    print "\n"
    bla = joensuu_experiment(k, trueks[k], writedata = writedata)
    

def grid_data_initialisation_test():
  """
  100 centers (0,0) ... (0,9) ... (9,0) ... (9,9). points 100 per cluster, center + N(0,I_2)
  """
  K = 100
  nperc = 100
  rK = int(np.sqrt(K))
  ndata = nperc*K
  X = np.empty(shape = (ndata, 2), dtype = np.float64)
  print nperc, rK
  X[:,0] = np.repeat(np.arange(rK), nperc*rK)
  X[:,1] = np.tile(np.arange(rK), nperc*rK)
  X += 1.0*npr.randn(ndata, 2)  
  time_per_alg = 14
  scores = initialisation_test(X, K, time_per_alg)
  plot_initialisation_test(scores)
  return scores
  
  








def plot_initialisation_test(scores, with_legend = False,  algs = ALGS):
  """
  Plots energies ranked from best run to worst run. Does not capture time.
  """
  fig = pl.gcf()
  ax = pl.gca()
  pl.cla()

  minscore = 10**55
  maxscore = 0
  for ik, k in enumerate(algs):
    minscore = min(minscore, scores[k]['mses'].min())
    maxscore = max(maxscore, scores[k]['mses'].max())

  lines = []
  labels = []
  
  for ik, k in enumerate(algs):
    color = colors[k]
    scores_sorted = scores[k]['mses'].copy()
    scores_sorted.sort()
    cum_counts = 1 - np.linspace(0, 1, scores_sorted.size)
    label = fullalgname[k]
    line = pl.plot(scores_sorted/minscore, cum_counts, marker = markers[k], linestyle =  linestyles[k], linewidth = linewidths[k], markersize = markersizes[k], label = label, color = color, alpha = 1.0)
    labels.append(label) 
    lines.append(line[0])
    
  if with_legend:
    pl.legend(loc = 'lower right')

  xlims = pl.xlim()
  xlim_diff = xlims[1] - xlims[0]
  pl.xticks([1, xlims[1]], [1, xlims[1]])#], [])
  
  pl.xlim([1 - 0.05*xlim_diff, xlims[1] + 0.05*xlim_diff])
  pl.ylim([-0.03,1.03])
  pl.yticks([])
  pl.subplots_adjust(top = 0.98)
  
  return lines, labels
    


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
'mnist']








  
def write_statistics():
  """
  Write the statistics : dataset, K, N, dim, TL. 
  """
  elapsed_dir = hardpaths.elapsed_time_kmeanspp
  basedir = hardpaths.clarans_paper_dir
  fn = os.path.join(basedir, "dataset_stats.txt")
  filly = open(fn, "w")
  datasets = trueks.keys()
  datasets.sort(key = str.lower)
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
      
    filly.write(r"""%s  &  %d &  %d &  %d &   %.2f \\
"""%(get_abridged(k), ndata, dimension, K, allocated_time))
  filly.close()



def plot_all_experiments(save = False, filename = "initplots.pdf"):
  resultsdir = hardpaths.initialisation_result_pickles
  pl.figure(num = 8, figsize = (6, 13))
  pl.clf()
  datasets = trueks.keys()
  datasets.sort(key = str.lower)
  for ik, dataset in enumerate(datasets):
    filly = open(os.path.join(resultsdir, "%s.pkl"%(dataset,)), "r")
    scores = cPickle.load(filly)
    filly.close()
    pl.subplot(7,3,ik+1)
      
    lines, labels = plot_initialisation_test(scores, with_legend = False)

    if ik%3 == 0:
      pl.yticks([0,1], [0,1])


    ax = pl.gca()
    text = dataset
    
    text = get_abridged(dataset)
    
    pl.text(0.96, 0.97, text, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    
    if ik in [18, 19]:
      pl.xlabel('energy', verticalalignment = 'bottom')
  
    if ik %3 == 0:
      pl.ylabel('fraction', verticalalignment = 'top')
    
    
  pl.subplot(7,3,21)
  
  lines = [copy.copy(l) for l in lines]
  for l in lines:
    l.set_linewidth(2)
    
  #Tracer()()
  pl.legend(lines, labels, loc = 'center', frameon = False, fontsize = 'medium')
  ax = pl.gca()
  ax.axis('off')
  
  
  if save: 
    ##pl.subplots_adjust(hspace = 0.14, wspace = 0.05)
    fname = os.path.join(hardpaths.clarans_paper_dir, filename)
    pl.savefig(fname)
    import commands
    commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))
  
  
  
  

def get_all_events(threshold = 1.01):
  """
  Get "The Champions league" data.
  """
  resultsdir = hardpaths.initialisation_result_pickles
  results = {}
  for ik, dataset in enumerate(trueks.keys()):
    filly = open(os.path.join(resultsdir, "%s.pkl"%(dataset,)), "r")
    scores = cPickle.load(filly)
    filly.close()
    results[dataset] = scores

  
  all_events = {}
  for ik, dataset in enumerate(trueks.keys()):
    all_mses = []
    all_times = []
    all_algs = []
    
    maxtime = NRUNS_KMEANSPP*results[dataset]['t_elapsed_1_kmeanspp'] + 0.
    
    for alg in ALGS:#results[dataset].keys():
      #print dataset, alg
      all_mses.append(results[dataset][alg]['mses'])
      all_times.append(results[dataset][alg]['times'])
      all_algs.append([alg]*results[dataset][alg]['times'].size)
    
    all_mses = np.concatenate(all_mses)
    all_times = np.concatenate(all_times)
    all_algs = np.concatenate(all_algs)
    
    order = np.argsort(all_times)
    all_mses = all_mses[order]
    all_times = all_times[order]
    all_algs = all_algs[order]
    
    events = [(9.999**19.999, 0, ['none'])] # (lowest energy, time, algs within threshold of energy)
    
    alg_lowest = dict.fromkeys(ALGS)
    for a in ALGS:
      alg_lowest[a] = 9.999**19.999
      
    for i in range(all_mses.size):
      if all_mses[i] < alg_lowest[all_algs[i]]:
        alg_lowest[all_algs[i]] = all_mses[i]
      
      if all_mses[i] <= threshold*events[-1][0]:
        event_time = all_times[i] /  maxtime
        lowest_mse = min(all_mses[i], events[-1][0])
        current_bests = [a for a in ALGS if alg_lowest[a] <= threshold*lowest_mse]
        events.append((lowest_mse, event_time, current_bests))
    
    events.append((-1, 1., ['none']))
    all_events[dataset] = events
  
  return all_events
  
 
def colorblob(threshold = 1.01, filename = 'colorblob.pdf'):
 
  pl.close('all')
  all_events = get_all_events(threshold)
  pl.figure(198, figsize = (8, 10))
  pl.clf()
  ax = pl.gca()
  #pl.clf()
  
  artists = []
  
  xibob = 1.6
  datasets = trueks.keys()
  datasets.sort(key = str.lower)
  datasets = datasets[-1::-1]
  for ik, dataset in enumerate(datasets):
    pl.plot([0, 1], [xibob*ik, xibob*ik], 'k', linewidth = .5)
    pl.plot([0, 1], [xibob*(ik)+1, xibob*(ik)+1], color = 'k', linewidth = 0.5) #(0.5, 0.5, 0.5)
    events = all_events[dataset]
    n_events = len(events) - 1
    for i in range(n_events):
      n_best = len(events[i][2])
      

      for best_i in range(n_best):
        trect = pl.Rectangle([events[i][1], xibob*ik + best_i/(n_best + 0.)], events[i+1][1] - events[i][1], 1./n_best, facecolor = colors[events[i][2][best_i]], edgecolor = 'none')
        artists.append(trect)
  
    for artsy_type in artists:
      ax.add_patch(artsy_type)

  pl.yticks(0.5 + xibob*np.arange(len(datasets)), [get_abridged(x) for x in datasets])
  pl.xlim(0, 1)
  pl.ylim(0,xibob*len(datasets) - (xibob - 1))
  pl.subplots_adjust(left = 0.43, bottom = 0.15, top = 0.95)

  pl.xlabel('fraction of time limit (TL) elapsed')
  fname = os.path.join(hardpaths.clarans_paper_dir, filename)
  pl.savefig(fname)
  
  import commands
  commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))



def colorblob_cheat(dataset = 'yeast', threshold = 1.01, filename = 'colorblob_example.pdf'):

  resultsdir = hardpaths.initialisation_result_pickles
  results = {}
  for ik, ds in enumerate(trueks.keys()):
    filly = open(os.path.join(resultsdir, "%s.pkl"%(ds,)), "r")
    scores = cPickle.load(filly)
    filly.close()
    results[ds] = scores
    

  all_events = get_all_events(threshold)

  min_E = 10**30
  for a in ALGS:
    min_E = min(min_E, results[dataset][a]['mses'].min())
  
  
  pl.figure(1011)
  pl.clf()

  gs = gridspec.GridSpec(4,1)
  
  print dataset, min_E
  pl.subplot(gs[0:3, 0:1])  
  for a in ALGS:
    pl.plot(results[dataset][a]['times'], results[dataset][a]['mses']/min_E, color = colors[a], marker = markers[a], linestyle = '.', markersize = 7, markeredgewidth = 2)
    
  ylim = pl.ylim()
  yrange = ylim[1] - ylim[0]
  pl.ylim(ylim[0] - yrange*0.05, ylim[1] + yrange*0.05)
  
  pl.xticks([])

      
  pl.subplot(gs[3:4, 0:1])
  ax = pl.gca()
  events = all_events[dataset]
  n_events = len(events) - 1

  artists = []
  for i in range(n_events):
    n_best = len(events[i][2])
    

    for best_i in range(n_best):
      trect = pl.Rectangle([events[i][1], 0 + best_i/(n_best + 0.)], events[i+1][1] - events[i][1], 1./n_best, facecolor = colors[events[i][2][best_i]], edgecolor = 'none')
      artists.append(trect)

  for artsy_type in artists:
    ax.add_patch(artsy_type)

  #pl.yticks(0.5 + np.arange(len(datasets)), [get_abridged(x) for x in datasets])
  pl.yticks([])
  pl.xlim(0, 1)
  pl.ylim(0,1)#len(datasets))

  #good_xlims = pl.xlim()
  
  pl.subplot(gs[0:3, 0:1])
  maxtime = NRUNS_KMEANSPP*results[dataset]['t_elapsed_1_kmeanspp'] + 0.
  pl.xlim(xmax = maxtime)

  pl.subplots_adjust(bottom = 0.6)

  fname = os.path.join(hardpaths.clarans_paper_dir, filename)
  pl.savefig(fname)
  
  import commands
  commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))



def clean_table(ALGSused =  ['cl_s1.0', 'cl_s3.0', 'pp']):

  datasets = trueks.keys()
  datasets.sort(key = str.lower)
  resultsdir = hardpaths.initialisation_result_pickles
  results = {}
  
  #Get the mimumum of ALGSused
  MINS = {}
  for ikop, dataset in enumerate(datasets):
    filly = open(os.path.join(resultsdir, "%s.pkl"%(dataset,)), "r")
    scores = cPickle.load(filly)
    #print dataset, scores
    filly.close()
    results[dataset] = scores

  for ikop, dataset in enumerate(datasets):
    print dataset
    nfactor = results[dataset]['pp']['mses'].mean()
    for alg in ALGSused:
      print alg,
      print results[dataset][alg]['mses'].size, 
      print results[dataset][alg]['mses'].mean()/nfactor
      #print results[dataset][alg]['mses'].std()/nfactor,
      #print results[dataset][alg]['mses'].min()/nfactor
      
  Tracer()()

#def mintracker(logy = True, ALGSused =  ['cl_s1.0', 'cl_s3.0', 'pp'], savefig = False):#, "un", "BF"]):

def mintracker(logy = True, ALGSused =  ['cl_s1.0-mc2', 'cl_s1.0', 'cl_s3.0', 'pp'], savefig = False):#, "un", "BF"]):
    
  figname = 'mintracker_%dalgs.pdf'%(len(ALGSused),)
 
 
  pl.close('all')
  pl.figure(1983, figsize = (11, 10.5)) #width, height.
  pl.clf()
  ax = pl.gca()
  
  xibob = 1.5
  datasets = trueks.keys()
  datasets.sort(key = str.lower)
  
  resultsdir = "/idiap/user/jnewling/zentasoutput/icml/pickles"
  #elapsed_dir = os.path.join(resultsdir, "elapsed_times_pp")
  #pickles_dir = os.path.join(resultsdir, "pickles")


  #resultsdir = #hardpaths.initialisation_result_pickles
  results = {}
  
  #Get the mimumum of ALGSused
  MINS = {}
  for ikop, dataset in enumerate(datasets):
    filly = open(os.path.join(resultsdir, "%s.pkl"%(dataset,)), "r")
    scores = cPickle.load(filly)
    #print dataset, scores
    #Tracer()()
    filly.close()
    results[dataset] = scores
    MINS[dataset] = 10**30
    for algi, alg in enumerate(ALGSused):
      MINS[dataset] = min(MINS[dataset], results[dataset][alg]['mses'].min())

  kmeanspp_final_relative = {}
  for iid, dataset in enumerate(datasets):
        
    pl.subplot(7, 3, iid + 1)
    
    figns_xlab = [18, 19, 20]      
    if iid + 1 in figns_xlab:
      pl.xlabel('fraction of TL elapsed')
      pl.xticks([0,1], [0,1])
    
    else:
      pl.xticks([])
    
    if iid == 9:
      pl.ylabel('$\log_{1.01} (E/E_{min})$', fontsize = 'large')

    
    lines = []
    MAXE = 0
    for algi, alg in enumerate(ALGSused):
      mins = []
      min_times = []
      
      #get the mimumum curve data 
      if True: #len(results[dataset][alg]['mses'].tolist()) > 2:
      #######################  
      
        print "---",len(results[dataset][alg]['mses'].tolist()), len(mins)
        for evi, evmse in enumerate(results[dataset][alg]['mses'].tolist()):
          if evi == 0 or evmse < mins[-1]:
            if results[dataset][alg]['times'][evi] < results[dataset]['t_elapsed_1_kmeanspp']*NRUNS_KMEANSPP:
              mins.append(evmse)
              min_times.append(results[dataset][alg]['times'][evi])
        
        
        yvalstoplot = np.log10(results[dataset][alg]['mses']/MINS[dataset])/np.log10(1.01)
        
        if alg in ALGSused: #   'un',  # ["BF"] #  ["BF"] # # ["BF"]# 'cl_s1.0', , 'cl_s4.0'
          MAXE = max(MAXE, yvalstoplot.max())
          
        pl.plot(results[dataset][alg]['times']/(results[dataset]['t_elapsed_1_kmeanspp']*NRUNS_KMEANSPP), yvalstoplot, color = colors[alg], marker = markers[alg], linestyle = '.', markersize = 4)
        
        if len (mins) == 0:
          raise RuntimeError ("No mins, is this an algorithm with no completions ? ")
        
        
        mins = np.array(mins).repeat(2)
        min_times = np.array(min_times).repeat(2)
        min_times = min_times[1::].tolist()
        min_times.append(results[dataset]['t_elapsed_1_kmeanspp']*NRUNS_KMEANSPP)
        min_times = np.array(min_times)/min_times[-1]
        label = get_label(alg)
        
        yvalstoplot = np.log10(mins/MINS[dataset])/np.log10(1.01)
          
        line = pl.plot(min_times, yvalstoplot, linestyle = '-', linewidth = 2, label = label, color = colors[alg], alpha = 0.6)
        lines.append(line[0])
      
        if 'pp' in alg:
          kmeanspp_final_relative[dataset] = mins[-1]/MINS[dataset]
      
      #######################
        
    ydiff = MAXE
    pl.ylim( - ydiff*0.05, MAXE + ydiff*0.15)
    
    pl.xlim([0,1])
    text = get_abridged(dataset)
    ax = pl.gca()

    pl.text(0.97, 0.94, text, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, bbox=dict(edgecolor = 'w', facecolor='white', alpha=0.1))

    ab = pl.yticks()
    a_a = copy.copy(ab[0])
    b_b = copy.copy(ab[1])
      
      
    if a_a.size > 0:
      if a_a[-3] == int(a_a[-3]):
        pl.yticks([0, int(a_a[-3])], [0, int(a_a[-3])])
      else:
        pl.yticks([0, a_a[-3]], [0, a_a[-3]])


  
  lines = [copy.copy(l) for l in lines]
  for l in lines:
    l.set_linewidth(2)
  
  labels = [get_label(alg) for alg in ALGSused]
  print labels
    
  pl.subplot(7, 3, 21)  
  
  if len(labels) > 3:
    nlegcols = 2
  else:
    nlegcols = 1
    
  leg = pl.legend(lines, labels, loc = (-0.05,-0.2), frameon = False, fontsize = 'large', ncol = nlegcols, columnspacing = 0.5)
  
  
  ax = pl.gca()
  ax.axis('off')

  pl.subplots_adjust(wspace = 0.18, hspace = 0.1, bottom = 0.1, top = 0.98)

  if savefig == True:
    fname = os.path.join(hardpaths.clarans_paper_dir, figname)
    pl.savefig(fname)
    
    import commands
    commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))


  
  sumofrels = 0
  for k in datasets:
    print '%s\t & %.4f \\\\'%(get_abridged(k), kmeanspp_final_relative[k])
    sumofrels += kmeanspp_final_relative[k]
  
  
  print "\nMEAN : ", sumofrels/len(kmeanspp_final_relative.keys())
  

  
  return leg


def make_all_figures():
  #mintracker(logy = True, ALGSused =  ['cl_s1.0', 'cl_s3.0', 'pp',"BF", "un"])
  leg = mintracker(logy = True, ALGSused =  ['cl_s1.0', 'cl_s3.0', 'pp'])
  #colorblob()
  #colorblob_cheat()
  #write_statistics()
  #plot_all_experiments()
