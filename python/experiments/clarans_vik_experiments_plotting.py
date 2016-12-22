# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

import hardpaths
reload(hardpaths)

import numpy as np
import matplotlib.pyplot as pl
pl.ion()
from IPython.core.debugger import Tracer
import os
import clarans_vik_experiments
reload(clarans_vik_experiments)



savedir = hardpaths.clarans_paper_dir


def get_experiment_results_from_file(experiment_name):
  """
  return a dictionary with keys : [algorithm_level][run_number]['E' or 'ttime']
  """
  results_dir = os.path.join(clarans_vik_experiments.clarans_vik_results_base, experiment_name)
  results = {}
  for alg_level in os.listdir(results_dir):
    results[alg_level] = {}
    dir_alg_level = os.path.join(results_dir, alg_level)
    for runx in os.listdir(dir_alg_level):
      run_number = int(runx.split("run")[-1])
      results[alg_level][run_number] = {'ttime':[], 'E':[]}
      filly = open(os.path.join(dir_alg_level, runx, "results.txt"), "r")
      lines = filly.readlines()
      filly.close()
      for l in lines:
        t, E = l.split()
        t = float(t)
        E = float(E)
        results[alg_level][run_number]['ttime'].append(t)
        results[alg_level][run_number]['E'].append(E)
  return results



def plot_experiment(experiment_name, relative = False):

  results = get_experiment_results_from_file(experiment_name)

  if relative == True:
    nruns = max([kp for kp in results.values()[0].keys()]) + 1
    min_E = 10**53
    for k in results.keys(): # names of algorithms used
      for run in range(nruns):
        min_E = min(min_E, min(results[k][run]['E']))
    

  nruns = max([kp for kp in results.values()[0].keys()]) + 1
  for run in range(nruns):   
    
    for k in results.keys(): # names of algorithms used
      print "--> ", k, run
      #energy_time = get_energy_time(results[k]['output'])

      label = None
      marker = None
      markersize = 1
      if k == 'clarans_3':
        label = 'cl(+)'
        linestyle = '-.'
        linewidth = 2
        color = 'r'
      elif k == 'clarans_0':
        label = 'cl(-)'
        linestyle = ':'
        linewidth = 2
        color = 'b'
      elif k == 'voronoi_0':
        label = 'vik'
        linestyle = '-'
        linewidth = 1
        marker = '.'
        markersize = 4
        color = 'g'
      else:
        raise RuntimeError("Unrecognised label (algorithm_level?), " + k)
      
      
      
      if run != 0:
        label = None
      
      if relative == True:
        results[k][run]['E'] = [x/min_E for x in results[k][run]['E']]
      
      pl.plot(np.array(results[k][run]['ttime']), results[k][run]['E'], label = label, linestyle = linestyle, linewidth = linewidth, color = color, alpha = 0.8, marker = marker, markersize = markersize)
  
  pl.xscale('log', basex = 2)

  #Tracer()()
  





def plot_sim_experiments(maxtime = 64.0): 
  """
  plot saved sim results
  """

  #pl.figure(156, figsize = (6,6))
  #pl.clf()
  #ax = pl.gca()  
  
  bottom_left = (0.05,0.05)
  #pl.close('all')  
  pl.figure(101, figsize = (6,6))
  pl.clf()
  ax = pl.gca()
  def_xlim = [0.01, maxtime]
  
  pl.subplots_adjust(hspace = 0.08)
  
  pl.subplot(2,2,1)
  ax = pl.gca()
  plot_experiment("binary_strings")
  pl.xlim(def_xlim)
  pl.xticks([])
  pl.ylabel('energy')
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "syn-1", transform=ax.transAxes)
  yticks = pl.yticks()
  
  pl.subplot(2,2,2)
  ax = pl.gca()
  plot_experiment("sparse_vectors")
  pl.xlim(def_xlim)
  pl.xticks([])
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "syn-2", transform=ax.transAxes)
  
  pl.subplot(2,2,3)
  ax = pl.gca()  
  plot_experiment("dense_vector_1")
  pl.xlim(def_xlim)
  pl.ylabel('energy')
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "syn-3", transform=ax.transAxes)
  pl.xlabel('time [s]')
  pl.xticks([2**-5, 2**-2, 2, 2**4])
  
  
  pl.subplot(2,2,4)
  ax = pl.gca()  
  plot_experiment("dense_vector_2")
  pl.xlim(def_xlim)
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "syn-4", transform=ax.transAxes)
  pl.xlabel('time [s]')
  pl.yticks([0, 0.1,0.2, 0.3, 0.4])
  pl.xticks([2**-5, 2**-2, 2, 2**4])
  
  pl.subplots_adjust(wspace = 0.18, bottom = 0.2)
  
  pl.legend(frameon = False)
  fname = os.path.join(os.path.expanduser('~'), savedir, 'vik_clarans_sim.pdf')
  pl.savefig(fname)
  
  import commands
  commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))



def get_yticks(yticks):
  vs, ts = yticks
  #new_vs, new_ts = [],[]
  #for i in range(len(vs)):
    ##Tracer()()
    #if ts[i].get_text()[-1] == 0:
      #new_vs.append(vs[i])
      #new_ts.append(vs[i])

  return vs, ts

def plot_real_experiments(maxtime = 2000.0): 
  """
  plot saved sim results
  """
  
  bottom_left = (0.05,0.05)
  pl.close('all')  
  pl.figure(102, figsize = (6,6))
  pl.clf()
  ax = pl.gca()
  def_xlim_left = [1., maxtime]
  def_xlim_right = [1.1*2**6, maxtime]
  pl.subplots_adjust(hspace = 0.08)

  
  pl.subplot(2,2,1)
  ax = pl.gca()
  plot_experiment("rcv1_parallel", relative = True)
  pl.xlim(def_xlim_left)
  pl.xticks([])
  pl.ylabel('relative energy')
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "rcv1", transform=ax.transAxes)
  pl.ylim(ymin = 0.99*pl.ylim()[0])
  pl.yticks([1, 1.1, 1.2, 1.3])
  
  pl.subplot(2,2,2)
  ax = pl.gca()
  plot_experiment("genome_parallel", relative = True)
  pl.xlim(def_xlim_right)
  pl.xticks([])
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "genome", transform=ax.transAxes)
  pl.ylim(ymin = 0.99*pl.ylim()[0])
  pl.yticks([1, 1.1])
  pl.xlim([2**6, 2**10.4])
    
  pl.subplot(2,2,3)
  ax = pl.gca()  
  plot_experiment("mnist_parallel", relative = True)
  pl.xlim(def_xlim_left)
  pl.ylabel('relative energy')
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "mnist", transform=ax.transAxes)
  pl.xlabel('time [s]')
  pl.xticks([2, 2**4, 2**7, 2**10])
  pl.ylim(ymin = 0.99*pl.ylim()[0])
  pl.legend(loc = 'upper right', handletextpad = 0.1, frameon = False)
  pl.yticks([1, 1.1, 1.2])
    
  pl.subplot(2,2,4)
  ax = pl.gca()  
  plot_experiment("english_words_parallel", relative = True)
  pl.xlim(def_xlim_right)
  pl.text(x = bottom_left[0], y = bottom_left[1], s = "words", transform=ax.transAxes)
  pl.xlabel('time [s]')
  pl.ylim(ymin = 0.99*pl.ylim()[0])
  pl.yticks([1, 1.1, 1.2])
  pl.xlim([2**6, 2**10.4])
  
  
  pl.subplots_adjust(wspace = 0.18, bottom = 0.2)
  fname = os.path.join(os.path.expanduser('~'), savedir, 'vik_clarans_real.pdf')
  pl.savefig(fname)
  
  
  import commands
  commands.getstatusoutput('pdfcrop %s %s'%(fname, fname))
