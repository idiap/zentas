# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

import sys
import os
sys.path.append("../../build/python")
import pyzentas
import numpy as np
import numpy.random as npr
import random 

import matplotlib 
import matplotlib.pyplot as pl
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


always_run_from_scratch = False

seed = 1011
ndata = int(0.5*10**6)
K = int(500)
dimension = 4
data = np.array(npr.randn(ndata, dimension), dtype = np.float32)
levels = [0,2,3]

def get_results():
  """
  cluster dense data points.
  """

  if always_run_from_scratch == False and os.path.exists(datapaths.datapaths["comparing_levels_fn"]):
    X = np.load(datapaths.datapaths["comparing_levels_fn"])
    results = {}
    for l in levels:
      results[l] = X[()][l]
  
  else:
    results = {}
    for level in levels:
      z = pyzentas.pyzen(init = "uniform", K = K, metric = 'l2', energy = 'quadratic', exponent_coeff = 0,  max_rounds = 20000, max_time = 60, seed = seed, nthreads = 4, patient = True, with_tests = False, algorithm = "clarans", level = level, capture_output = True)
    
      tangerine = z.den(data, do_vdimap = False, do_refinement = False)
      results[level] = pyzentas.get_processed_output(tangerine['output'])

  
  return results
  
if "results" not in locals() or always_run_from_scratch:
  results = get_results()

pl.figure(figsize = (4.5, 3.5))


colors = ["#900C3F", "#C70039", "#33FF93"] #, "#FF5733"
linestyles = ["-", "-", "-"]
linewidths = [1.3, 1.3, 1.5] #1.3, 

max_time = 10**20
max_ncalcs = 10**20
for li, l in enumerate(levels):
  
  label_map = {0:1, 2:2, 3:3}
  times = results[l]['Tt']/1000.
  ncalcs = 2**results[l]['lg2nc']
  
  max_time = min(times[-1], max_time)
  max_ncalcs = min(ncalcs[-1], max_ncalcs)
  
  pl.subplot(2,1,1)
  pl.plot(times, results[l]['mE'], label = label_map[l], color = colors[li], linestyle = linestyles[li], linewidth = linewidths[li], alpha = 0.5)
  
  pl.subplot(2,1,2)
  pl.plot(ncalcs, results[l]['mE'], label = None, color = colors[li], linestyle = linestyles[li], linewidth = linewidths[li], alpha = 0.5)  
  
  #pl.plot(results[l]['lg2nc'], results[l]['mE'], label = 1+l, color = colors[li], linestyle = linestyles[li])

pl.subplot(2,1,1)
pl.xlim(xmax = max_time*1.01)
pl.legend(labelspacing = 0.001, frameon = False, ncol = 4, columnspacing = 0.5, handletextpad = 0.01,borderaxespad = 0.01)
#, , borderpad = 0.04
pl.xlabel("time [s]")
pl.ylabel("$E$")
#pl.ylim(ymax = 0.44)

yticks = np.array(pl.yticks()[0])
yticks = yticks[0::2]
pl.yticks(yticks)


pl.subplot(2,1,2)
pl.xlim(xmax = max_ncalcs*1.01)
pl.xlabel("number of distance calculations")
pl.ylabel("$E$")
#pl.ylim(ymax = 0.44)

#pl.xticks([200000, 400000, 600000, 800000])
yticks = np.array(pl.yticks()[0])
yticks = yticks[0::2]
pl.yticks(yticks)

pl.subplots_adjust(bottom = 0.25, hspace = 0.7)

fn = datapaths.datapaths["comparing_levels_fig_fn"]
pl.savefig(fn)
import commands
commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))

#dense_data_example()

