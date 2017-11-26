# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

import matplotlib.pyplot as pl
import sys
import random
import numpy as np
import numpy.random as npr
import datapaths
import os

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../../build/python")
import pyzentas
from IPython.core.debugger import Tracer 
import time

import matplotlib
import matplotlib.pyplot as pl

always_run_from_scratch = False

# What label to use for Lloyd's algorithm,
LLOYDS = "LLOYD"

# Font solution, based on stack-overflow page
# https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
matplotlib.rcParams['text.latex.preamble']= r"\usepackage{cmbright} \renewcommand{\familydefault}{fos} \renewcommand{\seriesdefault}{l} \renewcommand{\bfdefault}{sb} "

# colors don't work 
# \usepackage{color} \usepackage{xcolor} :(
# https://stackoverflow.com/questions/24173250/colored-latex-labels-in-plots

#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : "sans-serif",
          'text.latex.unicode': True,
}

    
matplotlib.rcParams.update(params) 

arrow_color = "#8FD8D8"
init_line_kwargs = {"color" : arrow_color, 'linewidth':0.9, 'linestyle' : '-'}
width = 0.8
height = 0.6
eps = 0.2
text_x_eps = 0.06
col_kmeans_points = "red"
col_init_points = "green"
col_points = "k"

# make data. 
K = 12**2
ndata = K*25
n_per_row = int(np.sqrt(K))

delta_c = 1*(height + 0.)/n_per_row

if "X" not in locals().keys() or always_run_from_scratch:
  npr.seed(1014)
  X = np.zeros((ndata, 2))
  X[:,0] = np.repeat(np.arange(n_per_row), ndata/n_per_row)
  X[:,1] = np.tile(np.arange(n_per_row), ndata/n_per_row)  
  sigma = 0.13
  X += sigma*npr.randn(ndata, 2)


def go(z, do_rf):
    
  tangerine =  z.den(X = X,  do_vdimap = False, do_refinement = do_rf, rf_alg = "exponion")
    
  if do_rf == False:
    C = X[tangerine['indices_final']]
  else:
    C = pyzentas.get_vcenters(X, K, tangerine["labels"])


  energy = None
  lines = tangerine["output"].split("\n")
  for l in lines:
    if "mE=" in l:
      energy = float(l.split("mE=")[1].split()[0])

  if not energy:
    raise RuntimeError("no energy detected")
  print energy  
  return {"C":C, "d" : tangerine, "E":energy}


def get_on_grid(C, dims_before, dims_after):
  Ct = C.copy()
  for x in range(2):
    Ct[:, x] = C[:, x] - dims_before[x][0]
    Ct[:, x] *= (0. + dims_after[x][1] - dims_after[x][0]) / (dims_before[x][1]  - dims_before[x][0])
    Ct[:, x] += dims_after[x][0]
  
  return Ct

def get_g_results():

  #TODO :  a cache check, if cached then just load. 
  g_results = {}
  
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False, algorithm = "clarans", level = 0)
  g_results["uniform"] = {}
  for llo in [True, False]:
    g_results["uniform"][llo] = go(z,llo)
  
  g_results["voronoi"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 10, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False, algorithm = "voronoi", level = 0)
  for llo in [True, False]:
    g_results["voronoi"][llo] = go(z, llo)
  
  g_results["km_clarans"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 30, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "kmeans++-4", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    g_results["km_clarans"][llo] = go(z, llo)
  
  g_results["kmeans++"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "kmeans++-1", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    g_results["kmeans++"][llo] = go(z, llo)
  
  
  g_results["uni_clarans"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 30, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    g_results["uni_clarans"][llo] = go(z, llo)

  return g_results

if "g_results" not in locals().keys() or always_run_from_scratch:
  g_results = get_g_results()



def plot_X():
  X_on_grid = get_on_grid(X, [[0, n_per_row], [0, n_per_row]], [[2.02,2.02 + width], [3,3 + height]])
  pl.plot(X_on_grid[:,0], X_on_grid[:,1], linestyle = "none", marker = '.', color = col_points, markersize = 0.6, alpha = 0.8)



def jam(COOD, col, alg, kmeans, with_E = True, ghost = False):
  """
  """
  alpha = 0.7
  if ghost:
    col = 'k'
    alpha = 0.0
    
  C_on_grid = get_on_grid(g_results[alg][kmeans]["C"], [[0, n_per_row], [0, n_per_row]], [[COOD[0], COOD[0]+ width], [COOD[1],COOD[1] + height]])
  
  pl.plot(C_on_grid[:,0], C_on_grid[:,1], linestyle = "none", marker = '.', color = col, markersize = 2.5, alpha = alpha)
  
  if with_E:
    pl.text(COOD[0] -0.01, COOD[1] - 0.12, "$E=%.3f$"%(g_results[alg][kmeans]["E"]))


def jamarrow(x, y, dx, dy, **kwargs):
  print kwargs["linewidth"]
  pl.arrow(x, y, dx, dy, head_width = 0.04, head_length = 0.02, **kwargs)#, **init_line_kwargs)

def jamarrow2(X, Y, **kwargs):
  pl.plot([X[-3], X[-2]], [Y[-3], Y[-2]], **kwargs) 
  jamarrow(X[-2], Y[-2], X[-1] - X[-2], Y[-1] - Y[-2], **kwargs)
  

def plot_seeding(with_E):

  jam([1.5, 2], col_init_points, "uniform", False, with_E = with_E)

  jam([2.5, 2], col_init_points, "kmeans++", False, with_E = with_E)

  #arrow to uniform seeding
  jamarrow(1.5 + 3*width/4., 3. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
  pl.text(x = -0.9*text_x_eps + 1.5 + 3*width/4., y = 3 - (1 - height)/2, s = "uniform", verticalalignment = "top", horizontalalignment = "right")
  
  #arrow to k-means++ seeding
  jamarrow(2.5 + width/4., 3. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
  pl.text(x = text_x_eps + 2.5 + width/4., y = 3 - (1 - height)/2, s = "$K$-means++", verticalalignment = "top")


def poster_flow():

  with_E = True
  pl.figure(0, figsize = (7.5,7.5))
  pl.clf()
  
  plot_X()
  plot_seeding(with_E)
  
  VOR00 = [1 - eps, 1]
  jam(VOR00, col_init_points, "voronoi", False, with_E = with_E)
  
  UNICLA00 = [2 - eps, 1]
  jam(UNICLA00, col_init_points, "uni_clarans", False, with_E = with_E)
  
  UNIKM00 = [3 - eps, 1]
  jam(UNIKM00, col_init_points, "km_clarans", False, with_E = with_E)
  
  UNILL00 = [0, 0]
  jam(UNILL00, col_kmeans_points, "uniform", True, with_E = with_E)
  
  jam([1,0], col_kmeans_points, "voronoi", True, with_E = with_E)
  jam([2,0], col_kmeans_points, "uni_clarans", True, with_E = with_E)
  jam([3,0], col_kmeans_points, "km_clarans", True, with_E = with_E)
  
  KMPPLL00 = [4, 0]
  jam(KMPPLL00, col_kmeans_points, "kmeans++", True, with_E = with_E)
  
  
  #################### connecting lines ###########################
  
  #MEDLLOYD
  refinement_line_kwargs = {"color" : arrow_color, 'linewidth':0.9, 'linestyle' : '-'}
  jamarrow2(
  [1.5 - 3*delta_c, VOR00[0] + width/3., VOR00[0] + width/3.], 
  [2 + height/2., 2 + height/2., VOR00[1] + height + delta_c],
  **refinement_line_kwargs)
  pl.text(x = text_x_eps + VOR00[0] + width/3., y = 1 + height + (1 - height)/2, s = "MEDLLOYD", verticalalignment = "top")
  
  
  #CLARANS
  jamarrow(1.5 + 4.*width/5, 2 - 3*delta_c, 0, -(1 - height) + 4*delta_c, **refinement_line_kwargs)
  pl.text(x = text_x_eps + 1.5 + 4.*width/5, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")
  
  jamarrow(2.5 + 2.*width/3, 2 - 3*delta_c, 0, -(1 - height) + 4*delta_c, **refinement_line_kwargs)
  pl.text(x = text_x_eps + 2.5 + 2.*width/3, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")
  
  #LLOYDS
  lloyds_line_kwargs = {"color" : arrow_color, 'linewidth':0.9}
  jamarrow2([1.5 - 3*delta_c, UNILL00[0] + 2.*width/3, UNILL00[0] + 2.*width/3], [2 + 2.*height/3, 2 + 2.*height/3, height + delta_c], **lloyds_line_kwargs)
  
  pl.text(x = text_x_eps + UNILL00[0] + 2.*width/3, y = height + (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
  jamarrow(VOR00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
  pl.text(x = text_x_eps + VOR00[0] + 2.*width/3, y = height + (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
  jamarrow(UNICLA00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
  pl.text(x = text_x_eps + UNICLA00[0] + 2.*width/3, y = height + (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
  jamarrow(UNIKM00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
  pl.text(x = text_x_eps + UNIKM00[0] + 2.*width/3, y = height + (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
  jamarrow2(
  [2.5 + width + delta_c, KMPPLL00[0] + 1.*width/9, KMPPLL00[0] + 1.*width/9], 
  [2 + 2.*height/3, 2 + 2.*height/3, 0 + height  + delta_c],
  **lloyds_line_kwargs)
  
  pl.text(x = text_x_eps + KMPPLL00[0] + 1.*width/9, y = height + (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
  ax = pl.gca()
  ax.axis('off')
  
  
  import commands
  fn = datapaths.datapaths["nipsflow_poster"] 
  pl.savefig(fn)
  print commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))
  


def slide_flow(with_clarans = False):

  
  pl.figure(0, figsize = (8.3, 5.4))
  pl.clf()
  
  plot_X()
  plot_seeding(with_E = False)
  
  jam([1.5,1], col_kmeans_points, "uniform", True, with_E = True)  
  jam([2.5,1], col_kmeans_points, "kmeans++", True, with_E = True)


  pl.text(x = 2 + width + 0.1, y = 3.0 + height/4., s = "simulated data \n$K = 12^2, N = 25K$", horizontalalignment = "left")

  #arrow to uniform seeding
  jamarrow(1.5 + 3*width/4., 
          2. - 1*delta_c, 
          0, 
          -(1 - height) + 2*delta_c, **init_line_kwargs)
  pl.text(x = -0.9*text_x_eps + 1.5 + 3*width/4., y = 2 - (1 - height)/2, s = LLOYDS, verticalalignment = "top", horizontalalignment = "right")
  
  #arrow to k-means++ seeding
  jamarrow(2.5 + width/4., 2. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
  pl.text(x = text_x_eps + 2.5 + width/4., y = 2 - (1 - height)/2, s = LLOYDS, verticalalignment = "top")

  if (with_clarans == True):
    ghost = False
  
  else:
    ghost = True
    

  vertidown = 0
  horibar = 0.75
  
  if ghost == False:

    jam([0.5 - horibar,2 - vertidown], col_init_points, "uni_clarans", False, with_E = False, ghost = ghost)  
    jam([3.5 + horibar,2 - vertidown], col_init_points, "km_clarans", False, with_E = False, ghost = ghost)
  
    jam([0.5 - horibar,1 - vertidown], col_kmeans_points, "uni_clarans", True, with_E = True, ghost = ghost)  
    jam([3.5 + horibar,1 - vertidown], col_kmeans_points, "km_clarans", True, with_E = True, ghost = ghost)
  
  
    
    #CLARANS TO THE LEFT
    clarans_arrow_length = +(1 - width) - 2*delta_c + horibar
    jamarrow(1.5 - delta_c, 2.0 + height/2., -clarans_arrow_length,0, **init_line_kwargs)
    pl.text(x = 0.5 - horibar + width + 2*delta_c, y = 2 - vertidown + height/2., s = "CLARANS", verticalalignment = "bottom", horizontalalignment = "left")
  
    jamarrow(0.5 - horibar + width/2., 2. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
    pl.text(x = 0.5 - horibar + width/2. + delta_c, y = 2 - (1 - height)/2, s = LLOYDS, verticalalignment = "top")
  
    #CLARANS TO THE RIGHT
    jamarrow(2.5 + width + delta_c, 2.0 + height/2., clarans_arrow_length - delta_c,0, **init_line_kwargs)
    pl.text(x = 2.5 + width + delta_c, y = 2 - vertidown + height/2., s = "CLARANS", verticalalignment = "bottom", horizontalalignment = "left")
  
    jamarrow(2.5 + 3/2.*width + delta_c + horibar, 2. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
    pl.text(x = 2.5 + 3/2.*width + horibar, y = 2 - (1 - height)/2, s = LLOYDS, verticalalignment = "top", horizontalalignment = "right")

  x_bottom = 0.83
  pl.plot([-0.5, 5.2], [x_bottom, x_bottom], color = "0.1", alpha = 0.1)
  pl.plot([-0.5, 5.2], [3.6, 3.6], color = "0.1", alpha = 0.1)
  pl.plot([-0.5, -0.5], [x_bottom, 3.6], color = "0.1", alpha = 0.1)
  pl.plot([5.2, 5.2], [x_bottom, 3.6], color = "0.1", alpha = 0.1)
  #pl.xlim(-0.5, 5.2)
  #pl.ylim(0.9, 3.6)
  ax = pl.gca()
  ax.axis('off')

def slide_flow1():
  slide_flow(with_clarans = False)
  import commands
  fn = datapaths.datapaths["nipsflow_slide1"] 
  pl.savefig(fn, bbox_inches='tight')
  print "saving as ", fn
  print commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))


def slide_flow2():
  slide_flow(with_clarans = True)
  import commands
  fn = datapaths.datapaths["nipsflow_slide2"] 
  pl.savefig(fn, bbox_inches='tight')
  print "saving as ", fn
  print commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))



def impl_eval(with_stoppers = True):
  
  pl.figure(figsize = (6.5, 1.2))
  colors = {"uni_clarans":"r"} #, "km_clarans":"b"}
  
  for alg in colors.keys():
    lines = [l for l in g_results[alg][False]['d']['output'].split("\n")[0:-1] if "nprops=" in l]
    nprops = [int(l.split("nprops=")[-1]) - 3 for l in lines]
    pl.plot(nprops, marker = ".", linestyle = "none", label = alg.split("_")[0], color = "k", markersize = 2)
    pl.yscale("log", basey =10)
  
  pl.xlabel(r"implement (cumulative)", color = 'k')
  pl.ylabel(r"evaluate", color = 'k')
  
  
  pl.subplots_adjust(bottom = 0.33, left = 0.25)
  #pl.legend()
  pl.ylim(ymax = 0.8*(10**4))
  pl.xlim(xmax = 320)
  
  if with_stoppers:
    pl.plot([200, 200], [0, 10**4], linewidth = 2, color = "#119A51", alpha = 0.6)
    pl.plot([0, 320], [100, 100], linewidth = 2, color = "#8DD3AE", alpha = 0.6)
    pl.yticks([1, 10, 100, 1000], ["$10^0$", "$10^1$", "$R$", "$10^3$"])
    pl.xticks([0, 50, 100, 150, 200, 250, 300], [0, 50, 100, 150, "$S$",250, 300 ])
  
  fn = datapaths.datapaths["smld_impl_vs_eval_fn_ws%d"%(with_stoppers,)]
  pl.savefig(fn)
  import commands
  commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))



