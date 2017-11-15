# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

import matplotlib.pyplot as pl
import sys
import random
import numpy as np
import numpy.random as npr
import datapaths

#where is pyzentas.so ? Make sure this is correct.
sys.path.append("../../build/python")
import pyzentas
from IPython.core.debugger import Tracer 
import time

import matplotlib
import matplotlib.pyplot as pl



# Font solution, based on stack-overflow page
# https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
matplotlib.rcParams['text.latex.preamble']= r"\usepackage{cmbright} \renewcommand{\familydefault}{fos} \renewcommand{\seriesdefault}{l} \renewcommand{\bfdefault}{sb}"

#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : "sans-serif",
          'text.latex.unicode': True,
}
matplotlib.rcParams.update(params) 

# make data. 

npr.seed(1014)
K = 12**2
n_per_row = int(np.sqrt(K))
ndata = K*25
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


ncols = 10.
nrows = 6.
def get_on_grid(C, dims_before, dims_after):
  Ct = C.copy()
  for x in range(2):
    Ct[:, x] = C[:, x] - dims_before[x][0]
    Ct[:, x] *= (0. + dims_after[x][1] - dims_after[x][0]) / (dims_before[x][1]  - dims_before[x][0])
    Ct[:, x] += dims_after[x][0]
  
  return Ct

if False:    
  results = {}
  
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False, algorithm = "clarans", level = 0)
  results["uniform"] = {}
  for llo in [True, False]:
    results["uniform"][llo] = go(z,llo)
  
  results["voronoi"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 10, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False, algorithm = "voronoi", level = 0)
  for llo in [True, False]:
    results["voronoi"][llo] = go(z, llo)
  
  results["km_clarans"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 30, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "kmeans++-4", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    results["km_clarans"][llo] = go(z, llo)
  
  results["kmeans++"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 0, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "kmeans++-1", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    results["kmeans++"][llo] = go(z, llo)
  
  
  results["uni_clarans"] = {}
  z = pyzentas.pyzen(K = K, metric = 'l2', energy = 'quadratic', max_itok = 30, max_time = 100000, max_proposals = K**2, seed = 1011, patient = True, nthreads = 4, init = "uniform", with_tests = False, capture_output = True, rooted = False)
  for llo in [True, False]:
    results["uni_clarans"][llo] = go(z, llo)
  
width = 0.8
height = 0.6
eps = 0.2
text_x_eps = 0.06

pl.figure(0, figsize = (7.5,7.5))
pl.clf()

col_kmeans_points = "red"
col_init_points = "green"
col_points = "k"

X_on_grid = get_on_grid(X, [[0, n_per_row], [0, n_per_row]], [[2,2 + width], [3,3 + height]])
pl.plot(X_on_grid[:,0], X_on_grid[:,1], linestyle = "none", marker = '.', color = col_points, markersize = 1.3, alpha = 0.8)

def jam(COOD, col, alg, kmeans):
  C_on_grid = get_on_grid(results[alg][kmeans]["C"], [[0, n_per_row], [0, n_per_row]], [[COOD[0], COOD[0]+ width], [COOD[1],COOD[1] + height]])
  pl.plot(C_on_grid[:,0], C_on_grid[:,1], linestyle = "none", marker = '.', color = col, markersize = 1.4, alpha = 0.7)
  pl.text(COOD[0] -0.01, COOD[1] - 0.12, "$E=%.3f$"%(results[alg][kmeans]["E"]))

UNI00 = [1.5, 2]
jam(UNI00, col_init_points, "uniform", False)

KMPP00 = [2.5, 2]
jam(KMPP00, col_init_points, "kmeans++", False)

VOR00 = [1 - eps, 1]
jam(VOR00, col_init_points, "voronoi", False)

UNICLA00 = [2 - eps, 1]
jam(UNICLA00, col_init_points, "uni_clarans", False)

UNIKM00 = [3 - eps, 1]
jam(UNIKM00, col_init_points, "km_clarans", False)

UNILL00 = [0, 0]
jam(UNILL00, col_kmeans_points, "uniform", True)

jam([1,0], col_kmeans_points, "voronoi", True)
jam([2,0], col_kmeans_points, "uni_clarans", True)
jam([3,0], col_kmeans_points, "km_clarans", True)

KMPPLL00 = [4, 0]
jam(KMPPLL00, col_kmeans_points, "kmeans++", True)


#################### connecting lines ###########################

#KMEANS++
init_line_kwargs = {"color" : "#8b0000", 'linewidth':1.1, 'linestyle' : '-'}


#UNIFORM AT RANDOM
#pl.plot([UNI00[0] + 3*width/4., UNI00[0] + 3*width/4.], [UNI00[1] + height + 0.8*sigma, 3 - 0.8*sigma], **init_line_kwargs)

def jamarrow(x, y, dx, dy, **kwargs):
  pl.arrow(x, y, dx, dy, head_width = 0.04, head_length = 0.02, **kwargs)#, **init_line_kwargs)

  #, width = 0.008

def jamarrow2(X, Y, **kwargs):
  pl.plot([X[-3], X[-2]], [Y[-3], Y[-2]], **kwargs) 
  jamarrow(X[-2], Y[-2], X[-1] - X[-2], Y[-1] - Y[-2], **kwargs)
  
  
delta_c = 1*(height + 0.)/n_per_row

jamarrow(UNI00[0] + 3*width/4., 3. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
pl.text(x = -1.5*text_x_eps + UNI00[0] + 3*width/4., y = 3 - (1 - height)/2, s = "uniform", verticalalignment = "top", horizontalalignment = "right")

jamarrow(KMPP00[0] + width/4., 3. - 1*delta_c, 0, -(1 - height) + 2*delta_c, **init_line_kwargs)
pl.text(x = text_x_eps + KMPP00[0] + width/4., y = 3 - (1 - height)/2, s = "$K$-means++", verticalalignment = "top")


#MEDLLOYD
refinement_line_kwargs = {"color" : "#8b0000", 'linewidth':1.1, 'linestyle' : '-'}
jamarrow2(
[UNI00[0] - 3*delta_c, VOR00[0] + width/3., VOR00[0] + width/3.], 
[UNI00[1] + height/2., UNI00[1] + height/2., VOR00[1] + height + delta_c],
**refinement_line_kwargs)
pl.text(x = text_x_eps + VOR00[0] + width/3., y = 1 + height + (1 - height)/2, s = "MEDLLOYD", verticalalignment = "top")


#CLARANS
jamarrow(UNI00[0] + 4.*width/5, UNI00[1] - 3*delta_c, 0, -(1 - height) + 4*delta_c, **refinement_line_kwargs)
pl.text(x = text_x_eps + UNI00[0] + 4.*width/5, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")

jamarrow(KMPP00[0] + 2.*width/3, UNI00[1] - 3*delta_c, 0, -(1 - height) + 4*delta_c, **refinement_line_kwargs)
pl.text(x = text_x_eps + KMPP00[0] + 2.*width/3, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")

#LLOYDS
lloyds_line_kwargs = {"color" : "#8b0000", 'linewidth':1.1}
jamarrow2([UNI00[0] - 3*delta_c, UNILL00[0] + 2.*width/3, UNILL00[0] + 2.*width/3], [2 + 2.*height/3, 2 + 2.*height/3, height + delta_c], **lloyds_line_kwargs)

pl.text(x = text_x_eps + UNILL00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

jamarrow(VOR00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
pl.text(x = text_x_eps + VOR00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

jamarrow(UNICLA00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
pl.text(x = text_x_eps + UNICLA00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

jamarrow(UNIKM00[0] + 2.*width/3, 1  - 3*delta_c, 0, -(1 - height) + 4*delta_c, **lloyds_line_kwargs)
pl.text(x = text_x_eps + UNIKM00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

jamarrow2(
[KMPP00[0] + width + delta_c, KMPPLL00[0] + 1.*width/5, KMPPLL00[0] + 1.*width/5], 
[2 + 2.*height/3, 2 + 2.*height/3, 0 + height  + delta_c],
**lloyds_line_kwargs)


#pl.plot([KMPPLL00[0] + 1.*width/5, KMPP00[0] + width], [2 + 2.*height/3, 2 + 2.*height/3], **lloyds_line_kwargs)
#pl.plot([KMPPLL00[0] + 1.*width/5, KMPPLL00[0] + 1.*width/5], [2 + 2.*height/3, 0 + height], **lloyds_line_kwargs)
pl.text(x = text_x_eps + KMPPLL00[0] + 1.*width/5, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

ax = pl.gca()
ax.axis('off')


import commands
fn = datapaths.datapaths["nipsflow"] 
pl.savefig(fn)
print commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))

