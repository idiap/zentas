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

import matplotlib.pyplot as pl


# make data. 

npr.seed(1014)
K = 17**2
n_per_row = int(np.sqrt(K))
ndata = K*28
X = np.zeros((ndata, 2))
X[:,0] = np.repeat(np.arange(n_per_row), ndata/n_per_row)
X[:,1] = np.tile(np.arange(n_per_row), ndata/n_per_row)  
sigma = 0.1
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
text_x_eps = 0.02

pl.figure(0, figsize = (9, 11))
pl.clf()

col_kmeans_points = "red"
col_init_points = "green"
col_points = "k"

X_on_grid = get_on_grid(X, [[0, n_per_row], [0, n_per_row]], [[1.5,2.5 + width], [3,4 + height]])
pl.plot(X_on_grid[:,0], X_on_grid[:,1], linestyle = "none", marker = '.', color = col_points, markersize = 1)

def jam(COOD, col, alg, kmeans):
  C_on_grid = get_on_grid(results[alg][kmeans]["C"], [[0, n_per_row], [0, n_per_row]], [[COOD[0], COOD[0]+ width], [COOD[1],COOD[1] + height]])
  pl.plot(C_on_grid[:,0], C_on_grid[:,1], linestyle = "none", marker = '.', color = col, markersize = 2)
  pl.text(COOD[0] -0.01, COOD[1] - 0.1, "%.3f"%(results[alg][kmeans]["E"]))

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

#UNIFORM AT RANDOM
pl.plot([UNI00[0] + width/3., UNI00[0] + width/3.], [UNI00[1] + height, 3 - sigma/2], color = 'k')
pl.text(x = text_x_eps + UNI00[0] + width/3., y = 3 - (1 - height)/2, s = "uniform", verticalalignment = "top")

#KMEANS++
pl.plot([KMPP00[0] + width/3., KMPP00[0] + width/3.], [KMPP00[1] + height, 3 - sigma/2], color = 'k')
pl.text(x = text_x_eps + KMPP00[0] + width/3., y = 3 - (1 - height)/2, s = "$K$-means++", verticalalignment = "top")


#MEDLLOYD
pl.plot([UNI00[0], VOR00[0] + width/2.], [UNI00[1] + height/2., UNI00[1] + height/2.], color = 'k')
pl.plot([VOR00[0] + width/2., VOR00[0] + width/2.], [UNI00[1] + height/2., VOR00[1] + height], color = 'k')
pl.text(x = text_x_eps + VOR00[0] + width/2., y = 1 + height + (1 - height)/2, s = "MEDLLOYD", verticalalignment = "top")


#CLARANS
pl.plot([UNI00[0] + 2.*width/3, UNI00[0] + 2.*width/3], [UNI00[1], VOR00[1] + height], color = 'k')
pl.text(x = text_x_eps + UNI00[0] + 2.*width/3, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")

pl.plot([KMPP00[0] + 2.*width/3, KMPP00[0] + 2.*width/3], [UNI00[1], VOR00[1] + height], color = 'k')
pl.text(x = text_x_eps + KMPP00[0] + 2.*width/3, y = 1 + height + (1 - height)/2, s = "CLARANS", verticalalignment = "top")

#LLOYDS
pl.plot([UNILL00[0] + 2.*width/3, UNILL00[0] + 2.*width/3], [2 + 2.*height/3, 0 + height], color = 'k')
pl.plot([UNILL00[0] + 2.*width/3, UNI00[0]], [2 + 2.*height/3, 2 + 2.*height/3], color = 'k')
pl.text(x = text_x_eps + UNILL00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

pl.plot([VOR00[0] + 2.*width/3, VOR00[0] + 2.*width/3], [1, 0 + height], color = 'k')
pl.text(x = text_x_eps + VOR00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

pl.plot([UNICLA00[0] + 2.*width/3, UNICLA00[0] + 2.*width/3], [1, 0 + height], color = 'k')
pl.text(x = text_x_eps + UNICLA00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

pl.plot([UNIKM00[0] + 2.*width/3, UNIKM00[0] + 2.*width/3], [1, 0 + height], color = 'k')
pl.text(x = text_x_eps + UNIKM00[0] + 2.*width/3, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

pl.plot([KMPPLL00[0] + 1.*width/5, KMPP00[0] + width], [2 + 2.*height/3, 2 + 2.*height/3], color = 'k')
pl.plot([KMPPLL00[0] + 1.*width/5, KMPPLL00[0] + 1.*width/5], [2 + 2.*height/3, 0 + height], color = 'k')
pl.text(x = text_x_eps + KMPPLL00[0] + 1.*width/5, y = height + (1 - height)/2, s = "$K$-means", verticalalignment = "top")

ax = pl.gca()
ax.axis('off')


import commands
fn = pl.savefig(datapaths.datapaths["nipsflow"])
commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))

