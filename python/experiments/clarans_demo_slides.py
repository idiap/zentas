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




def clarans_demo():
  def get_at(radians):  
    X0 = np.array([0,0])
    Cs = np.array([X0 + 1.*np.array([np.cos(radians + 2*np.pi*j/3.), np.sin(radians +  + 2*np.pi*j/3.)]) for j in [0,1,2]])
    return Cs
    
  
  def demo_plot(iteration):
    
    X0 = 0.3*get_at(np.pi*2*(10./360)) + np.array([0,0])
    X1 = 0.3*get_at(np.pi*2*(30./360)) + np.array([3,0.8])
    X2 = 0.3*get_at(np.pi*2*(50./360)) + np.array([5,0])  
    
    points = np.vstack([X0, X1, X2])
    
    if iteration in [0,1,2,3,4,5,6]:
      center_indices = [1,2,4]
    
    else:
      center_indices = [1,4,8]
    
    centers = points[center_indices]
    
    X_lims = [points[:,0].min(), points[:,0].max()]
    Y_lims = [points[:,1].min(), points[:,1].max()]
  
    epso = 0.3  
    Xl_frame = [X_lims[0] - epso, X_lims[1] + epso]
    Yl_frame = [Y_lims[0] - epso, Y_lims[1] + epso]
    
    aspect = (X_lims[1] - X_lims[0] + 0.)/(Y_lims[1] - Y_lims[0])
    pl.figure(figsize = (10,1.5*10./aspect))
    ax = pl.gca()
  
    for i in [0,1]:
      pl.plot([Xl_frame[i], Xl_frame[i]], Yl_frame, color = "0.5", alpha = 0.5)
      pl.plot(Xl_frame, [Yl_frame[i], Yl_frame[i]], color = "0.5", alpha = 0.5)
  
  
    
    ncenter_kwargs = {"marker":".", "linestyle":"none", "markersize":7, "color":"k"}
    pl.plot(points[:,0], points[:,1], label = "non-center", **ncenter_kwargs) 
    
    center_kwargs = {"marker":".", "linestyle":"none", "markersize":20, "color":"g"}
    pl.plot(centers[:,0], centers[:,1], label = "center", **center_kwargs) 
    
    if iteration in [0,1]:
      pep_edge = 0.1
      pep_right = 0.3
      pep_down  = 0.2
      point_across = 0.25
      
      pl.plot([Xl_frame[0] + point_across], [Yl_frame[1] - pep_edge], **center_kwargs) 
      pl.text(Xl_frame[0] + pep_edge + pep_right, Yl_frame[1] - pep_edge, "center", fontsize = 17, verticalalignment = "center", horizontalalignment = "left")
      
      pl.plot([Xl_frame[0] + point_across], [Yl_frame[1] - pep_edge - pep_down], **ncenter_kwargs) 
      pl.text(Xl_frame[0] + pep_edge + pep_right, Yl_frame[1] - pep_edge - pep_down, "non-center", fontsize = 17, verticalalignment = "center", horizontalalignment = "left")
      
  
    if iteration in [0,1,4]:
      Lines = {2:[0], 4:[3,5,6,7,8]}
    
    elif iteration in [2,3,4]:
      Lines = {0:[3,4,5,6,7,8]}
    
    elif iteration in [5,6,7]:
      Lines = {1:[0,2], 4:[3,5], 8:[6,7]}
    
  
       
    for key in Lines.keys():
      center = points[key]
      
      for pindex in Lines[key]:
        point = points[pindex]
        eps = 0.16
        delta_x = center[0] - point[0]
        delta_y = center[1] - point[1]
        
        direction = np.array([delta_x, delta_y])
        direction /= np.sqrt(np.sum(direction*direction))
        
        Xs = [point[0] + eps*direction[0], center[0] - eps*direction[0]]
        Ys = [point[1] + eps*direction[1], center[1] - eps*direction[1]]
        pl.plot(Xs, Ys, color = "0.5")
    
    if iteration in [1,2,3]:
      sel_c = pl.Circle(points[4], radius = 0.12, facecolor = "none", edgecolor = "r")
      sel_nc = pl.Circle(points[0], radius = 0.12, facecolor = "none", edgecolor = "r")  
      ax.add_artist(sel_c)
      ax.add_artist(sel_nc)
      
    if iteration in [4,5,6]:
      sel_c = pl.Circle(points[2], radius = 0.12, facecolor = "none", edgecolor = "r")
      sel_nc = pl.Circle(points[8], radius = 0.12, facecolor = "none", edgecolor = "r")  
      ax.add_artist(sel_c)
      ax.add_artist(sel_nc)
    
    if iteration == 3:
      pl.text(0, 0.7, "reject the swap", fontsize = 17)
    
    if iteration == 6:
      pl.text(0, 0.7, "accept the swap", fontsize = 17)
    
    if iteration == 7:
      pl.text(0, 0.7, "implement the swap", fontsize = 17)
      
    #pl.ylim(ymax = Y_lims[1] + 0.2)
    ax.axis('off')
    
    
    
  for iteration in [0,1,2,3,4,5,6,7]:
    demo_plot(iteration)
    fn = os.path.join(datapaths.datapaths["smld_clarans_demo_dir"], "iter%d.pdf"%(iteration,))
    print fn
    pl.savefig(fn)
    import commands
    print commands.getstatusoutput("pdfcrop %s %s"%(fn, fn))
    
  
