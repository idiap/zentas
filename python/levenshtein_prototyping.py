# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas.
# zentas is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License version 3 as published by the Free Software Foundation.
# zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details. You should have received
# a copy of the GNU General Public License along with zentas. If not, see
# <http://www.gnu.org/licenses/>.

import hardpaths
reload(hardpaths)

import numpy as np
import numpy.random as npr

def levenshtein(v_vertical = np.array([1,1,1,1,1,0,1,1,1,1,1]), v_horizontal = np.array([1,1,1,1,1,1,1,1,1,1]), c_indel = 1.1, c_switch = 1.5, threshold = 3):  
  nrows = v_vertical.size
  ncols = v_horizontal.size


  threshold = min(threshold, c_indel*(nrows + ncols));

  
  #the number of positions next to the center which might be above threshold. 
  #TODO : specialise for threshold = infty.
  half_width = int(np.floor(threshold / c_indel))
   # if half_width = 1, then the values will be [threshold, x0, x1, x2, threshold]. The center is at 1 + half_width
  width = 3 + 2*(half_width)
  #A_prev = threshold*np.ones(width)
  #A_acti = threshold*np.ones(width)
  #A_acti[1+half_width::] = c_indel*np.arange(half_width +2)

  A_acti = np.abs(np.arange(width) - (half_width + 1))
  A_prev = np.abs(np.arange(width) - (half_width + 1))
    
  
  print "threshold : ", threshold
  
  #For debugging:
  SCORES = threshold*np.ones((nrows + 1, ncols + 1))

  row = 0
  min_distance = 0
  while row < nrows and min_distance < threshold:
    A_acti, A_prev = A_prev, A_acti
    row += 1
    
    min_distance = threshold
    for j in range(max(1, 1 + half_width - row), min(ncols - row + 2 + half_width, width-1)):      
      
      ########################################    
      #         equivalent form              #
      #for j in range(1, width-1):           #
      #  column = row + j - (1 + half_width) #
      #  if column >= 0 and column <= ncols: #
      ########################################

      column = row + j - (1 + half_width)
      A_acti[j] = min(
      c_indel + A_acti[j - 1], 
      c_indel + A_prev[j + 1], 
      A_prev[j] + c_switch*(v_vertical[row - 1] != v_horizontal[column - 1])
      )
      
      #print A_acti[j], " ",
      
      min_distance = min(min_distance, A_acti[j])
      #For debugging:
      SCORES[row, column] = min(threshold, A_acti[j])
    
    #print " "
    #print "---> ", min_distance
  #For debugging:  
  
  print SCORES[1::]
  
  
  #if min_distance < threshold:
  if np.abs(nrows - ncols) <= half_width + 1:
    print min(threshold, A_acti[1 + half_width + ncols - nrows])
  
  else:
    print threshold
  #else:
    #print threshold  
      


def levenshtein_stupid(v_vertical = np.array([1,1,1,1,1,0,1,1,1,1,1]), v_horizontal = np.array([1,1,1,1,1,1,1,1,1,1]), c_indel = 1.1, c_switch = 1.5, threshold = 3):
  nrows = v_vertical.size
  ncols = v_horizontal.size
  scores = np.zeros((nrows + 1, ncols + 1))
  scores[0,:] = c_indel*np.arange(ncols+1)
  scores[:,0] = c_indel*np.arange(nrows+1)
  for r in range(1, nrows + 1):
    for c in range(1, ncols + 1):
      scores[r,c] = min(c_indel + min(scores[r,c-1], scores[r-1,c]), scores[r-1, c-1] + c_switch*(v_vertical[r-1] != v_horizontal[c-1]))
  
  scores = scores*(scores < threshold) + threshold*(scores >= threshold) 
  print scores[1::,:]
  #print scores[-1,:]
    
def test(v_vertical = np.array([1,1,1,1,1,0,1,1,1,1,1], dtype = np.int), v_horizontal = np.array([2,2,2,2,1,1,1,1,1,1,1,1,1,1], dtype = np.int), c_indel = 1.1, c_switch = 1.5, threshold = 3):
  print v_vertical.reshape(-1,1)
  print v_horizontal
  levenshtein(v_vertical, v_horizontal, c_indel, c_switch, threshold)
  
  print "\n\n"
  levenshtein_stupid(v_vertical, v_horizontal, c_indel, c_switch, threshold)
 


def randomtest():
  seed = npr.randint(1000)
  npr.seed(seed)
  print seed
  #test(v_vertical = npr.randint(3, size = 8), v_horizontal = npr.randint(3, size = 12),    c_indel = 0.8, c_switch = 0.6, threshold = 4.7)
  v_vertical = np.ones(8)
  v_horizontal = np.ones(4)
  #v_horizontal[7::] = 0
  test(v_vertical, v_horizontal ,c_indel = 1.0, c_switch = 1.0, threshold = 30)

randomtest()
