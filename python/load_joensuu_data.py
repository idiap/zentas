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
import socket
import os
import sys

if (socket.gethostname() == 'idbean'):
  raise RuntimeError("Need to get this data...")

else:
  data_base_dir = hardpaths.joensuudata

def load_data(name = "KDDCUP04Bio"):
  txt_filename = os.path.join(data_base_dir, "txts", "%s.txt"%(name,))
  filly = open(txt_filename, "r")
  line1 = filly.readline()
  filly.close()
  
  if len(line1.split() ) == 2:
    nrows_skip = 1
  else:
    nrows_skip = 0
  
  print "loading data..."
  X  = np.loadtxt(txt_filename, skiprows = nrows_skip)
  
  return X
