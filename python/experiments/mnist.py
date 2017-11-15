# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

"""

Making and loading all things MNIST (without labels) for clustering tasks.

* Standard mnist dataset
* Enlarged mnist datasets
* 00 -> 99 mnist datasets (concatenated)
* Random projections of all above (projected) (TODO for concatenated still)

"""


import os, struct
from array import array
from cvxopt.base import matrix
import numpy as np
import socket
import commands

from IPython.core.debugger import Tracer
import numpy as np
import numpy.random as npr

import datapaths
reload(datapaths)


paths = {
'original' : os.path.join(datapaths.datapaths['infidpath'], 'generateddata/original'), 
'projected' : os.path.join(datapaths.datapaths['infidpath'], 'generateddata/projected'),
'concatenated' : os.path.join(datapaths.datapaths['infidpath'], 'generateddata/concatenated')}




def make_MNIST(dataset = "original", ndata = 100000, dimension = None):

  for p in paths.values():
    if not os.path.exists(p):
      print "making dir (if not exist) ", p
      os.makedirs(p)

  print "in make MNIST ( ", dataset, ", ", ndata, ", ", dimension, " ) "
  if dataset == "original":
    cwd = os.getcwd()
    os.chdir(datapaths.datapaths['infiexec'])
    
    
    datafn = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
    command = './infimnist pat 0 %d > %s'%(ndata, datafn)
    print "will execute command `" + command + "`"
    print commands.getstatusoutput(command)
    os.chdir(cwd)
    
  if dataset == "projected":
    datafn = os.path.join(paths['projected'], 'projected_%d_%d.npy'%(ndata,dimension))
    datafn_preprojected = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
    if not os.path.exists(datafn_preprojected):
      make_MNIST(dataset = "original", ndata = ndata)
      X = read_MNIST("original", ndata)
      os.remove(datafn_preprojected)
    else:
      make_MNIST(dataset = "original", ndata = ndata)
      print "reading original"
      X = read_MNIST("original", ndata)
    
    npr.seed(1011) 
    projection_matrix = npr.randn(784, dimension)
    X_proj = np.dot(X, projection_matrix).astype(np.float32)
    np.save(arr = X_proj, file=  datafn)      
        
def get_two_digs(X, ndata):
  n_orig, dimension = X.shape
  twod = np.empty((ndata, 2*dimension), dtype = np.float32)
  twod_left = twod[:, 0:dimension]
  twod_right = twod[:, dimension:2*dimension]
  for i in range(n_orig):
    if i*n_orig < ndata:
      twod_left[i*n_orig:min((i + 1)*n_orig, ndata)] = X[i].reshape(1, -1)
      nright = min((i + 1)*n_orig, ndata) - i*n_orig
      twod_right[i*n_orig:i*n_orig + nright, :] = X[0:nright, :]

  #Tracer()()
  return twod
 
  
def read_MNIST(dataset = "original", ndata = 100000, dimension = None):
    """
    Python function for importing the MNIST data set.
    """

    if dataset == "original":
      fname_img = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
      if not os.path.exists(fname_img):
        make_MNIST("original", ndata)
      
      
      
      
      fimg = open(fname_img, 'rb')
      magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
      img = array("B", fimg.read())
      fimg.close()
  
      images =  matrix(0, (ndata, rows*cols))
      for i in range(ndata):
        images[i, :] = img[ i*rows*cols : (i+1)*rows*cols ]

    
      return np.array(images, dtype = np.float32)
      
    elif dataset == "projected":
      datafn = os.path.join(paths['projected'], 'projected_%d_%d.npy'%(ndata,dimension))
      if not os.path.exists(datafn):
        make_MNIST("projected", ndata, dimension)
      
      return np.load(datafn)
    
    elif not isinstance(dataset, str):
      raise RuntimeError("in read_MNIST, parameter dataset should a string")
    
    #concatenated  (never actually saved as such)
    else:
      source_ndata = int(np.sqrt(ndata)) + 1
      sourcefn = os.path.join(paths['original'], 'original-%d-ubyte'%(source_ndata,))
      if not os.path.exists(sourcefn):
        make_MNIST("original", source_ndata)
      
      X = read_MNIST(dataset = "original", ndata = source_ndata)
      return get_two_digs(X, ndata)
      

def write_concatenated_train_test_split(N_train = 200000, N_test = 10000):
  Ns = {"train":N_train, "test":N_test}
  
  X = read_MNIST(dataset = "concatenated", ndata = N_train + N_test, dimension = None)
  import random
  test_indices = random.sample(xrange(N_train + N_test), N_test)
  train_indices = [x for x in range(N_train + N_test) if x not in test_indices]
  X_split = {"test":X[test_indices,:], "train":X[train_indices,:]}
  
  rootdir = hardpaths.trainandtestdir
  for tau in ["train", "test"]:
    np.savetxt(fname = "%s/infimnist_%s_withdims.txt"%(rootdir, tau), X = X_split[tau], fmt = '%d', header = "%d\t%d"%(Ns[tau], 2*(28**2)), comments = '')
    np.savetxt(fname = "%s/infimnist_%s_dimless.txt"%(rootdir, tau), X = X_split[tau],  fmt = '%d')
    

def write_mnist28by28_train_test_split(N_train = 400000, N_test = 20000):
  Ns = {"train":N_train, "test":N_test}
  
  X = read_MNIST(dataset = "original", ndata = N_train + N_test, dimension = None)
  import random
  test_indices = random.sample(xrange(N_train + N_test), N_test)
  train_indices = [x for x in range(N_train + N_test) if x not in test_indices]
  X_split = {"test":X[test_indices,:], "train":X[train_indices,:]}
  
  rootdir = trainandtestdir
  for tau in ["train", "test"]:
    np.savetxt(fname = "%s/infimnist28by28_%s_withdims.txt"%(rootdir, tau), X = X_split[tau], fmt = '%d', header = "%d\t%d"%(Ns[tau], (28**2)), comments = '')
    np.savetxt(fname = "%s/infimnist28by28_%s_dimless.txt"%(rootdir, tau), X = X_split[tau],  fmt = '%d')
    
 
    
      

