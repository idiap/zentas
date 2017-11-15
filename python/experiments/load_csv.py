# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>

import os
import numpy as np
import sys
import datapaths
reload(datapaths)

def get_csv_from_raw(dataset):
  txtfn = datapaths.datapath_filenames["%s_fn"%(dataset,)]
  fullfn = os.path.join(datapaths.datapaths["%s_raw"%(dataset,)], txtfn)

  if dataset == "epileptic":
    filly = open(fullfn, "r")
    filly.readline()
    n_per_line = len(filly.readline().split(","))
    print n_per_line
    filly.close()
    X = np.loadtxt(fullfn, delimiter=",", skiprows = 1, usecols = range(1, n_per_line))

  elif dataset == "mopac":
    filly = open(fullfn, "r")
    lines = filly.readlines()
    data = []
    for l in lines[1::]:
      data.append([float(x) for x in l.replace("?", "0").split(",")])
    X = np.array(data)
    filly.close()
    
  else:
    X = np.loadtxt(fullfn, delimiter=",")


  return X

def get_csv(dataset):
  """
  """  
  txtfn = datapaths.datapath_filenames["%s_fn"%(dataset,)]
  writepath = os.path.join(datapaths.datapaths["%s_write"%(dataset,)], "%s.npy"%(dataset,))
  if os.path.exists(writepath):
    print "loading npy file"
    X = np.load(writepath)
    print "data of shape ", X.shape
    return X
    
  else:
    print "loading from text file"
    
    if txtfn not in os.listdir(datapaths.datapaths["%s_raw"%(dataset,)]):
      raise RuntimeError("Expected file %s, not found, needs to be downloaded perhaps"%(txtfn,))

    # All but zeroth column.
    X = get_csv_from_raw(dataset)
    if not os.path.exists(datapaths.datapaths["%s_write"%(dataset,)]):
      os.makedirs(datapaths.datapaths["%s_write"%(dataset,)])

    filly = open(writepath, "w")
    np.save(filly, X)
    return X


def get_htru2():
  return get_csv("htru2")

def get_yearpredictionmsd():
  return get_csv("yearpredictionmsd")

def get_epileptic():
  return get_csv("epileptic")

def get_mopac():
  return get_csv("mopac")
