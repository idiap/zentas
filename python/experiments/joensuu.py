import os
import numpy as np
import sys
import datapaths
reload(datapaths)

def get_joensuu(dataset):
  """
  Raw data at:
  https://cs.joensuu.fi/sipu/datasets/
  
  (put it in datapaths["joensuuraw"])
  
  """
  
  txtfn = "%s.txt"%(dataset,)
  writepath = os.path.join(datapaths.datapaths["joensuuwrite"], "%s.npy"%(dataset,))
  if os.path.exists(writepath):
    print "loading npy file"
    X = np.load(writepath)
    print "data of shape ", X.shape
    return X
    
  else:
    print "loading from text file"    
    if txtfn not in os.listdir(datapaths.datapaths["joensuuraw"]):
      raise RuntimeError("Expected file %s, not found, needs to be downloaded perhaps"%(txtfn,))
        
    X = np.loadtxt(os.path.join(datapaths.datapaths["joensuuraw"], txtfn))
    if not os.path.exists(datapaths.datapaths["joensuuwrite"]):
      os.makedirs(datapaths.datapaths["joensuuwrite"])

    filly = open(writepath, "w")      
    np.save(filly, X)
    return X
