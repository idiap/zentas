import os
import numpy as np
import sys
import datapaths
reload(datapaths)

def get_rna():
  """
  Raw data at:
  https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#cod-rna
  
  (put it in datapaths["rnaraw"])
  
  """
  print datapaths.datapaths.keys()
  writepath = os.path.join(datapaths.datapaths["rnawrite"], "rna.npy")
  if os.path.exists(writepath):
    print "loading npy file"
    X = np.load(writepath)
    print "data of shape ", X.shape
    return X
    
  else:
    print "loading from text files"
    allvals = []
    
    n_goodfiles = 0
    for fn in os.listdir(datapaths.datapaths["rnaraw"]):
      if fn in ["cod-rna",  "cod-rna.r" ,  "cod-rna.t"]:
        n_goodfiles += 1
    
    if n_goodfiles != 3:
      raise RuntimeError("Expected files cod-rna, cod-rna.r and cod-rna, in rnaraw dir")
        
    
    for fn in os.listdir(datapaths.datapaths["rnaraw"]):
      if fn not in ["cod-rna",  "cod-rna.r" ,  "cod-rna.t"]:
        pass
      else:
        filly = open(os.path.join(datapaths.datapaths["rnaraw"], fn))
        lines = filly.readlines()
        for l in lines:
          vals = []
          for x in l.split(":")[1::]:
            vals.append(float(x.split()[0]))
          allvals.append(vals)
    
    X = np.array(allvals)
    if not os.path.exists(datapaths.datapaths["rnawrite"]):
      os.makedirs(datapaths.datapaths["rnawrite"])


    filly = open(writepath, "w")
      
    np.save(filly, X)
    return X
  #print lines[0], lines[-1]
