import socket

#I set the paths to data here, but you need to download data separately 
if socket.gethostname() == "goudurix12":
  import goudurix12paths
  reload(goudurix12paths)
  datapaths = goudurix12paths.datapaths


elif socket.gethostname() == "goudurix11":
  import goudurix11paths
  reload(goudurix11paths)
  datapaths = goudurix11paths.datapaths

elif socket.gethostname() == "idbean":
  import idbeanpaths
  reload(idbeanpaths)
  datapaths = idbeanpaths.datapaths
  
else:
  print "unknown host in datapaths.py, certain data paths may need to be set"

  tobeset = "/this/path/needs/to/be/set/in/datapaths.py"
  datapaths = {}
  
  # path to bin from http://leon.bottou.org/projects/infimnist
  # at http://leon.bottou.org/_media/projects/infimnist.tar.gz
  datapaths["infiexec"] = tobeset
  
  # path to where mnist data can be written
  datapaths["infidpath"] = tobeset
  

  # path to where the files cod-rna  cod-rna.r  cod-rna.t
  # from website 
  # https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#cod-rna
  #are.
  datapaths["rnaraw"] = tobeset
  datapaths["rnawrite"] = tobeset
