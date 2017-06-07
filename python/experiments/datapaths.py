import socket

#I set the paths to data here, but you need to download data separately 
if socket.gethostname() == "goudurix12":
  import goudurix12paths
  reload(goudurix12paths)
  datapaths = goudurix12paths.datapaths

else:
  print "unknown host in datapaths.py, certain data paths may need to be set"

  tobeset = "/this/path/needs/to/be/set/in/datapaths.py"
  datapaths = {}
  
  # path to bin from http://leon.bottou.org/projects/infimnist
  # at http://leon.bottou.org/_media/projects/infimnist.tar.gz
  datapaths["infiexec"] = tobeset
  
  # path to where mnist data can be written
  datapaths["infidpath"] = tobeset
  
