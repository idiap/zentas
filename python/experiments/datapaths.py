import socket

datapath_filenames = {}
datapath_filenames["yearpredictionmsd_fn"] = "YearPredictionMSD.txt"
datapath_filenames["htru2_fn"] = "HTRU_2.csv"
datapath_filenames["epileptic_fn"] = "data.csv"
datapath_filenames["mopac_fn"] = "allUsers.lcl.csv"

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
  
  #path to where to save a figure
  datapaths["nipsflow"] = tobeset
  datapaths["nips_plot1"] = tobeset
  datapaths["nips_plot2"] = tobeset
  datapaths["nips_plot3"] = tobeset

  #path to eakmeans install dir
  datapaths["eaklibdir"] = tobeset
  
  #path to where to save a figure
  datapaths["kscalingfigpath"] = tobeset
  
  # path to where the files cod-rna  cod-rna.r  cod-rna.t
  # from website https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#cod-rna
  datapaths["rnaraw"] = tobeset
  datapaths["rnawrite"] = tobeset
  
  #path to where files a1.txt, europediff.txt etc are
  # from website https://cs.joensuu.fi/sipu/datasets/
  datapaths["joensuuwrite"] = tobeset
  datapaths["joensuuraw"] = tobeset

  datapaths["yearpredictionmsd_write"] = tobeset
  #path to where YearPredictionMSD.txt is  
  #from https://archive.ics.uci.edu/ml/datasets/YearPredictionMSD
  datapaths["yearpredictionmsd_raw"] = tobeset

  datapaths["htru2_write"] = tobeset
  #path to where HTRU_2.csv is  
  #from https://archive.ics.uci.edu/ml/datasets/HTRU2
  datapaths["htru2_raw"] = tobeset


  datapaths["epileptic_write"] = tobeset
  #path to where data.csv is  
  #https://archive.ics.uci.edu/ml/datasets/Epileptic+Seizure+Recognition#
  datapaths["epileptic_raw"] = tobeset

  datapaths["mopac_write"] = tobeset
  #path to where allUsers.lcl.csv is  
  #https://archive.ics.uci.edu/ml/machine-learning-databases/00391/
  datapaths["mopac_raw"] = tobeset

  #name of cPickle file where to write results from experiment
  datapaths["pkl_results_dir"] = tobeset

  #name of directory where the sizes of datasets can be stored (so that they don't have to 
  #be loaded to check)
  datapaths["pkl_results_dir"] = tobeset


  # base direcory of font ttf files. 
  # The font used in the poster is downloadable at
  # https://fontlibrary.org/en/font/cmu-bright
  datapaths["font_dirs"] = tobesest

