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

import sys
import os
import socket

import load_joensuu_data
reload(load_joensuu_data)



zentas_dir = hardpaths.zentas_dir
infipath = hardpaths.infipath
infiexec = hardpaths.infiexec
zentas_output_dir = hardpaths.zentas_output_dir
zentas_data_dir = hardpaths.zentas_data_dir
clarans_vik_results_base = hardpaths.clarans_vik_results_base
english_words_dir = hardpaths.english_words_dir
fasta_data_dir = hardpaths.fasta_data_dir

sys.path.append(os.path.join(zentas_dir, "build/cleaninstall/lib"))

import pyzentas
import numpy as np
import numpy.random as npr
import random
import mnist
reload(mnist)
import pyzentas
pyzentas.pyhello()()


import genome_cutter
reload(genome_cutter)

def get_energy_time(output_string):
  """
  extract energy and time from output string.
  """
  E = []
  ctime = []
  ncalcs = []
  ttime = []
  for l in output_string.split("\n")[0:-1]:
    E.append(float(l.split("E: ")[1].split("\t")[0]))
    ctime.append(int(l.split("ctime: ")[1].split("\t")[0]))
    ttime.append(int(l.split("ttime: ")[1].split("\t")[0]))
    ncalcs.append(float(l.split("lg2 nc: ")[1].split("\t")[0]))
  
  return {'ttime' : np.array(ttime), 'ctime' : np.array(ctime), 'E': np.array(E), 'ncalcs' :np.array(ncalcs)}
  
  

def get_mut(base, n_muts):
  mut = [x for x in base]
  for m in range(n_muts):
    tau = npr.randint(3)
    if tau == 0: # flip
      index = npr.randint(len(mut))
      mut[index] = 1 - mut[index]
    elif tau == 1: # insert
      index = npr.randint(len(mut) + 1)
      mut.insert(index, npr.randint(2))
    else: # delete
      index = npr.randint(len(mut))
      mut.pop(index)
  return mut


def run_fromdata_experiment(experiment, algorithm, level, run, ndata, dimension, sizes, X, K, metric, energy, indices_s = None, maxtime = 10.0, critical_radius = 0.0, seed = None, rooted = False):
  """
  A single experiment, given the data (not from file)
  """
  
  print "ndata : ", ndata
  
  max_proposals = 1000000
  maxrounds = 100000
  exponent_coeff = 0.0
  if energy == 'exp':
    exponent_coeff = 1.0

  resultsdir = os.path.join(clarans_vik_results_base, "%s/%s_%d/run%d/"%(experiment, algorithm, level, run))
  if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)

  random.seed(seed)
  npr.seed(seed)
  indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
  indices_init.sort()

  results = pyzentas.pyzentas(ndata = ndata, dimension = dimension, sizes = sizes, X = X, K = K, indices_init = indices_init, algorithm = algorithm, level = level, max_proposals = max_proposals, capture_output = True, seed = seed, maxtime = maxtime, nthreads = 1, maxrounds = maxrounds, patient = True, metric = metric, rooted = rooted, energy = energy, with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None, indices_s = indices_s, exponent_coeff = exponent_coeff, critical_radius = critical_radius)

  filly = open("%s/results.txt"%(resultsdir,), "w")
  energy_time = get_energy_time(results['output'])
  
  for i in range(energy_time['ttime'].size):
    filly.write("%.5f\t%.5f\n"%(energy_time['ttime'][i]/1000.,  energy_time['E'][i])) #seconds
  filly.close()
  
  print "Complete!"

def run_fromdata_experiments_serial(ndata, dimension, sizes, X, K, metric, energy, indices_s = None, maxtime = 10.0, critical_radius = 0.0, nruns = 1, experiment = None, rooted = False):
  """
  """
  
  max_proposals = 1000000
  maxrounds = 100000
  exponent_coeff = 0.0
  if energy == 'exp':
    exponent_coeff = 1.0

  for run in range(nruns):
    print "run " , run
    print "***********"
    seed =  npr.randint(100000)
    random.seed(seed)
    npr.seed(seed)  
    indices_init = np.array(random.sample(xrange(ndata), K), dtype = np.uint64)
    indices_init.sort()   
    for algorithm, level in [["clarans", 3], ["clarans", 0], ["voronoi", 0]]:      
      print algorithm, " at " , level
      resultsdir = os.path.join(clarans_vik_results_base, "%s/%s_%d/run%d/"%(experiment, algorithm, level, run))
      if not os.path.exists(resultsdir):
        os.makedirs(resultsdir)
      filly = open("%s/results.txt"%(resultsdir,), "w")
      
      results = pyzentas.pyzentas(ndata = ndata, dimension = dimension, sizes = sizes, X = X, K = K, indices_init = indices_init, algorithm = algorithm, level = level, max_proposals = max_proposals, capture_output = True, seed = seed, maxtime = maxtime, nthreads = 1, maxrounds = maxrounds, patient = True, metric = metric, rooted = rooted, energy = energy, with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None, indices_s = indices_s, exponent_coeff = exponent_coeff, critical_radius = critical_radius)

      energy_time = get_energy_time(results['output'])
      for i in range(energy_time['ttime'].size):
        filly.write("%.5f\t%.5f\n"%(energy_time['ttime'][i]/1000.,  energy_time['E'][i])) #seconds
      filly.close()


def run_fromfile_experiment(algorithm, level, experiment, run, filenames_list, outfilename, costfilename, K, maxtime, metric, energy, seed):    
  
  max_proposals = 1000000
  maxrounds = 100000
  exponent_coeff = 0.0
  if energy == 'exp':
    exponent_coeff = 1.0

  print algorithm, " at " , level
  
  
  resultsdir = os.path.join(clarans_vik_results_base, "%s/%s_%d/run%d/"%(experiment, algorithm, level, run))
  if not os.path.exists(resultsdir):
    print "*", clarans_vik_results_base, "*"
    print resultsdir
    os.makedirs(resultsdir)

  filly = open("%s/results.txt"%(resultsdir,), "w")
  
  print "OUTFILENAME ", outfilename
  outfiledir = os.path.dirname(outfilename)
  if not os.path.isdir(outfiledir):
    os.makedirs(outfiledir)
  
  results = pyzentas.pyzentas(
  filenames_list = filenames_list, 
  outfilename = outfilename, 
  costfilename = costfilename, 
  K = K, 
  algorithm = algorithm, 
  level = level, 
  max_proposals = max_proposals, 
  capture_output = True, 
  maxtime = maxtime, 
  metric = metric, 
  nthreads = 1, 
  maxrounds = maxrounds, 
  patient = True, 
  energy = energy, 
  rooted = True, 
  seed = seed)
  
  energy_time = get_energy_time(results['output'])
  for i in range(energy_time['ttime'].size):
    filly.write("%.5f\t%.5f\n"%(energy_time['ttime'][i]/1000.,  energy_time['E'][i])) #seconds
  filly.close()
  
  print "completed!"
  

    
def run_fromfile_experiments_serial(filenames_list, outfilename, costfilename, K, maxtime, metric, energy, nruns = 1, experiment  = None):
  """
  """
  for run in range(nruns):
    print "run " , run
    print "***********"
    seed =  npr.randint(100000)
    random.seed(seed)
    npr.seed(seed)
    for algorithm, level in [["clarans", 3], ["clarans", 0], ["voronoi", 0]]:
      run_fromfile_experiment(algorithm, level, experiment, run, filenames_list, outfilename, costfilename, K, maxtime, metric,  energy,  seed)


def binary_strings_serial(base_length = 32, n_muts = 3, K = 200, N_per_K = 200, maxtime = 1.0, nruns = 1):
  """
  -- K 
  clusters, with "true" centers random binary strings of length 
  -- base_length. 
  Elements are 
  -- n_muts 
  mutations away from centers. There are exactly 
  -- N_per_K 
  elements per cluster
  """
  print "binary_strings_serial"
  ndata = K*N_per_K
  data = []
  sizes = []
  for k in range(K):
    base = [npr.randint(2) for i in range(base_length)]
    for n in range(N_per_K):
      mutated = get_mut(base, n_muts)
      data.extend(mutated)
      sizes.append(len(mutated))
      
  sizes = np.array(sizes, dtype = np.uint64)
  data = np.array(data, dtype = np.int8)
  
  run_fromdata_experiments_serial(ndata = ndata, dimension = None, sizes = sizes, X = data, K = K, metric = "levenshtein", energy = "identity", maxtime = maxtime, nruns = nruns, experiment = "binary_strings")
  

  
def sparse_vectors_serial(K = 200, average_per_cluster = 400, sparsity = 16, true_dimension = 10**6, maxtime = 1.0, nruns = 1):
  """
  sparse data. data from cluster k of {1, ..., K} : center(k) + Q * center(k') where Q is U[-1, 1] and k' is any center other than k.  
  center(k) is has 16 non-zero elements (of true_dimension), each N(0,1).
  """
  print "sparse_vectors_serial"
  ndata = K*average_per_cluster
  
  true_center_indices = [random.sample(xrange(true_dimension), sparsity/2)  for k in range(K)]
  true_center_values = [npr.randn(sparsity/2) for k in range(K)]

  sizes = sparsity*np.ones(ndata, dtype = np.uint64)
  indices_s = []
  data = []

  for d in range(ndata):
    
    k = npr.randint(K)    
    k_p = k
    k_p = npr.randint(K)
    while k_p == k:
      k_p = npr.randint(K)

    q = 0.5*(1 - 2*npr.rand())
    
    indices_local = np.concatenate([true_center_indices[k], true_center_indices[k]])    
    data_local = np.concatenate([true_center_values[k], q*true_center_values[k_p]])    
    indices_order = np.argsort(indices_local)
    indices_s.append(indices_local[indices_order])
    data.append(data_local[indices_order])
    
  indices_s = np.array(indices_s, dtype = np.uint64)
  data = np.array(data, dtype = np.float64)
  
  run_fromdata_experiments_serial(ndata = ndata, dimension = None, sizes = sizes, X = data, K = K, metric = "l2", energy = "quadratic", indices_s = indices_s, maxtime = maxtime, nruns = nruns, experiment = "sparse_vectors")



def dense_vector_1_serial(K = 121, N_per_K = 50, sigma = 1.0, maxtime = 1.0, nruns = 1):
  """
  On a grid
  """
  print "dense_vector_1_serial"
  rK = int(np.sqrt(K))
  ndata = N_per_K*K
  X = np.empty(shape = (ndata, 2), dtype = np.float64)
  X[:,0] = np.repeat(np.arange(rK), N_per_K*rK)
  X[:,1] = np.tile(np.arange(rK), N_per_K*rK)
  X += sigma*npr.randn(ndata, 2)

  run_fromdata_experiments_serial(ndata = ndata, dimension = 2, sizes = None, X = X, K = K, metric = "l1", energy = "exp", indices_s = None, maxtime = maxtime, nruns = nruns, experiment = "dense_vector_1")




def dense_vector_2_serial(K = 100, ndata = 10000, coverage = 1.1, maxtime = 1.0, nruns = 1):
  """
  coverage, If none of the blocks cross, what fraction of the grid is covered ?
  K*(2*r_square)**2 = coverage * 1.
  """
  print "dense_vector_2_serial"
  #data is in [0,1]^2
  data = npr.rand(ndata*2)
  r_square = np.sqrt(coverage / K) /2.
  print r_square
  run_fromdata_experiments_serial(ndata = ndata, dimension = 2, sizes = None, X = data, K = K, metric = "li", energy = "squarepotential", critical_radius = r_square, indices_s = None, maxtime = maxtime, nruns = nruns, experiment = "dense_vector_2")


def dense_high_d_serial(K = 10, ndata = 10000, dimension = 10, maxtime = 1.0, nruns = 1):
  print "dense_high_d_serial"
  data = np.array(npr.randn(ndata*dimension), dtype = np.float64)
  run_fromdata_experiments_serial(ndata = ndata, dimension = dimension, sizes = None, X = data, K = K, metric = "l2", energy = "identity", indices_s = None, maxtime = maxtime, nruns = nruns, experiment = "dense_high_d")


  

def do_sim_experiments_serial():
  maxtime = 64.0
  nruns = 4
  
  #syn-1
  N_per_K_syn1 = 50
  K_syn1 = 40
  N_syn1 = N_per_K_syn1*K_syn1
  binary_strings_serial(base_length = 16, n_muts = 2, K = K_syn1, N_per_K = N_per_K_syn1, maxtime = maxtime, nruns = nruns)
  
  #syn-2
  N_per_K_syn2 = 200
  K_syn2 = 100
  N_syn2 = N_per_K_syn2*K_syn2
  sparse_vectors_serial(K = K_syn2, average_per_cluster = N_per_K_syn2, sparsity = 16, true_dimension = 10**6, maxtime = maxtime, nruns = nruns)

  #syn-3
  N_per_K_syn3 = 200
  K_syn3 = 144
  N_syn3 = N_per_K_syn3*K_syn3  
  dense_vector_1_serial(K = K_syn3, N_per_K = N_per_K_syn3, sigma = 0.5, maxtime = maxtime, nruns = nruns)

  #syn-4
  K_syn4 = 100
  N_syn4 = K_syn4*200  
  dense_vector_2_serial(K = K_syn4, ndata = N_syn4, coverage = 1.1, maxtime = maxtime, nruns = nruns)

  counts_table = r"""
  
\begin{table}
\begin{center}
\begin{tabular}{|c|ccccc|}
\hline
& N & K & type & metric & $\energy(d)$ \\
\hline
syn-1 & %d & %d & string & Levenshtein & $d$ \\
syn-2 & %d & %d & sparse-v & $l_2$ & $d^2$ \\
syn-3 & %d & %d & dense-v & $l_1$ & $e^d$ \\
syn-4 & %d & %d & dense-v & $l_{\infty}$ & $I_{d > 0.05}$ \\
\hline
\end{tabular}
\end{center}

"""%(N_syn1, K_syn1, N_syn2, K_syn2, N_syn3, K_syn3, N_syn4, K_syn4)

  filly = open(path.join(hardpaths.clarans_paper_dir, "counts_table.txt"), "w")
  filly.write(counts_table)
  filly.close()

def english_words_serial(K = 600, maxtime = 2):
  experiment = "english_words_abridged"
  run_fromfile_experiments_serial(filenames_list = [os.path.join(english_words_dir, "words_abridged.txt")], outfilename = os.path.join(zentas_output_dir, "fromfile_output", "%s.txt"%(experiment,)), costfilename = os.path.join(zentas_data_dir, "costs.txt"), K = K, maxtime = maxtime, metric = "levenshtein", energy = "quadratic", nruns = 2, experiment = experiment)


def genome_serial(K = 100, n_snips = 20000, snip_range = [8, 10], maxtime = 20.):
  """
  cost matrix from http://www.hrbc-genomics.net/training/bcd/Curric/PrwAli/nodeD.html
  """
  
  experiment = "genome"
  genome_cutter.make_snippet_file(chr_38_fn = os.path.join(fasta_data_dir, "Homo_sapiens.GRCh38.dna.chromosome.10.fa"), n_snips = n_snips,snip_range = snip_range)

  root = os.path.join(zentas_output_dir, "data") #os.path.join(os.path.expanduser("~"), "libraries/zentas/data/")
  inputroot = fastadatadir
  
  run_fromfile_experiments_serial(filenames_list = [os.path.join(fasta_data_dir,"Homo_sapiens.GRCh38.dna.chromosome.10.fa.snips.fa")], outfilename = os.path.join(zentas_output_dir, "fromfile_output", "%s.txt"%(experiment,)), costfilename = os.path.join(zentas_data_dir, "ACTG_costs.txt"), K = K, maxtime = maxtime, metric = "normalised levenshtein", energy = "identity", nruns = 2, experiment = experiment)


def english_words_parallel(algorithm, level, run, K = 600, maxtime = 10.):  
  experiment = "english_words_parallel"
  filenames_list = [os.path.join(english_words_dir,"words.txt")]
  outfilename = os.path.join(zentas_output_dir, "fromfile_output", "%s_%s_%d_%d.txt"%(experiment,algorithm, level, run))
  costfilename = os.path.join(zentas_data_dir, "costs.txt")
  metric = "levenshtein"
  energy = "quadratic"
  seed =  npr.randint(100000)
  random.seed(seed)
  npr.seed(seed)
  run_fromfile_experiment(algorithm, level, experiment, run, filenames_list, outfilename, costfilename, K, maxtime, metric, energy, seed)


def genome_parallel(algorithm, level, run, K = 100, maxtime = 20.):  
  experiment = "genome_parallel"
  filenames_list = [os.path.join(fasta_data_dir,"Homo_sapiens.GRCh38.dna.chromosome.10.fa.snips.fa")]
  outfilename = os.path.join(zentas_output_dir, "fromfile_output", "%s_%s_%d_%d.txt"%(experiment,algorithm, level, run))
  costfilename = os.path.join(zentas_data_dir, "ACTG_costs.txt")
  metric = "normalised levenshtein"
  energy = "quadratic"
  seed =  npr.randint(100000)
  random.seed(seed)
  npr.seed(seed)
      
  run_fromfile_experiment(algorithm, level, experiment, run, filenames_list, outfilename, costfilename, K, maxtime, metric, energy, seed)    


def mnist_serial(K = 200, ndata = 60000, maxtime = 300., nruns = 2):
  X = mnist.read_MNIST(dataset = "original", ndata = ndata, dimension = None)
  ndata, dimension = X.shape
  run_fromdata_experiments_serial(ndata = ndata, dimension = dimension, sizes = None, X = X, K = K, metric = "l2", energy = "identity", indices_s = None, maxtime = maxtime, nruns = nruns, experiment = "mnist")


def mnist_parallel(algorithm, level, run, K = 400, maxtime = 20.): # ndata = 60000
  experiment = "mnist_parallel"
  
  #When grid is up again, I'll fix this:
  X = load_joensuu_data.load_data(name = "mnist")
  
   #npr.randn(ndata,50) #mnist.read_MNIST(dataset = "original", ndata = ndata, dimension = None)
  ndata, dimension = X.shape
  sizes = None
  metric = "l2"
  energy = "quadratic"
  seed = npr.randint(100000)
  run_fromdata_experiment(experiment, algorithm, level, run, ndata, dimension, sizes, X, K, metric, energy, indices_s = None, maxtime = maxtime, critical_radius = 0.0, seed = seed) 
  

def sparse_rcv1_serial(K = 200, maxtime = 80., nruns = 2):
  experiment = "rcv1"

  npzfiles = np.load(os.path.join(userdir, "a13-vector-files/lyrl2004_vectors_train.npz"))
  sizes = npzfiles['sizes']
  sizes = np.array(sizes, dtype = np.uint64)

  indices = npzfiles['indices']
  indices = np.array(indices, dtype = np.uint64)

  values = npzfiles['values']
  
  #from IPython.core.debugger import Tracer
  #Tracer()()
  ndata = len(sizes)
  
  print "ndata : ", ndata
  run_fromdata_experiments_serial(ndata = ndata, dimension = None, sizes = sizes, X = values, K = K, metric = "l2", energy = "quadratic", indices_s = indices, maxtime = maxtime, nruns = nruns, experiment = experiment, rooted = True)
    


def sparse_rcv1_parallel(algorithm, level, run, K = 50, maxtime = 80.):
  experiment = "rcv1_parallel"

  npzfiles = np.load(os.path.join(userdir, "a13-vector-files/lyrl2004_vectors_train.npz"))
  sizes = npzfiles['sizes']
  #sizes = np.array(sizes, dtype = np.uint64)
  ndata = len(sizes)
  dimension = None
  indices = npzfiles['indices']
  #indices = np.array(indices, dtype = np.uint64)

  values = npzfiles['values']
  metric = "l2"
  energy = "quadratic"
  seed = npr.randint(100000)
  
  print "ndata : ", ndata
  run_fromdata_experiment(experiment, algorithm, level, run, ndata, dimension, sizes, values, K, metric, energy, indices_s = indices, maxtime = maxtime, critical_radius = 0.0, seed = seed, rooted = True)
  


