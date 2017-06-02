import numpy as np
import random
import multiprocessing as mpr
from libcpp.string cimport string
from libcpp.vector cimport vector 
from libcpp cimport bool
cimport cython
cimport cython.floating

from numbers import Number
ctypedef fused char_or_int:
  cython.int
  cython.char

#important : every function with floating template needs its own fused type. Is this a cython bug?
ctypedef fused floating33:
  cython.float
  cython.double

ctypedef fused floating87:
  cython.float
  cython.double

cdef extern from "zentasinfo.hpp" namespace "nszen":
  string get_output_inf_string() except +;
  string get_python_paramater_string() except +;
  
  
cdef extern from "zentas.hpp" namespace "nszen":
  void hello() except +


  # dense vectors 
  void vzentas[T](size_t ndata, size_t dimension, const T * const ptr_datain, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff) except +;

  # sparse vectors 
  void sparse_vector_zentas[T](size_t ndata, const size_t * const sizes, const T * const ptr_datain, const size_t * const ptr_indices_s, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff) except +;

  # strings / sequences 
  void szentas[T](size_t ndata, const size_t * const sizes, const T * const ptr_datain, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, bool with_cost_matrices, size_t dict_size, double c_indel, double c_switch, const double * const c_indel_arr, const double * const c_switches_arr, double critical_radius, double exponent_coeff) except +;

  # sequences from text file
  void textfilezentas(vector[string] filenames, string outfilename, string costfilename, size_t K, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, string initialisation_method) except +;


def dangerwrap(f):
  """
  I assume f is a function which returns 
  an object and takes no parameters
  """
  event  = mpr.Event()
  q = mpr.Queue()
  
  def signalling_f():
    try:
      q.put(f())
    
    except Exception as e:
      print "Caught exception in dangerwrap:"
      print e
      q.put(e)
      event.set()
      return None
      
    event.set()
  
  f_process = mpr.Process(target = signalling_f)
  f_process.start()
  try:
    event.wait()
  
  except KeyboardInterrupt:
    f_process.terminate()
    f_process.join()
    raise KeyboardInterrupt("Caught KeyboardInterrupt in dangerwrap")
  
  return q.get()


def basehello():
  hello()
  return "goodbye ! :():"


def pyhello():
  """
  Example function
  """	
  return dangerwrap(lambda : basehello)
  
  
def pyzentas(ndata = None, dimension = None, sizes = None, X = None, K = 10, indices_init = None, initialisation_method = "from_indices_init", algorithm = "clarans", level = 3, max_proposals = 10**9, capture_output = False, seed = 1011, max_time = 3.0, min_mE = 0.0, metric = "l2", nthreads = 1, max_rounds = 10**9, patient = True, energy = "identity", with_tests = False, rooted = True, with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None, indices_s = None, critical_radius = 0, exponent_coeff = 0, filenames_list = None, outfilename = None, costfilename = None):
  """  
  
  
The function, as defined in python/pyzentas.pyx
----------------------------
def pyzentas(ndata = None, dimension = None, sizes = None, X = None, K = 10, indices_init = None, initialisation_method = "from_indices_init", algorithm = "clarans", level = 3, max_proposals = 10**9, capture_output = False, seed = 1011, max_time = 3.0, min_mE = 0.0, metric = "l2", nthreads = 1, max_rounds = 10**9, patient = True, energy = "identity", rooted = True, with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None, indices_s = None, critical_radius = 0, exponent_coeff = 0, filenames_list = None, outfilename = None, costfilename = None):



Description
----------------------------
Fast K-Medoids (CLARANS) for diverse datatypes and metrics, as described in xx-xx-xx.


Quick example (dense data)
----------------------------
>> import numpy as np
>> import numpy.random as npr
>> ndata = 10000
>> dimension = 4
>> X = npr.randn(ndata, dimension)
>> K = 100
>> indices_init = range(K)
>> results = pyzentas.pyzentas(X = X, K = K, indices_init = indices_init, max_rounds = 150, max_time = 10, capture_output = True)
>> results.keys()
      ['output', 'labels', 'indices_final']
      where:     
          'labels'        : cluster identities (vector of size ndata)
          'indices_final' : the indices of the final medoids (vector of size K)
          'output'        : a string containing analysis of times in each component of the algorithm, iteration-by-iteration (see below for details)




All the Input Options  
----------------------------
parameter (relevant when, * : always relevant) [type of input] {non-optional}
description
----------------------------
numtype : np.float32, np.float64
seqtype : np.int32, np.int8 ('|S1') 


ndata (data not from file) {non-optional}
  The number of data samples 
  
dimension (densevector) {non-optional}
  The dimension of dense vectors  

sizes (sparsevector, sequence) [np.array of np.uint64] {non-optional}
  Vector of length ndata containing 
  for sparsevectors : the number of non-zero elements
  for sequences : their lengths

X (*) {non-optional}
  Vector containing,
  for densevectors, ndata*dimension values, contiguous by vector [np.array of numtype]
  for sparsevectors, the sizes.sum() non-zero values, contiguous by vector [np.array of numtype]
  for sequences, the sizes.sum() seqtypes, contiguous by sequence [np.array of seqtype]

K (*)
  The number of clusters
  
indices_init (when filenames_list = None and initialisation_method = "from_indices_init") [np.array of np.uint64] {optional}
  A vector of K distinct integers in [0, ndata), the starting centers

initialisation_method (*) [string] {non-optional}
  "from_indices_init" : use indices_init
  "uniform" : initialised as random indces
  other methods still TODO

algorithm (*) [string]
  The algorithm to use. One of
  "clarans" - essentially the algorithm of Ng et al. (1994), see Newling et al. (??) for details. 
  "voronoi" - the algorithm of Park et al. (2009)
  "clarans" provides superior clustering to "voronoi", and should be one's first choice. see reference for details

level (*)
  The level of optimisation to use. 
  for "clarans", one of 0,1,2,3. 
  for "voronoi", this must be 0.
  We suggest that one uses "clarans" with level 3. see reference for details

max_proposals (clarans)
  The number of consecutively rejected proposals at which to halt

capture_output (*)
  if True: the statistics of the run are returned as a string (accesible through 'output' key is returned dict. 
  See entry in Output section for details.
  if False: the statistics of the run are output to terminal while the algorithm is running

seed (clarans)
  positive integer which determines the series of swap proposals in clarans

max_time (*)
  a new round will not be started after max_time seconds have passed

min_mE (*)
  a new round will not be started if mE is less than min_mE

metric (*)
  consider vectors <2,2> and <5,6> for sparsevectors and densevectors, this can be one of
  l0 :  a count of in how many dimensions two vectors differ (2.0)
  l1 :  sum of absolute difference across dimensions (7.0)
  l2 :  standard Euclidean norm (5.0)
  li :  infinity-norm : the largest absolute difference across dimensions (4.0)
  
  for sequences, this can be one of
  levenshtein : levenshtein distance between sequences, consider setting indel and mutation costs described below
  normalised levenshtein: a transformed version of levenshtein so that distance is relative to sequennce lengths. See Yujian (2007)

nthreads (*)
  number of threads to use
  see reference for details

max_rounds (*)
  the maximum number of rounds (center changes) the algorithm runs for.

patient (clarans)
  if False  a proposal is accepted as soon as it results in a decrease in total energy
  if True   a proposal is accepted only when the time spent in center updating exceeds the time spent in assignment updating
  For reproducibility with a fixed seed, patient should be False. However patient = True is generally a better choice, especially with threads > 1
  see reference for details. 
  

energy (*)
  Recall that we are trying to minimise : [sum over sample i] E (d(i)). where d(i) = min_k distance (x(i), x(c(k)))
  energy can be one of 
  identity :        E(d) = d. corresponds to vanilla k-medoids, minimise the sum of distances to nearest centers
  quadratic :       E(d) = d*d. energy corresponds to the same loss function as k-means, a good choice for k-means initialisation
  cubic :           E(d) = d*d*d.
  squarepotential : E(d) = { 0 if d <= critical_radius and 1 otherwise }   
  log :             E(d) = log (1 + d), logarithm base e.
  exp :             E(d) = exp ( exponent_coeff * d)
  critical_radius and exponent_coeff can be passed as parameters to this function


rooted (*)
  This parameter has no influence on final clustering. It relates to how algorithms are implemented. If
  rooted = False, data (sequences, vectors etc.) are stored contiguously by cluster. 
  This means moving data between clusters when there are reassignments, but has the advantage of better memory access.
  rooted = True means data is not moved, only pointers to data are stored in clusters. 
  rooted = False is generally faster but has a higher memory footprint.
  Note : the memory assigned to each sample is the same. So for sparsevectors and sequences, if rooted = True, very sparse sparsevectors and short sequences will need as much memory as the least sparse sparsevector / longest sequence. 

with_cost_matrices (normalised levenshtein, levenshtein)
  if True: c_indel and c_switch may be set, these are global insertion/deletion and mutation parameters, independant of characters involved.
  if False: the indel and switch costs are dependent on the seqtypes involved, and so c_indel_arr and c_switches_arr arrays should be provided.

dict_size (normalised levenshtein, levenshtein)
  only needs to be set if c_indel_arr and c_switches_arr are set : the largest seqtype in a sequence 

c_indel (normalised levenshtein, levenshtein)
  see with_cost_matrices
  
c_switch (normalised levenshtein, levenshtein)
  see with_cost_matrices

c_indel_arr (normalised levenshtein, levenshtein)
  array of size dict_size, where c_indel_arr[i] is the cost of inserting / deleting an i.

c_switches_arr (normalised levenshtein, levenshtein)
  dict_size x dict_size array, where c_switches_arr[i,j] is the cost of transforming i into j. Of course, c_switches_arr[i,i] should be zero. see with_cost_matrices

indices_s (sparsevector)
  the indices of the sizes.sum() non-zero values, contiguous by vector [np.array of numtype]

critical_radius (squarepotential)
  see energy

exponent_coeff (exp)
  see energy

filenames_list (normalised levenshtein, levenshtein)
  (to be set with outfilename and costfilename)
  a list of strings, names of files to read from. The files should either be
  FASTA formatted or  : although not restricted to nucleotide (ACTG) or amino acids.
  ordinary text files : each line is considered a sequence (could be ordinary words / sentences for example)
  
outfilename (normalised levenshtein, levenshtein)
  (to be set with filenames_list and costfilename)
  where to write the clustering results 
  
costfilename (normalised levenshtein, levenshtein)
  (to be set with filenames_list and outfilename)
  where to obtain the indel and switch costs. 
  file format:
  either - 
  * * v1
  * v2
  for a global indel cost if v2 and swap cost of v1
  or -
  X Y vs1
  .
  .
  .
  Q W vsn
  A vi1
  .
  .
  .
  S vin
  where vs1 is the cost of swapping X and Y,  vi1 is the indel cost for A, etc.

Output 
----------------------------
indices_final
  a K-element array, the indices of the samples which are the final centers. Specifically, indices_final[k] is an integer in [0, ndata) for 0 <= k < K.
labels
  the assigned cluster of every sample. Specifically, for 0 <= i < ndata, 0 <= labels[i] < K is the cluster of sample i
indices_init
  TODO

  


Examples
----------------------------
see ./examples.py


References
----------------------------
xx-xx-xx
  """	
  null_size_t = np.empty((1,), dtype = np.uint64)
  null_int = np.empty((1,), dtype = np.int32)
  null_float = np.empty((1,), dtype = np.float32)
  null_double = np.empty((1,), dtype = np.float64)  
  
  
  
  if outfilename == None and isinstance(X, np.ndarray) and len(X.shape) == 2:

    if ndata != None and ndata != X.shape[0]:
      raise RuntimeError("X is a 2-D array (hence assumed to be dense vector data), with shape[0] = " + str(X.shape[0]) + ". This is in disagreement with parameter `ndata', which is " + str(ndata) + ".")
      
    if dimension != None and dimension != X.shape[1]:
      raise RuntimeError("X is a 2-D array (hence assumed to be dense vector data), with shape[1] = " + str(X.shape[1]) + ". This is in disagreement with parameter `dimension', which is " + str(dimension) + ".")

    ndata, dimension = X.shape

  
  dimension_sizes_string = "Dimension should be set for dense vector data. Sizes should be set for string and sparse vector data."
  if filenames_list != None and dimension != None and sizes != None:
    raise RuntimeError("Both dimension and sizes have been set. Either one or the other should be set. " + dimension_sizes_string)

  if filenames_list == None and dimension == None and sizes == None:      
    raise RuntimeError("Neither dimension nor sizes is set. Either one or the other should be set. " + dimension_sizes_string)

  if (ndata == None and X != None):
    raise RuntimeError("With X provided, ndata must be provided as well. The only time ndata should not be provided is when reading sequence data from file")
  
  if (ndata != None and X== None):
    raise RuntimeError("With X not provided, we assume that data should be read from file. In which case, ndata should be None")
  
  
  if filenames_list != None:
    
    if (sizes != None or c_indel != 1 or c_switch != 1 or c_indel_arr != None or c_switches_arr != None or indices_init != None or ndata != None):
      raise RuntimeError("With filenames_list != None, it is expected that none of (sizes, c_indel, c_switch, c_indel_arr, c_switches_arr, indices_init, ndata) be set")
    
    if not (isinstance(filenames_list, list) and isinstance(outfilename, str) and isinstance(costfilename, str)):
      raise RuntimeError("Expected : filenames_list : list of strings, outfilename : string, costfilename : string")
    
    indices_init = null_size_t
    X = null_float
    ndata = 0
#    return dangerwrap(lambda : basezentas("f", ndata, dimension, null_size_t, null_int, X.ravel(), K, indices_init.ravel(), initialisation_method, algorithm, level, max_proposals, capture_output, seed, max_time, min_mE, metric, nthreads, max_rounds, patient, energy, rooted, False, 0, -1., -1., null_double, null_double, null_size_t, critical_radius, exponent_coeff, filenames_list, outfilename, costfilename))

  else:
    
    if not isinstance(initialisation_method, str):
      raise RuntimeError("initialisation_method should be a string")      
    
    if initialisation_method == "from_indices_init":
      if indices_init == None:
        raise RuntimeError("indices_init == None and initialisation_method is from_indices_init")

      if isinstance(indices_init, list):
        indices_init = np.array(indices_init)
    
      if not isinstance(indices_init, np.ndarray):
        raise RuntimeError("indices_init should be an np.ndarray")
        
      if indices_init.size != K:
        raise RuntimeError("indices_init should be of size K")
      
      if indices_init.max() > ndata -1:
        raise RuntimeError("max in indices_init should not exceed ndata - 1")

      if (indices_init.dtype != np.uint64):
        indices_init = np.array(indices_init, dtype =  np.uint64)


    else:
      if indices_init != None:
        raise RuntimeError("indices_init != None and initialisation_method is not from_indices_init")    
    
      if indices_init == None:
        indices_init = null_size_t
        
    # Dense vector data.
    if dimension != None:
      pass
#      return dangerwrap(lambda : basezentas("v", ndata, dimension, null_size_t, null_int, X.ravel(), K, indices_init.ravel(), initialisation_method, algorithm, level, max_proposals, capture_output, seed, max_time, min_mE, metric, nthreads, max_rounds, patient, energy, rooted, False, 0, -1., -1., null_double, null_double, null_size_t, critical_radius, exponent_coeff, [], "", ""))
    
    # Sparse vector data or string data.
    else:
      # String data.
      if indices_s == None:  
        if c_indel_arr == None:
          c_indel_arr = null_double
          c_switches_arr = null_double
        
        
#        return dangerwrap(lambda : basezentas("s", ndata, 0, sizes.ravel(), X.ravel(), null_float, K, indices_init.ravel(), initialisation_method, algorithm, level, max_proposals, capture_output, seed, max_time, min_mE, metric, nthreads, max_rounds, patient, energy, rooted, with_cost_matrices, dict_size, c_indel, c_switch, c_indel_arr.ravel(), c_switches_arr.ravel(), null_size_t, critical_radius, exponent_coeff, [], "", ""))
  
      #Sparse vector data.
      else:
        pass
#        return dangerwrap(lambda : basezentas("sv", ndata, 0, sizes.ravel(), null_int, X.ravel(), K, indices_init.ravel(), initialisation_method, algorithm, level, max_proposals, capture_output, seed, max_time, min_mE, metric, nthreads, max_rounds, patient, energy, rooted, False, 0, -1., -1., null_double, null_double, indices_s.ravel(), critical_radius, exponent_coeff, [], "", ""))




cdef class RetBundle:

  cdef size_t [:] labels
  cdef size_t [:] indices_final
  cdef string output_string
  def __init__(self, ndata, K):
    self.labels = np.empty((ndata,), dtype = np.uint64)
    self.indices_final = np.empty((K,), dtype = np.uint64) 

  def get_dict(self):
    return {"output": self.output_string, 'indices_final': np.array(self.indices_final), 'labels': np.array(self.labels)}



#def pyzentas(ndata = None, dimension = None, sizes = None, X = None, K = 10, indices_init = None, initialisation_method = "from_indices_init", algorithm = "clarans", level = 3, max_proposals = 10**9, capture_output = False, seed = 1011, max_time = 3.0, min_mE = 0.0, metric = "l2", nthreads = 1, max_rounds = 10**9, patient = True, energy = "identity", rooted = True, with_cost_matrices = False, dict_size = 0, c_indel = 1, c_switch = 1, c_indel_arr = None, c_switches_arr = None, indices_s = None, critical_radius = 0, exponent_coeff = 0, filenames_list = None, outfilename = None, costfilename = None):
  
class pyzen(object):
  "class string"

  
  def set_init_string(self):
    pass  

  def __init__(self, 
    # generated in zentasinfo.cpp 
    K = 0, 
    algorithm = 'clarans', 
    capture_output = False, 
    critical_radius = 0, 
    energy = 'quadratic', 
    exponent_coeff = 0, 
    init = 'kmeans++-20', 
    level = 3, 
    max_proposals = 1e6, 
    max_rounds = 1e5, 
    max_time = 10, 
    metric = 'l2', 
    min_mE = 0, 
    nthreads = 1, 
    patient = True, 
    rooted = False, 
    seed = 1011, 
    with_tests = False):
    #"init string"

    
    self.set_init_string()


    self.pms = {
    'K':K,
    'algorithm': algorithm, 
    'level': level,
    'max_proposals': max_proposals,
    'max_rounds': max_rounds,
    'max_time': max_time,
    'min_mE': min_mE,
    'patient': patient,
    'capture_output': capture_output,
    'nthreads': nthreads,
    'rooted': rooted,
    'metric': metric,
    'energy':energy,
    'with_tests':with_tests,
    'exponent_coeff': exponent_coeff,
    'seed': seed,
    'critical_radius': critical_radius}

    self.null = {
    'size_t': np.empty((1,), dtype = np.uint64),
    'int': np.empty((1,), dtype = np.int32),
    'float': np.empty((1,), dtype = np.float32), 
    'double': np.empty((1,), dtype = np.float64)}


    if isinstance(init, str):
      if init == "initialisation_from_indices":
        raise RuntimeError("initialisation_from_indices not a valid initialisation string for the Python lib ")
        
      self.pms['initialisation_method'] = init
      self.pms['indices_init'] = self.null['size_t']
    
    else:
      self.pms['initialisation_method'] = 'from_indices_init'
      self.pms['indices_init'] = np.array(init, dtype = np.uint64)
    
      if self.pms['indices_init'].size != self.pms['K']:
        raise RuntimeError("indices_init should be of size K")
          
#  __init__.__func__.__string__ = "101010101"



  ##################################################
  ################# dense vectors ##################
  ##################################################
  
  def base_vzentas(self, floating87 [:] X_v, size_t [:] indices_init, pms):

    cdef void (*cw_vzentas)(size_t, size_t, const floating87 * const, size_t, const size_t * const, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool, string &, size_t seed, double max_time, double min_mE, size_t * const i_f, size_t * const labs, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff) except +
  
    if floating87 is double:
      cw_vzentas=&vzentas[double]
      
    elif floating87 is float:
      cw_vzentas=&vzentas[float]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
    
    cw_vzentas(pms['ndata'], pms['dimension'], &X_v[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], &rb.indices_final[0], &rb.labels[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'])

    return rb.get_dict()

  def den(self, X):
    """
    make some comment
    """
    
    if isinstance(X, list):
      X = np.array(X, dtype = np.float32)
    
    if len(X.shape) != 2:
      raise RuntimeError("X should be 2-dimensional, i.e. X.shape = (ndata, dimension)")
    
    self.pms['ndata'], self.pms['dimension'] = X.shape    
    return dangerwrap(lambda : self.base_vzentas(X.ravel(), self.pms['indices_init'], self.pms))


  ####################################################
  ################## sparse vectors ##################
  ####################################################

  
  def base_sparse_vector_zentas(self, size_t [:] sizes, size_t [:] indices, floating87 [:] values, size_t [:] indices_init, pms):
  
    cdef void (*cw_sparse_vector_zentas)(size_t, const size_t * const, const floating87 * const, const size_t * const, size_t, const size_t * const, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool, string &, size_t seed, double max_time, double min_mE, size_t * const i_f, size_t * const labs, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff) except +

    if floating87 is double:
      cw_sparse_vector_zentas=&sparse_vector_zentas[double]
        
    elif floating87 is float:
      cw_sparse_vector_zentas=&sparse_vector_zentas[float]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
    
    cw_sparse_vector_zentas(pms['ndata'], &sizes[0], &values[0], &indices[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], &(rb.indices_final)[0], &(rb.labels)[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'])

    return rb.get_dict()
      
  def spa(self, sizes, indices, values):
    """
    make some comment
    """
    self.pms['ndata'] = sizes.size
    if (indices.size != sizes.sum()):
      raise RuntimeError("the sum of sizes is not the size of indices ")

    if (values.size != sizes.sum()):
      raise RuntimeError("the sum of sizes is not the size of values ")
            
    return dangerwrap(lambda : self.base_sparse_vector_zentas(sizes, indices, values, self.pms['indices_init'], self.pms))


  ##################################################
  ################# sequence data ##################
  ##################################################

  def base_szentas(self, size_t [:] sizes, char_or_int [:] values, size_t [:] indices_init, with_cost_matrices, dict_size, c_indel, c_switch, double [:] c_indel_arr, double [:] c_switches_arr, pms):

    cdef void (*cw_szentas)(size_t, const size_t * const, const char_or_int * const, size_t, const size_t * const, string initialisation_method, string, size_t, size_t, bool , string & , size_t, double, double min_mE, size_t * const , size_t * const, string, size_t, size_t, bool, string, bool, bool, bool, size_t, double, double, const double * const, const double * const, double critical_radius, double exponent_coeff) except +
    
    if char_or_int is int:
      cw_szentas = &szentas[int]
    
    elif char_or_int is char:
      cw_szentas = &szentas[char]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
        
    cw_szentas(pms['ndata'], &sizes[0], &values[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], &(rb.indices_final)[0], &(rb.labels)[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], with_cost_matrices, dict_size, c_indel, c_switch, &c_indel_arr[0], &c_switches_arr[0], pms['critical_radius'], pms['exponent_coeff'])
    
    return rb.get_dict()
  
  def seq(self, sizes, values, cost_indel, cost_switch):
    """
    a comment here
    """
    if isinstance(cost_indel, Number) and isinstance(cost_switch, Number):
      dict_size = 0
      c_indel_arr = self.null['double']
      c_switch_arr = self.null['double']
      c_indel = cost_indel
      c_switch = cost_switch
      with_cost_matrices = False
    
    elif isinstance(cost_indel, np.ndarray) and isinstance(cost_switch, np.ndarray):
      dict_size = cost_indel.size
      if (cost_switch.size != dict_size*dict_size):
        raise RuntimeError("(size of cost_switch array) = (size of cost_indel array)**2. ")
      
      cost_switch = cost_switch.reshape(dict_size, dict_size)
      c_indel_arr = np.array(cost_indel, np.float64)
      c_switch_arr = np.array(cost_switch, np.float64)
      c_indel = 0
      c_switch = 0
      with_cost_matrices = True 
    
      if (values.dtype == np.dtype('c')):
        raise RuntimeError("arrays for cost_indel and cost_switch are not supported if values are characters. Either (1) use txt_seq, or  (2) convert/map your data to integers, for which cost_indel and cost_switch can be array data.")
    
      
      if values.size != sizes.sum():
        raise RuntimeError("The size of `values' should be the sum of `sizes'")
      
      #confirm that dict size agrees with values. this could go in c++.
      for v in values:
        if v >= dict_size:
          raise RuntimeError("The size of the indel array (" + str(dict_size) + ") is too small, there is a value of (" + str(v) + ").  The size of the indel array should be larger than any value received for cost_indel[value] to be valid. ")
          
        if c_indel_arr[v] <= 0:
          raise RuntimeError("All indel (and switch) costs must be strictly positive. While checking this, we noticed that c_indel_arr[" + str(v) + "] =" + str(c_indel_arr[v]) )
    else:
      raise RuntimeError("(cost_indel, cost_switch) should be (float, float) or (float array, float array)")
      
    
    
    self.pms['ndata'] = sizes.size  
    return dangerwrap(lambda : self.base_szentas(sizes, values, self.pms['indices_init'], with_cost_matrices, dict_size, c_indel, c_switch, c_indel_arr.ravel(), c_switch_arr.ravel(), self.pms))
      
    
  ##################################################
  ################# from files #####################
  ##################################################    

  def base_fromfiles(self, filenames_list, outfilename, costfilename, pms):
    cdef vector[string] filenames_vec
    for fn in filenames_list:
      filenames_vec.push_back(fn)

    dummy_ndata = 0
    cdef RetBundle rb = RetBundle(dummy_ndata, self.pms['K'])

    textfilezentas(filenames_vec, outfilename, costfilename, pms['K'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'], pms['initialisation_method'])
    
    return rb.get_dict()

  def txt_seq(self, filenames_list, outfile, costfile):
    """
    make some comment
    """            
    return dangerwrap(lambda : self.base_fromfiles(filenames_list, outfile, costfile, self.pms))


pyzen.__init__.__func__.__doc__ = get_python_paramater_string()

