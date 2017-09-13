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

cdef extern from "zentas/zentasinfo.hpp" namespace "nszen":
  
  string get_python_init_string() except +;
  string get_python_txt_seq_string() except +;
  string get_python_seq_string() except +;
  string get_python_spa_string() except +;
  string get_python_den_string() except +;
  string get_output_verbose_string() except +;


cdef extern from "zentas/zentas.hpp" namespace "nszen":

  # dense vectors 
  void vzentas[T](size_t ndata, size_t dimension, const T * const ptr_datain, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, double max_itok, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, bool do_vdimap, bool do_refinement, string rf_alg,size_t rf_max_rounds, double rf_max_time, bool do_balance_labels) except +;

  # sparse vectors 
  void sparse_vector_zentas[T](size_t ndata, const size_t * const sizes, const T * const ptr_datain, const size_t * const ptr_indices_s, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, double max_itok, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, bool do_refinement, string rf_alg, size_t rf_max_rounds, double rf_max_time, bool do_balance_labels) except +;

  # strings / sequences 
  void szentas[T](size_t ndata, const size_t * const sizes, const T * const ptr_datain, size_t K, const size_t * const indices_init, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, double max_itok, size_t * const indices_final, size_t * const labels, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, bool with_cost_matrices, size_t dict_size, double c_indel, double c_switch, const double * const c_indel_arr, const double * const c_switches_arr, double critical_radius, double exponent_coeff, bool do_balance_labels) except +;

  # sequences from text file
  void textfilezentas(vector[string] filenames, string outfilename, string costfilename, size_t K, string algorithm, size_t level, size_t max_proposals, bool capture_output, string & text, size_t seed, double max_time, double min_mE, double max_itok, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, string initialisation_method, bool do_balance_labels) except +;


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


cdef class RetBundle:

  cdef size_t [:] labels
  cdef size_t [:] indices_final
  cdef string output_string
  def __init__(self, ndata, K):
    self.labels = np.empty((ndata,), dtype = np.uint64)
    self.indices_final = np.empty((K,), dtype = np.uint64) 

  def get_dict(self):
    return {"output": self.output_string, 'indices_final': np.array(self.indices_final), 'labels': np.array(self.labels)}

  
class pyzen(object):


  def set_init_string(self):
    pass  

  def __init__(self, 
    # generated in zentasinfo.cpp, to sync default params with info string *#14641#*
    K = 0,
    algorithm = 'clarans',
    capture_output = False,
    critical_radius = 0,
    do_balance_labels = False,
    energy = 'quadratic',
    exponent_coeff = 0,
    init = 'kmeans++-20',
    level = 3,
    max_itok = 1e7,
    max_proposals = 1e6,
    max_rounds = 1e5,
    max_time = 10,
    metric = 'l2',
    min_mE = 0,
    nthreads = 1,
    patient = True,
    rooted = False,
    seed = 1011,
    with_tests = False,
    # kwargs used just to give hints for the bad input argument message
    **kwargs):
      
    if kwargs.keys():
      argstring = ""
      for x in kwargs.keys():
        argstring += ("\'" + x + "\'")
      argstring += ". "
      for x in kwargs.keys():
        if "test" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'with_tests\'? "
        if "coeff" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'energy_coeff\'? "
        if "thread" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'nthreads\'? "
        if "proposal" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'max_proposals\'? "
        if "round" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'max_rounds\'? "
        if "time" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'max_time\'? "
        if "init" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'init\'? "
        if "verbose" in x or "output" in x or "capture" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'capture_output\'? "
        if "max_i" in x:
          argstring += "Maybe instead of \'" + x + "\', you mean \'max_itor\'? "

      
      
      raise TypeError("Got unexpected argument(s) in pyzen.__init__ : " + argstring)

    
    self.set_init_string()


    self.pms = {
    'K':K,
    'algorithm': algorithm, 
    'level': level,
    'max_proposals': max_proposals,
    'max_rounds': max_rounds,
    'max_time': max_time,
    'min_mE': min_mE,
    'max_itok':max_itok,
    'patient': patient,
    'capture_output': capture_output,
    'nthreads': nthreads,
    'rooted': rooted,
    'metric': metric,
    'energy':energy,
    'with_tests':with_tests,
    'exponent_coeff': exponent_coeff,
    'seed': seed,
    'critical_radius': critical_radius, 
    'do_balance_labels' : do_balance_labels
    }

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
          


  ##################################################
  ################# dense vectors ##################
  ##################################################
  
  
  def base_vzentas(self, floating87 [:] X_v, do_vdimap, do_refinement, rf_alg, rf_max_rounds, rf_max_time, size_t [:] indices_init, pms):

    cdef void (*cw_vzentas)(size_t, size_t, const floating87 * const, size_t, const size_t * const, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool, string &, size_t seed, double max_time, double min_mE, double max_itok, size_t * const i_f, size_t * const labs, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, bool do_vdimap, bool do_refinement, string rf_alg,  size_t rf_max_rounds, double rf_max_time, bool do_balance_labels) except +
  
    if floating87 is double:
      cw_vzentas=&vzentas[double]
      
    elif floating87 is float:
      cw_vzentas=&vzentas[float]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
    
    cw_vzentas(pms['ndata'], pms['dimension'], &X_v[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], pms['max_itok'], &rb.indices_final[0], &rb.labels[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'], do_vdimap, do_refinement, rf_alg, rf_max_rounds, rf_max_time, pms['do_balance_labels'])

    return rb.get_dict()

  def den(self, 
    X, 
    do_vdimap, 
    do_refinement = False, 
    rf_alg = "yinyang", 
    rf_max_rounds = 99999, 
    rf_max_time = 1e7
    ):
    """
    den(se) clustering
    """
    
    if isinstance(X, list):
      X = np.array(X, dtype = np.float32)
    
    if len(X.shape) != 2:
      raise RuntimeError("X should be 2-dimensional, i.e. X.shape = (ndata, dimension)")
    
    self.pms['ndata'], self.pms['dimension'] = X.shape    
    return dangerwrap(lambda : self.base_vzentas(X.ravel(), do_vdimap, do_refinement, rf_alg, rf_max_rounds, rf_max_time, self.pms['indices_init'], self.pms))

  ####################################################
  ################## sparse vectors ##################
  ####################################################

  
  def base_sparse_vector_zentas(self, size_t [:] sizes, size_t [:] indices, floating87 [:] values, do_refinement, rf_alg, rf_max_rounds, rf_max_time, size_t [:] indices_init, pms):
  
    cdef void (*cw_sparse_vector_zentas)(size_t, const size_t * const, const floating87 * const, const size_t * const, size_t, const size_t * const, string initialisation_method, string algorithm, size_t level, size_t max_proposals, bool, string &, size_t seed, double max_time, double min_mE, double max_itok, size_t * const i_f, size_t * const labs, string metric, size_t nthreads, size_t max_rounds, bool patient, string energy, bool with_tests, bool rooted, double critical_radius, double exponent_coeff, bool do_refinement, string rf_alg, size_t rf_max_rounds, double rf_max_time, bool do_balance_labels) except +

    if floating87 is double:
      cw_sparse_vector_zentas=&sparse_vector_zentas[double]
        
    elif floating87 is float:
      cw_sparse_vector_zentas=&sparse_vector_zentas[float]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
    
    cw_sparse_vector_zentas(pms['ndata'], &sizes[0], &values[0], &indices[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], pms['max_itok'], &(rb.indices_final)[0], &(rb.labels)[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'], do_refinement, rf_alg, rf_max_rounds, rf_max_time, pms['do_balance_labels'])

    return rb.get_dict()
      
  def spa(self, 
    sizes, 
    indices, 
    values, 
    do_refinement = False, 
    rf_alg = "yinyang", 
    rf_max_rounds = 99999, 
    rf_max_time = 1e7
    ):
    """
    spa(rse) clustering
    """
    self.pms['ndata'] = sizes.size
    if (indices.size != sizes.sum()):
      raise RuntimeError("the sum of sizes is not the size of indices ")

    if (values.size != sizes.sum()):
      raise RuntimeError("the sum of sizes is not the size of values ")
            
    return dangerwrap(lambda : self.base_sparse_vector_zentas(sizes, indices, values, do_refinement, rf_alg, rf_max_rounds, rf_max_time, self.pms['indices_init'], self.pms))


  ##################################################
  ################# sequence data ##################
  ##################################################

  def base_szentas(self, size_t [:] sizes, char_or_int [:] values, size_t [:] indices_init, with_cost_matrices, dict_size, c_indel, c_switch, double [:] c_indel_arr, double [:] c_switches_arr, pms):

    cdef void (*cw_szentas)(size_t, const size_t * const, const char_or_int * const, size_t, const size_t * const, string initialisation_method, string, size_t, size_t, bool , string & , size_t, double, double min_mE, double max_itok, size_t * const , size_t * const, string, size_t, size_t, bool, string, bool, bool, bool, size_t, double, double, const double * const, const double * const, double critical_radius, double exponent_coeff, bool do_balance_labels) except +
    
    if char_or_int is int:
      cw_szentas = &szentas[int]
    
    elif char_or_int is char:
      cw_szentas = &szentas[char]

    cdef RetBundle rb = RetBundle(self.pms['ndata'], self.pms['K'])
        
    cw_szentas(pms['ndata'], &sizes[0], &values[0], pms['K'], &indices_init[0], pms['initialisation_method'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], pms['max_itok'], &(rb.indices_final)[0], &(rb.labels)[0], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], with_cost_matrices, dict_size, c_indel, c_switch, &c_indel_arr[0], &c_switches_arr[0], pms['critical_radius'], pms['exponent_coeff'], pms['do_balance_labels'])
    
    return rb.get_dict()
  
  def seq(self, sizes, values, cost_indel, cost_switch):
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

    textfilezentas(filenames_vec, outfilename, costfilename, pms['K'], pms['algorithm'], pms['level'], pms['max_proposals'], pms['capture_output'], rb.output_string, pms['seed'], pms['max_time'], pms['min_mE'], pms['max_itok'], pms['metric'], pms['nthreads'], pms['max_rounds'], pms['patient'], pms['energy'], pms['with_tests'], pms['rooted'], pms['critical_radius'], pms['exponent_coeff'], pms['initialisation_method'], pms['do_balanace_centers'])
    
    return rb.get_dict()

  def txt_seq(self, filenames_list, outfile, costfile):
    return dangerwrap(lambda : self.base_fromfiles(filenames_list, outfile, costfile, self.pms))
    
  def get_output_verbose_string():
    return get_output_verbose_string()

def kmeans(X, K, seed = 1011, maxtime = 1e9, capture_output = True):
  """
  Quick-and-easy vanilla k-means for the impatient. For more options, such as different metrics and extended stopping criteria, consider creating a pyzen object and using the class function den for dense clustering. Etc.
  X:               2-D numpy array, ndata x dimension
  seed:            random seed
  maxtime:         maximum allotted time to run k-means (excluding initialisation)
  capture_output:  if True, output string is returned in dict.  
  """
  
  #max_itok will be the stopping criterion
  z = pyzen(init = "kmeans++-5", K = K, metric = 'l2', energy = 'quadratic', exponent_coeff = 0, max_rounds = 100000, max_time = 100000, max_itok = 10.0, seed = seed, nthreads = 1, patient = True, with_tests = False, algorithm = "clarans", level = 3, capture_output = capture_output)
  
  algorithm = "yinyang"
  if X.shape[1] < 5 :
    algorithm = "exponion"
  
  maxtime = 123123
  return z.den(X, True, True, algorithm, 1000000, maxtime)
  
  
pyzen.__init__.__func__.__doc__ = get_python_init_string()
pyzen.txt_seq.__func__.__doc__ = get_python_txt_seq_string()
pyzen.seq.__func__.__doc__ = get_python_seq_string()
pyzen.spa.__func__.__doc__ = get_python_spa_string()
pyzen.den.__func__.__doc__ = get_python_den_string()
pyzen.get_output_verbose_string.__func__.__doc__ = get_output_verbose_string()
