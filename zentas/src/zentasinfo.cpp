#include "zentasinfo.hpp"
#include <vector>
#include <algorithm>
#include <tuple>

#include <map>
namespace nszen{
  




/* R=16    mE=3.23235   Ti=0  Tb=27  Tc=4      Tu=5      Tr=0    Tt=38     lg2nc(c)=13.834   lg2nc=17.009   pc=0.80808   nprops=3 */


std::map<std::string, std::string> get_output_keys(){
  return {{"R",  "round (where a round is defined by 1 or several centers changing)"},
  {"mE",  "mean energy over samples"},
  {"Ti", "time [in milliseconds] taken for initialising the center indices, i.e. to run one of {kmeans++-INT, afk-mc2-INT, uniform, etc}."},
  {"Tb", "time placing samples in clusters, initialising necessary cluster statistics. This is the time after Ti, but before the main loop."},
  {"Tc", "time spent in center update. For clarans, this is the time spent evaluating proposals. For Voronoi, this is time in updating assignments"},
  {"Tu", "time spent in updating. For clarans, this is the time spent determining the nearest and second nearest centers following a swap. For voronoi, this is time spent determining medoids."},
  {"Tr", "time spent actually implementing moves. This cost involves moving data between clusters, while maintaining a random order within clusters. If rooted is true, this can be expected to be higher than if rooted is false, with a correspondlingly lower `Tc' when rooted is false"},
  {"Tt", "total time"},
  {"lg2nc(c)", "log_2 of the number of distance calculations made in center update (corresponding to Tc)"},
  {"lg2nc", "log_2 of the total number of distance calculations made"},
  {"pc", "distance calculations can be terminated early if it can be established that they are above some threshold. This is particularly effective when the Levenshtein or normalised Levenshtein distance is used. pc measures how effective the early stopping in a distance computation is. For vectors, it counts how many dimensions are actually looked at before halting, for sequences it measures the ratio of the computed cells in the dynamic alg. to the total number of cells (product of sequence lengths). "},
  {"nprops", "(for clarans) the number of rejected proposals before one is accepted."}};
}


const std::map<std::string, std::string> output_keys = get_output_keys();

std::string get_us(){
  return "Newling and Fleuret, 2017";
}

std::string get_basic_info(){
  std::stringstream ss;
  ss << "Partitional clustering around samples, also known as K-Medoids. The objective to minimise is: [sum over sample indices i] E(d(i)), where d(i) = [min over center indices k] distance (x(i), x(c(k))) and E is an energy function. Samples may be dense vectors (2-D numpy arrays), sparse vectors or sequences, with diverse metrics and energies supported as described below. For more information, see " << get_us() << ".";
  return ss.str();  
}


std::map<std::string, std::tuple<std::string, std::string> > get_parameter_info_map(){

  std::map<std::string, std::tuple<std::string, std::string>> pim;
  pim["K"] = std::make_tuple("The number of clusters", "0");

  std::stringstream ss_alg;
  ss_alg << "One of `clarans' (Ng and Han 1994, " << get_us() << "), and `voronoi' (Hastie et al. 2001, Park et al. 2009). Clarans is always the better choice.";
  pim["algorithm"] = std::make_tuple(ss_alg.str(), "'clarans'");
 
 
  std::stringstream ss_lev;
  ss_lev << "The level of optimisation. For `clarans' one of 0,1,2,3, for `voronoi' it must be 0. Clarans at level 3 is always fastest. See " << get_us() << " for details of the optimisations at each level.";
  pim["level"] = std::make_tuple(ss_lev.str(), "3");

  pim["max_proposals"] = std::make_tuple("The number of consecutively rejected swap proposals at which to clarans must halt.", "1e6");

  pim["max_rounds"] = std::make_tuple("The maximum number of rounds (center updates) after which clarans/voronoi must halt.", "1e5");

  pim["max_time"] = std::make_tuple("The maximum amount of time [s] at which clarans/voronoi will start a new round, after this the algorithm halts.", "10");
  
  pim["min_mE"] = std::make_tuple("The minimum energy for which clarans/voronoi will start a new round, below this the algorithm will halt.", "0");
  

  
  
  std::stringstream ss_pat;
  ss_pat << "clarans alternates between evaluating random swaps, and implementing swaps which reduce overall energy. the standard approach (patient is false) is to implement a swap as soon as an enegy reducing proposal has found/evaluated. But swap evaluations are much faster than swap implementations, and so it might be better to find several energy reducing swaps and only implement the best of these, something like PAM of Kaufman and Rousseeuw (1990). When patient is true, a swap is only implemented if the total time spent evaluating swaps exceeds the total time implementing swaps. patient should always be true if non-reproducible clustering is ok. See " << get_us() << " for more details.";
  
  pim["patient"] = std::make_tuple(ss_pat.str(), "True");

  pim["capture_output"] =  std::make_tuple("if capture_output is true, the statistics of the run are stored as a string which accessible to the user after the run, otherwise if capture_output is false, the statistics of the run are output to terminal while the algorithm is running", "False");

  pim["nthreads"] = std::make_tuple("the number of threads to use", "1");

  std::stringstream ss_roo;
  ss_roo << "This advanced setting relates to how data is stored, in particular if clusters store their underlying data contiguously, or just the pointers to data. If rooted is false, data (sequences, dense vectors, sparse vectors, etc.) are stored contiguously by cluster. This requires moving data between clusters when data is reassigned, but has the big advantage of better locality for memory access. if rooted is true, data is not moved, only pointers to data are stored in clusters. When rooted is false, algorithms are generally faster but have higher memory footprints. Note that, when rooted is false, the amount of memory assigned to each sample is the same. So, for sparse vectors and sequences, where memory required per sample varies, if rooted is false, very sparse sparse vectors and short sequences will need as much memory as the least sparse sparse vector or the longest sequence. For more details, see " << get_us(); 

  pim["rooted"] = std::make_tuple(ss_roo.str(), "False");

  std::stringstream ss_met;
  ss_met << "for sparse vectors and dense vectors, this is one of `l0', `l1', `l1', `li'. For sequences it is one of `levenshtein', and `normalised levenshtein' of Yujian (2007). To illustrate the vector norms, consider vectors in 2-D, {2,2} and {5,6}. The distance between them is, (l0) 2.0  (l1) 7.0  (l2) 5.0  li (4.0). In particular, l0 counts how many dimensions the vectors differ in, l1 is the sum of the absolute differences, l2 is the standard Euclidean distance, and li is the largest absolute difference across the dimensions. For more details, see " << get_us();

  pim["metric"] = std::make_tuple(ss_met.str(), "'l2'");

  pim["energy"] = std::make_tuple("Recall, we wish to minimise : [sum over samples i] E(d(i)), where d(i) = min_k distance (x(i), x(c(k))) and E is the energy function, which is one of the following 5 functions -- `identity' : E(d) = d. corresponds to vanilla k-medoids. `quadratic' : E(d) = d*d. as in k-means, a good choice for k-means initialisation. `cubic' : E(d) = d*d*d. `squarepotential'  : E(d) = { 0 if d <= critical_radius and 1 otherwise }  where critical_radius is a separate parameter. `log'  : E(d) = log (1 + d), logarithm base e. `exp'  : E(d) = exp (exponent_coeff * d) where exponent_coeff is a separate parameter", "'quadratic'");

  pim["with_tests"] = std::make_tuple("For debugging purposes, run tests to check that optimisations are reliable", "False");
        
  pim["exponent_coeff"] = std::make_tuple("see energy.", "0");

  pim["critical_radius"] = std::make_tuple("see energy.", "0");

  pim["seed"] = std::make_tuple("Positive integer, seeds all randomness in all algorithms. ", "1011");

  pim["init"] = std::make_tuple("One of `uniform', `afk-mc2-INT', `kmeans++-INT', or an array/list of K initialising indices. afk-mc2-INT is the AFK-MC^2 algorithm of Bachem et al (2016) with mixing chain length INT. kmeans++-INT is batched kmeans++ using a type of sub-sampling to improve memory access, where INT is the batch count, INT=1 is exact k-means++. `uniform' samples uniformly from the samples.", "'kmeans++-20'");

  pim["initialistion_method"] = std::make_tuple("(C++ specific) `uniform', `afk-mc2-INT', `kmeans++-INT' or `from_indices_init'", "\"kmeans++-20\"");
 
  pim["indices_init"] = std::make_tuple("(C++ specific) pointer to initialising indices (if from_indices_init) or nullptr", "nullptr");



  pim["(out) indices_final"] = std::make_tuple("A K-element array, the indices of the samples which are the final centers. Specifically, indices_final[k] is an integer in [0, ndata) for 0 <= k < K", "none");

  pim["(out) labels"] = std::make_tuple("The assigned cluster of every sample. Specifically, for 0 <= i < ndata, 0 <= labels[i] < K is the cluster of sample i", "none");

  pim["(out) output"] = std::make_tuple("A string containing the times, floating point operations, etc, performed in each part of the algorithm, iteration-by-iteration.", "none");

  //pim["(out) indices_init"] = std::make_tuple("The indices used to initialise voronoi/clarans. That is, the indices after running one the initialisers (kmeans-++-INT, afk-mc2-INT, uniform, etc) ", "none");
  


  return pim;
  

  
};

const std::map<std::string,std::tuple<std::string, std::string>> parameter_info_map = get_parameter_info_map();


std::vector<std::string> get_python_constructor_parms(){
  std::vector<std::string> X = {"K", "algorithm", "level", "max_proposals", "max_rounds", "max_time", "min_mE", "patient", "capture_output", "nthreads", "rooted", "metric", "energy", "with_tests", "exponent_coeff", "critical_radius", "seed", "init"};
  std::sort(X.begin(), X.end());
  return X;
}
 
 
const std::vector<std::string> python_constructor_parms =  get_python_constructor_parms();






std::string get_parameter_info_string(std::string parameter){
  if (parameter_info_map.count(parameter) != 0) {
    return std::get<0>(parameter_info_map.at(parameter));
  }
  
  else{
    return "none";
  }
}

std::string get_parameter_default_string(std::string parameter){
  if (parameter_info_map.count(parameter) != 0) {
    return std::get<1>(parameter_info_map.at(parameter));
  }
  
  else{
    return "none";
  }
}



std::string squish_shift(std::string s, size_t n_chars_per_line = 82, size_t indent_size = 4){
  size_t pos = 0;
  std::stringstream ss;
  while (pos < s.size()){
    size_t incr = std::min<size_t>(s.size() - pos, n_chars_per_line - indent_size);
    for (size_t i = 0; i < indent_size; ++i){
      ss << " ";
    }
    ss << s.substr(pos, n_chars_per_line - indent_size) << "\n";
    pos += incr;
  }
  return ss.str();  
}

std::string get_equals_line(size_t n){
  std::stringstream ss;
  for (size_t i = 0; i < n; ++i){
    ss << "=";
  }
  ss << "\n";
  return ss.str();
}


std::string get_output_info_string(){

  std::vector<std::string> oks = {"R", "mE", "Ti", "Tb", "Tc", "Tu", "Tr", "Tt", "lg2nc(c)", "lg2nc", "pc", "nprops"};

  std::stringstream ss;
  ss << get_equals_line(82);
  ss << "The output string contains the following statistics\n";
  ss << get_equals_line(82);     
  for (auto & ok : oks){
    std::string val = output_keys.at(ok);
    //std::string key = std::get<0>(x);
    //std::string val = std::get<1>(x);
    ss << ok << ":\n";
    ss << squish_shift(val, 82, 2);
    ss << "\n";
  }
  //ss << get_equals_line(82);
  return ss.str();

}


  
std::string get_python_paramater_string(){
  
  std::stringstream ss;
  ss << get_equals_line(82);
  ss << squish_shift(get_basic_info(), 82, 0);
  ss << get_equals_line(82);
  ss << "Constructor parameters (and default)\n";
  ss << get_equals_line(82);

  
  for (auto & x : python_constructor_parms){
    ss << x  << " (" << get_parameter_default_string(x) << ")\n";
    ss << squish_shift(get_parameter_info_string(x));
    ss << "\n";
  }
  ss << get_equals_line(82);
  ss << "Output dict, after calling a clustering function\n";
  ss << get_equals_line(82);


  std::vector<std::string> python_output_parms = {"(out) indices_final", "(out) labels", "(out) output"};
  for (auto & x : python_output_parms){
    ss << x << "\n";
    ss << squish_shift(get_parameter_info_string(x));
    ss << "\n";
  }
  
  
  ss << get_output_info_string();


    

  
  
  return ss.str();
  
}

std::string get_python_function_decl(){

  std::stringstream ss;
  for (auto & x : python_constructor_parms){
    ss << x << " = " << get_parameter_default_string(x) << ", ";
    ss << "\n";
  }
  return ss.str();


}


  

}
