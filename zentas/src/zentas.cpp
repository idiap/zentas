// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/claransl0.hpp>
#include <zentas/claransl1.hpp>
#include <zentas/claransl2.hpp>
#include <zentas/claransl3.hpp>
#include <zentas/dispatch.hpp>
#include <zentas/levenshtein.hpp>
#include <zentas/lpmetric.hpp>
#include <zentas/outputwriter.hpp>
#include <zentas/vdimap.hpp>
#include <zentas/voronoil0.hpp>
#include <zentas/zentas.hpp>
#include <zentas/zentaserror.hpp>

namespace nszen
{


void confirm_can_refine(std::string algorithm){
  if (algorithm == "voronoi"){
    throw zentas::zentas_error("Cannot currently perform refinement after Voronoi initialization");
  }
}
/* dense vectors */

// 10% -> 80% faster if unrooted (!) :)
template <typename T>
void vzentas(size_t              ndata,
             size_t              dimension,
             const T* const      ptr_datain,
             size_t              K,
             const size_t* const indices_init,
             std::string         initialisation_method,
             std::string         algorithm,
             size_t              level,
             size_t              max_proposals,
             bool                capture_output,
             std::string&        text,
             size_t              seed,
             double              max_time,
             double              min_mE,
             double              max_itok,
             size_t* const       indices_final,
             size_t* const       labels,
             std::string         metric,
             size_t              nthreads,
             size_t              max_rounds,
             bool                patient,
             std::string         energy,
             bool                with_tests,
             bool                rooted,
             double              critical_radius,
             double              exponent_coeff,
             bool                do_vdimap,
             bool                do_refinement,
             std::string         rf_alg,
             size_t              rf_max_rounds,
             double              rf_max_time,
             bool                do_balance_labels)
{

  auto bigbang = std::chrono::high_resolution_clock::now();

  std::vector<T> v_mapped;
  const T*       true_ptr_datain;
  size_t         true_dimension;

  if (do_vdimap == false)
  {
    true_ptr_datain = ptr_datain;
    true_dimension  = dimension;
  }

  else
  {
    vdimap::vdimap<T>(v_mapped, ptr_datain, ndata, dimension, seed);
    true_ptr_datain = v_mapped.data();
    true_dimension  = v_mapped.size() / ndata;
  }

  /* TODO here : ptr_datain: shuffle if requested. Another flag :) */

  LpMetricInitializer metric_initializer;
  

  if (do_refinement){
    confirm_can_refine(algorithm);
  }

  
  metric_initializer.reset(metric, do_refinement, rf_alg, rf_max_rounds, rf_max_time);

  EnergyInitialiser energy_initialiser(critical_radius, exponent_coeff);

  ConstLengthInitBundle<T> datain_ib(ndata, true_dimension, true_ptr_datain);
  if (rooted == true)
  {
    zentas_base<VDataRooted<DenseVectorDataRootedIn<T>>, LpMetric<DenseVectorDataRootedIn<T>>>(
      datain_ib,
      K,
      indices_init,
      initialisation_method,
      algorithm,
      level,
      max_proposals,
      capture_output,
      text,
      seed,
      max_time,
      min_mE,
      max_itok,
      indices_final,
      labels,
      nthreads,
      max_rounds,
      patient,
      energy,
      with_tests,
      metric_initializer,
      energy_initialiser,
      bigbang,
      do_balance_labels);
  }

  else
  {
    zentas_base<VData<DenseVectorDataUnrootedIn<T>>, LpMetric<DenseVectorDataUnrootedIn<T>>>(
      datain_ib,
      K,
      indices_init,
      initialisation_method,
      algorithm,
      level,
      max_proposals,
      capture_output,
      text,
      seed,
      max_time,
      min_mE,
      max_itok,
      indices_final,
      labels,
      nthreads,
      max_rounds,
      patient,
      energy,
      with_tests,
      metric_initializer,
      energy_initialiser,
      bigbang,
      do_balance_labels);
  }
}

template void vzentas(size_t              ndata,
                      size_t              dimension,
                      const double* const ptr_datain,
                      size_t              K,
                      const size_t* const indices_init,
                      std::string         initialisation_method,
                      std::string         algorithm,
                      size_t              level,
                      size_t              max_proposals,
                      bool                capture_output,
                      std::string&        text,
                      size_t              seed,
                      double              max_time,
                      double              min_mE,
                      double              max_itok,
                      size_t* const       indices_final,
                      size_t* const       labels,
                      std::string         metric,
                      size_t              nthreads,
                      size_t              max_rounds,
                      bool                patient,
                      std::string         energy,
                      bool                with_tests,
                      bool                rooted,
                      double              critical_radius,
                      double              exponent_coeff,
                      bool                do_vdimap,
                      bool                do_refinement,
                      std::string         rf_alg,
                      size_t              rf_max_rounds,
                      double              rf_max_time,
                      bool                do_balance_labels);

template void vzentas(size_t              ndata,
                      size_t              dimension,
                      const float* const  ptr_datain,
                      size_t              K,
                      const size_t* const indices_init,
                      std::string         initialisation_method,
                      std::string         algorithm,
                      size_t              level,
                      size_t              max_proposals,
                      bool                capture_output,
                      std::string&        text,
                      size_t              seed,
                      double              max_time,
                      double              min_mE,
                      double              max_itok,
                      size_t* const       indices_final,
                      size_t* const       labels,
                      std::string         metric,
                      size_t              nthreads,
                      size_t              max_rounds,
                      bool                patient,
                      std::string         energy,
                      bool                with_tests,
                      bool                rooted,
                      double              critical_radius,
                      double              exponent_coeff,
                      bool                do_vdimap,
                      bool                do_refinement,
                      std::string         rf_alg,
                      size_t              rf_max_rounds,
                      double              rf_max_time,
                      bool                do_balance_labels);

/* sparse vectors */

template <typename T>
void sparse_vector_zentas(size_t              ndata,
                          const size_t* const sizes,
                          const T* const      ptr_datain,
                          const size_t* const ptr_indices_s,
                          size_t              K,
                          const size_t* const indices_init,
                          std::string         initialisation_method,
                          std::string         algorithm,
                          size_t              level,
                          size_t              max_proposals,
                          bool                capture_output,
                          std::string&        text,
                          size_t              seed,
                          double              max_time,
                          double              min_mE,
                          double              max_itok,
                          size_t* const       indices_final,
                          size_t* const       labels,
                          std::string         metric,
                          size_t              nthreads,
                          size_t              max_rounds,
                          bool                patient,
                          std::string         energy,
                          bool                with_tests,
                          bool                rooted,
                          double              critical_radius,
                          double              exponent_coeff,
                          bool                do_refinement,
                          std::string         rf_alg,
                          size_t              rf_max_rounds,
                          double              rf_max_time,
                          bool                do_balance_labels)
{

  auto bigbang = std::chrono::high_resolution_clock::now();

  LpMetricInitializer metric_initializer;
  
  if (do_refinement){
    confirm_can_refine(algorithm);
  }
  
  metric_initializer.reset(metric, do_refinement, rf_alg, rf_max_rounds, rf_max_time);


  EnergyInitialiser energy_initialiser(critical_radius, exponent_coeff);

  SparseVectorDataInitBundle<T> datain_ib(ndata, sizes, ptr_datain, ptr_indices_s);

  if (rooted == true)
  {
    zentas_base<SparseVectorDataRooted<SparseVectorDataRootedIn<T>>,
                LpMetric<SparseVectorDataRootedIn<T>>>(datain_ib,
                                                       K,
                                                       indices_init,
                                                       initialisation_method,
                                                       algorithm,
                                                       level,
                                                       max_proposals,
                                                       capture_output,
                                                       text,
                                                       seed,
                                                       max_time,
                                                       min_mE,
                                                       max_itok,
                                                       indices_final,
                                                       labels,
                                                       nthreads,
                                                       max_rounds,
                                                       patient,
                                                       energy,
                                                       with_tests,
                                                       metric_initializer,
                                                       energy_initialiser,
                                                       bigbang,
                                                       do_balance_labels);
  }

  else
  {
    zentas_base<SparseVectorData<SparseVectorDataUnrootedIn<T>>,
                LpMetric<SparseVectorDataUnrootedIn<T>>>(datain_ib,
                                                         K,
                                                         indices_init,
                                                         initialisation_method,
                                                         algorithm,
                                                         level,
                                                         max_proposals,
                                                         capture_output,
                                                         text,
                                                         seed,
                                                         max_time,
                                                         min_mE,
                                                         max_itok,
                                                         indices_final,
                                                         labels,
                                                         nthreads,
                                                         max_rounds,
                                                         patient,
                                                         energy,
                                                         with_tests,
                                                         metric_initializer,
                                                         energy_initialiser,
                                                         bigbang,
                                                         do_balance_labels);
  }
}

template void sparse_vector_zentas(size_t              ndata,
                                   const size_t* const sizes,
                                   const double* const ptr_datain,
                                   const size_t* const ptr_indices_s,
                                   size_t              K,
                                   const size_t* const indices_init,
                                   std::string         initialisation_method,
                                   std::string         algorithm,
                                   size_t              level,
                                   size_t              max_proposals,
                                   bool                capture_output,
                                   std::string&        text,
                                   size_t              seed,
                                   double              max_time,
                                   double              min_mE,
                                   double              max_itok,
                                   size_t* const       indices_final,
                                   size_t* const       labels,
                                   std::string         metric,
                                   size_t              nthreads,
                                   size_t              max_rounds,
                                   bool                patient,
                                   std::string         energy,
                                   bool                with_tests,
                                   bool                rooted,
                                   double              critical_radius,
                                   double              exponent_coeff,
                                   bool                do_refinement,
                                   std::string         rf_alg,
                                   size_t              rf_max_rounds,
                                   double              rf_max_time,
                                   bool                do_balance_labels);

template void sparse_vector_zentas(size_t              ndata,
                                   const size_t* const sizes,
                                   const float* const  ptr_datain,
                                   const size_t* const ptr_indices_s,
                                   size_t              K,
                                   const size_t* const indices_init,
                                   std::string         initialisation_method,
                                   std::string         algorithm,
                                   size_t              level,
                                   size_t              max_proposals,
                                   bool                capture_output,
                                   std::string&        text,
                                   size_t              seed,
                                   double              max_time,
                                   double              min_mE,
                                   double              max_itok,
                                   size_t* const       indices_final,
                                   size_t* const       labels,
                                   std::string         metric,
                                   size_t              nthreads,
                                   size_t              max_rounds,
                                   bool                patient,
                                   std::string         energy,
                                   bool                with_tests,
                                   bool                rooted,
                                   double              critical_radius,
                                   double              exponent_coeff,
                                   bool                do_refinement,
                                   std::string         rf_alg,
                                   size_t              rf_max_rounds,
                                   double              rf_max_time,
                                   bool                do_balance_labels);

/* strings */

template <typename T>
void szentas(size_t              ndata,
             const size_t* const sizes,
             const T* const      ptr_datain,
             size_t              K,
             const size_t* const indices_init,
             std::string         initialisation_method,
             std::string         algorithm,
             size_t              level,
             size_t              max_proposals,
             bool                capture_output,
             std::string&        text,
             size_t              seed,
             double              max_time,
             double              min_mE,
             double              max_itok,
             size_t* const       indices_final,
             size_t* const       labels,
             std::string         metric,
             size_t              nthreads,
             size_t              max_rounds,
             bool                patient,
             std::string         energy,
             bool                with_tests,
             bool                rooted,
             bool                with_cost_matrices,
             size_t              dict_size,
             double              c_indel,
             double              c_switch,
             const double* const c_indel_arr,
             const double* const c_switches_arr,
             double              critical_radius,
             double              exponent_coeff,
             bool                do_balance_labels)
{

  auto bigbang = std::chrono::high_resolution_clock::now();

  EnergyInitialiser energy_initialiser(critical_radius, exponent_coeff);

  std::vector<std::string> possible_metrics       = {"levenshtein", "normalised levenshtein"};
  bool                     valid_metric           = false;
  std::string              possible_metric_string = " [ ";
  for (auto& possible_metric : possible_metrics)
  {
    possible_metric_string += " ";
    possible_metric_string += possible_metric;
    possible_metric_string += " ";

    if (metric.compare(possible_metric) == 0)
    {
      valid_metric = true;
      break;
    }
  }

  possible_metric_string += " ] ";

  if (valid_metric == true)
  {  // (metric.compare("levenshtein") == 0 || metric.compare("normalised levenshtein") == 0){

    bool normalised = false;
    if (metric.compare("normalised levenshtein") == 0)
    {
      normalised = true;
    }

    LevenshteinInitializer metric_initializer;

    if (with_cost_matrices == false)
    {
      if (c_indel <= 0)
      {
        throw zentas::zentas_error(
          "with_cost_matrices == false : c_indel should be a positive real number");
      }

      if (c_switch <= 0)
      {
        throw zentas::zentas_error(
          "with_cost_matrices == false : c_switch should be a positive real number");
      }

      metric_initializer = LevenshteinInitializer(c_indel, c_switch, normalised);  //
    }

    else
    {
      if (dict_size == 0)
      {
        throw zentas::zentas_error(
          "with_cost_matrics is true, and dict_size == 0. this is invalid.");
      }

      metric_initializer =
        LevenshteinInitializer(dict_size, c_indel_arr, c_switches_arr, normalised);  //
    }

    VariableLengthInitBundle<T> datain_ib(ndata, sizes, ptr_datain);
    if (rooted == true)
    {
      zentas_base<SDataRooted<StringDataRootedIn<T>>, LevenshteinMetric<StringDataRootedIn<T>>>(
        datain_ib,
        K,
        indices_init,
        initialisation_method,
        algorithm,
        level,
        max_proposals,
        capture_output,
        text,
        seed,
        max_time,
        min_mE,
        max_itok,
        indices_final,
        labels,
        nthreads,
        max_rounds,
        patient,
        energy,
        with_tests,
        metric_initializer,
        energy_initialiser,
        bigbang,
        do_balance_labels);
    }

    else
    {
      zentas_base<SData<StringDataUnrootedIn<T>>, LevenshteinMetric<StringDataUnrootedIn<T>>>(
        datain_ib,
        K,
        indices_init,
        initialisation_method,
        algorithm,
        level,
        max_proposals,
        capture_output,
        text,
        seed,
        max_time,
        min_mE,
        max_itok,
        indices_final,
        labels,
        nthreads,
        max_rounds,
        patient,
        energy,
        with_tests,
        metric_initializer,
        energy_initialiser,
        bigbang,
        do_balance_labels);
    }
  }

  else
  {

    std::string errm("Currently, only metrics ");
    errm = errm + possible_metric_string + " are implemented for string data";
    throw zentas::zentas_error(errm);
  }
}

template void szentas(size_t              ndata,
                      const size_t* const sizes,
                      const int* const    ptr_datain,
                      size_t              K,
                      const size_t* const indices_init,
                      std::string         initialisation_method,
                      std::string         algorithm,
                      size_t              level,
                      size_t              max_proposals,
                      bool                capture_output,
                      std::string&        text,
                      size_t              seed,
                      double              max_time,
                      double              min_mE,
                      double              max_itok,
                      size_t* const       indices_final,
                      size_t* const       labels,
                      std::string         metric,
                      size_t              nthreads,
                      size_t              max_rounds,
                      bool                patient,
                      std::string         energy,
                      bool                with_tests,
                      bool                rooted,
                      bool                with_cost_matrices,
                      size_t              dict_size,
                      double              c_indel,
                      double              c_switch,
                      const double* const c_indel_arr,
                      const double* const c_switches_arr,
                      double              critical_radius,
                      double              exponent_coeff,
                      bool                do_balance_labels);

template void szentas(size_t              ndata,
                      const size_t* const sizes,
                      const char* const   ptr_datain,
                      size_t              K,
                      const size_t* const indices_init,
                      std::string         initialisation_method,
                      std::string         algorithm,
                      size_t              level,
                      size_t              max_proposals,
                      bool                capture_output,
                      std::string&        text,
                      size_t              seed,
                      double              max_time,
                      double              min_mE,
                      double              max_itok,
                      size_t* const       indices_final,
                      size_t* const       labels,
                      std::string         metric,
                      size_t              nthreads,
                      size_t              max_rounds,
                      bool                patient,
                      std::string         energy,
                      bool                with_tests,
                      bool                rooted,
                      bool                with_cost_matrices,
                      size_t              dict_size,
                      double              c_indel,
                      double              c_switch,
                      const double* const c_indel_arr,
                      const double* const c_switches_arr,
                      double              critical_radius,
                      double              exponent_coeff,
                      bool                do_balance_labels);

}  // namespace nszen
