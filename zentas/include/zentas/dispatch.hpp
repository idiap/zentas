// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZEN_DISPATCH_HPP
#define ZEN_DISPATCH_HPP

#include <zentas/baseclusterer.hpp>
#include <zentas/energyinit.hpp>
#include <zentas/extrasbundle.hpp>
#include <zentas/voronoil0.hpp>
#include <zentas/claransl0.hpp>
#include <zentas/claransl1.hpp>
#include <zentas/claransl2.hpp>
#include <zentas/claransl3.hpp>

namespace nszen
{

                        
void scrutinize_input_1(const EnergyInitialiser& energy_initialiser,
                        std::string              energy,
                        size_t                   K,
                        std::string              algorithm,
                        size_t                   level,
                        size_t                   ndata);

template <typename TData, typename TMetric, typename TInitBundle>
void dispatch(std::string                                                 algorithm,
              size_t                                                      level,
              TInitBundle&                                                datain_ib,
              size_t                                                      K,
              const size_t* const                                         indices_init,
              std::string                                                 initialisation_method,
              size_t                                                      max_proposals,
              size_t                                                      seed,
              double                                                      max_time,
              double                                                      min_mE,
              double                                                      max_itok,
              size_t* const                                               indices_final,
              size_t* const                                               labels,
              size_t                                                      nthreads,
              size_t                                                      max_rounds,
              bool                                                        patient,
              std::string                                                 energy,
              bool                                                        with_tests,
              const typename TMetric::Initializer&                        metric_initializer,
              const EnergyInitialiser&                                    energy_initialiser,
              std::chrono::time_point<std::chrono::high_resolution_clock> bigbang,
              bool                                                        do_balance_labels)
{

  typedef typename TData::DataIn DataIn;

  DataIn datain(datain_ib);

  nszen::SkeletonClustererInitBundle sc(K,
                                        datain_ib.ndata,
                                        bigbang,
                                        indices_init,
                                        initialisation_method,
                                        max_time,
                                        min_mE,
                                        max_itok,
                                        max_rounds,
                                        nthreads,
                                        seed,
                                        energy,
                                        with_tests,
                                        indices_final,
                                        labels,
                                        &energy_initialiser,
                                        do_balance_labels);
  ExtrasBundle eb(max_proposals, patient);
  ClustererInitBundle<DataIn, TMetric> ib(sc, datain, metric_initializer, eb);

  //  BaseClaransInitBundle clib();

  if (algorithm.compare("clarans") == 0)
  {
    if (level == 0)
    {
      // throw zentas::zentas_error("clarans l0 not enabled, grep seow340jkosdm4 and uncomment here
      // to enable");
      nszen::Clusterer<TMetric, TData, ClaransL0> cc(ib);
      cc.go();
    }
    if (level == 1)
    {
      // throw zentas::zentas_error("clarans l1 not enabled, grep sdrfoweinsdima and uncomment here
      // to enable");
      nszen::Clusterer<TMetric, TData, ClaransL1> cc(ib);
      cc.go();
    }
    if (level == 2)
    {
      // throw zentas::zentas_error("clarans l2 not enabled, grep sdmi4sdfsdollll and uncomment
      // here to enable");
      nszen::Clusterer<TMetric, TData, ClaransL2> cc(ib);
      cc.go();
    }
    if (level == 3)
    {
      // throw zentas::zentas_error("clarans l3 not enabled, grep sdmi4sdfsdfollll and uncomment
      // here to enable");
      nszen::Clusterer<TMetric, TData, ClaransL3> cc(ib);
      cc.go();
    }
  }

  else if (algorithm.compare("voronoi") == 0)
  {
    // throw zentas::zentas_error("voronoi not enabled, grep dfseimmgrfiddiddidiid and uncomment
    // here to enable");
    if (level == 0)
    {
      nszen::Clusterer<TMetric, TData, VoronoiL0> cc(ib);
      cc.go();
    }
  }

  else
  {
    throw zentas::zentas_error("unrecognised algorithm in dispatch");
  }
}

/* This is the place to do all kinds of tests on the input: all user calls (R/Python/Terminal) will
 * pass through this function */
template <typename TData, typename TMetric, typename TInitBundle>
void zentas_base(const TInitBundle&                   datain_ib,
                 size_t                               K,
                 const size_t* const                  indices_init,
                 std::string                          initialisation_method,
                 std::string                          algorithm,
                 size_t                               level,
                 size_t                               max_proposals,
                 bool                                 capture_output,
                 std::string&                         text,
                 size_t                               seed,
                 double                               max_time,
                 double                               min_mE,
                 double                               max_itok,
                 size_t* const                        indices_final,
                 size_t* const                        labels,
                 size_t                               nthreads,
                 size_t                               max_rounds,
                 bool                                 patient,
                 std::string                          energy,
                 bool                                 with_tests,
                 const typename TMetric::Initializer& metric_initializer,
                 const EnergyInitialiser&             energy_initialiser,
                 const std::chrono::time_point<std::chrono::high_resolution_clock>& bigbang,
                 bool do_balance_labels)
{

/* used during experiments to see if openblas worth the effort. Decided not.
//openblas_set_num_threads(1);
*/

#ifndef COMPILE_FOR_R
  std::stringstream buffer;
  auto              cout_buff = std::cout.rdbuf();
  if (capture_output == true)
  {
    std::cout.rdbuf(buffer.rdbuf());
  }
  std::ofstream nowhere;
#endif

  scrutinize_input_1(energy_initialiser, energy, K, algorithm, level, datain_ib.ndata);

  if (initialisation_method == "from_init_indices" && indices_init == nullptr)
  {
    throw zentas::zentas_error(
      R"(initialisation_method == "from_init_indices" && indices_init == nullptr) is true)");
  }

  dispatch<TData, TMetric>(algorithm,
                           level,
                           datain_ib,
                           K,
                           indices_init,
                           initialisation_method,
                           max_proposals,
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

#ifndef COMPILE_FOR_R
  if (capture_output == true)
  {
    text = buffer.str();
    std::cout.rdbuf(cout_buff);
  }

  else
  {
    text = "capture_output was false, so nothing here";
  }
#endif
}


}

#endif
