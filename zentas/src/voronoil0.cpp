// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/voronoil0.hpp>

namespace nszen
{

VoronoiL0::VoronoiL0(const SkeletonClustererInitBundle& sb, const ExtrasBundle& eb)
  : SkeletonClusterer(sb)
{
  (void)eb;
}

void VoronoiL0::set_redistribute_order(std::vector<size_t>& redistribute_order)
{
  std::iota(redistribute_order.begin(), redistribute_order.end(), 0);
}

void VoronoiL0::custom_initialise_refinement()
{

  // Voronoi does not have second nearest information, so we need to build that here.

  std::vector<std::vector<XNearestInfo>> nearest_2_infos;
  nearest_2_infos.resize(K);

  std::unique_ptr<double[]> up_distances(new double[K]);
  double                    min_distance          = std::numeric_limits<double>::max();
  size_t                    nearest_center        = 0;
  size_t                    second_nearest_center = 0;
  double                    second_min_distance   = std::numeric_limits<double>::max();

  size_t total_visits = 0;
  for (size_t k2 = 0; k2 < K; ++k2)
  {
    for (size_t i = 0; i < nearest_1_infos[k2].size(); ++i)
    {
      ++total_visits;
      min_distance          = std::numeric_limits<double>::max();
      second_min_distance   = std::numeric_limits<double>::max();
      second_nearest_center = 0;
      nearest_center        = 0;

      for (size_t k = 0; k < K; ++k)
      {
        set_center_sample_distance(k, k2, i, second_min_distance, up_distances[k]);
        if (up_distances[k] <= second_min_distance)
        {
          if (up_distances[k] < min_distance)
          {
            second_min_distance   = min_distance;
            second_nearest_center = nearest_center;
            nearest_center        = k;
            min_distance          = up_distances[k];
          }
          else
          {
            second_nearest_center = k;
            second_min_distance   = up_distances[k];
          }
        }
      }
      if (nearest_center != k2 || nearest_1_infos[k2][i].a_x != k2)
      {
        std::stringstream ss;
        ss << " i = " << i << "  nearest_center = " << nearest_center << "   k2 = " << k2
           << " a_x = " << nearest_1_infos[k2][i].a_x << "  no way jose";
        throw zentas::zentas_error(ss.str());
      }
      else
      {
        nearest_2_infos[k2].emplace_back(
          second_nearest_center, second_min_distance, f_energy(second_min_distance));
      }
    }
  }

  prd->initialise_from_n1n2(nearest_1_infos, nearest_2_infos);
}

void VoronoiL0::put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1,
                                                                    size_t k2,
                                                                    double d2,
                                                                    double e2)
{

  (void)k1;
  (void)k2;
  (void)d2;
  (void)e2;
}

std::string VoronoiL0::get_round_summary()
{
  std::stringstream ss;
  ss << get_base_summary_string();
  return ss.str();
}

void VoronoiL0::update_sample_info()
{
  double adistance{0};
  double min_distance{0};
  size_t min_k{0};

  for (size_t k = 0; k < K; ++k)
  {
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      min_distance = std::numeric_limits<double>::max();
      for (size_t kp = 0; kp < K; ++kp)
      {
        set_center_sample_distance(kp, k, j, min_distance, adistance);
        if (adistance < min_distance)
        {
          min_distance = adistance;
          min_k        = kp;
        }
      }
      reset_nearest_info(k, j, min_k, min_distance, f_energy(min_distance));
    }
  }
}

bool VoronoiL0::update_centers()
{

  bool   modified = false;
  double E_old;
  double E_prop;
  double E_prop_best;
  size_t j_prop_best;

  for (size_t k = 0; k < K; ++k)
  {

    E_old = get_cluster_energy(k);

    // compute energies with all other centers
    E_prop_best = std::numeric_limits<double>::max();
    j_prop_best = 0;

    for (size_t j_prop = 0; j_prop < get_ndata(k); ++j_prop)
    {
      E_prop = get_e1(k, j_prop);
      for (size_t j = 0; j < get_ndata(k); ++j)
      {
        E_prop += f_energy(get_sample_sample_distance_nothreshold(k, j_prop, j));
      }
      if (E_prop < E_prop_best)
      {
        E_prop_best = E_prop;
        j_prop_best = j_prop;
      }
    }

    if (E_old > E_prop_best)
    {
      modified = true;
      swap_center_with_sample(k, j_prop_best);
    }
  }

  return modified;
}
}
