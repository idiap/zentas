// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <memory>
#include <tuple>
#include <zentas/skeletonclusterer.hpp>
namespace nszen
{

void SkeletonClusterer::rf_update_sample_info_yinyang_v1()
{

  std::vector<double> min_d1(prd->rf_n_groups, std::numeric_limits<double>::max());
  std::vector<size_t> min_k1(prd->rf_n_groups);

  std::vector<double> min_d2(prd->rf_n_groups, std::numeric_limits<double>::max());
  size_t              a_min_k2;

  size_t min_g_global;

  // will always have start_g first, the rest in increasing order.
  std::vector<size_t> groups_in_order(prd->rf_n_groups);
  std::iota(groups_in_order.begin(), groups_in_order.end(), 0);

  for (size_t start_g = 0; start_g < prd->rf_n_groups; ++start_g)
  {

    std::swap(groups_in_order[0], groups_in_order[start_g]);
    for (size_t k = prd->cum_in_group[start_g]; k < prd->cum_in_group[start_g + 1]; ++k)
    {
      for (size_t j = 0; j < get_ndata(k); ++j)
      {
        size_t ID             = prd->glt_ID[k][j];
        double global_l_bound = std::numeric_limits<double>::max();

        for (auto& g : groups_in_order)
        {
          prd->l_gps[ID * prd->rf_n_groups + g] -= prd->max_delta_group[g];
          global_l_bound = std::min(prd->l_gps[ID * prd->rf_n_groups + g], global_l_bound);
        }

        prd->upper_1[k][j] += prd->delta_C[k];
        if (global_l_bound < prd->upper_1[k][j])
        {

          set_rf_center_sample_distance_nothreshold(k, k, j, prd->upper_1[k][j]);
          if (global_l_bound < prd->upper_1[k][j])
          {

            min_d1[start_g] = prd->upper_1[k][j];
            min_k1[start_g] = k;
            min_d2[start_g] = -std::numeric_limits<double>::max();
            min_g_global    = start_g;

            for (auto& g : groups_in_order)
            {
              if (prd->l_gps[ID * prd->rf_n_groups + g] < prd->upper_1[k][j])
              {
                rf_set_2_smallest_group(k, j, g, min_k1[g], min_d1[g], a_min_k2, min_d2[g]);
                prd->l_gps[ID * prd->rf_n_groups + g] = min_d1[g];
                prd->t_gps[ID * prd->rf_n_groups + g] = prd->rf_round;
                if (min_d1[g] < prd->upper_1[k][j])
                {
                  prd->upper_1[k][j] = min_d1[g];
                  min_g_global       = g;
                }
              }
            }

            if (min_g_global != start_g)
            {
              prd->l_gps[ID * prd->rf_n_groups + start_g]      = min_d1[start_g];
              prd->l_gps[ID * prd->rf_n_groups + min_g_global] = min_d2[min_g_global];
            }

            // the group of the nearest is unchanged.
            // if we computed the distance to second nearest in this group, it should be the new
            // lower bound.
            // otherwise the lower bound should remain the same.
            else
            {
              prd->l_gps[ID * prd->rf_n_groups + start_g] =
                std::max(prd->l_gps[ID * prd->rf_n_groups + start_g], min_d2[start_g]);
            }

            if (min_k1[min_g_global] != k)
            {
              reset_nearest_info(
                k, j, min_k1[min_g_global], min_d1[min_g_global], f_energy(min_d1[min_g_global]));
            }
          }
        }
      }
    }
  }
}

// attempt at version in eakmeans
void SkeletonClusterer::rf_update_sample_info_yinyang_v2()
{

  double global_l_bound;

  size_t global_min_g;
  size_t global_min_k;

  double min_d1, min_d2;
  size_t min_k1, min_k2;
  size_t ID;

  for (size_t start_g = 0; start_g < prd->rf_n_groups; ++start_g)
  {
    for (size_t k = prd->cum_in_group[start_g]; k < prd->cum_in_group[start_g + 1]; ++k)
    {
      for (size_t j = 0; j < get_ndata(k); ++j)
      {

        ID             = prd->glt_ID[k][j];
        global_l_bound = std::numeric_limits<double>::max();
        for (size_t g = 0; g < prd->rf_n_groups; ++g)
        {
          prd->l_gps[ID * prd->rf_n_groups + g] -= prd->max_delta_group[g];
          global_l_bound = std::min(prd->l_gps[ID * prd->rf_n_groups + g], global_l_bound);
        }
        prd->upper_1[k][j] += prd->delta_C[k];

        if (global_l_bound < prd->upper_1[k][j])
        {
          set_rf_center_sample_distance_nothreshold(k, k, j, prd->upper_1[k][j]);

          if (global_l_bound < prd->upper_1[k][j])
          {
            global_min_g = start_g;
            global_min_k = k;

            for (size_t g = 0; g < prd->rf_n_groups; ++g)
            {
              rf_set_2_smallest_group(k, j, g, min_k1, min_d1, min_k2, min_d2);

              if (g != global_min_g)
              {
                if (min_d1 < prd->upper_1[k][j])
                {
                  if (g < global_min_g)
                  {
                    prd->l_gps[ID * prd->rf_n_groups + global_min_g] = std::min(
                      prd->l_gps[ID * prd->rf_n_groups + global_min_g], prd->upper_1[k][j]);
                  }
                  else
                  {
                    prd->l_gps[ID * prd->rf_n_groups + global_min_g] = prd->upper_1[k][j];
                  }
                  global_min_g                          = g;
                  global_min_k                          = min_k1;
                  prd->l_gps[ID * prd->rf_n_groups + g] = min_d2;
                  prd->upper_1[k][j]                    = min_d1;
                }
                else
                {
                  prd->l_gps[ID * prd->rf_n_groups + g] = min_d1;
                }
              }
              else
              {
                global_min_k                          = min_k1;
                prd->upper_1[k][j]                    = min_d1;
                prd->l_gps[ID * prd->rf_n_groups + g] = min_d2;
              }
            }
            reset_nearest_info(k, j, global_min_k, prd->upper_1[k][j], f_energy(global_l_bound));
          }
        }
      }
    }
  }
}
}
