// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <memory>
#include <tuple>
#include <zentas/skeletonclusterer.hpp>
namespace nszen
{

bool ExponionData::is_exponion() { return true; }

void ExponionData::specific_custom_append(size_t k_new, size_t k, size_t j)
{
  v_b[k_new].push_back(v_b[k][j]);
}

void ExponionData::specific_custom_replace_with_last(size_t k, size_t j)
{
  v_b[k][j] = v_b[k].back();
}

void ExponionData::specific_custom_remove_last(size_t k) { v_b[k].pop_back(); }

void ExponionData::specific_set_n_groups(size_t K, size_t ndata)
{
  (void)ndata;
  rf_n_groups = K / 10;
}

void ExponionData::specific_final_initialise_memory() {}

void ExponionData::specific_initialise_from_n1n2(const XNInfos& nearest_1_infos,
                                                 const XNInfos& nearest_2_infos)
{
  (void)nearest_1_infos;
  v_b.resize(K);
  for (size_t k = 0; k < K; ++k)
  {
    size_t ndata_k = nearest_2_infos[k].size();
    v_b[k].resize(ndata_k);
    for (size_t j = 0; j < nearest_2_infos[k].size(); ++j)
    {
      v_b[k][j] = nearest_2_infos[k][j].a_x;
    }
  }
  delta_C.resize(K, 0);
  u1_C.resize(K, std::numeric_limits<double>::max());
}

void ExponionData::specific_swap(size_t k1, size_t k2) { std::swap(v_b[k1], v_b[k2]); }

void ExponionData::specific_update_centers() {}

////////////////////////

bool YinyangData::is_exponion() { return false; }

void YinyangData::specific_custom_append(size_t k_new, size_t k, size_t j)
{
  glt_ID[k_new].push_back(glt_ID[k][j]);
}

void YinyangData::specific_custom_remove_last(size_t k) { glt_ID[k].pop_back(); }

void YinyangData::specific_custom_replace_with_last(size_t k, size_t j)
{
  glt_ID[k][j] = glt_ID[k].back();
}

void YinyangData::specific_initialise_from_n1n2(const XNInfos& nearest_1_infos,
                                                const XNInfos& nearest_2_infos)
{
  (void)nearest_1_infos;
  (void)nearest_2_infos;
}

void YinyangData::specific_set_n_groups(size_t K, size_t ndata)
{
  // K/10 if memory will not exceed 2 GB,
  double max_ngbs = 2.0;
  // otherwise (1e9 / (ndata*sizeof(double)))
  rf_n_groups = std::min<size_t>(
    K / 9, static_cast<size_t>(max_ngbs * 1e9 / ((sizeof(size_t) + sizeof(double)) * ndata)));
}

void YinyangData::specific_swap(size_t k1, size_t k2)
{
  (void)k1;
  (void)k2;
}

void YinyangData::specific_final_initialise_memory()
{
  // initialise group bounds (for data)
  l_gps.resize(ndata * rf_n_groups);

  // TODO : currently this is not updated or used
  t_gps.resize(ndata * rf_n_groups);

  size_t running_ID = 0;
  glt_ID.resize(K);
  for (size_t k = 0; k < K; ++k)
  {
    glt_ID[k].resize(get_ndata(k));
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      glt_ID[k][j] = running_ID;
      ++running_ID;
    }
  }
}

void YinyangData::specific_update_centers() {}

///////////////////

void RefinementData::custom_append(size_t k_new, size_t k, size_t j)
{
  lower_2[k_new].push_back(lower_2[k][j]);
  upper_1[k_new].push_back(upper_1[k][j]);
  specific_custom_append(k_new, k, j);
}

void RefinementData::custom_replace_with_last(size_t k, size_t j)
{
  lower_2[k][j] = lower_2[k].back();
  upper_1[k][j] = upper_1[k].back();
  specific_custom_replace_with_last(k, j);
}

void RefinementData::custom_remove_last(size_t k)
{
  lower_2[k].pop_back();
  upper_1[k].pop_back();
  specific_custom_remove_last(k);
}

void RefinementData::swap(size_t k1, size_t k2)
{
  std::swap(lower_2[k1], lower_2[k2]);
  std::swap(upper_1[k1], upper_1[k2]);
  specific_swap(k1, k2);
}

void RefinementData::initialise_from_n1n2(const XNInfos& nearest_1_infos,
                                          const XNInfos& nearest_2_infos)
{
  lower_2.resize(K);
  upper_1.resize(K);
  for (size_t k = 0; k < K; ++k)
  {
    size_t ndata_k = nearest_1_infos[k].size();
    lower_2[k].resize(ndata_k);
    upper_1[k].resize(ndata_k);
    for (size_t j = 0; j < nearest_2_infos[k].size(); ++j)
    {
      lower_2[k][j] = nearest_2_infos[k][j].d_x;
      upper_1[k][j] = nearest_1_infos[k][j].d_x;
    }
  }
  delta_C.resize(K, 0);
  specific_initialise_from_n1n2(nearest_1_infos, nearest_2_infos);
}

size_t RefinementData::get_ndata(size_t k)
{
  // will break for lower_2 not used.
  return lower_2[k].size();
}

void SkeletonClusterer::rf_swap(size_t k1, size_t k2)
{
  swap_data(k1, k2);
  std::swap(nearest_1_infos[k1], nearest_1_infos[k2]);
  std::swap(sample_IDs[k1], sample_IDs[k2]);
}

void SkeletonClusterer::initialise_refinement()
{

  if (rf_alg == "yinyang")
  {
    prd.reset(new YinyangData);
  }

  else if (rf_alg == "exponion")
  {
    prd.reset(new ExponionData);
  }

  else
  {
    throw zentas::zentas_error("either yinyang or exponion should be chosen as rf_alg");
  }

  prd->K     = K;
  prd->ndata = ndata;
  prd->specific_set_n_groups(K, ndata);

  /* put centers into clusters, initialise rf_center and old_rf_center */
  for (size_t k = 0; k < K; ++k)
  {
    put_sample_in_cluster(center_IDs[k]);
    append_zero_to_rf_center_data();
    append_zero_to_rf_sum_data();
    append_zero_to_old_rf_center_data();
  }

  custom_initialise_refinement();

  // cluster the centers
  prd->groups = get_subclustered_centers_labels(prd->rf_n_groups);

  if (prd->groups.size() != K)
  {
    throw zentas::zentas_error(
      "groups.size() != K just after center clustering, this is strange (wrong), logic error");
  }

  if (*std::max_element(prd->groups.begin(), prd->groups.end()) >= prd->rf_n_groups)
  {
    throw zentas::zentas_error("max group index too large, logic error");
  }

  prd->n_in_group.resize(prd->rf_n_groups, 0);
  for (size_t k = 0; k < K; ++k)
  {
    ++prd->n_in_group[prd->groups[k]];
  }

  prd->cum_in_group.resize(prd->rf_n_groups + 1, 0);
  for (size_t g = 0; g < prd->rf_n_groups; ++g)
  {
    prd->cum_in_group[g + 1] = prd->cum_in_group[g] + prd->n_in_group[g];
  }

  // sort the clusters by group
  size_t p_in  = 0;
  size_t p_out = 0;

  for (size_t g = 0; g < prd->rf_n_groups; ++g)
  {

    if (p_in != prd->cum_in_group[g])
    {
      throw zentas::zentas_error("This does not make sense, p_in should be auto here");
    }

    p_out = prd->cum_in_group[g + 1];

    while (p_in < prd->cum_in_group[g] || p_out < K)
    {
      if (prd->groups[p_in] == g)
      {
        ++p_in;
      }

      else if (prd->groups[p_out] != g)
      {
        ++p_out;
      }

      else
      {
        rf_swap(p_in, p_out);
        prd->swap(p_in, p_out);
        std::swap(prd->groups[p_in], prd->groups[p_out]);
      }
    }
  }

  /* also clear-up unnec. mem */
  custom_rf_clear_initmem();

  prd->max_delta_group.resize(prd->rf_n_groups);
  prd->l_gC.resize(K * prd->rf_n_groups, std::numeric_limits<double>::max());
  for (size_t k = 0; k < K; ++k)
  {
    // initialise sum
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      rf_increment_sum(k, j);
    }
  }

  for (size_t k = 0; k < K; ++k)
  {
    // initialise group bounds (for centers)
    double adist;
    for (size_t g = 0; g < prd->rf_n_groups; ++g)
    {
      for (size_t kp = prd->cum_in_group[g]; kp < prd->cum_in_group[g + 1]; ++kp)
      {
        set_rf_center_center_distance(k, kp, prd->l_gC[k * prd->rf_n_groups + g], adist);
        prd->l_gC[k * prd->rf_n_groups + g] = std::min(prd->l_gC[k * prd->rf_n_groups + g], adist);
      }
    }
  }

  prd->rf_round = 1;

  prd->specific_final_initialise_memory();
}

void SkeletonClusterer::rf_tighten_nearest(size_t k)
{
  double min_d1;
  for (size_t j = 0; j < get_ndata(k); ++j)
  {
    set_rf_center_sample_distance_nothreshold(k, k, j, min_d1);
    reset_nearest_info(k, j, k, min_d1, f_energy(min_d1));
    prd->upper_1[k][j] = min_d1;
  }
}

void SkeletonClusterer::run_refinement()
{

  while (halt_refinement() == false)
  {

    /* ************** *
    * UPDATE CENTERS *
    * ****************/
    set_t_ncalcs_update_centers_start();
    bool modified_centers = rf_update_centers();

    // no change, prepare to break.
    if (modified_centers == false)
    {
      update_t_ncalcs_center_end();
      break;
    }

    /* ************************* *
    * UPDATE CENTER CENTER INFO *
    * ***************************/
    rf_update_center_center_info();
    if (with_tests == true)
      rf_post_center_update_test();
    update_t_ncalcs_center_end();

    for (auto& x : to_leave_cluster)
    {
      x.resize(0);
    }

    /* ****************** *
    * UPDATE SAMPLE INFO *
    * ********************/
    rf_update_sample_info();
    if (with_tests == true)
      rf_post_sample_update_test();
    update_t_ncalcs_sample_update_end();

    /* ************ *
    * REDISTRIBUTE *
    * **************/
    rf_redistribute();
    if (with_tests == true)
      rf_post_redistribute_test();
    update_t_redistibute_end();

    /* ************************* *
    * UPDATE CLUSTER STATISTICS *
    * ***************************/
    if (is_rf_correct_d1_round() == true)
    {
      rf_update_energies();
      if (with_tests == true)
        rf_post_update_statistics_test();
    }
    else
    {
      E_total = -1;
    }
    update_t_update_all_cluster_stats_end();

    if (with_tests == true)
      mowri << zentas::Endl;
    mowri << rf_get_round_summary() << zentas::Endl;

    ++round;
    ++prd->rf_round;
  }

  for (size_t k = 0; k < K; ++k)
  {
    rf_tighten_nearest(k);
  }
  rf_update_energies();

  output_halt_refinement_reason();
}

bool SkeletonClusterer::is_rf_correct_d1_round()
{
  return (prd->rf_round < 5 || prd->rf_round % 26 == 0);
}

bool SkeletonClusterer::is_rf_tighten_cluster_radius_round()
{
  return (prd->rf_round < 5 || prd->rf_round % 15 == 0);
}

bool SkeletonClusterer::rf_update_centers()
{

  bool   has_changed = false;
  double sum_abs;

  for (size_t g = 0; g < prd->rf_n_groups; ++g)
  {
    prd->max_delta_group[g] = 0;
    for (size_t k = prd->cum_in_group[g]; k < prd->cum_in_group[g + 1]; ++k)
    {
      set_old_rf_center_data(k);
      set_rf_center_data(k);
      set_delta_rf_new_and_old(k, std::numeric_limits<double>::max(), prd->delta_C[k]);
      prd->max_delta_group[g] = std::max(prd->max_delta_group[g], prd->delta_C[k]);
      set_sum_abs_rf_new_and_old(k, sum_abs);
      cluster_has_changed[k] = (prd->delta_C[k] > 1e-7 * sum_abs);
      has_changed            = cluster_has_changed[k] == true ? true : has_changed;
    }
  }

  prd->specific_update_centers();

  return has_changed;
}

// An exponion-LIKE algorithm
void SkeletonClusterer::rf_update_sample_info_exponion()
{

  size_t min_k1, min_k2;
  double min_d1, min_d2;
  double worker_dist;

  // the mimimum distance to a center other than `self'.
  double min_cc;

  // sort in ascending order by first (double) argument
  auto tupsorter = [](std::tuple<double, size_t>& x, std::tuple<double, size_t>& y) {
    return std::get<0>(x) < std::get<0>(y);
  };

  //  not all inter-centroid distances are required.
  //  these vectors only store those which are required.
  //  ie there exists a sample with a sufficiently large compl_exp_rad.
  size_t              n_cc_required;
  std::vector<double> cc_all(K);
  std::vector<std::tuple<double, size_t>> cc_required_tuple(K);
  std::vector<size_t> cc_required_indices(K);
  std::vector<double> cc_required_distances(K);

  for (size_t k = 0; k < K; ++k)
  {

    // update upper bound on distance to furthest member
    if (is_rf_tighten_cluster_radius_round() == true)
    {
      prd->u1_C[k] = 0;
      for (size_t j = 0; j < get_ndata(k); ++j)
      {
        prd->u1_C[k] = std::max(prd->u1_C[k], prd->upper_1[k][j]);
      }
    }
    prd->u1_C[k] += prd->delta_C[k];

    // update lower bound on distance to nearest member of each group
    for (size_t g = 0; g < prd->rf_n_groups; ++g)
    {
      prd->l_gC[k * prd->rf_n_groups + g] -= prd->max_delta_group[g];
      prd->l_gC[k * prd->rf_n_groups + g] -= prd->delta_C[k];
    }

    min_cc = std::numeric_limits<double>::max();
    std::vector<size_t> groups_in_vicinity;

    // in own group, set distance to all centers other than self
    for (size_t kp = prd->cum_in_group[prd->groups[k]]; kp < prd->cum_in_group[prd->groups[k] + 1];
         ++kp)
    {
      if (kp != k)
      {
        set_rf_center_center_distance(k, kp, std::numeric_limits<double>::max(), cc_all[kp]);
        min_cc = std::min(min_cc, cc_all[kp]);
      }
    }
    groups_in_vicinity.push_back(prd->groups[k]);

    // we want to determine which centers are within the maximum exponion radius :
    // 2*(max distance to an element of k) + (min distance to another center from k) (**)
    // we have (1) a LOWER bound on the distance on the centers in each group to center k, and
    //         (2) a UPPER bound on (**) of 2*prd->u1_C[k] + (min_cc so far) +

    for (size_t g = 0; g < prd->rf_n_groups; ++g)
    {
      if ((g != prd->groups[k]) &&
          (prd->l_gC[k * prd->rf_n_groups + g] < 2 * prd->u1_C[k] + min_cc))
      {  // do not need to add s to make max_delta_possible_second correct.
        groups_in_vicinity.push_back(g);
        prd->l_gC[k * prd->rf_n_groups + g] = std::numeric_limits<double>::max();
        for (size_t kp = prd->cum_in_group[g]; kp < prd->cum_in_group[g + 1]; ++kp)
        {
          set_rf_center_center_distance(k, kp, std::numeric_limits<double>::max(), cc_all[kp]);
          min_cc = std::min(min_cc, cc_all[kp]);
          prd->l_gC[k * prd->rf_n_groups + g] =
            std::min(prd->l_gC[k * prd->rf_n_groups + g], cc_all[kp]);
        }
      }
    }

    n_cc_required                    = 0;
    double max_delta_possible_second = 0;
    for (auto& g : groups_in_vicinity)
    {
      for (size_t kp = prd->cum_in_group[g]; kp < prd->cum_in_group[g + 1]; ++kp)
      {
        if ((kp != k) && (cc_all[kp] < 2 * prd->u1_C[k] + min_cc))
        {
          cc_required_tuple[n_cc_required] = std::make_tuple(cc_all[kp], kp);
          ++n_cc_required;
          max_delta_possible_second = std::max(max_delta_possible_second, prd->delta_C[kp]);
        }
      }
    }

    cc_required_tuple[n_cc_required] =
      std::make_tuple(std::numeric_limits<double>::max(), std::numeric_limits<size_t>::max());
    ++n_cc_required;

    std::sort(cc_required_tuple.begin(), cc_required_tuple.begin() + n_cc_required, tupsorter);

    for (size_t jp = 0; jp < n_cc_required; ++jp)
    {
      cc_required_distances[jp] = std::get<0>(cc_required_tuple[jp]);
      cc_required_indices[jp]   = std::get<1>(cc_required_tuple[jp]);
    }

    double compl_exp_rad_x;
    double exprad;
    size_t kp;

    if (is_rf_correct_d1_round() == true)
    {
      rf_tighten_nearest(k);  //_and_u1
    }

    else
    {
      if (prd->delta_C[k] > 0)
      {
        for (size_t j = 0; j < get_ndata(k); ++j)
        {
          prd->upper_1[k][j] += prd->delta_C[k];
        }
      }
    }

    for (size_t j = 0; j < get_ndata(k); ++j)
    {

      prd->lower_2[k][j] -= max_delta_possible_second;  // this is better than max_prd->delta_C.
                                                        // but is it correct? Yes, but TODO the
                                                        // name is wrong
      if (prd->lower_2[k][j] < prd->upper_1[k][j])
      {

        set_rf_center_sample_distance_nothreshold(k, k, j, min_d1);
        reset_nearest_info(k, j, k, min_d1, f_energy(min_d1));

        if (prd->lower_2[k][j] < prd->upper_1[k][j])
        {

          set_rf_center_sample_distance_nothreshold(prd->v_b[k][j], k, j, min_d2);
          compl_exp_rad_x = (min_d2 < min_d1 ? 2 * min_d1 : min_d1 + min_d2);
          exprad          = (1 + 1e-6) * std::min(2 * min_d1 + min_cc, compl_exp_rad_x);
          min_k1          = k;
          min_k2          = prd->v_b[k][j];

          if (min_d1 > min_d2)
          {
            std::swap(min_d1, min_d2);
            std::swap(min_k1, min_k2);
          }

          // get insertion index:
          size_t insertion_index = 0;
          while (exprad > cc_required_distances[insertion_index])
          {
            ++insertion_index;
          }

          for (size_t kpi = 0; kpi < insertion_index; ++kpi)
          {
            kp = cc_required_indices[kpi];
            if (kp != prd->v_b[k][j])
            {

              set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);

              if (worker_dist < min_d1)
              {
                min_d2 = min_d1;
                min_k2 = min_k1;
                min_d1 = worker_dist;
                min_k1 = kp;
              }

              else if (worker_dist < min_d2)
              {
                min_d2 = worker_dist;
                min_k2 = kp;
              }
            }
          }

          reset_nearest_info(k, j, min_k1, min_d1, f_energy(min_d1));
          prd->lower_2[k][j] = min_d2;
          prd->v_b[k][j]     = min_k2;
          prd->upper_1[k][j] = min_d1;

          if (min_d1 > prd->u1_C[min_k1])
          {
            prd->u1_C[min_k1] = min_d1;
          }
        }
      }
    }
  }
}

void SkeletonClusterer::update_lt_pgs(size_t ID, size_t g, double lv)
{
  prd->l_gps[ID * prd->rf_n_groups + g] = lv;
}

double SkeletonClusterer::rf_bound(size_t ID, size_t g)
{
  return prd->l_gps[ID * prd->rf_n_groups + g];
}

// TODO : this is half as fast as ns- version of eakmeans. work on this. some ideas:
// (1) ns bounding
// (2) tests within group bound failure, using center-specific distance moved
// TODO : comment on the current code, how it works

void SkeletonClusterer::rf_update_sample_info_yinyang()
{

  std::vector<double> min_d1(prd->rf_n_groups, std::numeric_limits<double>::max());
  std::vector<size_t> min_k1(prd->rf_n_groups);
  std::vector<double> min_d2(prd->rf_n_groups, std::numeric_limits<double>::max());
  // size_t a_min_k2;
  size_t min_g_global;

  std::vector<size_t> other_groups(prd->rf_n_groups - 1);
  std::iota(other_groups.begin(), other_groups.end(), 1);

  auto tupsorter = [](std::tuple<double, size_t>& x, std::tuple<double, size_t>& y) {
    return std::get<0>(x) < std::get<0>(y);
  };
  std::vector<std::tuple<double, size_t>> v_dk(K);
  std::vector<std::vector<size_t>> group_ks(prd->rf_n_groups);

  for (size_t g = 0; g < prd->rf_n_groups; ++g)
  {
    group_ks[g].resize(prd->cum_in_group[g + 1] - prd->cum_in_group[g]);
  }

  for (size_t start_g = 0; start_g < prd->rf_n_groups; ++start_g)
  {
    for (size_t k = prd->cum_in_group[start_g]; k < prd->cum_in_group[start_g + 1]; ++k)
    {

      if (is_rf_correct_d1_round() == true)
      {
        rf_tighten_nearest(k);
      }

      for (size_t kp = prd->cum_in_group[start_g]; kp < prd->cum_in_group[start_g + 1]; ++kp)
      {
        set_rf_center_center_distance(
          k, kp, std::numeric_limits<double>::max(), std::get<0>(v_dk[kp]));
        std::get<1>(v_dk[kp]) = kp;
      }

      std::sort(v_dk.begin() + prd->cum_in_group[start_g],
                v_dk.begin() + prd->cum_in_group[start_g + 1],
                tupsorter);
      for (size_t ki = 0; ki < prd->cum_in_group[start_g + 1] - prd->cum_in_group[start_g]; ++ki)
      {
        group_ks[start_g][ki] = std::get<1>(v_dk[prd->cum_in_group[start_g] + ki]);
      }

      for (size_t j = 0; j < get_ndata(k); ++j)
      {

        size_t ID             = prd->glt_ID[k][j];
        double global_l_bound = std::numeric_limits<double>::max();

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

            min_d1[start_g] = prd->upper_1[k][j];
            min_k1[start_g] = k;
            min_d2[start_g] = -std::numeric_limits<double>::max();
            min_g_global    = start_g;

            if (rf_bound(ID, start_g) < prd->upper_1[k][j])
            {
              rf_set_2_smallest_ws(k,
                                   j,
                                   group_ks[start_g].data() + 1,
                                   group_ks[start_g].size() - 1,
                                   min_k1[start_g],
                                   min_d1[start_g],
                                   min_d2[start_g]);
              update_lt_pgs(ID, start_g, min_d1[start_g]);
              if (min_d1[start_g] < prd->upper_1[k][j])
              {
                prd->upper_1[k][j] = min_d1[start_g];
              }
            }

            for (auto& g : other_groups)
            {
              if (rf_bound(ID, g) < prd->upper_1[k][j])
              {
                // like  60% of the time is spent here :
                rf_set_2_smallest_group(k, j, g, min_k1[g], min_d1[g], min_d2[g]);

                update_lt_pgs(ID, g, min_d1[g]);
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
            // if we computed the distance to second nearest in this group,
            // it should be the new lower bound.
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
    if (start_g < prd->rf_n_groups - 1)
    {
      other_groups[start_g] -= 1;
    }
  }
}

void SkeletonClusterer::rf_update_sample_info()
{
  if (prd->is_exponion())
  {
    rf_update_sample_info_exponion();
  }
  else
  {
    rf_update_sample_info_yinyang();
  }
}

void SkeletonClusterer::rf_redistribute()
{

  size_t k_new;
  size_t j;

  for (size_t k = 0; k < K; ++k)
  {
    for (size_t ji = to_leave_cluster[k].size(); ji-- > 0;)
    {

      j     = to_leave_cluster[k][ji];
      k_new = nearest_1_infos[k][j].a_x;

      rf_cumulative_correction(k_new, k, j);  // correct sums or whatever.

      nearest_1_infos[k_new].push_back(nearest_1_infos[k][j]);
      sample_IDs[k_new].push_back(sample_IDs[k][j]);
      append_across(k_new, k, j);
      prd->custom_append(k_new, k, j);

      // this is `remove_with_tail_pull'
      if (j != get_ndata(k) - 1)
      {
        nearest_1_infos[k][j] = *(nearest_1_infos[k].end() - 1);
        sample_IDs[k][j]      = *(sample_IDs[k].end() - 1);
        replace_with_last_element(k, j);
        prd->custom_replace_with_last(k, j);
      }
      nearest_1_infos[k].pop_back();
      sample_IDs[k].pop_back();
      remove_last(k);
      prd->custom_remove_last(k);
    }
  }
}

void SkeletonClusterer::rf_update_center_center_info()
{
  // TODO
}

std::string SkeletonClusterer::rf_get_round_summary()
{
  return get_base_summary_string() + "---";  // get_round_summary();
}

void SkeletonClusterer::rf_post_center_update_test()
{
  std::cout << "TODO : rf_post_center_update_test" << std::endl;
}

void SkeletonClusterer::rf_post_sample_update_test()
{
  std::cout << "TODO : rf_post_sample_update_test" << std::endl;
}

void SkeletonClusterer::rf_post_redistribute_test()
{
  std::cout << "TODO : rf_post_redistribute_test" << std::endl;
}

void SkeletonClusterer::rf_post_update_statistics_test()
{
  std::cout << "TODO : rf_post_update_statistics_test" << std::endl;
}

// TODO : any custom cluster stats here ?
void SkeletonClusterer::rf_update_energies()
{
  old_E_total = E_total;
  E_total     = 0;
  for (size_t k = 0; k < K; ++k)
  {
    cluster_energies[k] = 0;
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      cluster_energies[k] += nearest_1_infos[k][j].e_x;
    }
    cluster_mean_energies[k] = cluster_energies[k] / static_cast<double>(get_ndata(k));
    E_total += cluster_energies[k];
  }
}

void SkeletonClusterer::rf_set_2_smallest(
  size_t k, size_t j, size_t& min_k1, double& min_d1, size_t& min_k2, double& min_d2)
{

  min_k1 = 0;
  min_k2 = 0;
  min_d1 = std::numeric_limits<double>::max();
  min_d2 = std::numeric_limits<double>::max();

  double worker_dist;
  for (size_t kp = 0; kp < K; ++kp)
  {
    set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
    if (worker_dist < min_d2)
    {
      if (worker_dist < min_d1)
      {
        min_d2 = min_d1;
        min_k2 = min_k1;
        min_d1 = worker_dist;
        min_k1 = kp;
      }
      else
      {
        min_d2 = worker_dist;
        min_k2 = kp;
      }
    }
  }
}

// assumes min_k1, min_d1 are already valid.
void SkeletonClusterer::rf_set_2_smallest_ws(
  size_t k, size_t j, const size_t* ks, size_t nks, size_t& min_k1, double& min_d1, double& min_d2)
{

  double worker_dist;
  min_d2 = std::numeric_limits<double>::max();
  size_t kp;

  for (size_t ki = 0; ki < nks; ++ki)
  {
    kp = ks[ki];
    set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
    if (worker_dist < min_d2)
    {
      if (worker_dist < min_d1)
      {
        min_d2 = min_d1;
        min_d1 = worker_dist;
        min_k1 = kp;
      }
      else
      {
        min_d2 = worker_dist;
      }
    }
  }
}

void SkeletonClusterer::rf_set_2_smallest_group(
  size_t k, size_t j, size_t g, size_t& min_k1, double& min_d1, size_t& min_k2, double& min_d2)
{

  min_k1 = 0;
  min_k2 = 0;
  min_d1 = std::numeric_limits<double>::max();
  min_d2 = std::numeric_limits<double>::max();

  double worker_dist;
  for (size_t kp = prd->cum_in_group[g]; kp < prd->cum_in_group[g + 1]; ++kp)
  {
    set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
    if (worker_dist < min_d2)
    {
      if (worker_dist < min_d1)
      {
        min_d2 = min_d1;
        min_k2 = min_k1;
        min_d1 = worker_dist;
        min_k1 = kp;
      }
      else
      {
        min_d2 = worker_dist;
        min_k2 = kp;
      }
    }
  }
}

void SkeletonClusterer::rf_set_2_smallest_group(
  size_t k, size_t j, size_t g, size_t& min_k1, double& min_d1, double& min_d2)
{

  min_k1 = 0;
  min_d1 = std::numeric_limits<double>::max();
  min_d2 = std::numeric_limits<double>::max();

  double worker_dist;
  for (size_t kp = prd->cum_in_group[g]; kp < prd->cum_in_group[g + 1]; ++kp)
  {
    set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
    if (worker_dist < min_d2)
    {
      if (worker_dist < min_d1)
      {
        min_d2 = min_d1;
        min_d1 = worker_dist;
        min_k1 = kp;
      }
      else
      {
        min_d2 = worker_dist;
      }
    }
  }
}

// Standard algorithm.
void SkeletonClusterer::rf_update_sample_info_standard()
{

  size_t min_k1, min_k2;
  double min_d1, min_d2;

  for (size_t k = 0; k < K; ++k)
  {
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      rf_set_2_smallest(k, j, min_k1, min_d1, min_k2, min_d2);
      reset_nearest_info(k, j, min_k1, min_d1, f_energy(min_d1));
    }
  }
}

void SkeletonClusterer::set_rf_center_sample_distance_nothreshold(size_t  k,
                                                                  size_t  k1,
                                                                  size_t  j1,
                                                                  double& distance)
{
  set_rf_center_sample_distance(k, k1, j1, std::numeric_limits<double>::max(), distance);
}

std::vector<size_t> SkeletonClusterer::get_subclustered_centers_labels(size_t sub_K)
{

  mowri << "clustering the centers ... " << zentas::Flush;
  auto sub_bigbang = std::chrono::high_resolution_clock::now();

  const size_t* sub_indices_init          = nullptr;
  std::string   sub_initialisation_method = "kmeans++-5";
  double        sub_max_time =
    SkeletonClusterer::time_total / (1000 * 1000 * 25.);  // spent 1/25th of time so far.
  double              sub_min_mE        = 0;
  double              sub_max_itok      = 1e8;
  size_t              sub_max_rounds    = 10000000;
  size_t              sub_max_proposals = 10000000;
  bool                sub_patient       = true;
  size_t              sub_nthreads      = 1;
  size_t              sub_seed          = 1011;
  std::string         sub_energy        = "cubic";
  bool                sub_with_tests    = SkeletonClusterer::with_tests;
  std::vector<size_t> sub_v_indices_final(sub_K);
  std::vector<size_t> sub_v_labels(K);

  std::string sub_algorithm("clarans");
  size_t      sub_level          = 3;
  bool        sub_capture_output = true;
  std::string sub_output_text;

  bool sub_do_balance_labels = true;

  perform_subclustering(sub_K,
                        sub_indices_init,
                        sub_initialisation_method,
                        sub_algorithm,
                        sub_level,
                        sub_max_proposals,
                        sub_capture_output,
                        sub_output_text,
                        sub_seed,
                        sub_max_time,
                        sub_min_mE,
                        sub_max_itok,
                        sub_v_indices_final.data(),
                        sub_v_labels.data(),
                        sub_nthreads,
                        sub_max_rounds,
                        sub_patient,
                        sub_energy,
                        sub_with_tests,
                        sub_bigbang,
                        sub_do_balance_labels);

  mowri << "done, the final line was:" << zentas::Endl;

  auto firstx   = sub_output_text.find_last_of("R");
  auto frag     = sub_output_text.substr(firstx);
  auto lastx    = frag.find("\n");
  auto lastline = frag.substr(0, lastx);
  mowri << '[' << lastline << ']';
  mowri << "\n";

  return sub_v_labels;
}

bool SkeletonClusterer::halt_refinement()
{
  auto   time_now = std::chrono::high_resolution_clock::now();
  size_t time_in_refine =
    std::chrono::duration_cast<std::chrono::microseconds>(time_now - refine_start).count();
  return ((rf_max_rounds < prd->rf_round) || (rf_max_time_micros <= time_in_refine));
}

void SkeletonClusterer::output_halt_refinement_reason()
{
  auto   time_now = std::chrono::high_resolution_clock::now();
  size_t time_in_refine =
    std::chrono::duration_cast<std::chrono::microseconds>(time_now - refine_start).count();
  unsigned n_reasons = 0;
  mowri << "halted refinement because: ";
  if (rf_max_time_micros <= time_in_refine)
  {
    mowri << "  [" << n_reasons + 1 << "] exceeded rf_max_time (" << rf_max_time_micros / 1000.
          << ")";
    ++n_reasons;
  }
  if (rf_max_rounds < prd->rf_round)
  {
    mowri << "  [" << n_reasons + 1 << "] exceeded rf_max_rounds (" << rf_max_rounds
          << ") of refinement";
    ++n_reasons;
  }
  if (n_reasons == 0)
  {
    mowri << "   round without any center update";
  }
  mowri << zentas::Endl;
}
}
