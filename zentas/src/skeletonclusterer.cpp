// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/skeletonclusterer.hpp>

namespace nszen
{

void P2Bundle::initialise_data(size_t n)
{
  ndata = n;
  d12_.reset(new std::array<double, 2>[n]);
  k_1_.reset(new size_t[n]);
  k_2_.reset(new size_t[n]);
  ori_.reset(new size_t[n]);

  for (unsigned i = 0; i < n; ++i)
  {
    d_1(i) = std::numeric_limits<double>::max();
    d_2(i) = std::numeric_limits<double>::max();
  }
  for (unsigned i = 0; i < n; ++i)
  {
    /* note that in the first round of kmeanspp, the lookup in cc requires this to be zero  */
    k_1(i) = 0;
  }
}

P2Bundle::P2Bundle(size_t n) { initialise_data(n); }

P2Bundle::P2Bundle(const std::vector<size_t>& ori)
{
  initialise_data(ori.size());
  for (size_t i = 0; i < ndata; ++i)
  {
    ori_[i] = ori[i];
  }
}

void XNearestInfo::reset(size_t new_a_x, double new_d_x, double new_e_x)
{
  a_x = new_a_x;
  d_x = new_d_x;
  e_x = new_e_x;
}

void XNearestInfo::reset(XNearestInfo& nearest_x_infos)
{
  a_x = nearest_x_infos.a_x;
  d_x = nearest_x_infos.d_x;
  e_x = nearest_x_infos.e_x;
}

XNearestInfo::XNearestInfo(size_t a_x_, double d_x_, double e_x_) : a_x(a_x_), d_x(d_x_), e_x(e_x_) 
{}

std::string XNearestInfo::get_string()
{
  return std::to_string(a_x) + "\t " + std::to_string(d_x) + "\t " + std::to_string(e_x);
}

SkeletonClustererInitBundle::SkeletonClustererInitBundle(
  size_t                                                      K_,
  size_t                                                      nd_,
  std::chrono::time_point<std::chrono::high_resolution_clock> bb_,
  const size_t* const                                         center_indices_init_predefined_,
  const std::string&                                          init_method_,
  double                                                      max_time_,
  double                                                      min_mE_,
  double                                                      max_itok_,
  size_t                                                      max_rounds_,
  size_t                                                      nthreads_,
  size_t                                                      seed_,
  const std::string&                                          energy_,
  bool                                                        with_tests_,
  size_t* const                                               indices_final_,
  size_t* const                                               labels_,
  const EnergyInitialiser*                                    ptr_energy_initialiser_,
  bool                                                        do_balance_labels_)
  : K(K_),
    ndata(nd_),
    bigbang(bb_),
    center_indices_init_predefined(center_indices_init_predefined_),
    init_method(init_method_),
    max_time(max_time_),
    min_mE(min_mE_),
    max_itok(max_itok_),
    max_rounds(max_rounds_),
    nthreads(nthreads_),
    seed(seed_),
    energy(energy_),
    with_tests(with_tests_),
    indices_final(indices_final_),
    labels(labels_),
    ptr_energy_initialiser(ptr_energy_initialiser_),
    do_balance_labels(do_balance_labels_)
{
}

SkeletonClusterer::SkeletonClusterer(const SkeletonClustererInitBundle& sb)
  : mowri(true, false, ""),
    K(sb.K),
    bigbang(sb.bigbang),
    ndata(sb.ndata),
    initialisation_method(sb.init_method),
    nearest_1_infos(K),
    sample_IDs(K),
    to_leave_cluster(K),
    cluster_has_changed(K, true),
    cluster_energies(K, 0),
    cluster_mean_energies(K),
    E_total(std::numeric_limits<double>::max()),
    old_E_total(0),
    round(0),
    v_center_indices_init(K),
    center_indices_init(v_center_indices_init.data()),
    max_time_micros(static_cast<size_t>(sb.max_time * 1000000.)),
    min_mE(sb.min_mE),
    max_itok(sb.max_itok),
    labels(sb.labels),
    nthreads(sb.nthreads),
    nthreads_fl(static_cast<double>(nthreads)),
    max_rounds(sb.max_rounds),
    energy(sb.energy),
    with_tests(sb.with_tests),
    gen(sb.seed),
    kmoo_cc(K * K),
    kmoo_p2bun(ndata),
    do_balance_labels(sb.do_balance_labels)

{

  if (energy.compare("identity") == 0)
  {
    f_energy = nszen::Identity();
  }

  else if (energy.compare("quadratic") == 0)
  {
    f_energy = nszen::Quadratic();
  }

  else if (energy.compare("cubic") == 0)
  {
    f_energy = nszen::Cubic();
  }

  else if (energy.compare("squarepotential") == 0)
  {
    f_energy = nszen::SquarePotential(sb.ptr_energy_initialiser->get_critical_radius());
  }

  else if (energy.compare("log") == 0)
  {
    f_energy = nszen::Log();
  }

  else if (energy.compare("exp") == 0)
  {
    f_energy = nszen::Exponential(sb.ptr_energy_initialiser->get_exponent_coeff());
  }

  else if (energy.compare("sqrt") == 0)
  {
    f_energy = nszen::SquareRoot();
  }

  else
  {
    std::stringstream errmss;
    errmss << "Unrecognised energy function, `" << energy << "'. ";
    if (energy == "linear")
    {
      errmss << "Perhaps instead of `linear' you meant `identity'?";
    }
    throw zentas::zentas_error(std::string() + energy);
  }

  /* confirm that f_energy(0) is 0 */
  if (f_energy(0) != 0)
  {
    std::stringstream ss;
    ss << "the energy function, f_energy, has f_energy(0) = ";
    ss << f_energy(0) << ". This is a problem for k-means++ initialisation.";
    ss << "If you're not using k-means++ initialisation, this should not cause any problems, ";
    ss << "but as k-means++ is the default, we are being cautious and throwing an error.";
    throw zentas::zentas_error(ss.str());
  }

  /* initialisation from indices. */
  if (initialisation_method == "from_indices_init")
  {
    if (sb.center_indices_init_predefined == nullptr)
    {
      throw zentas::zentas_error(
        "center_indices_init_predefined is nullptr, logic error in skeleton");
    }
    init::populate_from_indices_init(
      sb.center_indices_init_predefined, center_indices_init, K, ndata);
  }

  else
  {
    // will do in go ( )
  }

  // not sure about this.
  center_IDs = sb.indices_final;
}

size_t SkeletonClusterer::get_time_in_update_centers() { return time_in_update_centers; }

size_t SkeletonClusterer::get_time_in_update_sample_info() { return time_in_update_sample_info; }

double SkeletonClusterer::get_time_remaining()
{
  auto t1    = std::chrono::high_resolution_clock::now();
  time_total = std::chrono::duration_cast<std::chrono::microseconds>(t1 - bigbang).count();
  if (max_time_micros > time_total)
  {
    return max_time_micros - time_total;
  }
  else
  {
    return -1;
  }
}

size_t SkeletonClusterer::get_start(size_t ti, size_t nthreads_, size_t j_A, size_t j_Z)
{
  size_t n_js = j_Z - j_A;
  double t_fl = static_cast<double>(ti);
  size_t j_a  = j_A + (t_fl / static_cast<double>(nthreads_)) * n_js;
  return j_a;
}

size_t SkeletonClusterer::get_end(size_t ti, size_t nthreads_, size_t j_A, size_t j_Z)
{
  size_t n_js = j_Z - j_A;
  double t_fl = static_cast<double>(ti);
  size_t j_z  = j_A + ((t_fl + 1.) / static_cast<double>(nthreads_)) * n_js;
  return j_z;
}

size_t SkeletonClusterer::get_sample_from(std::vector<double>& v_cum_nearest_energies)
{

  /* kmeans++ is exhausted : everything sample has a center at distance 0 (obviously duplicated
   * data) however, this case should have been caught upstream, so throwing an error. */
  if (v_cum_nearest_energies.back() == 0)
  {
    throw zentas::zentas_error("exhausted in get_sample_from (kmeans++). in particular, the "
                               "cumulative energy is zero. this should have been caught "
                               "upstream: logic error in zentas");
  }

  return std::distance(v_cum_nearest_energies.begin(),
                       std::lower_bound(v_cum_nearest_energies.begin(),
                                        v_cum_nearest_energies.end(),
                                        dis_uni01(gen) * v_cum_nearest_energies.back())) -
         1;
}

std::string SkeletonClusterer::get_base_summary_string()
{

  auto t_now   = std::chrono::high_resolution_clock::now();
  time_total   = std::chrono::duration_cast<std::chrono::microseconds>(t_now - bigbang).count();
  ncalcs_total = get_ncalcs();

  std::string st_round = "R=" + std::to_string(round);
  st_round.resize(3 + 6, ' ');

  std::stringstream st_mE_ss;
  std::string       st_mE;
  if (E_total > 0)
  {
    st_mE_ss << "mE=" << std::setprecision(7) << E_total / static_cast<double>(ndata);
    st_mE = st_mE_ss.str();
  }
  else
  {
    st_mE = '-';
  }
  st_mE.resize(std::max<size_t>(st_mE.size() + 1, 3 + 14), ' ');

  std::string st_Tp = "Tp=" + std::to_string(time_prehistory / 1000);
  st_Tp += "  ";

  std::string st_Ti = "Ti=" + std::to_string(time_to_initialise_centers / 1000);
  st_Ti += "  ";

  std::string st_Tb = "Tb=" + std::to_string(time_initialising / 1000);
  st_Tb += "  ";

  std::string st_Tc = "Tc=" + std::to_string(time_in_update_centers / 1000);
  st_Tc.resize(3 + 8, ' ');

  std::string st_Tu = "Tu=" + std::to_string(time_in_update_sample_info / 1000);
  st_Tu.resize(3 + 7, ' ');

  std::string st_Tr = "Tr=" + std::to_string(time_in_redistribute / 1000);
  st_Tr.resize(3 + 5, ' ');

  std::string st_Tt = "Tt=" + std::to_string(time_total / 1000);
  st_Tt.resize(3 + 8, ' ');

  std::stringstream ncc_ss;
  ncc_ss << "lg2nc(c)=" << std::setprecision(5) << std::log2(ncalcs_in_update_centers);
  std::string ncc = ncc_ss.str();
  ncc.resize(9 + 9, ' ');

  std::stringstream nc_ss;
  nc_ss << "lg2nc=" << std::setprecision(5) << std::log2(ncalcs_total);
  std::string nc = nc_ss.str();
  nc.resize(6 + 9, ' ');

  std::stringstream pc_ss;
  pc_ss << "pc=" << std::setprecision(5) << static_cast<double>(get_rel_calccosts());
  std::string pc = pc_ss.str();
  pc.resize(3 + 10, ' ');

  std::ostringstream out;
  out << st_round << st_mE << st_Tp << st_Ti << st_Tb << st_Tc << st_Tu << st_Tr << st_Tt << ncc
      << nc << pc;
  return out.str();
}

/* rule : functions with suffix 'basic' will not touch to_leave_cluster */
void SkeletonClusterer::reset_nearest_info_basic(
  size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest)
{
  nearest_1_infos[k][j].reset(k_nearest, d_nearest, e_nearest);
  cluster_has_changed[k] = true;
}

void SkeletonClusterer::reset_nearest_info(
  size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest)
{
  reset_nearest_info_basic(k, j, k_nearest, d_nearest, e_nearest);
  if (k != k_nearest)
  {
    std::lock_guard<std::mutex> lock(mutex0);
    to_leave_cluster[k].push_back(j);
  }
}

void SkeletonClusterer::signal_cluster_change(size_t k) { cluster_has_changed[k] = true; }

size_t SkeletonClusterer::draw_j_uniform(size_t k) { return dis(gen) % get_ndata(k); }

void SkeletonClusterer::kmoo_prepare()
{
  kmoo_cc.resize(K * K);
  kmoo_p2bun = P2Bundle(ndata);
}

void SkeletonClusterer::print_ndatas()
{
  for (size_t k = 0; k < K; ++k)
  {
    mowri << get_ndata(k) << " ";
  }
  mowri << zentas::Endl;
}

void SkeletonClusterer::default_initialise_with_kmeanspp()
{
  unsigned n_bins;
  if (initialisation_method == "kmeans++")
  {
    n_bins = 1;
  }

  else
  {
    std::string prefix = "kmeans++-";
    n_bins             = init::extract_INT(initialisation_method, prefix.size());
    if (n_bins == 0)
    {
      throw zentas::zentas_error("n_bins passed to kmeans++ should be positive.");
    }
  }
  kmoo_prepare();
  km_is_exhausted = false;
  if (n_bins != 1)
  {
    triangular_kmeanspp_aq2(n_bins);

    if (km_is_exhausted)
    {
      mowri << "exhausted in " << initialisation_method << ", will revert to kmeans++-1"
            << zentas::Endl;
      kmoo_prepare();
    }
  }

  if (km_is_exhausted == true || n_bins == 1)
  {
    km_is_exhausted = false;

    triangular_kmeanspp();
    if (km_is_exhausted)
    {
      mowri << "exhausted in kmeans++-1, used uniform sampling to complete initialisation"
            << zentas::Endl;
    }
  }

  /* rearrange so that center_indices_init are in order */
  std::vector<std::array<size_t, 2>> vi(K);
  for (size_t k = 0; k < K; ++k)
  {
    std::get<0>(vi[k]) = center_indices_init[k];
    std::get<1>(vi[k]) = k;
  }

  auto fg = std::less<size_t>();
  std::sort(vi.begin(), vi.end(), [&fg](std::array<size_t, 2>& lhs, std::array<size_t, 2>& rhs) {
    return fg(std::get<0>(lhs), std::get<0>(rhs));
  });

  std::vector<size_t> old_to_new(K);
  for (unsigned k = 0; k < K; ++k)
  {
    center_indices_init[k]         = std::get<0>(vi[k]);
    old_to_new[std::get<1>(vi[k])] = k;
  }

  for (unsigned i = 0; i < ndata; ++i)
  {
    kmoo_p2bun.k_1(i) = old_to_new[kmoo_p2bun.k_1(i)];
    kmoo_p2bun.k_2(i) = old_to_new[kmoo_p2bun.k_2(i)];
  }
}

void SkeletonClusterer::triangular_kmeanspp_aq2(size_t n_bins)
{

  /* experiments so far show that multithreading does not help here, can hurt. what's weird is
   * that even if nthreads = 1 in the pll version, it's sig slower than the serial version.  */
  bool   multithread_kmpp = true;
  double a_distance;

  /* non_tail_k will be how many k's for which only 1/n_bins of the data is used.
   * it will be (n_bins - 1)/n_bins, rounded down to the nearest multiple of n_bins.
   * an analysis suggests that the tail should be sqrt(K/nbins), so this is conservative.
   * we also ensure that the tail is at least of length min(K, 50).
   * */
  size_t non_tail_k;
  non_tail_k =
    static_cast<size_t>(static_cast<double>(K) * (1. - 1. / static_cast<double>(n_bins)));
  size_t minimal_tail_k     = std::min<size_t>(K, 50);
  size_t maximal_non_tail_k = K - minimal_tail_k;
  non_tail_k                = std::min(non_tail_k, maximal_non_tail_k);

  /* tail_k is how many k's have use all the data (at the end) */
  size_t tail_k = K - non_tail_k;
  reset_p2buns_dt(n_bins);
  aq2p_original_indices = std::vector<std::vector<size_t>>(n_bins);
  aq2p_p2buns           = std::vector<P2Bundle>(n_bins);

  std::vector<size_t> bin_indices(n_bins);
  for (size_t bin = 0; bin < n_bins; ++bin)
  {
    bin_indices[bin] = bin;
  }
  for (size_t i = 0; i < ndata; ++i)
  {
    if (i % n_bins == 0)
    {
      /* shuffle the indices */
      for (size_t bin = 0; bin < n_bins; ++bin)
      {
        std::swap(bin_indices[bin], bin_indices[bin + dis(gen) % (n_bins - bin)]);
      }
    }
    size_t bin = bin_indices[i % n_bins];
    append_aq2p_p2buns(bin, i);
    aq2p_original_indices[bin].push_back(i);
  }

  for (size_t bin = 0; bin < n_bins; ++bin)
  {
    aq2p_p2buns[bin] = P2Bundle(aq2p_original_indices[bin]);
  }

  auto update_nearest_info =
    [&a_distance, this](size_t bin, size_t k_start, size_t k_end, size_t i0, size_t i1) {
      for (size_t i = i0; i < i1; ++i)
      {
        for (size_t k = k_start; k < k_end; ++k)
        {
          if (kmoo_cc[aq2p_p2buns[bin].k_1(i) * K + k] <
              aq2p_p2buns[bin].d_1(i) + aq2p_p2buns[bin].d_2(i))
          {
            set_center_sample_pp_distance(k, bin, i, a_distance);
            kmpp_inner(i, k, a_distance, aq2p_p2buns[bin]);
          }
        }
      }
    };

  for (size_t bin = 0; bin < n_bins; ++bin)
  {
    size_t k0 = (bin * non_tail_k) / n_bins;
    size_t k1 = ((bin + 1) * non_tail_k) / n_bins;

    /* update nearest info from 0 to k0 of this bin*/
    if (nthreads == 1 || multithread_kmpp == false)
    {
      update_nearest_info(bin, 0, k0, 0, aq2p_p2buns[bin].get_ndata());
    }

    else
    {
      std::vector<std::thread> threads;
      for (size_t ti = 0; ti < get_nthreads(); ++ti)
      {
        threads.emplace_back([this, ti, bin, k0, &update_nearest_info]() {
          update_nearest_info(bin,
                              0,
                              k0,
                              get_start(ti, get_nthreads(), 0, aq2p_p2buns[bin].get_ndata()),
                              get_end(ti, get_nthreads(), 0, aq2p_p2buns[bin].get_ndata()));
        });
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }

    /* we don't even try to parallelize to center by center kmeans++ part, too many syncs. */
    triangular_kmeanspp_after_initial(bin, k0, k1, false);
    /* km_is_exhausted is true, we will bail as handling this here (with n_bins > 1) is too
     * complex. Will try again with n_bins = 1 */
    if (km_is_exhausted)
    {
      return;
    }
  }

  /* update all until the tail k's (*&^*) */
  for (size_t bin = 0; bin < n_bins; ++bin)
  {
    size_t k1 = ((bin + 1) * non_tail_k) / n_bins;

    if (nthreads == 1 || multithread_kmpp == false)
    {
      update_nearest_info(bin, k1, non_tail_k, 0, aq2p_p2buns[bin].get_ndata());
    }

    else
    {
      std::vector<std::thread> threads;
      for (size_t ti = 0; ti < get_nthreads(); ++ti)
      {
        threads.emplace_back([this, ti, bin, k1, non_tail_k, &update_nearest_info]() {
          update_nearest_info(bin,
                              k1,
                              non_tail_k,
                              get_start(ti, get_nthreads(), 0, aq2p_p2buns[bin].get_ndata()),
                              get_end(ti, get_nthreads(), 0, aq2p_p2buns[bin].get_ndata()));
        });
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }
  }

  // get p2bun up to date.
  for (size_t bin = 0; bin < n_bins; ++bin)
  {
    for (size_t i = 0; i < aq2p_p2buns[bin].get_ndata(); ++i)
    {
      size_t index0          = aq2p_p2buns[bin].ori(i);
      kmoo_p2bun.d_1(index0) = aq2p_p2buns[bin].d_1(i);
      kmoo_p2bun.d_2(index0) = aq2p_p2buns[bin].d_2(i);
      kmoo_p2bun.k_1(index0) = aq2p_p2buns[bin].k_1(i);
      kmoo_p2bun.k_2(index0) = aq2p_p2buns[bin].k_2(i);
      kmoo_p2bun.ori(index0) = index0;
    }
  }

  // set the tail k's
  triangular_kmeanspp_after_initial(std::numeric_limits<size_t>::max(), K - tail_k, K, true);
  /* km_is_exhausted is true, we will bail as handling this here (with n_bins > 1) is too
   * complex. Will try again with n_bins = 1 */
  if (km_is_exhausted)
  {
    return;
  }
}

void SkeletonClusterer::triangular_kmeanspp()
{
  for (size_t i = 0; i < ndata; ++i)
  {
    kmoo_p2bun.ori(i) = i;
  }
  size_t k0 = 0;
  size_t k1 = K;
  triangular_kmeanspp_after_initial(std::numeric_limits<size_t>::max(), k0, k1, true);
}

void SkeletonClusterer::kmpp_inner(size_t i, size_t k, double a_distance, P2Bundle& p2bun)
{

  if (a_distance < p2bun.d_2(i))
  {
    if (a_distance < p2bun.d_1(i))
    {
      p2bun.d_2(i) = p2bun.d_1(i);
      p2bun.k_2(i) = p2bun.k_1(i);
      p2bun.d_1(i) = a_distance;
      p2bun.k_1(i) = k;
    }
    else
    {
      p2bun.d_2(i) = a_distance;
      p2bun.k_2(i) = k;
    }
  }
}

size_t SkeletonClusterer::get_nthreads() { return nthreads; }

double SkeletonClusterer::get_nthreads_fl() { return nthreads_fl; }

void SkeletonClusterer::triangular_kmeanspp_after_initial(size_t aq2p_bin,
                                                          size_t k0,
                                                          size_t k1,
                                                          bool   from_full_data)
{

  if (from_full_data != (aq2p_bin == std::numeric_limits<size_t>::max()))
  {
    throw zentas::zentas_error("logic error in triangular_kmeanspp_after_initial : "
                               "from_full_data != (aq2p_bin == "
                               "std::numeric_limits<size_t>::max()");
  }

  P2Bundle& p2bun          = from_full_data ? kmoo_p2bun : aq2p_p2buns[aq2p_bin];
  size_t    kmeanspp_ndata = from_full_data ? ndata : get_p2buns_dt_ndata(aq2p_bin);

  if (kmeanspp_ndata - 1 < k1 - k0)
  {
    std::stringstream ss;
    ss << "triangular_kmeanspp_after_initial, attempting to find more centers than data - 1 : ";
    ss << "kmeanspp_ndata = " << kmeanspp_ndata << "   and   k1 - k0 = " << k1 - k0;
    throw zentas::zentas_error(ss.str());
  }

  test_parameters_to_tkai(k0, k1, p2bun.get_ndata(), kmeanspp_ndata);

  double a_distance;

  /* the cumulative sampling distribution */
  std::vector<double> v_cum_nearest_energies(kmeanspp_ndata + 1);
  v_cum_nearest_energies[0] = 0.;

  if (k0 == 0)
  {
    for (size_t i = 0; i < kmeanspp_ndata; ++i)
    {
      v_cum_nearest_energies[i + 1] = v_cum_nearest_energies[i] + 1.;
    }
  }

  else
  {
    for (size_t i = 0; i < kmeanspp_ndata; ++i)
    {
      v_cum_nearest_energies[i + 1] = v_cum_nearest_energies[i] + f_energy(p2bun.d_1(i));
    }
  }

  /* depending on whether we're using all the data or the data in ptr_p2bun_dt, distance
   * calculations are performed differently. */
  std::function<void(size_t, size_t)> set_distance_kk;
  std::function<void(size_t, size_t)> set_distance_ik;
  std::function<void(size_t)> update_c_ind_c_dt;

  /* whether to revert to uniform sampling when all samples already have a center at distance 0
   * (obviously duplicated data)*/
  bool try_uniform_when_exhausted;

  if (from_full_data == true)
  {

    set_distance_kk = [this](size_t k, size_t kp) {
      set_sampleID_sampleID_distance_nothreshold(
        center_indices_init[k], center_indices_init[kp], kmoo_cc[k * K + kp]);
    };

    set_distance_ik = [this, &a_distance, &p2bun](size_t i, size_t k) {
      set_sampleID_sampleID_distance(center_indices_init[k], i, p2bun.d_2(i), a_distance);
    };

    update_c_ind_c_dt = [this, &v_cum_nearest_energies](size_t k) {
      center_indices_init[k] = get_sample_from(v_cum_nearest_energies);
      append_pp_from_ID(center_indices_init[k]);
    };

    /* we can handle the exhausted case here */
    try_uniform_when_exhausted = true;
  }

  else
  {

    set_distance_kk = [this](size_t k, size_t kp) { set_cc_pp_distance(k, kp); };

    set_distance_ik = [this, &a_distance, &p2bun, aq2p_bin](size_t i, size_t k) {
      set_center_sample_pp_distance(k, aq2p_bin, i, a_distance);
    };

    update_c_ind_c_dt = [this, &v_cum_nearest_energies, &p2bun, aq2p_bin](size_t k) {
      size_t sami            = get_sample_from(v_cum_nearest_energies);
      center_indices_init[k] = p2bun.ori(sami);
      append_pp_from_bin(aq2p_bin, sami);

    };

    /* handling the exhausted case here is too complex (cross batch indices needed, argh)
     * , rather just bail and try again from full data */
    try_uniform_when_exhausted = false;
  }

  /* k-means++, at last */
  for (size_t k = k0; k < k1; ++k)
  {

    if (v_cum_nearest_energies.back() != 0)
    {
      km_is_exhausted = false;
      update_c_ind_c_dt(k);
    }

    else
    {
      km_is_exhausted = true;
      if (try_uniform_when_exhausted == false)
      {
        /* bailing, because
         * (1) every sample has a center a distance 0 : k-means++ fails
         * (2) must not switch to uniform sampling.
         * */
        return;
      }

      else
      {
        if (from_full_data == false)
        {
          std::stringstream ss;
          ss << "try uniform when exhausted is true, but from full data is false, this looks "
                "like a logic error.";
          throw zentas::zentas_error(ss.str());
        }

        size_t new_index = dis(gen) % (K);
        while (std::find(center_indices_init, center_indices_init + k, new_index) !=
               center_indices_init + k)
        {
          new_index = dis(gen) % (K);
        }
        center_indices_init[k] = new_index;
        append_pp_from_ID(center_indices_init[k]);
      }
    }

    /* update cc */
    kmoo_cc[k * K + k] = 0.;
    for (size_t kp = 0; kp < k; ++kp)
    {
      set_distance_kk(k, kp);
      kmoo_cc[kp * K + k] = kmoo_cc[k * K + kp];
    }

    /* update nearest distances, second nearest distances, and centers */
    for (size_t i = 0; i < kmeanspp_ndata; ++i)
    {
      /* note that if k0 = 0, k_1(i) is 0 so this is fine. */
      if (kmoo_cc[p2bun.k_1(i) * K + k] < p2bun.d_1(i) + p2bun.d_2(i))
      {
        set_distance_ik(i, k);
        kmpp_inner(i, k, a_distance, p2bun);
      }
      v_cum_nearest_energies[i + 1] = v_cum_nearest_energies[i] + f_energy(p2bun.d_1(i));
    }
  }
}

void SkeletonClusterer::test_parameters_to_tkai(size_t k0,
                                                size_t k1,
                                                size_t ndata_1,
                                                size_t ndata_2)
{
  if (k0 >= K || k1 < k0 || k1 > K)
  {
    std::stringstream errm;
    errm << "invalid k0, k1 in kmeanspp : ";
    errm << "k0 = " << k0 << "     k1 = " << k1 << "     K = " << K << ".";
    throw zentas::zentas_error(errm.str());
  }

  if (ndata_1 != ndata_2)
  {
    throw zentas::zentas_error("ndata_1 != ndata_2 in test_parameters_to_tkai");
  }
}

// set the distance from center k to the j1'th element of cluster k1.
void SkeletonClusterer::set_center_sample_distance_nothreshold(size_t  k,
                                                               size_t  k1,
                                                               size_t  j1,
                                                               double& distance)
{
  set_center_sample_distance(k, k1, j1, std::numeric_limits<double>::max(), distance);
}

void SkeletonClusterer::set_sampleID_sampleID_distance_nothreshold(size_t  i1,
                                                                   size_t  i2,
                                                                   double& distance)
{
  set_sampleID_sampleID_distance(i1, i2, std::numeric_limits<double>::max(), distance);
}

double SkeletonClusterer::get_center_sample_distance_nothreshold(size_t k, size_t k1, size_t j1)
{
  double distance;
  set_center_sample_distance_nothreshold(k, k1, j1, distance);
  return distance;
}

double SkeletonClusterer::get_sample_sample_distance_nothreshold(size_t k, size_t j1, size_t j2)
{
  double adistance;
  set_sample_sample_distance(k, j1, k, j2, std::numeric_limits<double>::max(), adistance);
  return adistance;
}

void SkeletonClusterer::set_center_center_distance_nothreshold(size_t  k1,
                                                               size_t  k2,
                                                               double& adistance)
{
  set_center_center_distance(k1, k2, std::numeric_limits<double>::max(), adistance);
}

void SkeletonClusterer::output_halt_kmedoids_reason()
{

  unsigned n_reasons = 0;
  mowri << "halted k-medoids because: " << zentas::Flush;
  if (time_total >= max_time_micros)
  {
    mowri << "  [" << n_reasons + 1 << "] exceeded max_time (" << max_time_micros / 1000. << ")"
          << zentas::Flush;
    ++n_reasons;
  }
  if (round >= max_rounds)
  {
    mowri << "  [" << n_reasons + 1 << "] exceeded max_rounds (" << max_rounds << ")"
          << zentas::Flush;
    ++n_reasons;
  }
  if (E_total / static_cast<double>(ndata) < min_mE)
  {
    mowri << "  [" << n_reasons + 1 << "] mE below termination min_mE (" << min_mE << ")"
          << zentas::Flush;
    ;
    ++n_reasons;
  }

  if (time_total >= (max_itok + 1) * (time_to_initialise_centers + time_initialising))
  {
    mowri << "  [" << n_reasons + 1
          << "] exceeded ratio max_itok (itialize : k-medoids <= " << max_itok << ")"
          << zentas::Flush;
    ++n_reasons;
  }

  if (n_reasons == 0)
  {
    mowri << "   round without any center update" << zentas::Flush;
  }

  mowri << zentas::Endl;
}

bool SkeletonClusterer::halt_kmedoids()
{
  bool do_not_halt =
    (time_total < max_time_micros) && (round < max_rounds) &&
    ((E_total / static_cast<double>(ndata)) >= min_mE) &&
    time_total < (max_itok + 1) * (time_to_initialise_centers + time_initialising);

  return (do_not_halt == false);
}

void SkeletonClusterer::set_t_ncalcs_update_centers_start()
{
  t_update_centers_start      = std::chrono::high_resolution_clock::now();
  ncalcs_update_centers_start = get_ncalcs();
}

void SkeletonClusterer::update_t_ncalcs_center_end()
{
  t_update_centers_end = std::chrono::high_resolution_clock::now();
  time_in_update_centers += std::chrono::duration_cast<std::chrono::microseconds>(
                              t_update_centers_end - t_update_centers_start)
                              .count();
  ncalcs_update_centers_end = get_ncalcs();
  ncalcs_in_update_centers += ncalcs_update_centers_end - ncalcs_update_centers_start;
}

void SkeletonClusterer::update_t_ncalcs_sample_update_end()
{
  t_update_sample_info_end = std::chrono::high_resolution_clock::now();
  time_in_update_sample_info += std::chrono::duration_cast<std::chrono::microseconds>(
                                  t_update_sample_info_end - t_update_centers_end)
                                  .count();
  ncalcs_update_sample_info_end = get_ncalcs();
  ncalcs_in_update_sample_info += ncalcs_update_sample_info_end - ncalcs_update_centers_end;
}

void SkeletonClusterer::update_t_redistibute_end()
{
  t_redistribute_end = std::chrono::high_resolution_clock::now();
  time_in_redistribute += std::chrono::duration_cast<std::chrono::microseconds>(
                            t_redistribute_end - t_update_sample_info_end)
                            .count();
}

void SkeletonClusterer::update_t_update_all_cluster_stats_end()
{
  t_update_all_cluster_statistics_end = std::chrono::high_resolution_clock::now();
  time_in_update_all_cluster_statistics +=
    std::chrono::duration_cast<std::chrono::microseconds>(t_update_all_cluster_statistics_end -
                                                          t_redistribute_end)
      .count();
}

void SkeletonClusterer::run_kmedoids()
{
  core_kmedoids_loops();
  output_halt_kmedoids_reason();
}

void SkeletonClusterer::core_kmedoids_loops()
{

  while (halt_kmedoids() == false)
  {

    /* ************** *
    * UPDATE CENTERS *
    * ****************/
    set_t_ncalcs_update_centers_start();
    bool modified_centers = update_centers();

    if (modified_centers == false)
    {
      update_t_ncalcs_center_end();
      mowri << get_round_summary() << zentas::Endl;
      break;
    }

    /* ************************* *
    * UPDATE CENTER CENTER INFO *
    * ***************************/
    update_center_center_info();
    if (with_tests == true)
      post_center_update_test();
    update_t_ncalcs_center_end();

    for (auto& x : to_leave_cluster)
    {
      x.resize(0);
    }

    /* ****************** *
    * UPDATE SAMPLE INFO *
    * ********************/
    update_sample_info();
    if (with_tests == true)
      post_sample_update_test();
    update_t_ncalcs_sample_update_end();

    /* ************ *
    * REDISTRIBUTE *
    * **************/
    redistribute();
    if (with_tests == true)
      post_redistribute_test();
    update_t_redistibute_end();

    /* ************************* *
    * UPDATE CLUSTER STATISTICS *
    * ***************************/
    update_all_cluster_statistics();
    if (with_tests == true)
      post_update_statistics_test();
    update_t_update_all_cluster_stats_end();

    if (with_tests == true)
      mowri << zentas::Endl;
    mowri << get_round_summary() << zentas::Endl;
    ++round;
  }
}

void SkeletonClusterer::populate_labels()
{
  for (size_t k = 0; k < K; ++k)
  {
    labels[center_IDs[k]] = k;
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      labels[sample_IDs[k][j]] = k;
    }
  }
}

void SkeletonClusterer::go()
{

  if (with_tests == true)
  {
    mowri << "\n\nRUNNING WITH TESTS ENABLED : WILL BE SLOW" << zentas::Endl;
  }

  mowri <<
    R"((The prevent output to terminal, set capture_output to false)
(For a description of column statistics, consider function get_output_verbose_string())
)";

  std::string init1foo = "Will now initialize with method: " + initialisation_method + ".";
  mowri << '\n' << init1foo << '\n';
  // mowri << get_char_line(init1foo.size(), '-');
  initialise_all();

  std::string kmefact =
    "Will now run k-medoids a.k.a  k-centers with method: " + get_kmedoids_method_string() + ".";
  mowri << '\n' << kmefact << '\n';
  mowri << get_equals_line(get_round_summary().size());
  run_kmedoids();

  if (do_refinement == true)
  {
    refine_start = std::chrono::high_resolution_clock::now();

    std::string inifact = "Will now initialize for refinement with method: " + rf_alg + ".";
    mowri << '\n' << inifact << '\n';
    mowri << get_equals_line(get_round_summary().size());
    initialise_refinement();

    std::string runfact = "Will now run refinement.";
    mowri << '\n' << runfact << '\n';
    mowri << get_equals_line(rf_get_round_summary().size());
    run_refinement();

    std::string finalref = "Refinement complete, presenting final statistics.";
    mowri << '\n' << finalref << '\n';
    mowri << get_equals_line(rf_get_round_summary().size());
    auto final_line = rf_get_round_summary();
    mowri << final_line << '\n';
  }
  populate_labels();

  if (do_balance_labels == true && do_refinement == true)
  {
    throw zentas::zentas_error("balancing labels is only an option when there is no refinement");
  }

  if (do_balance_labels == true)
  {
    balance_the_labels();
  }
}

void SkeletonClusterer::balance_the_labels()
{

  // clusters with less than 1/ballyfactor * (ndata/K) will greedily steal members from nearest
  // clusters,
  // until all clusters have at least 1/ballyfactor * (ndata/K) members.

  const double ballyfactor = 1.5;

  if (ndata < 3 * K)
  {
    throw zentas::zentas_error("zentas does not support balancing when ndata < 3K");
  }

  size_t min_per_cluster = std::max<size_t>(2, static_cast<size_t>((ndata) / (ballyfactor * K)));
  std::vector<std::vector<size_t>> ids(K);
  for (size_t i = 0; i < ndata; ++i)
  {
    ids[labels[i]].push_back(i);
  }

  if (min_per_cluster * K > ndata)
  {
    throw zentas::zentas_error(" min_per_cluster*K > ndata : does not make sense");
  }

  double cc_dist;
  std::vector<std::tuple<double, size_t>> dk(K);
  for (size_t k = 0; k < K; ++k)
  {
    // if cluster k is too small,
    if (ids[k].size() < min_per_cluster)
    {

      // get the distances to the other clusters, and sort from nearest to farthest.
      for (size_t kp = 0; kp < K; ++kp)
      {
        set_center_center_distance_nothreshold(k, kp, cc_dist);
        dk[kp] = std::make_tuple(cc_dist, kp);
      }
      std::sort(
        dk.begin(), dk.end(), [](std::tuple<double, size_t>& x1, std::tuple<double, size_t>& x2) {
          return std::get<0>(x1) < std::get<0>(x2);
        });
      if (std::get<0>(dk[0]) > std::get<0>(dk[1]))
      {
        std::stringstream ss;
        ss << "sorting order seems incorrect in balance_the_labels:\n";
        for (size_t kpp = 0; kpp < K; ++kpp)
        {
          ss << std::get<0>(dk[kpp]) << " ";
        }
        throw zentas::zentas_error(ss.str());
      }

      // going through the `other' clusters,
      for (auto& x : dk)
      {
        size_t kp = std::get<1>(x);
        // while the other cluster can lose elements while remaining large enough and k is still
        // too small, make a switch.
        while (ids[kp].size() > min_per_cluster && ids[k].size() < min_per_cluster)
        {
          labels[ids[kp].back()] = k;
          ids[k].push_back(ids[kp].back());
          ids[kp].pop_back();
          if (ids[k].size() == min_per_cluster)
          {
            break;
          }
        }
        if (ids[k].size() == min_per_cluster)
        {
          break;
        }
      }
      if (ids[k].size() != min_per_cluster)
      {
        std::stringstream ss;
        ss << "ids[k] size should be min_per_cluster here.\n";
        ss << "ids[k] has size " << ids[k].size() << ", and min_per_cluster is " << min_per_cluster
           << ".\n";
        ss << "ndata is " << ndata << " and K is " << K << "\n";
        throw zentas::zentas_error(ss.str());
      }
    }
  }

  for (size_t k = 0; k < K; ++k)
  {
    if (ids[k].size() < min_per_cluster)
    {
      throw zentas::zentas_error("balancing appears to have failed, this is a logic error of jn.");
    }
  }
}

void SkeletonClusterer::post_initialisation_test()
{
  mowri << "((post_initialisation_test):" << zentas::Flush;
  mowri << "(injective)" << zentas::Flush;
  injective_ID_test();
  mowri << "(ndata)" << zentas::Flush;
  ndata_tests();
  mowri << "(info)" << zentas::Flush;
  info_tests();
  mowri << "(statistics)" << zentas::Flush;
  mowri << ")" << zentas::Endl;
  cluster_statistics_test();
}

void SkeletonClusterer::center_center_info_test() {}

void SkeletonClusterer::post_center_update_test()
{
  mowri << "((post_center_update_test):" << zentas::Flush;
  mowri << "(injective)" << zentas::Flush;
  injective_ID_test();
  mowri << "(ndata)" << zentas::Flush;
  ndata_tests();
  mowri << "(center_center_info)" << zentas::Flush;
  center_center_info_test();
  mowri << ")" << zentas::Endl;
}

void SkeletonClusterer::post_sample_update_test()
{
  mowri << "(post_sample_update_test)" << zentas::Flush;
  mowri << "(ndata)" << zentas::Flush;
  ndata_tests();
  mowri << "(info)" << zentas::Flush;
  info_tests();
  mowri << "(to_leave_cluster)" << zentas::Flush;
  to_leave_cluster_test();
  mowri << "(injective)" << zentas::Flush;
  injective_ID_test();
  mowri << ")" << zentas::Endl;
}

void SkeletonClusterer::post_redistribute_test()
{
  mowri << "(post_redistribute_test)" << zentas::Flush;
  mowri << "(ndata)" << zentas::Flush;
  ndata_tests();
  mowri << "(injective)" << zentas::Flush;
  injective_ID_test();
  mowri << "(as_assigned_test)" << zentas::Flush;
  as_assigned_test();
  mowri << ")" << zentas::Endl;
}

void SkeletonClusterer::post_update_statistics_test()
{
  mowri << "(injective)" << zentas::Flush;
  injective_ID_test();
  mowri << "(cluster_statistics)" << zentas::Flush;
  cluster_statistics_test();
  mowri << ")" << zentas::Endl;
}

void SkeletonClusterer::to_leave_cluster_test()
{
  std::string errm("error in leave_cluster_test");
  for (size_t k = 0; k < K; ++k)
  {
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      if (nearest_1_infos[k][j].a_x != k)
      {
        bool is_listed = false;
        for (auto& x : to_leave_cluster[k])
        {
          if (x == j)
          {
            is_listed = true;
          }
        }
        if (is_listed == false)
        {
          mowri << "\nis_listed == false" << zentas::Endl;
          throw zentas::zentas_error(errm);
        }
      }
    }

    for (auto& x : to_leave_cluster[k])
    {

      size_t nappears = 0;
      for (auto& xp : to_leave_cluster[k])
      {
        if (xp == x)
        {
          ++nappears;
        }
      }
      if (nappears != 1)
      {
        mowri << std::to_string(x) << " appears " << nappears
              << " times in the to leave list of cluster " << k << zentas::Endl;
      }

      if (nearest_1_infos[k][x].a_x == k)
      {
        mowri << "\ncluster and size ( " << k << " : " << get_ndata(k) << ") " << zentas::Endl;
        mowri << "j : " << x << zentas::Endl;
        mowri << "size of to_leave_cluster[k] " << to_leave_cluster[k].size() << zentas::Endl;
        mowri << "members of to_leave_cluster[k]" << zentas::Endl;
        for (auto& x_ : to_leave_cluster[k])
        {
          mowri << x_ << " " << zentas::Flush;
        }
        mowri << "\n(a,d,e) : " << nearest_1_infos[k][x].get_string() << zentas::Endl;
        mowri << "\nto leave cluster but k is 'k'" << zentas::Endl;
        throw zentas::zentas_error(errm);
      }
    }
  }
}

void SkeletonClusterer::custom_cluster_statistics_test() {}

void SkeletonClusterer::cluster_statistics_test()
{
  std::string errs("cluster statistics test failed.");
  // energies.
  double E__0 = 0;

  for (size_t k = 0; k < K; ++k)
  {
    double E__0k = 0;
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      E__0k += nearest_1_infos[k][j].e_x;
    }
    if (E__0k != cluster_energies[k])
    {
      throw zentas::zentas_error(errs);
    }
    E__0 += E__0k;
  }

  if (E__0 != get_E_total())
  {
    throw zentas::zentas_error(errs);
  }

  custom_cluster_statistics_test();
}

void SkeletonClusterer::fundamental_triangle_inequality_test()
{

  size_t n_tests = ndata * ndata;

  for (unsigned ti = 0; ti < n_tests; ++ti)
  {
    size_t index0 = dis(gen) % ndata;
    size_t index1 = dis(gen) % ndata;
    size_t index2 = dis(gen) % ndata;
    double d01, d02, d12;

    set_sampleID_sampleID_distance_nothreshold(index0, index1, d01);
    set_sampleID_sampleID_distance_nothreshold(index0, index2, d02);
    set_sampleID_sampleID_distance_nothreshold(index1, index2, d12);

    if ((1 - 1e-5) * d02 > d01 + d12)
    {
      std::stringstream ss;
      ss << "Triangle inequality failed :\n";
      ss << "Sample 0: " << string_from_ID(index0) << "\n";
      ss << "Sample 1: " << string_from_ID(index1) << "\n";
      ss << "Sample 2: " << string_from_ID(index2) << "\n";
      ss << "d02 = " << d02 << ", d01 = " << d01 << ", d12 = " << d12
         << ", d01 + d12 = " << (d01 + d12);
      throw zentas::zentas_error(ss.str());
    }
  }
}

void SkeletonClusterer::custom_info_test() {}

void SkeletonClusterer::info_tests()
{
  std::string errm("info test failed. ");
  for (size_t k = 0; k < K; ++k)
  {
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      std::vector<double> distances(K);
      size_t              k_first_nearest;
      double              d_first_nearest = std::numeric_limits<double>::max();
      for (size_t k2 = 0; k2 < K; ++k2)
      {
        set_center_sample_distance_nothreshold(k2, k, j, distances[k2]);
        if (distances[k2] < d_first_nearest)
        {
          d_first_nearest = distances[k2];
          k_first_nearest = k2;
        }
      }
      double e_first_nearest = f_energy(d_first_nearest);
      if (d_first_nearest != nearest_1_infos[k][j].d_x)
      {

        std::stringstream strerrm;

        strerrm << "error detected : d_first_nearest != d_first_nearest\n";
        strerrm << "k=" << k << "\n";
        strerrm << "j=" << j << "\n";
        strerrm << "the " << j << "'th sample in cluster " << k << " is " << string_for_sample(k, j)
             << "\n";
        strerrm << "cluster size is " << nearest_1_infos[k].size() << "\n";
        strerrm << std::setprecision(20);
        strerrm << "get_a1(k,j)=" << get_a1(k, j) << "\n";
        strerrm << "just computed first nearest center index: " << k_first_nearest << "\n";
        strerrm << "the " << k_first_nearest
             << "'th center is: " << string_for_center(k_first_nearest) << "\n";
        strerrm << "just computed distance to this center is: " << d_first_nearest << "\n";
        strerrm << "the recorded first nearest center index: " << nearest_1_infos[k][j].a_x << "\n";
        strerrm << "the " << nearest_1_infos[k][j].a_x << "'th center is "
             << string_for_center(nearest_1_infos[k][j].a_x) << "\n";
        strerrm << "the recorded distance to this center is " << nearest_1_infos[k][j].d_x << "\n";
        throw zentas::zentas_error(strerrm.str());
      }

      if (e_first_nearest != nearest_1_infos[k][j].e_x)
      {
        throw zentas::zentas_error(errm + "e_first_nearest != nearest_1_infos[k][j].e_x");
      }
    }
  }
  custom_info_test();
}

void SkeletonClusterer::custom_ndata_test() {}

void SkeletonClusterer::ndata_tests()
{
  std::string errm("ndata test failure");
  size_t      ndata_internal = 0;
  for (size_t k = 0; k < K; ++k)
  {
    if (get_ndata(k) != nearest_1_infos[k].size())
    {
      throw zentas::zentas_error(errm);
    }
    ndata_internal += get_ndata(k);
  }

  custom_ndata_test();

  size_t ndata_target;
  ndata_target = ndata - K;

  if (ndata_internal != ndata_target)
  {
    throw zentas::zentas_error(errm);
  }
}

void SkeletonClusterer::as_assigned_test()
{
  for (size_t k = 0; k < K; ++k)
  {
    for (size_t j = 0; j < get_ndata(k); ++j)
    {
      if (nearest_1_infos[k][j].a_x != k)
      {
        std::string errstring = "A sample in cluster " + std::to_string(k) + " has a_x " +
                                std::to_string(nearest_1_infos[k][j].a_x);
        throw zentas::zentas_error(errstring);
      }
    }
  }
}

void SkeletonClusterer::injective_ID_test()
{
  for (size_t i = 0; i < ndata; ++i)
  {
    size_t              count_i = 0;
    std::vector<size_t> as_center_ID;
    std::vector<size_t> in_cluster;
    std::vector<size_t> at_index;

    for (size_t k = 0; k < K; ++k)
    {

      if (i == center_IDs[k])
      {
        ++count_i;
        as_center_ID.push_back(k);
      }

      for (size_t j = 0; j < get_ndata(k); ++j)
      {
        if (i == sample_IDs[k][j])
        {
          ++count_i;
          in_cluster.push_back(k);
          at_index.push_back(j);
        }
      }
    }
    if (count_i != 1)
    {
      std::stringstream errm;
      errm << "Error caught. ID " << i << " appears not exactly once. ";
      errm << "It appears as centers of [";
      for (auto& x : as_center_ID)
      {
        errm << " `" << x << "' ";
      }
      errm << " ], as well as member (k,j) {size of k}) in [";
      for (size_t q = 0; q < in_cluster.size(); ++q)
      {
        errm << " (" << in_cluster[q] << "," << at_index[q] << " {" << get_ndata(in_cluster[q])
             << "}) ";
      }

      errm << "].";
      throw zentas::zentas_error(errm.str());
    }
  }
}

void SkeletonClusterer::post_initialise_centers_test()
{

  for (size_t k = 0; k < K; ++k)
  {
    if (center_indices_init[k] == center_indices_init[(k + 1) % ndata])
    {
      std::stringstream errm_ss;
      errm_ss << "initialising indices must be unique, received " << center_indices_init[k]
              << " at least twice. The initialising indices are \n";
      for (size_t k2 = 0; k2 < K; ++k2)
      {
        errm_ss << " " << center_indices_init[k2] << " ";
      }
      errm_ss << "\n";
      throw zentas::zentas_error(errm_ss.str());
    }
  }
}

void SkeletonClusterer::populate_afk_mc2()
{

  size_t chain_length = init::extract_INT(initialisation_method, 8);
  mowri << "afk_mc2 chain length : " << chain_length << zentas::Endl;

  double adistance;

  std::vector<bool> is_center(ndata, false);
  /* the first center is chosen at random */
  size_t index0          = dis(gen) % ndata;
  center_indices_init[0] = index0;
  is_center[index0]      = true;

  /* energies */
  std::vector<double> e0(ndata);
  /* total unnormalised energy */
  double e0_sum = 0;
  /* cumulative normalised energies */
  std::vector<double> E0(ndata);

  /* set raw energies and Z (E0_sum) */
  for (size_t i = 0; i < ndata; ++i)
  {
    set_sampleID_sampleID_distance_nothreshold(index0, i, adistance);
    e0[i] = f_energy(adistance);
    e0_sum += e0[i];
  }

  E0[0] = e0[0];
  for (size_t i = 1; i < ndata; ++i)
  {
    E0[i] = E0[i - 1] + e0[i];
  }
  for (size_t i = 0; i < ndata; ++i)
  {
    E0[i] /= e0_sum;
  }

  /* will now sample 2*ndata samples from (almost) q of Bachem et al. */
  size_t              n_pseudo_samples = 2 * ndata;
  std::vector<size_t> pseudo_sample(n_pseudo_samples);

  double                                 d_ndata = static_cast<double>(ndata);
  std::uniform_real_distribution<double> dis_uni01_2(0.0, 1.0);
  double                                 rand_offset = dis_uni01_2(gen) / d_ndata;

  unsigned running_index = 0;
  for (size_t i = 0; i < ndata; ++i)
  {
    /* the uniform distribution component (with replacement, not what Bachem et al do but, but
     * good enough) */
    pseudo_sample[2 * i] = i;
    /* the non-uniform distribution component. Again, the sampling is not independent as in
     * Bachem et al, but good enough  */
    while (E0[running_index] < (static_cast<double>(i) + rand_offset) / d_ndata)
    {
      ++running_index;
    }
    pseudo_sample[2 * i + 1] = running_index;
  }

  /* shuffle the samples */
  size_t swap_i;
  for (size_t i = 0; i < ndata; ++i)
  {
    swap_i = dis(gen) % (ndata - i);
    std::swap(pseudo_sample[swap_i], pseudo_sample[ndata - 1 - i]);
  }

  size_t q_index = 0;
  /* x of Bachem et al */
  size_t sample_index_current;
  /* d_x of Bachem et al */
  double e_current;

  /* y of Bachem et al */
  size_t sample_index_proposal;
  /* d_y of Bachem et al */
  double e_proposal;

  for (size_t k = 1; k < K; ++k)
  {

    do
    {
      q_index += 1;
      q_index %= n_pseudo_samples;
      sample_index_current = pseudo_sample[q_index];
    } while (is_center[sample_index_current] == true);

    e_current = std::numeric_limits<double>::max();
    for (size_t kp = 0; kp < k; ++kp)
    {
      set_sampleID_sampleID_distance_nothreshold(
        sample_index_current, center_indices_init[kp], adistance);
      e_current = std::min<double>(e_current, f_energy(adistance));
    }

    for (size_t m = 1; m < chain_length; ++m)
    {
      do
      {
        q_index += 1;
        q_index %= n_pseudo_samples;
        sample_index_proposal = pseudo_sample[q_index];
      } while (is_center[sample_index_proposal] == true);

      e_proposal = std::numeric_limits<double>::max();
      for (size_t kp = 0; kp < k; ++kp)
      {
        set_sampleID_sampleID_distance_nothreshold(
          sample_index_proposal, center_indices_init[kp], adistance);
        e_proposal = std::min<double>(e_proposal, f_energy(adistance));
      }

      if ((e_proposal / e_current) * ((e0[sample_index_current] * 2 * ndata + e0_sum) /
                                      (e0[sample_index_proposal] * 2 * ndata + e0_sum)) >
          dis_uni01_2(gen))
      {
        e_current            = e_proposal;
        sample_index_current = sample_index_proposal;
      }
    }
    is_center[sample_index_current] = true;
    center_indices_init[k]          = sample_index_current;
  }
}

// single threaded redistribution. mulithreading here will be very tricky.
// requires that to_leave_cluster be reliably set.
// centers are unchanged by this function. this function simply moves samples between clusters.
// additional complexity required to maintain randomness (immigrant samples inserted into random
// indices)
void SkeletonClusterer::redistribute()
{

  size_t k_new;
  size_t j;
  size_t j_new;

  std::vector<size_t> redistribute_order(K, 0);
  set_redistribute_order(redistribute_order);

  std::vector<size_t> redistribute_rank(K, 0);
  for (size_t k = 0; k < K; ++k)
  {
    redistribute_rank[redistribute_order[k]] = k;
  }

  bool target_is_a_mover;
  for (auto& k : redistribute_order)
  {
    if (nthreads > 1)
    {
      std::sort(to_leave_cluster[k].begin(), to_leave_cluster[k].end());
    }
    // thought : it would be nice if one could auto over a vector in reverse
    for (size_t ji = to_leave_cluster[k].size(); ji-- > 0;)
    {

      j = to_leave_cluster[k][ji];

      /* size_t max indicates that the mover has moved to the tail, because of an earlier
       * immigration */
      if (j != std::numeric_limits<size_t>::max())
      {
        k_new = nearest_1_infos[k][j].a_x;

        cluster_has_changed[k]     = true;
        cluster_has_changed[k_new] = true;

        // now for the redistribution
        // (1) make space somewhere in k_new (k,j) will be inserted at k_new, j_new:
        j_new = dis(gen) % (get_ndata(k_new) + 1);

        /* try to swap with a non-mover, as it makes swapping quicker. but if after 10 attempts
         * no luck,
         * can still swap with movers */
        size_t insert_attempts = 0;
        while (j_new != get_ndata(k_new) && nearest_1_infos[k_new][j_new].a_x != k_new &&
               insert_attempts < 10)
        {
          j_new = dis(gen) % (get_ndata(k_new) + 1);
          ++insert_attempts;
        }

        /* case 1 : k,j goes on the tail */
        if (j_new == get_ndata(k_new))
        {
          nearest_1_infos[k_new].push_back(nearest_1_infos[k][j]);
          sample_IDs[k_new].push_back(sample_IDs[k][j]);
          append_across(k_new, k, j);
          custom_append(k_new, k, j);
        }

        /* case 2 : k_new, j_new goes on the tail, then k,j goes in at k_new, j_new */
        else
        {

          if (nearest_1_infos[k_new][j_new].a_x != k_new)
          {
            target_is_a_mover = true;
          }

          else
          {
            target_is_a_mover = false;
          }

          /* putting k_new, j_new on the tail to make space for k, j */
          nearest_1_infos[k_new].push_back(nearest_1_infos[k_new][j_new]);
          sample_IDs[k_new].push_back(sample_IDs[k_new][j_new]);

          append_across(k_new, k_new, j_new);
          custom_append(k_new, k_new, j_new);

          /* putting k, j in where space has been made for it */
          nearest_1_infos[k_new][j_new] = nearest_1_infos[k][j];
          sample_IDs[k_new][j_new]      = sample_IDs[k][j];
          replace_with(k_new, j_new, k, j);

          custom_replace_with(k_new, j_new, k, j);

          /* k_new, j_new is still going to move when k_new immigration starts.
           * Record that k_new, j_new has moved to tail, and leave a notive (max size_t)
           * so that when k_new immigration starts we know */
          if (target_is_a_mover && redistribute_rank[k_new] > redistribute_rank[k])
          {
            /* if to_leave_cluster is sorted we can improve this search. */
            for (auto& x : to_leave_cluster[k_new])
            {
              if (x == j_new)
              {
                x = std::numeric_limits<size_t>::max();
                to_leave_cluster[k_new].push_back(get_ndata(k_new) - 1);
                break;
              }
            }
          }
        }

        // swap with tail.
        remove_with_tail_pull(k, j);

      }  // j !=  max size_t
    }
  }
}

void SkeletonClusterer::set_all_cluster_statistics()
{
  E_total = 0;
  for (size_t k = 0; k < K; ++k)
  {
    set_cluster_statistics(k);
    cluster_has_changed[k] = false;
    E_total += cluster_energies[k];
  }
}

void SkeletonClusterer::update_all_cluster_statistics()
{
  old_E_total = E_total;
  E_total     = 0;
  for (size_t k = 0; k < K; ++k)
  {
    if (cluster_has_changed[k] == true)
    {

      set_cluster_statistics(k);
      cluster_has_changed[k] = false;
    }
    E_total += cluster_energies[k];
  }
}

/* set energy of cluster k, and set as custom cluster statistics */
void SkeletonClusterer::set_cluster_statistics(size_t k)
{
  cluster_energies[k] = 0;
  set_to_zero_custom_cluster_statistics(k);
  for (size_t j = 0; j < get_ndata(k); ++j)
  {
    cluster_energies[k] += nearest_1_infos[k][j].e_x;
    increment_custom_cluster_statistics(k, j);  // includes max
  }
  cluster_mean_energies[k] = cluster_energies[k] / static_cast<double>(get_ndata(k));
  set_normalised_custom_cluster_statistics(k);
}

// remove the j'th sample from cluster k, and (if it is not the last element) fill the hole with
// the tail (size drops by 1)
void SkeletonClusterer::remove_with_tail_pull(size_t k, size_t j)
{
  if (j != get_ndata(k) - 1)
  {
    nearest_1_infos[k][j] = *(nearest_1_infos[k].end() - 1);
    sample_IDs[k][j]      = *(sample_IDs[k].end() - 1);
    replace_with_last_element(k, j);
    custom_replace_with_last(k, j);
  }

  nearest_1_infos[k].pop_back();
  sample_IDs[k].pop_back();
  remove_last(k);
  custom_remove_last(k);
  cluster_has_changed[k] = true;
}

void SkeletonClusterer::initialise_all()
{

  if (with_tests == true)
  {
    fundamental_triangle_inequality_test();
  }

  tstart_initialise_centers = std::chrono::high_resolution_clock::now();
  time_prehistory =
    std::chrono::duration_cast<std::chrono::microseconds>(tstart_initialise_centers - bigbang)
      .count();

  initialise_center_indices();

  for (size_t k = 0; k < K; ++k)
  {
    append_to_centers_from_ID(center_indices_init[k]);
    /* here we make an empty cluster, using datain to determine nec. `shape' parameters */
    add_empty_cluster();
  }

  /* if the indices are from the user or in debug mode, we run a test that initialising indices
   * look fine. */
  if (with_tests || initialisation_method == "from_indices_init")
  {
    post_initialise_centers_test();
  }

  auto t_endit = std::chrono::high_resolution_clock::now();
  time_to_initialise_centers =
    std::chrono::duration_cast<std::chrono::microseconds>(t_endit - tstart_initialise_centers)
      .count();
  tstart = std::chrono::high_resolution_clock::now();

  // initialisation
  set_center_center_info();
  put_samples_in_clusters();
  set_all_cluster_statistics();

  /* checkpoint: all assignments and cluster statistics must be correct. Tests to pass : all */
  if (with_tests == true)
  {
    post_initialisation_test();
  }

  auto t1             = std::chrono::high_resolution_clock::now();
  time_initialising   = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
  ncalcs_initialising = get_ncalcs();
  mowri << get_equals_line(get_round_summary().size());
  mowri << get_round_summary() << zentas::Endl;
}

void SkeletonClusterer::put_samples_in_clusters()
{

  // get sorted center IDs.
  std::unique_ptr<size_t[]> up_sorted_center_IDs(new size_t[K]);
  auto                      sorted_center_IDs = up_sorted_center_IDs.get();
  std::copy(center_IDs, center_IDs + K, sorted_center_IDs);
  std::sort(sorted_center_IDs, sorted_center_IDs + K);

  // make array with non center IDs
  std::vector<size_t> non_center_IDs(ndata - K);
  size_t              non_center_number = 0;
  size_t              center_number     = 0;

  for (size_t ID = 0; ID < ndata; ++ID)
  {
    if (center_number == K || ID != sorted_center_IDs[center_number])
    {
      non_center_IDs[non_center_number] = ID;
      ++non_center_number;
    }

    else
    {
      ++center_number;
    }
  }
  // shuffle non center IDs
  for (size_t i = ndata - K - 1; i > 0; --i)
  {
    size_t swap_with = dis(gen) % (i + 1);
    std::swap(non_center_IDs[i], non_center_IDs[swap_with]);
  }

  std::vector<std::thread> threads;
  for (size_t ti = 0; ti < nthreads; ++ti)
  {
    threads.emplace_back([this, ti, &non_center_IDs]() {
      pll_put_samples_in_cluster(get_start(ti, get_nthreads(), 0, ndata - K),
                                 get_end(ti, get_nthreads(), 0, ndata - K),
                                 non_center_IDs);
    });
  }
  for (auto& t : threads)
  {
    t.join();
  }
  kmoo_finish_with();
}

void SkeletonClusterer::pll_put_samples_in_cluster(size_t               i_a,
                                                   size_t               i_z,
                                                   std::vector<size_t>& non_center_IDs)
{

  std::string prefix = "kmeans++";
  if (initialisation_method.substr(0, prefix.size()) == prefix)
  {
    for (size_t i = i_a; i < i_z; ++i)
    {
      size_t nc_ID = non_center_IDs[i];
      size_t k1    = kmoo_p2bun.k_1(nc_ID);
      size_t k2    = kmoo_p2bun.k_2(nc_ID);
      double d1    = kmoo_p2bun.d_1(nc_ID);
      double d2    = kmoo_p2bun.d_2(nc_ID);
      final_push_into_cluster_post_kmeanspp(nc_ID, k1, k2, d1, d2);
    }
  }

  else
  {
    for (size_t i = i_a; i < i_z; ++i)
    {
      put_sample_in_cluster(non_center_IDs[i]);
    }
  }
}

void SkeletonClusterer::final_push_into_cluster_post_kmeanspp(
  size_t i, size_t k1, size_t k2, double d1, double d2)
{
  std::lock_guard<std::mutex> lockraii(mutex0);
  final_push_into_cluster_basic(i, k1, d1);
  put_nearest_2_infos_margin_in_cluster_post_kmeanspp(k1, k2, d2, f_energy(d2));
}

void SkeletonClusterer::final_push_into_cluster_basic(size_t i,
                                                      size_t nearest_center,
                                                      double min_distance)
{
  append_from_ID(nearest_center, i);
  nearest_1_infos[nearest_center].emplace_back(
    nearest_center, min_distance, f_energy(min_distance));
  sample_IDs[nearest_center].push_back(i);
}

void SkeletonClusterer::reset_multiple_sample_infos(size_t k_to, size_t j_a, size_t j_z)
{
  for (size_t j = j_a; j < j_z; ++j)
  {
    reset_sample_infos(k_to, j);
  }
}

void SkeletonClusterer::base_put_sample_in_cluster(size_t i)
{
  std::unique_ptr<double[]> up_distances(new double[K]);
  double                    min_distance        = std::numeric_limits<double>::max();
  size_t                    nearest_center      = 0;
  double                    second_min_distance = std::numeric_limits<double>::max();

  for (size_t k = 0; k < K; ++k)
  {
    set_center_sampleID_distance(k, i, second_min_distance, up_distances[k]);
    if (up_distances[k] <= second_min_distance)
    {
      if (up_distances[k] < min_distance)
      {
        second_min_distance = min_distance;
        nearest_center      = k;
        min_distance        = up_distances[k];
      }
      else
      {
        second_min_distance = up_distances[k];
      }
    }
  }
  final_push_into_cluster(i, nearest_center, min_distance, up_distances.get());
}

void SkeletonClusterer::triangular_put_sample_in_cluster(size_t i, const double* const cc)
{
  std::unique_ptr<double[]> up_distances(new double[K]);
  double                    min_distance        = std::numeric_limits<double>::max();
  double                    second_min_distance = std::numeric_limits<double>::max();

  size_t nearest_center = 0;

  for (size_t k = 0; k < K; ++k)
  {
    up_distances[k] = std::numeric_limits<double>::max();
    if (cc[nearest_center * K + k] < second_min_distance + min_distance)
    {
      set_center_sampleID_distance(k, i, second_min_distance, up_distances[k]);
      if (up_distances[k] < second_min_distance)
      {
        if (up_distances[k] < min_distance)
        {
          second_min_distance = min_distance;
          nearest_center      = k;
          min_distance        = up_distances[k];
        }
        else
        {
          second_min_distance = up_distances[k];
        }
      }
    }
  }
  final_push_into_cluster(i, nearest_center, min_distance, up_distances.get());
}

void SkeletonClusterer::swap_center_with_sample(size_t k, size_t j)
{
  cluster_has_changed[k] = true;
  size_t c_id_k          = center_IDs[k];
  center_IDs[k]          = sample_IDs[k][j];
  /* jn : bug fix, voronoi 11 june 2017. the following was missing: */
  sample_IDs[k][j] = c_id_k;
  swap_center_data(k, j);
  reset_sample_infos_basic(k, j);
}

void SkeletonClusterer::move_center_into_its_own_cluster(size_t k)
{
  // move it from the original data (more likely to have a cache miss, but for now less code =
  // good)
  // note that we use the base_version of the function, at the overloaded version (for clarans
  // l23) uses cc which might not
  // be reliable at this point.
  base_put_sample_in_cluster(center_IDs[k]);
  center_IDs[k]          = 911911;
  cluster_has_changed[k] = true;
}

size_t SkeletonClusterer::draw_k_uniform() { return dis(gen) % K; }

size_t SkeletonClusterer::draw_k_prop_ndata()
{
  std::unique_ptr<size_t[]> up_cum_ndatas(new size_t[K]);
  size_t                    cum_ndata = 0;
  for (size_t k = 0; k < K; ++k)
  {
    cum_ndata += get_ndata(k);
    up_cum_ndatas[k] = cum_ndata;
  }
  if (up_cum_ndatas[K - 1] != (ndata - K))
  {
    mowri << "Weird : cum_ndatas[K-1] != ndata : up_cum_ndatas[K-1] = " << up_cum_ndatas[K - 1]
          << " and ndata = " << ndata << zentas::Endl;
    throw zentas::zentas_error("(see above)");
  }
  size_t i       = dis(gen) % (ndata - K);
  size_t k_below = 0;
  while (i >= up_cum_ndatas[k_below])
  {
    ++k_below;
  }

  return k_below;
}

/* basic, because has no effect on to_leave_cluster */
void SkeletonClusterer::reset_sample_infos_basic(size_t k, size_t j)
{

  std::unique_ptr<double[]> up_distances(new double[K]);
  double                    min_distance        = std::numeric_limits<double>::max();
  size_t                    nearest_center      = 0;
  double                    second_min_distance = std::numeric_limits<double>::max();

  for (size_t kp = 0; kp < K; ++kp)
  {

    set_center_sample_distance(kp, k, j, second_min_distance, up_distances[kp]);

    if (up_distances[kp] < second_min_distance)
    {
      if (up_distances[kp] < min_distance)
      {
        nearest_center = kp;
        min_distance   = up_distances[kp];
      }
      else
      {
        second_min_distance = up_distances[kp];
      }
    }
  }

  reset_nearest_info_basic(k, j, nearest_center, min_distance, f_energy(min_distance));
  reset_sample_custom(k, j, nearest_center, up_distances.get());
}

void SkeletonClusterer::reset_sample_infos(size_t k, size_t j)
{

  reset_sample_infos_basic(k, j);
  if (k != nearest_1_infos[k][j].a_x)
  {
    std::lock_guard<std::mutex> lock(mutex0);
    to_leave_cluster[k].push_back(j);
  }
}

void SkeletonClusterer::final_push_into_cluster(size_t              i,
                                                size_t              nearest_center,
                                                double              min_distance,
                                                const double* const distances)
{
  // get your lock on, time for polyphonics.
  std::lock_guard<std::mutex> lockraii(mutex0);

  /* the common part of pushing into a cluster */
  final_push_into_cluster_basic(i, nearest_center, min_distance);

  /* as usual: the final parameter up_distances.get() guarantees correctness
   * only of lowest and second lowest, other values may exceed */
  /* TODO : there may be computation in here which should not be under the lock */
  put_sample_custom_in_cluster(i, nearest_center, distances);
}

void SkeletonClusterer::overwrite_center_with_sample(size_t k1, size_t k2, size_t j2)
{
  centers_replace_with_sample(k1, k2, j2);
  center_IDs[k1]          = sample_IDs[k2][j2];
  cluster_has_changed[k1] = true;
}

void SkeletonClusterer::initialise_center_indices()
{
  /* initialisation from indices. */
  if (initialisation_method == "from_indices_init")
  {
    // already done in constructor
  }

  /* initialisation uniformly. */
  else if (initialisation_method == "uniform")
  {
    init::populate_uniformly(center_indices_init, K, ndata, dis, gen);
  }

  /* initialisation uniformly according to Bachem et al. */
  else if (initialisation_method.substr(0, 8) == "afk-mc2-")
  {
    populate_afk_mc2();
  }

  /* initialisation with kmeans++ */
  else if (initialisation_method == "kmeans++" ||
           initialisation_method.substr(0, 9) == "kmeans++-")
  {
    initialise_with_kmeanspp();
  }

  else
  {
    std::stringstream vims_ss;
    vims_ss << "The valid strings for initialisation_method are [from_indices_init, uniform, "
               "kmeans++-INT, afk-mc2-INT (for some positive INT)]";
    throw zentas::zentas_error(vims_ss.str());
  }

  std::sort(v_center_indices_init.begin(), v_center_indices_init.end());
  std::copy(center_indices_init, center_indices_init + K, center_IDs);
}

}  // samenpace
