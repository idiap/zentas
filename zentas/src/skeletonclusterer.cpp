#include "skeletonclusterer.hpp"


namespace nszen{

void P2Bundle::initialise_data(size_t n){
  ndata = n;
  d12_.reset(new std::array<double, 2> [n]);
  k_1_.reset(new size_t[n]);
  k_2_.reset(new size_t[n]);
  ori_.reset(new size_t[n]);
  
  for (unsigned i = 0; i < n; ++i){
    d_1(i) = std::numeric_limits<double>::max();
    d_2(i) = std::numeric_limits<double>::max();
  }
  for (unsigned i = 0; i < n; ++i){
    /* note that in the first round of kmeanspp, the lookup in cc requires this to be zero  */ 
    k_1(i) = 0; 
  }
}

P2Bundle::P2Bundle(size_t n){
  initialise_data(n);
}

P2Bundle::P2Bundle(const std::vector<size_t> & ori){
  initialise_data(ori.size());
  for (size_t i = 0; i < ndata; ++i){
    ori_[i] = ori[i];
  }
}


void XNearestInfo::reset(size_t new_a_x, double new_d_x, double new_e_x){
  a_x = new_a_x;
  d_x = new_d_x;
  e_x = new_e_x;
}

void XNearestInfo::reset(XNearestInfo & nearest_x_infos){
  a_x = nearest_x_infos.a_x;
  d_x = nearest_x_infos.d_x;
  e_x = nearest_x_infos.e_x;
}

XNearestInfo::XNearestInfo(size_t a_x, double d_x, double e_x):a_x(a_x), d_x(d_x), e_x(e_x){}

std::string XNearestInfo::get_string(){
  return  std::to_string(a_x) + "\t " + std::to_string(d_x) + "\t " + std::to_string(e_x);
}


SkeletonClusterer::SkeletonClusterer(size_t K, std::chrono::time_point<std::chrono::high_resolution_clock> bigbang, const size_t * const center_indices_init_, size_t ndata_in, std::string init_method, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t max_rounds, std::string energy, bool with_tests, size_t random_sd, const EnergyInitialiser & energy_initialiser):
    
  mowri(true, false, ""), K(K), bigbang(bigbang), ndata(ndata_in), initialisation_method(init_method), nearest_1_infos(K), sample_IDs(K), to_leave_cluster(K), cluster_has_changed(K, true), cluster_energies(K,0), cluster_mean_energies(K), E_total(std::numeric_limits<double>::max()), old_E_total(0), round(0), v_center_indices_init(K), center_indices_init(v_center_indices_init.data()), /*bigbang(bigbang),*/ max_time_micros(static_cast<size_t>(max_time*1000000.)), min_mE(min_mE), labels(labels), nthreads(nthreads), nthreads_fl(static_cast<double> (nthreads)), max_rounds(max_rounds), energy(energy), with_tests(with_tests), gen(random_sd) {
  
  /* TODO here: set-up the center refiner function. 
   * Takes in center by reference and samples by const reference? */ 
  
  if (energy.compare("identity") == 0){
    f_energy = nszen::Identity(); 
  } 
   
  else if (energy.compare("quadratic") == 0){
    f_energy = nszen::Quadratic();
  }
  
  else if (energy.compare("cubic") == 0){
    f_energy = nszen::Cubic();
  }
  
  else if (energy.compare("squarepotential") == 0){
    f_energy = nszen::SquarePotential(energy_initialiser.get_critical_radius());
  }
  
  else if (energy.compare("log") == 0){
    f_energy = nszen::Log();
  }
  
  else if (energy.compare("exp") == 0){
    f_energy = nszen::Exponential(energy_initialiser.get_exponent_coeff());
  }
  
  else if (energy.compare("sqrt") == 0){
    f_energy = nszen::SquareRoot();
  }
  
  else{
    std::stringstream errmss;
    errmss << "Unrecognised energy function, `" << energy << "'. ";
    if (energy == "linear"){
      errmss << "Perhaps instead of `linear' you meant `identity'?";
    }
    throw zentas::zentas_error(std::string() + energy);
  }
  
  
  /* confirm that f_energy(0) is 0 */
  if (f_energy(0) != 0){
    std::stringstream ss;
    ss << "the energy function, f_energy, has f_energy(0) = ";
    ss << f_energy(0) << ". This is a problem for k-means++ initialisation.";
    ss << "If you're not using k-means++ initialisation, this should not cause any problems, ";
    ss << "but as k-means++ is the default, we are being cautious and throwing an error.";
    throw zentas::zentas_error(ss.str());
    
  }

  /* initialisation from indices. */
  if (initialisation_method == "from_indices_init"){
    populate_from_indices_init(center_indices_init_, center_indices_init, K, ndata);
  }
  
  else{
    //will do in go ( )
  }

  center_IDs = indices_final;

}









    size_t SkeletonClusterer::get_time_in_update_centers(){
      return time_in_update_centers;
    }
    
    size_t SkeletonClusterer::get_time_in_update_sample_info(){
      return time_in_update_sample_info;
    }
    
    double SkeletonClusterer::get_time_remaining(){
      auto t1 = std::chrono::high_resolution_clock::now();
      time_total = std::chrono::duration_cast<std::chrono::microseconds>(t1 - bigbang).count();        
      
      if (max_time_micros > time_total){
        return max_time_micros - time_total;
      }
      else{
        return -1;
      }
    }


    size_t SkeletonClusterer::get_start(size_t ti, size_t nthreads, size_t j_A, size_t j_Z){
      size_t n_js = j_Z - j_A;
      double t_fl = static_cast<double>(ti);
      size_t j_a = j_A + (t_fl / static_cast<double> (nthreads)) * n_js;
      return j_a;
    }

    size_t SkeletonClusterer::get_end(size_t ti, size_t nthreads, size_t j_A, size_t j_Z){
      size_t n_js = j_Z - j_A;
      double t_fl = static_cast<double>(ti);
      size_t j_z = j_A + ( ( t_fl + 1.) / static_cast<double> (nthreads)) * n_js;
      return j_z;
    }

    size_t SkeletonClusterer::get_sample_from(std::vector<double> & v_cum_nearest_energies){
      
      /* kmeans++ is exhausted : everything sample has a center at distance 0 (obviously duplicated data)
       * however, this case should have been caught upstream, so throwing an error. */
      if (v_cum_nearest_energies.back() == 0){
        throw zentas::zentas_error("exhausted in get_sample_from (kmeans++). in particular, the cumulative energy is zero. this should have been caught upstream: logic error in zentas");  
      }
      
      return std::distance(
      v_cum_nearest_energies.begin(), 
      std::lower_bound(v_cum_nearest_energies.begin(), v_cum_nearest_energies.end(), dis_uni01(gen)*v_cum_nearest_energies.back())) - 1;
    }

    std::string SkeletonClusterer::get_base_summary_string(){
      
      auto t_now = std::chrono::high_resolution_clock::now();
      time_total = std::chrono::duration_cast<std::chrono::microseconds>(t_now - bigbang).count();
      ncalcs_total = get_ncalcs();
                
      
      std::string st_round = "R=" + std::to_string(round);
      st_round.resize(3 + 6, ' ');
      
      std::stringstream st_mE_ss;
      st_mE_ss << "mE=" << std::setprecision(7) << E_total / static_cast<double>(ndata);
      std::string st_mE = st_mE_ss.str();
      st_mE.resize(std::max<size_t>(st_mE.size() + 1,3 + 11), ' ');

      std::string st_Tp = "Tp=" + std::to_string(time_prehistory/1000);
      st_Tp += "  ";

      std::string st_Ti = "Ti=" + std::to_string(time_to_initialise_centers/1000);
      st_Ti += "  ";
      
      std::string st_Tb = "Tb=" + std::to_string(time_initialising/1000);
      st_Tb += "  ";
      
      std::string st_Tc = "Tc=" + std::to_string(time_in_update_centers/1000);
      st_Tc.resize(3 + 8, ' ');
      
      std::string st_Tu = "Tu=" + std::to_string(time_in_update_sample_info/1000);
      st_Tu.resize(3 + 7, ' ');
      
      std::string st_Tr = "Tr=" + std::to_string(time_in_redistribute/1000);
      st_Tr.resize(3 + 5, ' ');

      
      std::string st_Tt = "Tt="  + std::to_string(time_total/1000);
      st_Tt.resize(3 + 8, ' ');


      std::stringstream ncc_ss;
      ncc_ss <<  "lg2nc(c)="  << std::setprecision(5) << std::log2(ncalcs_in_update_centers);
      std::string ncc = ncc_ss.str();
      ncc.resize(9+9, ' ');


      std::stringstream nc_ss;
      nc_ss <<  "lg2nc="  << std::setprecision(5) << std::log2(ncalcs_total);
      std::string nc = nc_ss.str();
      nc.resize(6+9, ' ');


      std::stringstream pc_ss;
      pc_ss <<  "pc="  << std::setprecision(5) << static_cast<double> (get_rel_calccosts());
      std::string pc = pc_ss.str();
      pc.resize(3+10, ' ');

      
      std::ostringstream out;
      out << st_round << st_mE << st_Tp << st_Ti << st_Tb << st_Tc << st_Tu << st_Tr << st_Tt << ncc << nc << pc;
      return out.str();
    

    }
    


    /* rule : functions with suffix 'basic' will not touch to_leave_cluster */
    void SkeletonClusterer::reset_nearest_info_basic(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest){
      nearest_1_infos[k][j].reset(k_nearest, d_nearest, e_nearest);      
      cluster_has_changed[k] = true;
    }
    
    void SkeletonClusterer::reset_nearest_info(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest){
      reset_nearest_info_basic(k, j, k_nearest, d_nearest, e_nearest);      
      if (k != k_nearest){
        std::lock_guard<std::mutex> lock(mutex0);
        to_leave_cluster[k].push_back(j);
      }
    }
    

}
