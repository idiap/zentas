// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/vdimap.hpp>
#include <zentas/zentaserror.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

namespace nszen
{
namespace vdimap
{

/* make the e'th eigenvector normalised and perpendicular to all previous ones. */
template <typename T>
void make_one_nperp(T* const evs, size_t dim, size_t e)
{

  /* normalise */
  T l2_norm = 0;
  for (size_t dp = 0; dp < dim; ++dp)
  {
    l2_norm += evs[e * dim + dp] * evs[e * dim + dp];
  }
  l2_norm = std::sqrt(l2_norm);

  for (size_t dp = 0; dp < dim; ++dp)
  {
    evs[e * dim + dp] /= l2_norm;
  }

  /* make perpendicular to preceding vectors */
  for (size_t ep = 0; ep < e; ++ep)
  {
    T inner_eigens = 0;
    for (size_t dp = 0; dp < dim; ++dp)
    {
      inner_eigens += evs[e * dim + dp] * evs[ep * dim + dp];
    }
    T norm2_ = 0;
    for (size_t dp = 0; dp < dim; ++dp)
    {
      evs[e * dim + dp] -= inner_eigens * evs[ep * dim + dp];
      norm2_ += evs[e * dim + dp] * evs[e * dim + dp];
    }
    norm2_ = std::sqrt(norm2_);
    for (size_t dp = 0; dp < dim; ++dp)
    {
      evs[e * dim + dp] /= norm2_;
    }
  }
}

template <typename T>
void make_all_nperp(T* const evs, size_t dim, size_t neigs)
{
  for (size_t e = 0; e < neigs; ++e)
  {
    make_one_nperp(evs, dim, e);
  }
}

template void make_all_nperp<float>(float* const evs, size_t dim, size_t neigs);
template void make_all_nperp<double>(double* const evs, size_t dim, size_t neigs);

template <typename T>
void confirm_all_nperp(T* const evs, size_t dim, size_t neigs)
{
  for (size_t e = 0; e < neigs; ++e)
  {

    double l2_norm = 0;
    for (size_t dp = 0; dp < dim; ++dp)
    {
      l2_norm += evs[e * dim + dp] * evs[e * dim + dp];
    }

    if (std::abs(l2_norm - 1.) > 1e-4)
    {
      std::stringstream ss;
      ss << "norm2 is not 1 confirm_all_nperp. ";
      ss << l2_norm;
      throw zentas::zentas_error(ss.str());
    }

    for (size_t ep = 0; ep < e; ++ep)
    {
      double inner = 0;
      for (size_t dp = 0; dp < dim; ++dp)
      {
        inner += evs[e * dim + dp] * evs[ep * dim + dp];
      }

      if (std::abs(inner) > 1e-4)
      {
        std::stringstream ss;
        ss << "is not zero in confirm_all_nperp. ";
        ss << std::abs(inner);
        throw zentas::zentas_error(ss.str());
      }
    }
  }
}

template <typename T>
void set_mean_variance(std::vector<T>& mean,
                       std::vector<T>& variance,
                       size_t          ndata,
                       size_t          dimension,
                       const T* const  ptr_datain)
{
  mean.resize(dimension, 0);
  variance.resize(dimension, 0);

  for (size_t i = 0; i < ndata; ++i)
  {
    for (size_t d = 0; d < dimension; ++d)
    {
      mean[d] += ptr_datain[i * dimension + d];
      variance[d] += ptr_datain[i * dimension + d] * ptr_datain[i * dimension + d];
    }
  }

  for (size_t d = 0; d < dimension; ++d)
  {
    variance[d] -= mean[d] * mean[d] / ndata;
    variance[d] /= ndata;
    mean[d] /= static_cast<T>(ndata);
  }
}

template <typename T>
std::vector<size_t> get_indices_decreasing(size_t dimension, const T* const variance)
{
  std::vector<size_t> indices_decreasing;
  std::vector<std::tuple<T, size_t>> v_i(dimension);
  for (size_t d = 0; d < dimension; ++d)
  {
    v_i[d] = std::make_tuple(variance[d], d);
  }

  std::sort(v_i.begin(), v_i.end(), [](std::tuple<T, size_t> a, std::tuple<T, size_t> b) {
    return std::get<0>(a) > std::get<0>(b);
  });
  if (std::get<0>(v_i.at(0)) < std::get<0>(v_i.back()))
  {
    throw zentas::zentas_error("logic error : the variance is highest variance dimension is lower "
                               "than the variance of the lowest variance dimension ");
  }

  /* the data is constant in all dimensions */
  if (std::get<0>(v_i[0]) == 0)
  {
    throw zentas::zentas_error("the data appears to be constant. this can be handled, but it "
                               "seems sufficiently weird that i'm bailing ");
  }

  else
  {
    for (size_t dp = 0; dp < dimension; ++dp)
    {
      /* if the variance in the dimension is not zero, plop it on the back. */
      if (std::get<0>(v_i[dp]) > T(1e-16) * std::get<0>(v_i[0]))
      {
        indices_decreasing.push_back(std::get<1>(v_i[dp]));
      }
    }
  }
  return indices_decreasing;
}

template <typename T>
void set_unimat(size_t                                 neigs,
                size_t                                 dim,
                T* const                               evs,
                std::uniform_int_distribution<size_t>& dis,
                std::default_random_engine&            gen)
{
  for (size_t e = 0; e < neigs; ++e)
  {
    for (size_t dp = 0; dp < dim; ++dp)
    {
      T rv         = T(1e-3) * T(dis(gen) % 10000 - 5000);
      evs[e * dim + dp] = rv;
    }
  }
  make_all_nperp(evs, dim, neigs);
}

template <typename T>
void power_method(T* const       eigenvectors,
                  T* const       eigenvector_norms,
                  size_t         neigs,
                  size_t         dim,
                  size_t         n_samples_for_covariance,
                  const T* const cov_samples,
                  size_t         n_power_iterations)
{

  std::vector<T> sigma_v(dim, 0);
  T              inner_product;

  std::vector<T> sigma_matrix;

  /* if you perform sigma v implicitly, looping through the covariance samples : */
  size_t flops_by_product = neigs * dim * n_samples_for_covariance * n_power_iterations;

  /* if you precompute the sigma matrix from the samples, and then explicity use then */
  size_t flops_by_matrix =
    dim * dim * n_samples_for_covariance + dim * dim * neigs * n_power_iterations;

  bool product_by_samples = (flops_by_product < flops_by_matrix);

  /* explicit sigma matrix */
  if (product_by_samples == false)
  {
    sigma_matrix.resize(dim * dim, 0);

    for (size_t si = 0; si < n_samples_for_covariance; ++si)
    {
      for (size_t d1 = 0; d1 < dim; ++d1)
      {
        for (size_t d2 = 0; d2 <= d1; ++d2)
        {
          sigma_matrix[d1 * dim + d2] += cov_samples[si * dim + d1] * cov_samples[si * dim + d2];
        }
      }
    }
    for (size_t d1 = 0; d1 < dim; ++d1)
    {
      for (size_t d2 = 0; d2 <= d1; ++d2)
      {
        sigma_matrix[d1 * dim + d2] *= (1. / n_samples_for_covariance);
        sigma_matrix[d2 * dim + d1] = sigma_matrix[d1 * dim + d2];
      }
    }
  }

  /* we compute eigenvectors in parallel.
   * this means better data reuse, fewer cache misses.
   * to do the standard serial way, just switch the
   * first two for loops here.
   * */
  for (size_t pi = 0; pi < n_power_iterations; ++pi)
  {
    for (size_t e = 0; e < neigs; ++e)
    {

      make_one_nperp(eigenvectors, dim, e);
      for (auto& x : sigma_v)
      {
        x = 0;
      }

      /* compute sigma_hat eigen_e, implicitly */
      if (product_by_samples == true)
      {
        for (size_t si = 0; si < n_samples_for_covariance; ++si)
        {
          inner_product = 0;
          for (size_t dp = 0; dp < dim; ++dp)
          {
            inner_product += eigenvectors[e * dim + dp] * cov_samples[si * dim + dp];
          }

          for (size_t dp = 0; dp < dim; ++dp)
          {
            sigma_v[dp] +=
              T(T(1) / n_samples_for_covariance) * inner_product * cov_samples[si * dim + dp];
          }
        }
      }

      /* compute sigma_hat eigen_e, explicitly */
      else
      {
        for (size_t d1 = 0; d1 < dim; ++d1)
        {
          sigma_v[d1] = 0;
          for (size_t d2 = 0; d2 < dim; ++d2)
          {
            sigma_v[d1] += eigenvectors[e * dim + d2] * sigma_matrix[d1 * dim + d2];
          }
        }
      }

      /* update eigenvalue and eigenvector */
      eigenvector_norms[e] = 0;
      for (size_t dp = 0; dp < dim; ++dp)
      {
        eigenvector_norms[e] += sigma_v[dp] * sigma_v[dp];
      }
      eigenvector_norms[e] = std::sqrt(eigenvector_norms[e]);

      for (size_t dp = 0; dp < dim; ++dp)
      {
        eigenvectors[e * dim + dp] = sigma_v[dp];
      }
      make_one_nperp(eigenvectors, dim, e);
    }
  }
}

template <typename T>
void vdimap(
  std::vector<T>& v_mapped, const T* const ptr_datain, size_t ndata, size_t dimension, size_t seed)
{

  bool print_times = false;

  /* timing  */
  auto t0 = std::chrono::high_resolution_clock::now();

  std::uniform_int_distribution<size_t> dis;
  std::default_random_engine            gen(seed);

  /* (1) set the mean and variance of all the variables */
  std::vector<T> mean;
  std::vector<T> variance;
  set_mean_variance(mean, variance, ndata, dimension, ptr_datain);

  if (print_times)
  {
    auto t00             = std::chrono::high_resolution_clock::now();
    auto time_mean_var00 = std::chrono::duration_cast<std::chrono::microseconds>(t00 - t0).count();
    std::cout << "Time first set mean variance : " << std::to_string(time_mean_var00 / 1000)
              << std::endl;
  }

  /* (2) get the order of the indices, in decreasing order of variance.
   * note that all zero or very small with respect to largest variance
   * are not included here */
  std::vector<size_t> indices_decreasing = get_indices_decreasing<T>(dimension, variance.data());

  /* the number of non-zero indices */
  size_t dim_nz = indices_decreasing.size();

  /* TODO : how to set ? */
  size_t neigs = std::min<size_t>(dim_nz / 5, 5);

  /* TODO : how to set ?
   * we use 20*neigs to approximate the covariance
   * matrix. For now, matrix multiply is done implicitly,
   * although when d < 10*neigs it should be explicit (TODO)  */
  size_t n_samples_for_covariance = std::min<size_t>(100 * neigs, ndata);

  /* TODO : how to set ? */
  size_t n_power_iterations = 30;

  /* samples to use for covariance, will be the first n_samples_for_covariance from this */
  std::vector<size_t> sample_indices(ndata);
  std::iota(sample_indices.begin(), sample_indices.end(), 0);
  for (size_t i = 0; i < n_samples_for_covariance; ++i)
  {
    std::swap(sample_indices[i], sample_indices[i + dis(gen) % (ndata - i)]);
  }
  std::sort(sample_indices.begin(), sample_indices.begin() + n_samples_for_covariance);

  /* The dimension of v_mapped */
  size_t dim_m = neigs + dim_nz;

  /* v_mapped will be of size dim_m x ndata.
   * the first neigs columns will be the innerproducts with eigenvectors,
   * the remaining will be the residual of the vectors, included so as to make
   * distances unchanged. */
  v_mapped.resize(dim_m * ndata);

  std::vector<T> cov_samples(n_samples_for_covariance * dim_nz);
  for (size_t i = 0; i < n_samples_for_covariance; ++i)
  {
    for (size_t dp = 0; dp < dim_nz; ++dp)
    {
      /* we implicitly perform the first transformation (by variance pre-eigen) */
      cov_samples[i * dim_nz + dp] =
        ptr_datain[sample_indices[i] * dimension + indices_decreasing[dp]] -
        mean[indices_decreasing[dp]];
    }
  }

  /* we initialise the eigenvectors and eigenvalues (to be obtained by power iteration) to uniform
   * noise */
  std::vector<T> eigenvectors(neigs * dim_nz, 0);
  std::vector<T> eigenvector_norms(neigs, 0);

  /* initialise the eigenvectors */
  set_unimat(neigs, dim_nz, eigenvectors.data(), dis, gen);

  if (print_times)
  {
    auto t1             = std::chrono::high_resolution_clock::now();
    auto time_pre_power = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    std::cout << "Time pre power : " << std::to_string(time_pre_power / 1000) << std::endl;
  }

  /* run the power method */
  power_method(eigenvectors.data(),
               eigenvector_norms.data(),
               neigs,
               dim_nz,
               n_samples_for_covariance,
               cov_samples.data(),
               n_power_iterations);

  if (print_times)
  {
    auto t2              = std::chrono::high_resolution_clock::now();
    auto time_post_power = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t0).count();
    std::cout << "Time post power : " << std::to_string(time_post_power / 1000) << std::endl;
  }

  /* we have the eigenvectors, now set first neigs dimensions of v_mapped, and remove from the
   * remaining dims */
  T inner_product;

  std::vector<T> final_variance(dim_m, 0);
  std::vector<T> final_mean(dim_m, 0);

  for (size_t i = 0; i < ndata; ++i)
  {

    /* the original transformation, delayed til here so that only a single pass through v_mapped.
     * this is the pre-eigen by variance shuffle.  */
    for (size_t dp = 0; dp < dim_nz; ++dp)
    {
      v_mapped[neigs + dim_m * i + dp] =
        ptr_datain[i * dimension + indices_decreasing[dp]] - mean[indices_decreasing[dp]];
    }

    for (unsigned e = 0; e < neigs; ++e)
    {
      inner_product = 0;
      for (unsigned dp = 0; dp < dim_nz; ++dp)
      {
        inner_product += eigenvectors[e * dim_nz + dp] * v_mapped[neigs + dim_m * i + dp];
      }
      v_mapped[dim_m * i + e] = inner_product;
      /* make just residual, can do here as eigenvectors are perp */
      for (unsigned dp = 0; dp < dim_nz; ++dp)
      {
        v_mapped[neigs + (dim_nz + neigs) * i + dp] -=
          inner_product * eigenvectors[e * dim_nz + dp];
      }
    }

    for (unsigned d = 0; d < dim_m; ++d)
    {
      final_mean[d] += v_mapped[i * dim_m + d];
      final_variance[d] += v_mapped[i * dim_m + d] * v_mapped[i * dim_m + d];
    }
  }

  for (unsigned d = 0; d < dim_m; ++d)
  {
    final_variance[d] -= (final_mean[d] / ndata) * final_mean[d];
  }

  if (print_times)
  {
    auto t3 = std::chrono::high_resolution_clock::now();
    auto time_post_residual =
      std::chrono::duration_cast<std::chrono::microseconds>(t3 - t0).count();
    std::cout << "Time post residual : " << std::to_string(time_post_residual / 1000) << std::endl;
  }

  /* do a final sort by variance */
  auto           final_order = get_indices_decreasing<T>(dim_m, final_variance.data());
  std::vector<T> v_temp(dim_m);

  for (size_t i = 0; i < ndata; ++i)
  {
    for (size_t d = 0; d < dim_m; ++d)
    {
      v_temp[d] = v_mapped[i * dim_m + final_order[d]];
    }
    for (size_t d = 0; d < dim_m; ++d)
    {
      v_mapped[i * dim_m + d] = v_temp[d];
    }
  }

  if (print_times)
  {
    auto t4 = std::chrono::high_resolution_clock::now();
    auto time_post_final_sort =
      std::chrono::duration_cast<std::chrono::microseconds>(t4 - t0).count();
    std::cout << "Time post final sort : " << std::to_string(time_post_final_sort / 1000)
              << std::endl;
  }

  /* another pass through data, really nec ?
   * If so, at least compute updated variance in loop above
   * TODO */

  /* done ! */
  /* confirmation */
  for (size_t p = 0; p < 7; ++p)
  {
    size_t i = dis(gen) % ndata;
    size_t j = dis(gen) % ndata;

    T d1 = 0;
    for (size_t d = 0; d < dimension; ++d)
    {
      T diff = ptr_datain[i * dimension + d] - ptr_datain[j * dimension + d];
      d1 += diff * diff;
    }
    d1 = std::sqrt(d1);

    T d2 = 0;
    for (size_t d = 0; d < dim_m; ++d)
    {
      T diff = v_mapped[i * dim_m + d] - v_mapped[j * dim_m + d];
      d2 += diff * diff;
    }
    d2 = std::sqrt(d2);
    if (std::isnan(d2) || d1 == std::isinf(d2))
    {
      throw zentas::zentas_error("inf or nan detected in vdimap");
    }
    
    // on getting the floating point of literals correct:
    // https://stackoverflow.com/questions/22380778/what-is-the-best-way-to-express-a-templated-numeric-literal
    T score = std::abs(d1 - d2) / std::abs(d1 + d2 + T(1e-7));
    T toll = T(1e-5);
    if (score > toll)
    {
      std::stringstream ss;
      ss << "numerical issue in vdimap. ";
      ss << "data, indices " << i << "  and " << j << ".  ";
      ss << "distance using original data : " << d1 << ".  ";
      ss << "distance using transformed data : " << d2 << ".   ";
      ss << "this is a difference of " << std::abs(d1 - d2) << ".  ";
      ss << "and std::abs(d1 - d2) / std::abs(d1 + d2 + 1e-7) is " << score << ".  ";
      ss << "this is greater than the tolerance, " << toll << ".  ";
      throw zentas::zentas_error(ss.str());
    }
  }

  auto t5         = std::chrono::high_resolution_clock::now();
  auto time_total = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t0).count();
  (void)time_total;
}

template void vdimap<float>(std::vector<float>& v_mapped,
                            const float* const  ptr_datain,
                            size_t              ndata,
                            size_t              dimension,
                            size_t              seed);
template void vdimap<double>(std::vector<double>& v_mapped,
                             const double* const  ptr_datain,
                             size_t               ndata,
                             size_t               dimension,
                             size_t               seed);
}
}
