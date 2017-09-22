// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <random>
#include <zentas/vdimap.hpp>

int main()
{

  std::vector<float> v_mapped;
  std::vector<float> v_data;

  size_t ndata     = 100000;
  size_t dimension = 25;
  size_t seed      = 3;

  std::vector<float> eigenvalues        = {11, 12, 13, 14, 15};
  size_t             n_true_eigenvalues = eigenvalues.size();
  std::vector<float> eigenvectors(n_true_eigenvalues * dimension);

  std::uniform_int_distribution<size_t> dis;
  std::default_random_engine            gen(seed);

  /* populate the eigenvectors */
  for (unsigned e = 0; e < n_true_eigenvalues; ++e)
  {
    for (unsigned d = 0; d < dimension; ++d)
    {
      int rv                          = dis(gen) % 10000 - 5000;
      eigenvectors[e * dimension + d] = rv;
    }
  }

  nszen::vdimap::make_all_nperp(eigenvectors.data(), dimension, n_true_eigenvalues);

  auto get_mp1 = [&dis, &gen]() { return (2. / 100) * (dis(gen) % 100) - 1.; };

  v_data.resize(ndata * dimension);
  for (unsigned i = 0; i < ndata; ++i)
  {
    for (unsigned d = 0; d < dimension; ++d)
    {
      v_data[i * dimension + d] = 0.001f * static_cast<float>(get_mp1());
    }

    for (unsigned e = 0; e < n_true_eigenvalues; ++e)
    {
      auto coeff = std::sqrt(3) * get_mp1();
      for (unsigned d = 0; d < dimension; ++d)
      {
        v_data[i * dimension + d] +=
          coeff * eigenvectors[e * dimension + d] * std::sqrt(eigenvalues[e]);
      }
    }
  }

  std::cout << "entering nszen::vdimap::vdimap<float> " << std::endl;
  nszen::vdimap::vdimap<float>(v_mapped, v_data.data(), ndata, dimension, seed);
}
