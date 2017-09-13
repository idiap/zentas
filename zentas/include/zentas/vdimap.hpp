// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_VDIMAP_HPP
#define ZENTAS_VDIMAP_HPP

#include <iostream>
#include <vector>

namespace nszen
{
namespace vdimap
{

template <typename T>
void make_all_nperp(T* const evs, size_t dim, size_t neigs);

template <typename T>
void vdimap(std::vector<T>& v_mapped,
            const T* const  ptr_datain,
            size_t          ndata,
            size_t          dimension,
            size_t          seed);
}
}

#endif
