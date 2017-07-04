#ifndef ZENTAS_SPARSEVECTORRFCENTER_HPP
#define ZENTAS_SPARSEVECTORRFCENTER_HPP

#include <unordered_map>

// TODO : currently computing distances involiving SparseVectorRfCenters is slow
// it can probably be optimised is a hand-crafted data structure is used to store sparse centers

template <typename TAtomic>
using SparseVectorRfCenter = std::unordered_map<size_t, TAtomic>;


template <typename TAtomic>
struct SparseVectorSample{
  public:
    size_t size;
    const TAtomic * const values;
    const size_t * const indices;
    SparseVectorSample(size_t size, const TAtomic * const values, const size_t * const indices): size(size), values(values), indices(indices) {}
};



#endif
