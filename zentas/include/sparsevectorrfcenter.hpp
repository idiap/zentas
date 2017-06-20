#ifndef ZENTAS_SPARSEVECTORRFCENTER_HPP
#define ZENTAS_SPARSEVECTORRFCENTER_HPP

template <typename TAtomic>
struct SparseVectorRfCenter{ /* Sparse vector sample */
  // TODO : a map of some sort.

};


template <typename TAtomic>
struct SparseVectorSample{
  public:
    size_t size;
    const TAtomic * const values;
    const size_t * const indices;
    SparseVectorSample(size_t size, const TAtomic * const values, const size_t * const indices): size(size), values(values), indices(indices) {}
};



#endif
