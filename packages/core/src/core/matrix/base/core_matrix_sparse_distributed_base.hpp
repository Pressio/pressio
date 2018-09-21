
#ifndef CORE_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixSparseDistributedBase
  : private core::details::CrtpBase<
  MatrixSparseDistributedBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==0,
  "OOPS: non-distributed matrix inheriting from sparse distributed base!");

  using traits_t = details::traits<derived_type>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  
public:
  void insertGlobalValues(GO_t row,
			  LO_t numEntries,
			  const sc_t * values,
			  const GO_t * indices){
    this->underlying().insertGlobalValuesImpl(row,
					      numEntries,
					      values, indices);
  }

private:  
  friend derived_type;
  friend core::details::CrtpBase<
    MatrixSparseDistributedBase<derived_type>>;

  MatrixSparseDistributedBase() = default;
  ~MatrixSparseDistributedBase() = default; 

};//end class  
} // end namespace core
#endif
