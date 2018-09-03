
#ifndef CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixDenseSharedMemBase
  : private core::details::CrtpBase<MatrixDenseSharedMemBase<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==0 &&
		 details::traits<derived_type>::isSharedMem==1,
		 "OOPS: distributed matrix inheriting \
from dense sharedMem base!");
  
// private:
//   using sc_t = typename details::traits<derived_type>::scalar_t;
//   using der_t = typename details::traits<derived_type>::derived_t;
//   using wrap_t = typename details::traits<derived_type>::wrapped_t;
//   using ord_t = typename details::traits<derived_type>::ordinal_t;  

// public:
//   // the () operators below are placed here becasue they serve the
//   // purpose of scripting operator for matrices, but this is just meant
//   // to be for DENSE SHAREDMEM matrices. So we leave them here.
//   // Every dense sharedMem matrix should define these.
//   sc_t & operator() (ord_t row, ord_t col);
//   sc_t const & operator() (ord_t row, ord_t col) const;

private:  
  friend derived_type;
  friend core::details::CrtpBase<MatrixDenseSharedMemBase<derived_type>>;

  MatrixDenseSharedMemBase() = default;
  ~MatrixDenseSharedMemBase() = default;
  
};//end class
  
} // end namespace core
#endif
