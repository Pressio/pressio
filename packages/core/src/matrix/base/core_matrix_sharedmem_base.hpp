
#ifndef CORE_MATRIX_BASE_MATRIX_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixSharedMemBase
  : private core::details::CrtpBase<MatrixSharedMemBase<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==0 &&
		 details::traits<derived_type>::isSharedMem==1,
		 "OOPS: distributed matrix inheriting \
from sharedMem base!");
  
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;  

public:
  ord_t rows() const{
    return this->underlying().rowsImpl();
  }

  ord_t cols() const{
    return this->underlying().colsImpl();
  }

  void resize(ord_t nrows, ord_t ncols){
    this->underlying().resizeImpl(nrows, ncols);
  }

  void addToDiagonal(sc_t value) {
    return this->underlying().addToDiagonalImpl(value);
  };
  
private:  
  friend derived_type;
  friend core::details::CrtpBase<MatrixSharedMemBase<derived_type>>;
  MatrixSharedMemBase() = default;
  ~MatrixSharedMemBase() = default;
  
};//end class
  
} // end namespace core
#endif
