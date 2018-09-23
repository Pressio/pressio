
#ifndef CORE_MATRIX_BASE_MATRIX_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type>
class MatrixSharedMemBase
  : private core::details::CrtpBase<
  MatrixSharedMemBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed matrix inheriting from sharedMem base!");
  
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;  

public:
  ord_t rows() const{
    return this->underlying().rowsImpl();}

  ord_t cols() const{
    return this->underlying().colsImpl();}
  
private:  
  friend derived_type;
  friend core::details::CrtpBase<MatrixSharedMemBase<derived_type>>;
  MatrixSharedMemBase() = default;
  ~MatrixSharedMemBase() = default;
  
};//end class
  
} // end namespace core
}//end namespace rompp
#endif
