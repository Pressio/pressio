
#ifndef CORE_MATRIX_BASE_MATRIX_SPARSE_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SPARSE_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type>
class MatrixSparseBase
  : private core::details::CrtpBase<
  MatrixSparseBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_sparse==1,
  "OOPS: dense matrix inheriting from sparse base!");

  using traits_t = details::traits<derived_type>;

public:
  size_t nonZerosCount()const{
    return this->underlying().nonZerosCountImpl();}
    
private:  
  friend derived_type;
  friend core::details::CrtpBase<
    MatrixSparseBase<derived_type>>;

  MatrixSparseBase() = default;
  ~MatrixSparseBase() = default;
 
};//end class
  
} // end namespace core
}//end namespace rompp
#endif
