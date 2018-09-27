
#ifndef CORE_MATRIX_MATRIX_SUBTRACT_OPERATOR_HPP_
#define CORE_MATRIX_MATRIX_SUBTRACT_OPERATOR_HPP_

#include "core_matrix_traits.hpp"
#include "core_matrix_meta.hpp"

namespace rompp{
namespace core{

  
template <typename T1, 
	  core::meta::enable_if_t<
	    core::meta::is_eigen_dense_matrix_wrapper<T1>::value or
	    core::meta::is_eigen_sparse_matrix_wrapper<T1>::value
	    > * = nullptr>
auto operator-(const T1 & A, const T1 & B) {
  assert( A.rows() == B.rows() );
  assert( A.cols() == B.cols() );
  T1 C(A.rows(), A.cols());
  *C.data() = *A.data_ - *B.data();
  return C;
}
  

}//end namespace core
}//end namespace rompp
#endif
