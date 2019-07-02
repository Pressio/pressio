
#ifndef CONTAINERS_MATRIX_MATRIX_PLUS_OPERATOR_HPP_
#define CONTAINERS_MATRIX_MATRIX_PLUS_OPERATOR_HPP_

#include "containers_matrix_traits.hpp"
#include "containers_matrix_meta.hpp"

namespace pressio{ namespace containers{
  
template <typename T1, 
	  ::pressio::mpl::enable_if_t<
	    containers::meta::is_dense_matrix_wrapper_eigen<T1>::value or 
	    containers::meta::is_sparse_matrix_wrapper_eigen<T1>::value
	    > * = nullptr>
T1 operator+(const T1 & A, const T1 & B) {
  assert( A.rows() == B.rows() );
  assert( A.cols() == B.cols() );
  T1 C(A.rows(), A.cols());
  *C.data() = *A.data() + *B.data();
  return C;
}
  
}}//end namespace containers::pressio
#endif
