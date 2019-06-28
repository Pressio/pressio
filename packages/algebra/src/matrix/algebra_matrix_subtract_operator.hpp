
#ifndef ALGEBRA_MATRIX_MATRIX_SUBTRACT_OPERATOR_HPP_
#define ALGEBRA_MATRIX_MATRIX_SUBTRACT_OPERATOR_HPP_

#include "algebra_matrix_traits.hpp"
#include "algebra_matrix_meta.hpp"

namespace rompp{ namespace algebra{

template <typename T1, 
	  ::rompp::mpl::enable_if_t<
	    algebra::meta::is_dense_matrix_wrapper_eigen<T1>::value or
	    algebra::meta::is_sparse_matrix_wrapper_eigen<T1>::value
	    > * = nullptr>
T1 operator-(const T1 & A, const T1 & B) {
  assert( A.rows() == B.rows() );
  assert( A.cols() == B.cols() );
  T1 C(A.rows(), A.cols());
  *C.data() = *A.data_ - *B.data();
  return C;
}
  

}}//end namespace algebra::rompp
#endif
