
#ifndef CORE_MATRIX_MATRIX_TIMES_OPERATOR_HPP_
#define CORE_MATRIX_MATRIX_TIMES_OPERATOR_HPP_

#include "core_matrix_traits.hpp"
#include "core_matrix_meta.hpp"

namespace rompp{
namespace core{

  
template <typename T1, 
	  ::rompp::mpl::enable_if_t<
	    core::meta::is_dense_matrix_wrapper_eigen<T1>::value or
	    core::meta::is_sparse_matrix_wrapper_eigen<T1>::value
	    > * = nullptr>
T1 operator*(const T1 & A, const T1 & B) {
  assert( A.rows() == B.rows() );
  assert( A.cols() == B.cols() );
  T1 C(A.rows(), A.cols());
  *C.data() = (*A.data_) * (*B.data());
  return C;
}


}//end namespace core
}//end namespace rompp
#endif
