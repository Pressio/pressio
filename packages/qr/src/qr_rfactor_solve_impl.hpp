
#ifndef QR_RFACTOR_SOLVE_HPP_
#define QR_RFACTOR_SOLVE_HPP_

#include "qr_ConfigDefs.hpp"
#include "../../CORE_VECTOR"
#include "../../CORE_MATRIX"
#ifdef HAVE_TRILINOS
#include "Teuchos_SerialDenseMatrix.hpp"
#endif

namespace rompp{ namespace qr{ namespace impl{

// #ifdef HAVE_TRILINOS
// template <typename T, typename enable = void>
// struct is_teuchos_serial_dense_matrix : std::false_type{};

// template <typename T>
// struct is_teuchos_serial_dense_matrix<
//   T,
//   typename std::enable_if<
//     std::is_same<T,
// 		 Teuchos::SerialDenseMatrix<
// 		   typename T::ordinalType,
// 		   typename T::scalarType>
// 		 >::value
//     >::type
//   > : std::true_type{};
// #endif


template<typename vector_type,
	 typename R_type,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_eigen_vector_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_eigen_dense_matrix_wrapper<R_type>::value
	   > * = nullptr>
void solve(const vector_type & rhs,
	   const R_type & Rmatrix,
	   vector_type & y){
  auto nat_rhs = rhs.data();
  auto nat_y = y.data();
  auto nat_R = Rmatrix.data();

  *nat_y = (*nat_R).template triangularView<Eigen::Upper>().solve(*nat_rhs);
}


}}}//end namespace rompp::qr::impl
#endif
