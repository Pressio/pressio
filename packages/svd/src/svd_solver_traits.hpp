
#ifndef SVD_SOLVER_TRAITS_HPP_
#define SVD_SOLVER_TRAITS_HPP_

#include "svd_forward_declarations.hpp"
#include "matrix/core_matrix_traits.hpp"
//#include meta/core_matrix_meta.hpp

namespace svd{
namespace details{

template<typename T, typename enable = void>
struct svd_traits{
  using derived_t = void;
  using native_matrix_t = void;
  using scalar_t = void;
  using lsv_t = void;
  using rsv_t = void;
  
};

template<typename T> 
struct svd_traits<const T> : svd_traits<T> {};

  
template <typename matrix_type,
	  template <typename...> class lsv_type,
	  template <typename...> class rsv_type>
struct svd_traits<Solver<matrix_type,
			 lsv_type,
			 rsv_type,
			 typename
			 std::enable_if<
			   core::meta::is_core_matrix_wrapper<
			     matrix_type>::value==true &&
			   core::details::traits<matrix_type>::isEpetra==1 
			   >::type
			 >
		  >
{
  using derived_t = Solver<matrix_type,
			   lsv_type,
			   rsv_type>;

  using matrix_t = matrix_type;
  using native_matrix_t = typename core::details::traits<matrix_type>::wrapped_t;
  using scalar_t = typename core::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<Epetra_MultiVector>;
  using rsv_t = rsv_type<Epetra_MultiVector>;
};

  
}//end namespace details
}//end namespace svd

#endif



// using wrapped_solver_t = typename std::conditional<which_impl==svdKind::EigenJacobi,
// 						   Eigen::JacobiSVD<native_matrix_t>,
// 						   Eigen::BDCSVD<native_matrix_t>
// 						   >::type;
