
#ifndef SVD_SOLVER_TRAITS_HPP_
#define SVD_SOLVER_TRAITS_HPP_

#include "svd_forward_declarations.hpp"
#include "core_meta.hpp"
#include "matrix/core_matrix_traits.hpp"

namespace svd{
namespace details{

template <typename matrix_type>
struct traits<Solver<matrix_type,
		     precond_type,
		     typename
		     std::enable_if<
		       core::meta::is_core_matrix_wrapper<matrix_type>::value==true &&
		       core::details::traits<matrix_type>::isEpetra==1 
		       >::type
		     >
{
  using derived_t = Solver<matrix_type, precond_type>;

  using native_matrix_t = typename core::details::traits<matrix_type>::wrapped_t;
  using scalar_t = typename core::details::traits<matrix_type>::scalar_t;
  using lsv_t = core::Matrix<Epetra_MultiVector>;
  using rsv_t = core::Matrix<Epetra_MultiVector>;
};
    
}//end namespace details
}//end namespace svd

#endif



// using wrapped_solver_t = typename std::conditional<which_impl==svdKind::EigenJacobi,
// 						   Eigen::JacobiSVD<native_matrix_t>,
// 						   Eigen::BDCSVD<native_matrix_t>
// 						   >::type;
