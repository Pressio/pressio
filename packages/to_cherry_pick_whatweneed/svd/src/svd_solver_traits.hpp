
#ifndef SVD_SOLVER_TRAITS_HPP_
#define SVD_SOLVER_TRAITS_HPP_

#include "svd_forward_declarations.hpp"
#include "core_meta.hpp"
#include "matrix/core_matrix_traits.hpp"
#include <Eigen/SVD>


namespace svd{
namespace details{

  //*******************************
  // eigen svd
  //*******************************   
  template <typename matrix_type,
	    svdKind which_impl>
  struct traits< solver<matrix_type, which_impl,
			typename
			std::enable_if< core::details::traits<matrix_type>::isEigen==1 &&
					(which_impl == svdKind::EigenJacobi ||
					 which_impl == svdKind::EigenBDCSVD)
					>::type
			>
		 >
  {
    using derived_t = solver<matrix_type, which_impl>;

    using native_matrix_t = typename core::details::traits<matrix_type>::wrapped_t;

    using wrapped_solver_t = typename std::conditional<which_impl==svdKind::EigenJacobi,
						Eigen::JacobiSVD<native_matrix_t>,
						Eigen::BDCSVD<native_matrix_t>
    					     >::type;
    using scalar_t = typename core::details::traits<matrix_type>::scalar_t;
    using u_matrix_t = core::Matrix<native_matrix_t>;
    using v_matrix_t = core::Matrix<native_matrix_t>;
    
    enum {
      isSerial = 1,
      isDistributed = !isSerial,
      isEigen = 1
    };

    svdKind method = which_impl;
  };
  
  
}//end namespace details
}//end namespace svd

#endif
