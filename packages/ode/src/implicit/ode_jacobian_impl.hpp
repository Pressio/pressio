
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../core/src/matrix/core_matrix_meta.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template <typename jacobian_type, typename scalar_type,
	  core::meta::enable_if_t<
	    core::meta::is_eigen_sparse_matrix_wrapper<jacobian_type>::value or
	    core::meta::is_epetra_sparse_matrix_wrapper<jacobian_type>::value or
	    core::meta::is_eigen_dense_matrix_wrapper<jacobian_type>::value
	    > * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   scalar_type dt){

  jac.scale(-dt);
  jac.addToDiagonal(static_cast<scalar_type>(1));
}
//---------------------------------------------------------------

template <typename jacobian_type, typename scalar_type,
	  core::meta::enable_if_t<
	    core::meta::is_eigen_sparse_matrix_wrapper<jacobian_type>::value
	    > * = nullptr
	  >
void implicit_bdf2_time_discrete_jacobian(jacobian_type & jac,
					   scalar_type dt){
  
  using namespace ::rompp::ode::impl::coeffs;
  jac.scale(-bdf2<scalar_type>::c3*dt);
  jac.addToDiagonal(static_cast<scalar_type>(1));
}


//--------------------------------------------
//need to add also other overloads
//--------------------------------------------


}}}//end namespace rompp::ode::impl
#endif
