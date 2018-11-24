
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../core/src/matrix/core_matrix_meta.hpp"
#include "../../../core/src/multi_vector/core_multi_vector_meta.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template <typename jacobian_type, typename scalar_type,
	  core::meta::enable_if_t<
	    core::meta::is_eigen_sparse_matrix_wrapper<jacobian_type>::value or
	    #ifdef HAVE_TRILINOS
	    	core::meta::is_epetra_sparse_matrix_wrapper<jacobian_type>::value or
	    #endif
	    core::meta::is_eigen_dense_matrix_wrapper<jacobian_type>::value
	    > * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   scalar_type dt){

  jac.scale(-dt);
  jac.addToDiagonal(static_cast<scalar_type>(1));
}
//---------------------------------------------------------------


#ifdef HAVE_TRILINOS
template <typename jacobian_type, typename scalar_type,
	  typename basis_type,
	  core::meta::enable_if_t<
    core::meta::is_epetra_multi_vector_wrapper<jacobian_type>::value and
    core::meta::is_epetra_multi_vector_wrapper<basis_type>::value
	    > * = nullptr
	  >
void implicit_euler_time_discrete_jacobian(jacobian_type & jac,
					   scalar_type dt,
					   const basis_type & phi){

  jac.scale(-dt);
  jac += phi;
}
#endif
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
