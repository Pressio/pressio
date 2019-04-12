
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../core/src/matrix/core_matrix_meta.hpp"
#include "../../../core/src/multi_vector/core_multi_vector_meta.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    (odeMethod == ::rompp::ode::ImplicitEnum::Euler) and
    (core::meta::is_sparse_matrix_wrapper_eigen<jacobian_type>::value or
#ifdef HAVE_TRILINOS
    core::meta::is_sparse_matrix_wrapper_epetra<jacobian_type>::value or
#endif
    core::meta::is_dense_matrix_wrapper_eigen<jacobian_type>::value)
    > * = nullptr
  >
  void time_discrete_jacobian(jacobian_type & jac,
			      scalar_type dt){
  jac.scale(-dt);
  jac.addToDiagonal(static_cast<scalar_type>(1));
}

template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    (odeMethod == ::rompp::ode::ImplicitEnum::BDF2) and
    core::meta::is_sparse_matrix_wrapper_eigen<jacobian_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt){
  using namespace ::rompp::ode::coeffs;
  jac.scale(-bdf2<scalar_type>::c3_*dt);
  jac.addToDiagonal(static_cast<scalar_type>(1));
}

}}}//end namespace rompp::ode::impl
#endif
