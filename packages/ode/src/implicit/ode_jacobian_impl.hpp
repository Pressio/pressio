
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
  constexpr auto one = ::rompp::core::constants::one<scalar_type>();
  jac.scale(-dt);
  jac.addToDiagonal(one);
}

#ifdef HAVE_PYBIND11
template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    (odeMethod == ::rompp::ode::ImplicitEnum::Euler) and
    core::meta::is_cstyle_array_pybind11<jacobian_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt){
  using namespace ::rompp::ode::coeffs;
  constexpr auto one = ::rompp::core::constants::one<scalar_type>();

  if (jac.ndim() != 2)
    throw std::runtime_error("Tensors with dim>2 not supported");

  double *ptr = jac.mutable_data();
  const size_t rows = jac.shape()[0];
  const size_t cols = jac.shape()[1];

  for (size_t irow = 0; irow < rows; irow++){
    for (size_t icol = 0; icol < cols; icol++){
      ptr[irow*cols + icol] *= -dt;
      if (irow == icol and irow == 2)
	ptr[irow*cols + icol] += one;
    }
  }

  // no need to reshape array since it has already right shape
}
#endif


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
  constexpr auto one = ::rompp::core::constants::one<scalar_type>();
  jac.scale(-bdf2<scalar_type>::c3_*dt);
  jac.addToDiagonal(one);
}


}}}//end namespace rompp::ode::impl
#endif
