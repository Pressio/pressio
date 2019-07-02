
#ifndef ODE_JACOBIAN_IMPL_HPP_
#define ODE_JACOBIAN_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/matrix/containers_matrix_meta.hpp"
#include "../../../containers/src/multi_vector/containers_multi_vector_meta.hpp"
#include "ode_implicit_constants.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    (odeMethod == ::pressio::ode::ImplicitEnum::Euler) and
    (containers::meta::is_sparse_matrix_wrapper_eigen<jacobian_type>::value or
#ifdef HAVE_TRILINOS
     containers::meta::is_sparse_matrix_wrapper_epetra<jacobian_type>::value or
#endif
    containers::meta::is_dense_matrix_wrapper_eigen<jacobian_type>::value)
    > * = nullptr
  >
  void time_discrete_jacobian(jacobian_type & jac,
			      scalar_type dt){
  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  jac.scale(-dt);
  jac.addToDiagonal(one);
}

#ifdef HAVE_PYBIND11
template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    (odeMethod == ::pressio::ode::ImplicitEnum::Euler) and
    containers::meta::is_cstyle_array_pybind11<jacobian_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt){
  using namespace ::pressio::ode::coeffs;
  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();

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
  ::pressio::mpl::enable_if_t<
    (odeMethod == ::pressio::ode::ImplicitEnum::BDF2) and
    containers::meta::is_sparse_matrix_wrapper_eigen<jacobian_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt){
  using namespace ::pressio::ode::coeffs;
  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  jac.scale(-bdf2<scalar_type>::c3_*dt);
  jac.addToDiagonal(one);
}


}}}//end namespace pressio::ode::impl
#endif
