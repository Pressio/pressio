
#ifndef ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_
#define ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_

#include "../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace rompp{ namespace rom{ namespace impl{

#ifdef HAVE_TRILINOS
template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  typename decoder_jac_type,
  core::meta::enable_if_t<
    odeMethod == ode::ImplicitEnum::Euler and
    core::meta::is_multi_vector_wrapper_epetra<jacobian_type>::value and
    core::meta::is_multi_vector_wrapper_epetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac, scalar_type dt,
			    const decoder_jac_type & phi){
  jac.scale(-dt);
  jac += phi;
}

template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  typename decoder_jac_type,
  core::meta::enable_if_t<
    odeMethod == ode::ImplicitEnum::Euler and
    core::meta::is_multi_vector_wrapper_tpetra<jacobian_type>::value and
    core::meta::is_multi_vector_wrapper_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
  void time_discrete_jacobian(jacobian_type & jac,
			      scalar_type dt,
			      const decoder_jac_type & phi){
  jac.scale(-dt);
  jac += phi;
}


template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  typename decoder_jac_type,
  core::meta::enable_if_t<
    (odeMethod == ode::ImplicitEnum::BDF2) and
    core::meta::is_multi_vector_wrapper_epetra<jacobian_type>::value and
    core::meta::is_multi_vector_wrapper_epetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt,
			    const decoder_jac_type & phi){
  jac.scale(-ode::coeffs::bdf2<scalar_type>::c3_*dt);
  jac += phi;
}


template <
  ode::ImplicitEnum odeMethod,
  typename jacobian_type,
  typename scalar_type,
  typename decoder_jac_type,
  core::meta::enable_if_t<
    (odeMethod == ode::ImplicitEnum::BDF2) and
    core::meta::is_multi_vector_wrapper_tpetra<jacobian_type>::value and
    core::meta::is_multi_vector_wrapper_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(jacobian_type & jac,
			    scalar_type dt,
			    const decoder_jac_type & phi){
  jac.scale(-ode::coeffs::bdf2<scalar_type>::c3_*dt);
  jac += phi;
}
#endif

}}}//end namespace rompp::rom::impl
#endif
