
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type,
	 typename scalar_type,
	 core::meta::enable_if_t<
      ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value == false
	   > * = nullptr
	 >
void implicit_euler_time_discrete_residual(const state_type & yn,
					   const state_type & ynm1,
					   state_type & R,
					   scalar_type dt){
  // On input: R should contain the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // so that on output, it contains the time discrete residual
  R = yn - ynm1 - dt*R;
}
//-------------------------------------------------------


template<typename state_type,
	 typename scalar_type,
	 core::meta::enable_if_t<
	   ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
	   > * = nullptr
	 >
void implicit_euler_time_discrete_residual(const state_type & yn,
					   const state_type & ynm1,
					   state_type & R,
					   scalar_type dt){
  // On input: R should contain the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // so that on output, it contains the time discrete residual
  R.data()->update(1.0, *yn.data(), -1.0, *ynm1.data(), -dt);
}
//-------------------------------------------------------


template<typename state_type,
	 typename scalar_type,
	 core::meta::enable_if_t<
	   !::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
	   > * = nullptr
	 >
void implicit_bdf2_time_discrete_residual(const state_type & yn,
					  const state_type & ynm1,
					  const state_type & ynm2,
					  state_type & R,
					  scalar_type dt){
  // On input: R should contain the application RHS
  using namespace ::rompp::ode::coeffs;
  R = yn - bdf2<scalar_type>::c1_*ynm1
      + bdf2<scalar_type>::c2_*ynm2
      - bdf2<scalar_type>::c3_*dt*R;
}
//-------------------------------------------------------

template<typename state_type,
	 typename scalar_type,
	 core::meta::enable_if_t<
	   ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
	   > * = nullptr
	 >
void implicit_bdf2_time_discrete_residual(const state_type & yn,
					  const state_type & ynm1,
					  const state_type & ynm2,
					  state_type & R,
					  scalar_type dt){
  // On input: R should contain the application RHS
  using namespace ::rompp::ode::coeffs;

  const scalar_type oneSc = static_cast<scalar_type>(1);
  const scalar_type c2 = bdf2<scalar_type>::c2_;
  R.data()->update(c2, *ynm2.data(), -bdf2<scalar_type>::c3_*dt);
  R.data()->update(oneSc, *yn.data(), -bdf2<scalar_type>::c1_, *ynm1.data(), oneSc);
}


}}}//end namespace rompp::ode::impl
#endif
