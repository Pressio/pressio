
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type,
	 typename scalar_type,
	 typename std::enable_if<
	   ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value == false
	   >::type * = nullptr
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
	 typename std::enable_if<
	   ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
	   >::type * = nullptr
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

template<typename state_type, typename scalar_t>
void implicit_bdf2_time_discrete_residual(const state_type & yn,
					  const state_type & ynm1,
					  const state_type & ynm2,
					  state_type & R,
					  scalar_t dt){
  using namespace ::rompp::ode::impl::coeffs;
  R = yn - bdf2<scalar_t>::c1*ynm1 +
      bdf2<scalar_t>::c2*ynm2 -
      bdf2<scalar_t>::c3*dt*R;
}

}}}//end namespace rompp::ode::impl
#endif
