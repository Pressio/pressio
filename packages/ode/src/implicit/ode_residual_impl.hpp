
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<::rompp::ode::ImplicitEnum odeMethod,
	  int numStates,
	  typename state_type,
	  typename scalar_type,
	  core::meta::enable_if_t<
  (odeMethod == ::rompp::ode::ImplicitEnum::Euler)
#ifdef HAVE_TRILINOS  
   and ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value == false
#endif	
  > * = nullptr
	  >
void implicit_time_discrete_residual(const state_type & yn,
				     const std::array<state_type,numStates> & ynm,
				     state_type & R,
				     scalar_type dt){
  // On input: R should contain the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  R = yn - ynm[0] - dt*R;
}
//-------------------------------------------------------

template<::rompp::ode::ImplicitEnum odeMethod,
	  int numStates,
	  typename state_type,
	  typename scalar_type,
	  core::meta::enable_if_t<
  (odeMethod == ::rompp::ode::ImplicitEnum::BDF2) 
#ifdef HAVE_TRILINOS
 and ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value == false
#endif 
  > * = nullptr
	 >
void implicit_time_discrete_residual(const state_type & yn,
				     const std::array<state_type,numStates> & ynm,
				     state_type & R,
				     scalar_type dt){
  // On input: R should contain the application RHS
  using namespace ::rompp::ode::coeffs;
  R = yn - bdf2<scalar_type>::c1_*ynm[1]
      + bdf2<scalar_type>::c2_*ynm[0]
      - bdf2<scalar_type>::c3_*dt*R;
}
//-------------------------------------------------------



#ifdef HAVE_TRILINOS 
template<::rompp::ode::ImplicitEnum odeMethod,
	  int numStates,
	  typename state_type,
	  typename scalar_type,
	  core::meta::enable_if_t<
  (odeMethod == ::rompp::ode::ImplicitEnum::Euler)
  and ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
  > * = nullptr
	  >
void implicit_time_discrete_residual(const state_type & yn,
				     const std::array<state_type,numStates> & ynm,
				     state_type & R,
				     scalar_type dt){
  // On input: R should contain the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  R.data()->update(1.0, *yn.data(), -1.0, *ynm[0].data(), -dt);
}
//-------------------------------------------------------


template< ::rompp::ode::ImplicitEnum odeMethod,
	  int numStates,
	  typename state_type,
	  typename scalar_type,
	  core::meta::enable_if_t<
  (odeMethod == ::rompp::ode::ImplicitEnum::BDF2) 
 and ::rompp::core::meta::is_tpetra_vector_wrapper<state_type>::value
  > * = nullptr
	 >
void implicit_time_discrete_residual(const state_type & yn,
				     const std::array<state_type,numStates> & ynm,
				     state_type & R,
				     scalar_type dt){
  // // On input: R should contain the application RHS
  using namespace ::rompp::ode::coeffs;
  // R.data()->update(bdf2<scalar_type>::c2_,
		//    *ynm[1].data(),
		//    -bdf2<scalar_type>::c3_*dt);
  const scalar_type oneSc = static_cast<scalar_type>(1);
  const scalar_type c2 = bdf2<scalar_type>::c2_;
  R.data()->update(c2, *ynm[0].data(), -bdf2<scalar_type>::c3_*dt);
  R.data()->update(oneSc, *yn.data(), -bdf2<scalar_type>::c1_, *ynm[1].data(), oneSc);
}
#endif

}}}//end namespace rompp::ode::impl
#endif
