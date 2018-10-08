
#ifndef ODE_RESIDUAL_IMPL_HPP_
#define ODE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type, typename scalar_type>
void implicit_euler_time_discrete_residual(const state_type & yn,
					   const state_type & ynm1,
					   state_type & R,
					   scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // On output, it contains the residual
  R = yn - ynm1 - dt*R;
  // for (decltype(R.size()) i=0; i < R.size(); i++){
  //   R[i] = yn[i] - ynm1[i] - dt*R[i];
  // }
}
//-------------------------------------------------------


template<typename state_type, typename scalar_type>
void implicit_bdf2_time_discrete_residual(const state_type & ynp1,
					  const state_type & yn,
					  const state_type & ynm1,
					  state_type & R,
					  scalar_type dt){
  
  constexpr scalar_type c1 = static_cast<scalar_type>(4.)/3.;
  constexpr scalar_type c2 = static_cast<scalar_type>(1.)/3.;
  constexpr scalar_type c3 = static_cast<scalar_type>(2.)/3.;
  
  // On input: R contains the application RHS, i.e. if
  //           dudt = f(x,u,...), R contains f(...)
  // On output, it contains the residual
  R = ynp1 - c1*yn +c2*ynm1 - c3*dt*R;
}

      
}}}//end namespace rompp::ode::impl
#endif 
